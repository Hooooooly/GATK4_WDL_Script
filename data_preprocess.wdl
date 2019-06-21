
## This is for clean fastq reads processing
## in oder to get analysis-ready bam files



# TASK DEFINITION PART

## Generate a file array for trim
task ArrayCleanData {
  String CleanDataDir
  
  command <<<
    set -e
    set -o pipefail

    ls ${CleanDataDir}/*.fq.gz | \
    xargs -n2 | tr ' ' '\t' > data4convertionlist
    samplebasename=$(basename $(head -1 data4convertionlist | cut -f1) _1.fq.gz)
    echo $samplebasename | cut -d "-" -f1
    echo $samplebasename | cut -d "_" -f1-2
  >>>
  output {
    File outlist = "data4convertionlist"
    Array[String] sampleinfo = read_lines(stdout())
  }
}

## Convert clean fastq reads to ubam
task Fastq2ubam {
  File cleanFq_R1
  File cleanFq_R2
  String picard_path
  String? platform_override
  String platform=select_first([platform_override, "COMPLETE"])
  Array[String] sampleinfo
  String pattern1=sampleinfo[0]+"-"
  String pattern2=sampleinfo[1]+"_"
  String out_ubam_basename=basename(cleanFq_R1, "_1.fq.gz")
  String readgroup=sub(out_ubam_basename, pattern1, "")
  String library=sub(sub(out_ubam_basename, pattern2, ""), "-"+readgroup, "")
  String samplename

  command <<<
    ${picard_path} -Xms3g -Xmx8g FastqToSam \
    F1=${cleanFq_R1} \
    F2=${cleanFq_R2} \
    O=${out_ubam_basename}_unmapped.bam \
    PL=${platform} \
    SM=${samplename} \
    LB=${library} \
    RG=${readgroup} \
    SEQUENCING_CENTER=BGI \
    RUN_DATE=`date +"%y-%m-%dT%H:%M:%S"`

  >>>
  output {
    File out_ubam = "${out_ubam_basename}_unmapped.bam"
  }
}

# Collect sequencing yield quality metrics
task CollectQualityYieldMetrics {
  File in_ubam
  String metrics_filename
  String picard_path

  command <<<
   ${picard_path} -Xmx2g \
      CollectQualityYieldMetrics \
      INPUT=${in_ubam} \
      OQ=true \
      OUTPUT=${metrics_filename}
  >>>
  output {
    File metrics = "${metrics_filename}"
  }
}


task SamToFastqAndBwaMem {
  File in_ubam
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  # -K process K input bases in each batch regardless of nThreads
  # -v verbosity level: 3 message
  # -t threads
  # -p 如果参数 –p 被设定，那么， mem 命令会认为 read.fq 中的 第 2i-th 和
  # 第 (2i + 1)-th 的 reads 组成一个 read 对 （a read pair），这种方式也
  # 被成为交错式的（interleaved paired-end)。
  # -Y Use soft clipping CIGAR operation for supplementary alignments.
  # By default, BWA-MEM uses soft clipping for the primary alignment
  # and hard clipping for supplementary alignments.
  String bwa_path  
  String bwa_commandline=bwa_path + " mem -K 100000000 -p -v 3 -t 10 -Y $bash_ref_fasta"
  String output_bam_basename=basename(in_ubam, "_unmapped.bam")
  String picard_path
  String samtools_path
 

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit), 
  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy 
  # references such as b37 and hg19.
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}

    ${picard_path} -Xms3g -Xmx8g \
    SamToFastq \
    INPUT=${in_ubam} \
    FASTQ=/dev/stdout \
    INTERLEAVE=true \
    NON_PF=true | \
    ${bwa_commandline} /dev/stdin -  2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) | \
    ${samtools_path} view -@ 6 -1 - > ${output_bam_basename}.bam

  >>>
  output {
    File output_bam = "${output_bam_basename}.bam"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
  }
}

# Merge original input uBAM file with BWA-aligned BAM file
# PROGRAM_RECORD_ID, PROGRAM_GROUP_VERSION, PROGRAM_GROUP_COMMAND_LINE,
# PROGRAM_GROUP_NAME must all be included or none of them
task MergeBamAlignment {
  File unmapped_bam
  File aligned_bam
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  String bwa_path
  String gatk4_path
  String out_bam_basename=basename(unmapped_bam, "_unmapped.bam")
  String bwa_commandline=bwa_path + " mem -K 100000000 -p -v 3 -t 10 -Y $bash_ref_fasta"

  command <<<
   # set the bash variable needed for the command-line
   bash_ref_fa=${ref_fasta}
   ${gatk4_path} --java-options "-Dsamjdk.compression_level=4 -Xms3g" \
     MergeBamAlignment \
     --VALIDATION_STRINGENCY=SILENT \
     --EXPECTED_ORIENTATIONS=FR \
     --ATTRIBUTES_TO_RETAIN=X0 \
     --ALIGNED_BAM=${aligned_bam} \
     --UNMAPPED_BAM=${unmapped_bam} \
     --OUTPUT=${out_bam_basename}.bam \
     --REFERENCE_SEQUENCE=${ref_fasta} \
     --PAIRED_RUN=true \
     --SORT_ORDER="unsorted" \
     --IS_BISULFITE_SEQUENCE=false \
     --ALIGNED_READS_ONLY=false \
     --CLIP_ADAPTERS=false \
     --MAX_RECORDS_IN_RAM=2000000 \
     --ADD_MATE_CIGAR=true \
     --MAX_INSERTIONS_OR_DELETIONS=-1 \
     --PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
     --PROGRAM_RECORD_ID "bwamem" \
     --PROGRAM_GROUP_VERSION "0.7.17-r1188" \
     --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
     --PROGRAM_GROUP_NAME "bwamem" \
     --UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
     --ALIGNER_PROPER_PAIR_FLAGS=true \
     --UNMAP_CONTAMINANT_READS=true

  >>>
  output {
    File out_bam = "${out_bam_basename}.bam"
  }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  Array[File] input_bams
  Array[String] sampleinfo

  String output_bam_basename = sampleinfo[0] + ".aligned.unsorted.duplicates_marked"
  String metrics_filename = sampleinfo[0] + ".duplicate_metrics"
  String gatk4_path

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command <<<
    ${gatk4_path} --java-options "-Xms5g" \
      MarkDuplicates \
      --INPUT ${sep=' --INPUT ' input_bams} \
      --OUTPUT ${output_bam_basename}.bam \
      --METRICS_FILE ${metrics_filename} \
      --VALIDATION_STRINGENCY SILENT \
      --READ_NAME_REGEX=null \
      --ASSUME_SORT_ORDER "queryname" \
      --CREATE_MD5_FILE true
   >>>
   output {
     File output_bam = "${output_bam_basename}.bam"
     File duplicate_metrics = "${metrics_filename}"
   }
}



# WORKFLOW DEFINITION 

workflow PreProcessing4VariantDiscovery_GATK4 {
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
    
  String CleanDataDir
  String bwa_path
  String picard_path
  String samtools_path
  String gatk4_path

  call ArrayCleanData {
    input:
      CleanDataDir = CleanDataDir
  }

  Array[Array[File]] inputSamples = read_tsv(ArrayCleanData.outlist)
  scatter (fastq in inputSamples) {
    call Fastq2ubam {
      input:
        cleanFq_R1 = fastq[0],
        cleanFq_R2 = fastq[1],
        sampleinfo = ArrayCleanData.sampleinfo,
        samplename = "AT2",
        picard_path = picard_path
    }

    call SamToFastqAndBwaMem {
      input:
        bwa_path = bwa_path,
        samtools_path = samtools_path,
        picard_path = picard_path,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_bwt = ref_bwt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        in_ubam = Fastq2ubam.out_ubam

    }

    # QC the unmapped BAM
    call CollectQualityYieldMetrics {
      input:
        picard_path=picard_path,
        in_ubam = Fastq2ubam.out_ubam,
        metrics_filename = basename(Fastq2ubam.out_ubam, ".bam") + ".quality_yield_metrics"
    }

    call MergeBamAlignment {
      input:
        bwa_path = bwa_path,
        unmapped_bam = Fastq2ubam.out_ubam,
        aligned_bam = SamToFastqAndBwaMem.output_bam,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        gatk4_path = gatk4_path
    }
  }

  call MarkDuplicates {
    input:
      input_bams = MergeBamAlignment.out_bam,
      sampleinfo = ArrayCleanData.sampleinfo,
      gatk4_path = gatk4_path
  }

}




