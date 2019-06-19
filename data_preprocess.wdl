
## This is for clean fastq reads processing
## in oder to get analysis-ready bam files



# TASK DEFINITION PART

## Generate a file array for trim
task ArrayCleanData {
  String CleanDataDir
  
  command <<<
    ls ${CleanDataDir}/*.fq.gz | \
    xargs -n2 | tr ' ' '\t' > data4convertionlist
  >>>
  output {
    File outlist = "data4convertionlist"
  }
}

## Convert clean fastq reads to ubam
task Fastq2ubam {
  File cleanFq_R1
  File cleanFq_R2
  String picard_path
  String? platform_override
  String platform=select_first([platform_override, "COMPLETE"])
  String library=basename(cleanFq_R1, "")
  String samplename
  String readgroup
  String out_ubam_basename=basename(cleanFq_R1, "_1.fq.gz")
  
  command <<<
    ${picard_path} --javaOptions "-Xms5g -Xmx8g" FastqToSam \
    F1=${cleanFq_R1} \
    F2=${cleanFq_R2} \
    O=${out_ubam_basename}.ubam \
    PL=${platform} \
    SM=${samplename} \
    LB=${library} \
    RG=${readgroup} \
    SEQUENCING_CENTER=BGI \
    RUN_DATE=`date +"%y-%m-%d_%H:%M:%S"`

  >>>
  output {
    File out_ubam = "${out_ubam_basename}.ubam"
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
  String output_bam_basename=basename(in_ubam, ".ubam")
  String picard_path
  String samtools_path
 

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit), 
  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy 
  # references such as b37 and hg19.
  File? ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fa=${ref_fasta}

    ${picard_path} --javaOptions "-Xms3g -Xmx8g" \
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
    
    File stringtie = '/services/tools/stringtie/1.3.3b/stringtie'
    File annotation = '/home/projects/pr_99006/data/rna/genome/hg38_tran/hg38_ucsc.annotated.filtered.gtf'
    
    String CleanDataDir = "/home/projects/pr_99006/data/rna/hiPSC-derived_NPC_Neuron/data4Trim"
    String bwa_path
    String picard_path
    String samtools_path

    call ArrayCleanData {
      input:
        CleanDataDir = CleanDataDir
    }

    Array[Array[File]] inputSamples = read_tsv(ArrayCleanData.outlist)
    scatter (fastq in inputSamples) {
      call Fastq2ubam {
        input:
          cleanFq_R1 = fastq[0],
          cleanFq_R2 = fastq[1]
      }

    }

    Array[File] ubam_files = Fastq2ubam.out_ubam
    scatter (ubam in ubam_files) {
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
          in_ubam = ubam

      }
    }
}













