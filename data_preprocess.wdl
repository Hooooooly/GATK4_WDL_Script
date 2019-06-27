## This is for clean fastq reads processing
## in oder to get analysis-ready bam files


##########################
## TASK DEFINITION PART ##
##########################

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
    Array[Array[File]] outlist = read_tsv("data4convertionlist")
    Array[String] sampleinfo = read_lines(stdout())
  }
}

## Convert clean fastq reads to ubam
task Fastq2ubam {
  File cleanFq_R1
  File cleanFq_R2
  String picard_path
  String? platform_override
  String platform=select_first([platform_override, "Complete"])
  Array[String] sampleinfo
  String pattern1=sampleinfo[0]+"-"
  String pattern2=sampleinfo[1]+"_"
  String out_ubam_basename=basename(cleanFq_R1, "_1.fq.gz")
  String readgroup=sub(out_ubam_basename, pattern1, "")
  String library=sub(sub(out_ubam_basename, pattern2, ""), "-"+readgroup, "")
  String samplename

  command <<<
    ${picard_path} -Xms3g -Xmx16g FastqToSam \
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
   ${picard_path} -Xmx6g \
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

    ${picard_path} -Xms3g -Xmx16g \
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

# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
  File input_bam
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  String output_bam_basename
  String gatk4_path

  command {
    set -o pipefail

    ${gatk4_path} --java-options "-Xms5g" \
      SortSam \
      --INPUT ${input_bam} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
    | \
    ${gatk4_path} --java-options "-Xms3g" \
      SetNmMdAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT ${output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ${ref_fasta}
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  File ref_dict  

  String python="/home/tang-lab/biotools/anaconda3/bin/python"

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter. 
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    ${python} <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  File input_bam
  File input_bam_index
  File dbSNP_vcf
  File dbSNP_vcf_index
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  Array[String] sequence_group_interval
  String gatk4_path
  String recalibration_report_filename

  command {
    ${gatk4_path} --java-options "-Xms5g" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${recalibration_report_filename} \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
      -L ${sep=" -L " sequence_group_interval}
  }
  output {
    File recalibration_report = "${recalibration_report_filename}"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
# Note that when run from GATK 3.x the tool is not a walker and is invoked differently.
task GatherBqsrReports {
  Array[File] input_bqsr_reports
  String output_report_filename
  String gatk4_path

  command {
    ${gatk4_path} --java-options "-Xms5g" \
      GatherBQSRReports \
      -I ${sep=' -I ' input_bqsr_reports} \
      -O ${output_report_filename}
  }
  output {
    File output_bqsr_report = "${output_report_filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  File input_bam
  File input_bam_index
  File recalibration_report
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  Array[String] sequence_group_interval
  String output_bam_basename
  String gatk4_path

  command {
    ${gatk4_path} --java-options "-Xms4g" \
      ApplyBQSR \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${output_bam_basename}.bam \
      -L ${sep=" -L " sequence_group_interval} \
      -bqsr ${recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities
  }
  output {
    File recalibrated_bam = "${output_bam_basename}.bam"
  }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  Array[File] input_bams
  String output_bam_basename
  String gatk4_path

  command {
    ${gatk4_path} --java-options "-Xms5g" \
      GatherBamFiles \
      --INPUT ${sep=' --INPUT ' input_bams} \
      --OUTPUT ${output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

# Collect base quality and insert size metrics
task CollectUnsortedReadgroupBamQualityMetrics {
  File in_bam
  String out_bam_prefix
  String picard_path

  command {
    ${picard_path} -Xmx16g \
      CollectMultipleMetrics \
      INPUT=${in_bam} \
      OUTPUT=${out_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM="null" \
      PROGRAM="CollectBaseDistributionByCycle" \
      PROGRAM="CollectInsertSizeMetrics" \
      PROGRAM="MeanQualityByCycle" \
      PROGRAM="QualityScoreDistribution" \
      METRIC_ACCUMULATION_LEVEL="null" \
      METRIC_ACCUMULATION_LEVEL="ALL_READS"
  }
  output {
    File base_distribution_by_cycle_pdf = "${out_bam_prefix}.base_distribution_by_cycle.pdf"
    File base_distribution_by_cycle_metrics = "${out_bam_prefix}.base_distribution_by_cycle_metrics"
    File insert_size_histogram_pdf = "${out_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "${out_bam_prefix}.insert_size_metrics"
    File quality_by_cycle_pdf = "${out_bam_prefix}.quality_by_cycle.pdf"
    File quality_by_cycle_metrics = "${out_bam_prefix}.quality_by_cycle_metrics"
    File quality_distribution_pdf = "${out_bam_prefix}.quality_distribution.pdf"
    File quality_distribution_metrics = "${out_bam_prefix}.quality_distribution_metrics"
  }
}

# Collect alignment summary and GC bias quality metrics
task CollectReadgroupBamQualityMetrics {
  File in_bam
  File in_bai
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  String picard_path
  String out_bam_prefix

  command {
    ${picard_path} -Xmx16g \
      CollectMultipleMetrics \
      INPUT=${in_bam} \
      REFERENCE_SEQUENCE=${ref_fasta} \
      OUTPUT=${out_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM="null" \
      PROGRAM="CollectAlignmentSummaryMetrics" \
      PROGRAM="CollectGcBiasMetrics" \
      METRIC_ACCUMULATION_LEVEL="null" \
      METRIC_ACCUMULATION_LEVEL="READ_GROUP"      
  }
  output {
    File alignment_summary_metrics = "${out_bam_prefix}.alignment_summary_metrics"
    File gc_bias_detail_metrics = "${out_bam_prefix}.gc_bias.detail_metrics"
    File gc_bias_pdf = "${out_bam_prefix}.gc_bias.pdf"
    File gc_bias_summary_metrics = "${out_bam_prefix}.gc_bias.summary_metrics"
  }
}

# Validate the output bam file
task ValidateSamFile {
  File in_bam
  File in_bai
  String report_filename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int? max_output
  Array[String]? ignore
  String picard_path

  command {
    ${picard_path} -Xmx16g \
      ValidateSamFile \
      INPUT=${in_bam} \
      OUTPUT=${report_filename} \
      REFERENCE_SEQUENCE=${ref_fasta} \
      ${"MAX_OUTPUT=" + max_output} \
      IGNORE=${default="null" sep=" IGNORE=" ignore} \
      MODE=VERBOSE \
      IS_BISULFITE_SEQUENCED=false
  }
  output {
    File report = "${report_filename}"
  }
}


##############################
## WORKFLOW DEFINITION PART ##
##############################

workflow PreProcessing4VariantDiscovery_GATK4 {
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
    
#  String CleanDataDir
#  String samplename
  String bwa_path
  String picard_path
  String samtools_path
  String gatk4_path


  call ArrayCleanData

#  Array[Array[File]] inputSamples = read_tsv(ArrayCleanData.outlist)
  scatter (fastq in ArrayCleanData.outlist) {
    call Fastq2ubam {
      input:
        cleanFq_R1 = fastq[0],
        cleanFq_R2 = fastq[1],
        sampleinfo = ArrayCleanData.sampleinfo,
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
        picard_path = picard_path,
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

    # QC the aligned but unsorted readgroup BAM
    # No reference needed as the input here is unsorted; providing a reference would cause an error
    call CollectUnsortedReadgroupBamQualityMetrics {
      input:
        picard_path = picard_path,
        in_bam = MergeBamAlignment.out_bam,
        out_bam_prefix = basename(Fastq2ubam.out_ubam, "_unmapped.bam") + ".readgroup"
    }

    # Sort and fix tags in the merged BAM
    call SortAndFixTags as SortAndFixReadGroupBam {
      input:
        gatk4_path = gatk4_path,
        input_bam = MergeBamAlignment.out_bam,
        output_bam_basename = basename(Fastq2ubam.out_ubam, "_unmapped.bam") + ".sorted",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }

    # Validate the aligned and sorted readgroup BAM
    # This is called to help in finding problems early.
    # If considered too time consuming and not helpful, can be removed.
    call ValidateSamFile as ValidateReadGroupSamFile {
      input:
        picard_path = picard_path,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        in_bam = SortAndFixReadGroupBam.output_bam,
        in_bai = SortAndFixReadGroupBam.output_bam_index,
        report_filename = basename(Fastq2ubam.out_ubam, "_unmapped.bam") + ".validation_report"
    }
  }

  call MarkDuplicates {
    input:
      input_bams = MergeBamAlignment.out_bam,
      sampleinfo = ArrayCleanData.sampleinfo,
      gatk4_path = gatk4_path
  }

  # Sort aggregated+deduped BAM file and fix tags
  call SortAndFixTags as SortAndFixSampleBam {
    input:
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = basename(MarkDuplicates.output_bam, ".aligned.unsorted.duplicates_marked.bam") + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      gatk4_path = gatk4_path
  }

  # Create list of sequences for scatter-gather parallelization 
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict
  }

  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        input_bam = SortAndFixSampleBam.output_bam,
        input_bam_index = SortAndFixSampleBam.output_bam_index,
        recalibration_report_filename = basename(MarkDuplicates.output_bam, ".aligned.unsorted.duplicates_marked.bam") + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        gatk4_path = gatk4_path
    }  
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports {
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = basename(MarkDuplicates.output_bam, ".aligned.unsorted.duplicates_marked.bam") + ".recal_data.csv",
      gatk4_path = gatk4_path
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {

    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        input_bam = SortAndFixSampleBam.output_bam,
        input_bam_index = SortAndFixSampleBam.output_bam_index,
        output_bam_basename = basename(GatherBqsrReports.output_bqsr_report, ".recal_data.csv") + ".aligned.duplicates_marked.recalibrated",
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        gatk4_path = gatk4_path
    }
  } 

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = basename(GatherBqsrReports.output_bqsr_report, ".recal_data.csv"),
      gatk4_path = gatk4_path
  }

 # Validate the final BAM
 call ValidateSamFile as ValidateAggregatedSamFile {
   input:
     picard_path = picard_path,
     in_bam = GatherBamFiles.output_bam,
     in_bai = GatherBamFiles.output_bam_index,
     report_filename = basename(GatherBqsrReports.output_bqsr_report, ".recal_data.csv") + ".validation_report",
     ref_dict = ref_dict,
     ref_fasta = ref_fasta,
     ref_fasta_index = ref_fasta_index
 }

 # QC the final BAM (consolidated after scattered BQSR)
 call CollectReadgroupBamQualityMetrics {
   input:
     picard_path = picard_path,
     in_bam = GatherBamFiles.output_bam,
     in_bai = GatherBamFiles.output_bam_index,
     out_bam_prefix = basename(GatherBqsrReports.output_bqsr_report, ".recal_data.csv") + ".readgroup",
     ref_dict = ref_dict,
     ref_fasta = ref_fasta,
     ref_fasta_index = ref_fasta_index
 }


  # Outputs that will be retained when execution is complete  
  output {
    File duplication_metrics = MarkDuplicates.duplicate_metrics
    File bqsr_report = GatherBqsrReports.output_bqsr_report
    File analysis_ready_bam = GatherBamFiles.output_bam
    File analysis_ready_bam_index = GatherBamFiles.output_bam_index
    File analysis_ready_bam_md5 = GatherBamFiles.output_bam_md5
  } 
}



