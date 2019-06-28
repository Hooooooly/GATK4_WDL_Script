## This WDL workflow runs GATK4 Mutect 2 on a single tumor-normal pair,
## and performs additional filtering and functional annotation tasks.
##
## Main requirements/expectations :
## - One analysis-ready BAM file (and its index) for each sample


##########################
## TASK DEFINITION PART ##
##########################

task Mutect2 {
  File? intervals
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File tumor_bam
  File tumor_bam_index
  File? normal_bam
  File? normal_bam_index
  File? pon
  File? pon_idx
  File? gnomad
  File? gnomad_idx
  String? m2_extra_args
  Boolean? make_bamout
  Boolean? run_ob_filter
  Boolean analyzeSoftClippedBase
  Boolean compress

  File? variants_for_contamination
  File? variants_for_contamination_idx

  String gatk4_path
  String output_vcf = "output" + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"
  String output_stats = output_vcf + ".stats"

  command <<<
    set -e

    # We need to create these files regardless, even if they stay empty
    touch bamout.bam
    touch f1r2.tar.gz
    echo "" > normal_name.txt

    ${gatk4_path} --java-options "-Xmx4g" GetSampleName \
    -R ${ref_fasta} -I ${tumor_bam} \
    -O tumor_name.txt -encode
    tumor_command_line="-I ${tumor_bam} -tumor `cat tumor_name.txt`"

    if [[ ! -z "${normal_bam}" ]]; then
        ${gatk4_path} --java-options "-Xmx4g" GetSampleName \
        -R ${ref_fasta} -I ${normal_bam} \
        -O normal_name.txt -encode
        normal_command_line="-I ${normal_bam} -normal `cat normal_name.txt`"
    fi

    ${gatk4_path} --java-options "-Xmx12g" Mutect2 \
        -R ${ref_fasta} \
        $tumor_command_line \
        $normal_command_line \
        ${"--germline-resource " + gnomad} \
        ${"-pon " + pon} \
        ${"-L " + intervals} \
        -O "${output_vcf}" \
        ${true='--bam-output bamout.bam' false='' make_bamout} \
        ${true='--f1r2-tar-gz f1r2.tar.gz' false='' run_ob_filter} \
        ${true='' false='--dont-use-soft-clipped-bases true' analyzeSoftClippedBase} \
        ${m2_extra_args}

    ### GetPileupSummaries
    # These must be created, even if they remain empty, as cromwell doesn't support optional output
    touch tumor-pileups.table
    touch normal-pileups.table

    if [[ ! -z "${variants_for_contamination}" ]]; then
        ${gatk4_path} --java-options "-Xmx8g" GetPileupSummaries \
        -R ${ref_fasta} \
        -I ${tumor_bam} \
        ${"--interval-set-rule INTERSECTION -L " + intervals} \
        -V ${variants_for_contamination} \
        -L ${variants_for_contamination} \
        -O tumor-pileups.table

        if [[ ! -z "${normal_bam}" ]]; then
            ${gatk4_path} --java-options "-Xmx8g" GetPileupSummaries \
            -R ${ref_fasta} \
            -I ${normal_bam} \
            ${"--interval-set-rule INTERSECTION -L " + intervals} \
            -V ${variants_for_contamination} \
            -L ${variants_for_contamination} \
            -O normal-pileups.table
        fi
    fi
  >>>
  output {
    File unfiltered_vcf = "${output_vcf}"
    File unfiltered_vcf_idx = "${output_vcf_idx}"
    File output_bamOut = "bamout.bam"
    String tumor_sample = read_string("tumor_name.txt")
    String normal_sample = read_string("normal_name.txt")
    File stats = "${output_stats}"
    File f1r2_counts = "f1r2.tar.gz"
    File tumor_pileups = "tumor-pileups.table"
    File normal_pileups = "normal-pileups.table"
  }
}

task MergeVCFs {

  Array[File] input_vcfs
  Array[File] input_vcf_indices
  Boolean compress
  String gatk4_path
  String output_name
  String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

  # using MergeVcfs instead of GatherVcfs so we can create indices
  # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
  command {
    ${gatk4_path} --java-options "-Xmx12g" MergeVcfs \
      -I ${sep=' -I ' input_vcfs} \
      -O ${output_vcf}
  }
  output {
    File merged_vcf = "${output_vcf}"
    File merged_vcf_idx = "${output_vcf_idx}"
  }
}

task MergeBamOuts {
  # inputs
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  Array[File]+ bam_outs
  String gatk4_path
  String output_vcf_name

  command <<<
    # This command block assumes that there is at least one file in bam_outs.
    #  Do not call this task if len(bam_outs) == 0
    ${gatk4_path} --java-options "-Xmx12g" GatherBamFiles \
      -I ${sep=" -I " bam_outs} \
      -O unsorted.out.bam \
      -R ${ref_fasta}

    # We must sort because adjacent scatters may have overlapping (padded) assembly regions, hence
    # overlapping bamouts

    ${gatk4_path} --java-options "-Xmx16g" SortSam \
      -I unsorted.out.bam \
      -O ${output_vcf_name}.out.bam \
      --SORT_ORDER coordinate \
      -VALIDATION_STRINGENCY LENIENT

    ${gatk4_path} --java-options "-Xmx12g" BuildBamIndex \
      -I ${output_vcf_name}.out.bam \
      -VALIDATION_STRINGENCY LENIENT
  >>>
  output {
    File merged_bam_out = "${output_vcf_name}.out.bam"
    File merged_bam_out_index = "${output_vcf_name}.out.bai"
  }
}


task MergeStats {
  # inputs
  Array[File]+ stats

  String gatk4_path

  command {
    ${gatk4_path} --java-options "-Xmx12g" MergeMutectStats \
      -stats ${sep=" -stats " stats} \
      -O merged.stats
  }
  output {
    File merged_stats = "merged.stats"
  }
}

task MergePileupSummaries {
  # input_tables needs to be optional because GetPileupSummaries is in an if-block
  Array[File?] input_tables
  File ref_dict
  String output_name
  String gatk4_path

  command {
    ${gatk4_path} --java-options "-Xmx12g" GatherPileupSummaries \
      --sequence-dictionary ${ref_dict} \
      -I ${sep=' -I ' input_tables} \
      -O ${output_name}.tsv
  }
  output {
    File merged_table = "${output_name}.tsv"
  }
}

# Calculates sum of a list of floats
task SumFloats {
  Array[Float] sizes

  command <<<
    python -c 'print(${sep="+" sizes})'
  >>>

  output {
    Float total_size = read_float(stdout())
  }
}

# Learning step of the orientation bias mixture model, which is the recommended orientation bias filter as of September 2018
task LearnReadOrientationModel {
  Array[File] f1r2_tar_gz
  String gatk4_path

  command {
    ${gatk4_path} --java-options "-Xmx16g" LearnReadOrientationModel \
      -I ${sep=" -I " f1r2_tar_gz} \
      -O "artifact-priors.tar.gz"
  }
  output {
    File artifact_prior_table = "artifact-priors.tar.gz"
  }
}

task CalculateContamination {
  # inputs
  String? intervals
  String gatk4_path
  File tumor_pileups
  File? normal_pileups

  command {
    ${gatk4_path} --java-options "-Xmx20g" CalculateContamination \
    -I ${tumor_pileups} \
    -O contamination.table \
    --tumor-segmentation segments.table \
    ${"-matched " + normal_pileups}
  }
  output {
    File contamination_table = "contamination.table"
    File maf_segments = "segments.table"
  }
}






##############################
## WORKFLOW DEFINITION PART ##
##############################

workflow Mutect2_GATK4 {

  File ref_fasta
  File ref_fasta_index
  File ref_dict
  
  File tumor_bam
  File tumor_bam_index
  File normal_bam
  File normal_bam_index
  
  File? pon
  File? pon_idx

  File variants_for_contamination
  File variants_for_contamination_idx
  File gnomad
  File gnomad_idx

  Array[File] scattered_calling_intervals

  Boolean? run_ob_filter
  Boolean? make_bamout
  Boolean is_ob_filter_ran = select_first([run_ob_filter, false])
  Boolean make_bamout_or_default = select_first([make_bamout, false])

  String gatk4_path
  String? m2_extra_args

  # logic about output file names -- these are the names *without* .vcf extensions
  String output_basename
  String unfiltered_name = output_basename + "-unfiltered"
  String filtered_name = output_basename + "-filtered"
  String funcotated_name = output_basename + "-funcotated"


  # Somatic variant calling with MuTect2 and summarization
  # of read support for a set number of known variants
  scatter (subintervals in scattered_calling_intervals) {
    call Mutect2 {
      input:
        intervals = subintervals,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        tumor_bam = tumor_bam,
        tumor_bam_index = tumor_bam_index,
        normal_bam = normal_bam,
        normal_bam_index = normal_bam_index,
        pon = pon,
        pon_idx = pon_idx,
        gnomad = gnomad,
        gnomad_idx = gnomad_idx,
        m2_extra_args = m2_extra_args,
        run_ob_filter = run_ob_filter,
        make_bamout = make_bamout,
        variants_for_contamination = variants_for_contamination,
        variants_for_contamination_idx = variants_for_contamination_idx,
        gatk4_path = gatk4_path
    }

    Float sub_vcf_size = size(Mutect2.unfiltered_vcf, "GB")
    Float sub_bamout_size = size(Mutect2.output_bamOut, "GB")
  }

  call SumFloats as SumSubVcfs {
    input:
      sizes = sub_vcf_size
  }

  if (is_ob_filter_ran) {
    call LearnReadOrientationModel {
      input:
        f1r2_tar_gz = Mutect2.f1r2_counts,
        gatk4_path = gatk4_path
    }
  }

  call MergeVCFs {
    input:
      input_vcfs = Mutect2.unfiltered_vcf,
      input_vcf_indices = Mutect2.unfiltered_vcf_idx,
      output_name = unfiltered_name,
      gatk4_path = gatk4_path
}


  if (make_bamout_or_default) {
    call SumFloats as SumSubBamouts {
      input:
        sizes = sub_bamout_size
    }

    call MergeBamOuts {
      input:
        ref_fasta = ref_fasta,
        ref_fasta_index =ref_fasta_index,
        ref_dict = ref_dict,
        bam_outs = Mutect2.output_bamOut,
        output_vcf_name = basename(MergeVCFs.merged_vcf, ".vcf.gz"),
        gatk4_path = gatk4_path
    }
}

  call MergeStats {
    input:
      stats = Mutect2.stats,
      gatk4_path = gatk4_path
}

  if (defined(variants_for_contamination)) {
    call MergePileupSummaries as MergeTumorPileups {
      input:
        input_tables = Mutect2.tumor_pileups,
        output_name = output_basename,
        ref_dict = ref_dict,
        gatk4_path = gatk4_path
    }

    if (defined(normal_bam)){
      call MergePileupSummaries as MergeNormalPileups {
        input:
          input_tables = Mutect2.normal_pileups,
          output_name = output_basename,
          ref_dict = ref_dict,
          gatk4_path = gatk4_path
      }
    }

    call CalculateContamination {
      input:
        tumor_pileups = MergeTumorPileups.merged_table,
        normal_pileups = MergeNormalPileups.merged_table,
        gatk4_path = gatk4_path
    }
}

  # Outputs that will be retained when execution is complete  
  output {

  }

}
