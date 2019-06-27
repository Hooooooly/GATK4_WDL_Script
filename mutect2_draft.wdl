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
  File ref_fai
  File ref_dict
  File tumor_bam
  File tumor_bai
  File? normal_bam
  File? normal_bai
  File? pon
  File? pon_idx
  File? gnomad
  File? gnomad_idx
  String? m2_extra_args
  Boolean? make_bamout
  Boolean? run_ob_filter
  Boolean compress
  File? gga_vcf
  File? gga_vcf_idx
  File? variants_for_contamination
  File? variants_for_contamination_idx

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
        ${"--alleles " + gga_vcf} \
        -O "${output_vcf}" \
        ${true='--bam-output bamout.bam' false='' make_bamout} \
        ${true='--f1r2-tar-gz f1r2.tar.gz' false='' run_ob_filter} \
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












##############################
## WORKFLOW DEFINITION PART ##
##############################

workflow Mutect2_GATK4 {

   Array[File] scattered_calling_intervals



  scatter (subintervals in scattered_calling_intervals) {
    call Mutect2 {
      input:
        intervals = subintervals,
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        tumor_bam = tumor_bam,
        tumor_bai = tumor_bai,
        normal_bam = normal_bam,
        normal_bai = normal_bai,
        pon = pon,
        pon_idx = pon_idx,
        gnomad = gnomad,
        gnomad_idx = gnomad_idx,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries,
        m2_extra_args = m2_extra_args,
        variants_for_contamination = variants_for_contamination,
        variants_for_contamination_idx = variants_for_contamination_idx,
        make_bamout = make_bamout_or_default,
        run_ob_filter = run_ob_filter,
        compress = compress,
        gga_vcf = gga_vcf,
        gga_vcf_idx = gga_vcf_idx
    }

    Float sub_vcf_size = size(M2.unfiltered_vcf, "GB")
    Float sub_bamout_size = size(M2.output_bamOut, "GB")
  }

  # Somatic variant calling with MuTect2
  call MuTect2 {
    input:
      GATK4_LAUNCH=gatk4_launch,
      contamination = CheckContamination_tumor.contamination,
      in_bam_tumor = GatherBamFiles_tumor.out_bam,
      in_bai_tumor = GatherBamFiles_tumor.out_bai,
      in_bam_normal = GatherBamFiles_normal.out_bam,
      in_bai_normal = GatherBamFiles_normal.out_bai,
      sample_name = sample_name,
      sample_name_tumor = sample_name+'_tumor',
      sample_name_normal = sample_name+'_normal',
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx,
      transcript_intervals = transcript_intervals,
      gnomad_exome_vcf = gnomad_exome_vcf,
      gnomad_exome_vcf_idx = gnomad_exome_vcf_idx,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_idx = dbSNP_vcf_idx
  }

}
