
## Whole Pipeline GATK4 Somatic Variants Discovery


import "data_preprocess.wdl" as PreProcess
import "mutect2.wdl" as M2


##############################
## WORKFLOW DEFINITION PART ##
##############################


workflow FullSomaticPipeline {
  
  ## PreProcess parameters
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
    
  String bwa_path
  String picard_path
  String samtools_path
  String gatk4_path


  call PreProcess.PreProcessing4VariantDiscovery_GATK4 as PreProcess_Tumor {
    input:
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,

      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,

      bwa_path = bwa_path,
      picard_path = picard_path,
      samtools_path = samtools_path,
      gatk4_path = gatk4_path
  }

  call PreProcess.PreProcessing4VariantDiscovery_GATK4 as PreProcess_Normal {
    input:
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,

      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,

      bwa_path = bwa_path,
      picard_path = picard_path,
      samtools_path = samtools_path,
      gatk4_path = gatk4_path
  }

  call M2.Mutect2_GATK4 {
    input:
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      tumor_bam = PreProcess_Tumor.analysis_ready_bam,
      tumor_bam_index = PreProcess_Tumor.analysis_ready_bam_index,
      normal_bam = PreProcess_Normal.analysis_ready_bam,
      normal_bam_index = PreProcess_Normal.analysis_ready_bam_index,
      gatk4_path = gatk4_path
  }

}


