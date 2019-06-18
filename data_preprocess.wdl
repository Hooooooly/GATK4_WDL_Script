## This is for clean fastq reads processing
## in oder to get analysis-ready bam files





## Generate a file array for trim
task ArrayRawData {
  String rawDatadir
  
  command <<<
    #cd ${rawDatadir}
    ls ${rawDatadir}/*.fastq.gz | \
    xargs -n2 | tr ' ' '\t' > data4trimlist
  >>>
  output {
    File outlist = "data4trimlist"
  }
}

## Convert clean fastq reads to ubam
task Fastq2ubam {
  File cleanFq_R1
  File cleanFq_R2
  File PICARD
  String ?Platform
  String basenameR1=basename(cleanFq_R1, ".fastq.gz")
  String basenameR1=basename(cleanFq_R2, ".fastq.gz")
  
  command <<<
    ${PICARD} --javaOptions "-Xms5g -Xmx8g" FastqToSam \
    F1=${cleanFq_R1} \
    F2=${cleanFq_R2} \
    PL=COMPLETE \
    SM=${samplename} \
    LB=${library} \
    RG=${readgroup} \
    O=${basenameR1}.ubam

  >>>
  output {
    File out_ubam = "${basenameR1}.ubam"
  }
}





workflow Raw_Reads_Trimming_Using_Trimmomatic_wf {
    File genome_fa = '/home/projects/pr_99006/data/rna/genome/hg38/hg38.fa'
    File bedtools = '/services/tools/bedtools/2.26.0/bin/bedtools'
    File samtools = '/services/tools/samtools/1.4.1/bin/samtools'
    File stringtie = '/services/tools/stringtie/1.3.3b/stringtie'
    File annotation = '/home/projects/pr_99006/data/rna/genome/hg38_tran/hg38_ucsc.annotated.filtered.gtf'
    String rawDatadir = "/home/projects/pr_99006/data/rna/hiPSC-derived_NPC_Neuron/data4Trim"

    call ArrayRawData4Trim {
      input:
        rawDatadir = rawDatadir
    }

    Array[Array[File]] inputSamples = read_tsv(ArrayRawData4Trim.outlist)
    scatter (fastq in inputSamples) {
      call TrimReads {
        input:
          in_fastqR1 = fastq[0],
          in_fastqR2 = fastq[1]
      }
    
      call FastQC {
        input:
          in_R1 = TrimReads.out_R1,
          in_R2 = TrimReads.out_R2
     }

      call Mapping as Mapping4stringTie {
        input:
          in_R1 = TrimReads.out_R1,
          in_R2 = TrimReads.out_R2,
          tag = "--dta",
          SAMTOOLS = samtools
      }

      call Mapping as Mapping4cufflinks {
        input:
          in_R1 = TrimReads.out_R1,
          in_R2 = TrimReads.out_R2,
          tag = "--dta-cufflinks",
          SAMTOOLS = samtools
      }

      call extractUniqreadsAndSort as extractUniqreadsAndSort4stringTie {
        input:
          in_bam = Mapping4stringTie.out_bam,
          SAMTOOLS = samtools
      }

      call extractUniqreadsAndSort as extractUniqreadsAndSort4cufflinks {
        input:
          in_bam = Mapping4cufflinks.out_bam,
          SAMTOOLS = samtools
      }

      call Assembly {
        input:
          STRINGTIE = stringtie,
          in_uniq_bam = extractUniqreadsAndSort4stringTie.out_uniq_bam,
          annotation = annotation
      }

      call DepthCalculation {
        input:
          in_uniq_bam = extractUniqreadsAndSort4stringTie.out_uniq_bam,
          in_bam_idx = extractUniqreadsAndSort4stringTie.out_bam_idx,
          FwdExon1Range = ExtractFirst2exononPTransAndLast2exononNTrans.out_FwdExon1Range,
          RevExon1Range = ExtractFirst2exononPTransAndLast2exononNTrans.out_RevExon1Range,
          SAMTOOLS = samtools
      }

      call ConvertDepthToFPKM as ConvertDepthToFPKM_FwdExon {
        input:
          TotalDep = DepthCalculation.out_totalDep,
          Exon1Dep = DepthCalculation.out_FwdExon1Depth
      }

      call ConvertDepthToFPKM as ConvertDepthToFPKM_RevExon {
        input:
          TotalDep = DepthCalculation.out_totalDep,
          Exon1Dep = DepthCalculation.out_RevExon1Depth
      }

      call GetSplicedReads {
        input:
          in_uniq_bam = extractUniqreadsAndSort4stringTie.out_uniq_bam
      }

      call QuantificationPerSample {
        input:
          genome_fa = genome_fa,
          in_bam = extractUniqreadsAndSort4cufflinks.out_uniq_bam,
          annotation = MergeGTF.out_merged_gtf
      }
    }

      call MergeGTF {
        input:
          STRINGTIE = stringtie,
          annotation = annotation,
          in_gtf = Assembly.out_gtf
      }

      call ExpressionLevelCal {
        input:
          in_cxb = QuantificationPerSample.out_abundance_cxb,
          mergedGTF = MergeGTF.out_merged_gtf
      }

      call TranscriptMT200bp {
        input:
          mergedGTF = MergeGTF.out_merged_gtf
      }

      call GenesMT1FPKMand1Transcript {
        input:
          geneFPKM = ExpressionLevelCal.out_genefpkm,
          geneMT1Tran = TranscriptMT200bp.out_geneID_genes_MT1transcript,
          tranMT200bpGTF = TranscriptMT200bp.out_gtf_tranMT200bp_furtherfiltered_mulT
      }

      call SingleExonTranscriptsRemoval {
        input:
          genesMT1TranAnd1FPKM = GenesMT1FPKMand1Transcript.out_gtf_Genes_moreThan1transANDmoreThan1FPKM
      }

      call ExtractFirst2exononPTransAndLast2exononNTrans {
        input:
          GenesMT1transandMT1FPKMwithSingleExonRemoved = SingleExonTranscriptsRemoval.out_gtf_GenesMT1transandMT1FPKMwithSingleExonRemoved
      }

      call ScreenBasedonExon1Cov as ScreenBasedonExon1Cov_FwdExon {
        input:
          in_fpkm = ConvertDepthToFPKM_FwdExon.out_fpkm
      }

      call ScreenBasedonExon1Cov as ScreenBasedonExon1Cov_RevExon {
        input:
          in_fpkm = ConvertDepthToFPKM_RevExon.out_fpkm
      }

      call First1ExonandLast1ExonGTF {
        input:
          in_fwd_gtf = ExtractFirst2exononPTransAndLast2exononNTrans.out_gtf_First2ExonFwd,
          in_rev_gtf = ExtractFirst2exononPTransAndLast2exononNTrans.out_gtf_Last2ExonRev
      }

      call CorrespondingTrans as CorrespondingTrans_Fwd {
        input:
          in_list = ScreenBasedonExon1Cov_FwdExon.out_Exon1RangeFPKM_screened,
          in_gtf = First1ExonandLast1ExonGTF.out_gtf_First1ExonFwd, 
          in_GTF = ExtractFirst2exononPTransAndLast2exononNTrans.out_gtf_First2ExonFwd
      }

      call CorrespondingTrans as CorrespondingTrans_Rev {
        input:
          in_list = ScreenBasedonExon1Cov_RevExon.out_Exon1RangeFPKM_screened,
          in_gtf = First1ExonandLast1ExonGTF.out_gtf_Last1ExonRev, 
          in_GTF = ExtractFirst2exononPTransAndLast2exononNTrans.out_gtf_Last2ExonRev
      }

      call SATSTrans {
      input:
        in_Fwd_ID = CorrespondingTrans_Fwd.out_ID_correspondingTrans,
        in_Rev_ID = CorrespondingTrans_Rev.out_ID_correspondingTrans,
        GenesMT1transandMT1FPKMwithSingleExonRemoved = SingleExonTranscriptsRemoval.out_gtf_GenesMT1transandMT1FPKMwithSingleExonRemoved
      }

      call bam2tdf {
        input:
          in_bam = Mapping4stringTie.out_bam,
          BEDTOOLS = bedtools
      }
}


