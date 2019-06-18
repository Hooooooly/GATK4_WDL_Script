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






