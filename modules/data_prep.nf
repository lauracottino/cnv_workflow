#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def check_files(file_list) {
    file_list.each { myfile ->
        if (!file(myfile).exists() && !file(myfile).isFile()) exit 1, "|-- ERROR: File ${myfile} not found. Please check your config file."
    }
}

process splitSampleBams {
    tag { sample }
    
    input:
    tuple val(sample), file(bam), val(chr)
    
    output:
    tuple val(sample), file("${sample}_chr${chr}.bam"), file("${sample}_chr${chr}.bam.bai")

    """
    samtools view -bh ${sample}.bam chr${chr} > ${sample}_chr${chr}.bam
    samtools index ${sample}_chr${chr}.bam
    """
}