#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process lumpy_preProcess {
    tag { sample }
    cpus 4

    input:
    tuple val(sample), file(bam), file(bai)

    output:
    tuple val(sample), file(bam), file(bai), 
        file("${bam.baseName}.discordants.bam"), file("${bam.baseName}.splitters.bam")
    
    """
    samtools view -b -F 1294 ${bam} > ${bam.baseName}.discordants.unsorted.bam
    samtools view -h ${bam} | \
        /opt/exp_soft/bioinf/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | \
         samtools view -Sb - > ${bam.baseName}.splitters.unsorted.bam
    samtools sort ${bam.baseName}.discordants.unsorted.bam --threads 4 -o ${bam.baseName}.discordants.bam
    samtools sort ${bam.baseName}.splitters.unsorted.bam --threads 4 -o ${bam.baseName}.splitters.bam
    """
}

process run_lumpy {
    tag { sample }
    memory '5 GB'

    input:
    set sample, file(bam), file(bai), file(disc), file(spl) from sampleBAMsChrSet_prelumpy

    output:
    set sample, file("${bam.baseName}_tmp.vcf") into sampleVCFsChr_lumpy
    
    """
    /opt/exp_soft/bioinf/lumpy-sv/bin/lumpyexpress \
        -B ${bam} \
        -S ${spl} \
        -D ${disc} \
        -K /home/phelelani/nf-workflows/nf-sve/bin/lumpyexpress.config \
        -o ${bam.baseName}_tmp.vcf
    """
}

process lumpy_postProcess {
    tag { sample }

    input:
    set sample, file(vcf) from sampleVCFsChr_lumpy

    output:
    set sample, file("${vcf.baseName.replaceAll("_tmp","")}.vcf") into sampleVCFsChr_postlumpy
    
    """
    cp ${vcf} tmp.vcf

    while read line
    do
        sed -i '/^#CHROM/i '"\$line"'' tmp.vcf
    done < ${contigs}

    bcftools sort tmp.vcf -Ov -o ${vcf.baseName.replaceAll("_tmp","")}.vcf
    """
}

sampleVCFsChr_postlumpy
    .groupTuple(by: 0)
    .set { sampleVCFsChr_lumpy_grp }

process lumpy_mergeSampleChr {
    tag { sample }
    publishDir "${out_dir}/lumpy", overwrite: true, mode: 'copy'

    input:
    set sample, file(vcf) from sampleVCFsChr_lumpy_grp

    output:
    set sample, file("${sample}_lumpy.vcf.gz") into sampleVCFs_lumpy

    """
    bcftools concat \
        ${vcf.find { it =~ '_chr1.vcf$' } } \
        ${vcf.find { it =~ '_chr2.vcf$' } } \
        ${vcf.find { it =~ '_chr3.vcf$' } } \
        ${vcf.find { it =~ '_chr4.vcf$' } } \
        ${vcf.find { it =~ '_chr5.vcf$' } } \
        ${vcf.find { it =~ '_chr6.vcf$' } } \
        ${vcf.find { it =~ '_chr7.vcf$' } } \
        ${vcf.find { it =~ '_chr8.vcf$' } } \
        ${vcf.find { it =~ '_chr9.vcf$' } } \
        ${vcf.find { it =~ '_chr10.vcf$' } } \
        ${vcf.find { it =~ '_chr11.vcf$' } } \
        ${vcf.find { it =~ '_chr12.vcf$' } } \
        ${vcf.find { it =~ '_chr13.vcf$' } } \
        ${vcf.find { it =~ '_chr14.vcf$' } } \
        ${vcf.find { it =~ '_chr15.vcf$' } } \
        ${vcf.find { it =~ '_chr16.vcf$' } } \
        ${vcf.find { it =~ '_chr17.vcf$' } } \
        ${vcf.find { it =~ '_chr18.vcf$' } } \
        ${vcf.find { it =~ '_chr19.vcf$' } } \
        ${vcf.find { it =~ '_chr20.vcf$' } } \
        ${vcf.find { it =~ '_chr21.vcf$' } } \
        ${vcf.find { it =~ '_chr22.vcf$' } } \
        -Oz -o ${sample}_lumpy_unsorted.vcf.gz
    bcftools sort ${sample}_lumpy_unsorted.vcf.gz -Oz -o ${sample}_lumpy.vcf.gz
    tabix -p vcf ${sample}_lumpy.vcf.gz
    """
}