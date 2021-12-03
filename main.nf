#!/usr/bin/env nextflow

// INPUTS
// ref_hg19  = file("/home/phelelani/projects/laura/data/indexes/hg19/human_g1k_v37_decoy.fa", type:'file')
ref_hg38  = file("/home/phelelani/projects/laura/data/indexes/hg38/Homo_sapiens_assembly38.fa", type:'file')
bams = Channel.fromFilePairs("/external/diskC/build38/datasets/bam/1000G/NA1848*{.bam,.bam.bai}", size:-2)
contigs = file("/home/phelelani/nf-workflows/nf-sve/data/contigs.txt", type:'file')

// DELLY PARAMS
chromosomes = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
chroms = ( 1..22 ).collect { "chr$it" }

// OUTPUT DIR
out_dir = file("${PWD}/test_out", type: 'dir')
out_dir.mkdir()

// =============================================== SPLIT ============================================== //
process splitSampleBams {
    tag { sample }
    
    input:
    each chr from chromosomes
    set sample, file(bam) from bams
    
    output:
    set sample, file("${sample}_chr${chr}.bam"), file("${sample}_chr${chr}.bam.bai") into sampleBAMsChr

    """
    samtools view -bh ${sample}.bam chr${chr} > ${sample}_chr${chr}.bam
    samtools index ${sample}_chr${chr}.bam
    """
}

// CREATE MUTLIPLE CHANNELS FOR EACH PROCESS
sampleBAMsChr.into { sampleBAMsChr_delly; sampleBAMsChr_gatk; sampleBAMsChr_prelumpy; sampleBAMsChr_manta }
// =============================================== SPLIT ============================================== //
//
//
// =============================================== DELLY ============================================== //
process delly {
    tag { sample }
    memory '5 GB'

    input:
    set sample, file(bam), file(bai) from sampleBAMsChr_delly

    output:
    set sample, file("${bam.baseName}.bcf"), file("${bam.baseName}.bcf.csi") into sampleBCFsChr_delly
    
    """
    /opt/exp_soft/bioinf/bin/delly call \
        --outfile ${bam.baseName}.bcf \
        --genome ${ref_hg38} \
        ${bam}
    """
}

sampleBCFsChr_delly
    .groupTuple(by: 0)
    .set { sampleBCFsChr_delly_grp }

process delly_mergeSampleChr {
    tag { sample }
    publishDir "${out_dir}/delly", mode: 'copy'

    input:
    set sample, file(bcf), file(csi) from sampleBCFsChr_delly_grp

    output:
    set sample, file("${sample}_delly.vcf.gz"), file("${sample}_delly.vcf.gz.tbi") into sampleVCFs_delly

    """
    bcftools concat \
        ${bcf.find { it =~ '_chr1.bcf$' } } \
        ${bcf.find { it =~ '_chr2.bcf$' } } \
        ${bcf.find { it =~ '_chr3.bcf$' } } \
        ${bcf.find { it =~ '_chr4.bcf$' } } \
        ${bcf.find { it =~ '_chr5.bcf$' } } \
        ${bcf.find { it =~ '_chr6.bcf$' } } \
        ${bcf.find { it =~ '_chr7.bcf$' } } \
        ${bcf.find { it =~ '_chr8.bcf$' } } \
        ${bcf.find { it =~ '_chr9.bcf$' } } \
        ${bcf.find { it =~ '_chr10.bcf$' } } \
        ${bcf.find { it =~ '_chr11.bcf$' } } \
        ${bcf.find { it =~ '_chr12.bcf$' } } \
        ${bcf.find { it =~ '_chr13.bcf$' } } \
        ${bcf.find { it =~ '_chr14.bcf$' } } \
        ${bcf.find { it =~ '_chr15.bcf$' } } \
        ${bcf.find { it =~ '_chr16.bcf$' } } \
        ${bcf.find { it =~ '_chr17.bcf$' } } \
        ${bcf.find { it =~ '_chr18.bcf$' } } \
        ${bcf.find { it =~ '_chr19.bcf$' } } \
        ${bcf.find { it =~ '_chr20.bcf$' } } \
        ${bcf.find { it =~ '_chr21.bcf$' } } \
        ${bcf.find { it =~ '_chr22.bcf$' } } \
        -Oz -o ${sample}_delly_unsorted.vcf.gz
    bcftools sort ${sample}_delly_unsorted.vcf.gz -Oz -o ${sample}_delly.vcf.gz
    tabix ${sample}_delly.vcf.gz
    """
}
// =============================================== DELLY ============================================== //
//
//
// =============================================== GATK =============================================== //
process gatk {
    tag { sample }
    cpus 5
    memory '3 GB'
    clusterOptions = '--exclude=n12,n01,n02'
    
    input:
    set sample, file(bam), file(bai) from sampleBAMsChr_gatk

    output:
    set sample, file("${bam.baseName}.vcf") into sampleVCFsChr_gatk
    
    """
    /bin/hostname
    gatk --java-options "-Xmx3g -XX:+UseParallelGC -XX:ParallelGCThreads=1" HaplotypeCaller \
        --pairHMM AVX_LOGLESS_CACHING_OMP --native-pair-hmm-threads 4 \
        --reference ${ref_hg38} \
        --input ${bam} \
        --output ${bam.baseName}.vcf
    """
}

sampleVCFsChr_gatk
    .groupTuple(by: 0)
    .set { sampleVCFsChr_gatk_grp }

process gatk_mergeSampleChr {
    tag { sample }
    publishDir "${out_dir}/gatk_haplo", mode: 'copy'
    
    input:
    set sample, file(vcf) from sampleVCFsChr_gatk_grp

    output:
    set sample, file("${sample}_gatk_haplo.vcf.gz"), file("${sample}_gatk_haplo.vcf.gz.tbi") into sampleVCFs_gatk
    
    """
    gatk --java-options "-Xmx2g" MergeVcfs \
        --INPUT ${vcf.find { it =~ '_chr1.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr2.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr3.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr4.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr5.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr6.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr7.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr8.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr9.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr10.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr11.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr12.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr13.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr14.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr15.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr16.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr17.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr18.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr19.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr20.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr21.vcf$' } } \
        --INPUT ${vcf.find { it =~ '_chr22.vcf$' } } \
        -OUTPUT ${sample}_gatk_haplo.vcf.gz \
        --COMPRESSION_LEVEL 1
    """
}
// =============================================== GATK =============================================== //
//
//
// =============================================== LUMPY ============================================== //
process lumpy_preProcess {
    tag { sample }
    cpus 4

    input:
    set sample, file(bam), file(bai) from sampleBAMsChr_prelumpy

    output:
    set sample, file(bam), file(bai), 
        file("${bam.baseName}.discordants.bam"), file("${bam.baseName}.splitters.bam") into sampleBAMsChrSet_prelumpy
    
    """
    samtools view -b -F 1294 ${bam} > ${bam.baseName}.discordants.unsorted.bam
    samtools view -h ${bam} | \
        /opt/exp_soft/bioinf/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | \
         samtools view -Sb - > ${bam.baseName}.splitters.unsorted.bam
    samtools sort ${bam.baseName}.discordants.unsorted.bam --threads 4 -o ${bam.baseName}.discordants.bam
    samtools sort ${bam.baseName}.splitters.unsorted.bam --threads 4 -o ${bam.baseName}.splitters.bam
    """
}

process lumpy {
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
// =============================================== LUMPY ============================================== //
//
//
// =============================================== MANTA ============================================== //
process manta {
	tag { sample }
    cpus 5
    memory '5 GB'

    input:
    set sample, file(bam), file(bai) from sampleBAMsChr_manta
    
    output:
    set val("${bam.baseName.replaceAll("${sample}_","")}"),
        file("${bam.baseName}.vcf.gz"), file("${bam.baseName}.vcf.gz.tbi") into chrVCFs_manta
	set val("${bam.baseName.replaceAll("${sample}_","")}"), file(bam), file(bai) into chrBAMs_graphtyper

	"""
	/opt/exp_soft/bioinf/manta/bin/configManta.py \
	    --bam ${bam} \
	        --referenceFasta ${ref_hg38} \
		    --runDir .

		    ./runWorkflow.py -j 5 -g 5
    cp results/variants/diploidSV.vcf.gz ${bam.baseName}.vcf.gz
    cp results/variants/diploidSV.vcf.gz.tbi ${bam.baseName}.vcf.gz.tbi
    """
}

chrVCFs_manta
    .groupTuple(by: 0)
    .set { chrVCFs_manta_grp }

process svimmer {
    tag { chr }
    cpus 4

    input:
    set chr, file(vcf), file(tbi) from chrVCFs_manta_grp

    output:
    set chr, file("svimmer_merged_${chr}.vcf.gz"), file("svimmer_merged_${chr}.vcf.gz.tbi") into chrVCFs_svimmer
    
    """
    ls *_${chr}.vcf.gz > input_vcfs
    /opt/exp_soft/bioinf/svimmer/svimmer --threads 4 input_vcfs ${chr} --output svimmer_merged_${chr}.vcf
    sed -i "/chr6.*BND/d" svimmer_merged_${chr}.vcf
    bgzip svimmer_merged_${chr}.vcf
    tabix svimmer_merged_${chr}.vcf.gz
    """
}

chrBAMs_graphtyper
    .groupTuple(by: 0)
    .join(chrVCFs_svimmer,by: 0)
    .set { chrBAMsVCFs_graphtyper }

process graphTyper {
   tag { chr }
    errorStrategy 'finish'
    memory '50 GB'
    cpus 4

    input:
    set chr, file(bam), file(bai), file(vcf), file(tbi) from chrBAMsVCFs_graphtyper

    output:
    set chr, file("sv_results/${chr}") into chrVCFs_graphtyper

    """
    #!/bin/bash
    ls *_${chr}.bam > bamlist
    graphtyper genotype_sv ${ref_hg38} ${vcf} --sams=bamlist --region=${chr} --threads=4
    """
}

process graphTyper_mergeByChr {
    tag { chr }

    input:
    set chr, file(chrom) from chrVCFs_graphtyper

    output:
    set file("${chr}.vcf.gz"), file("${chr}.vcf.gz.tbi") into chrVCFsMerged_graphtyper
    
    """
    bcftools concat $chrom/*vcf.gz -Oz -o ${chr}_unsorted.vcf.gz
    bcftools sort ${chr}_unsorted.vcf.gz -Oz -o ${chr}.vcf.gz
    tabix ${chr}.vcf.gz
    """
}

process graphTyper_mergeAllChroms {
    publishDir "${out_dir}/manta", mode: 'copy'

    input:
    file(chroms) from chrVCFsMerged_graphtyper.collect()

    output:
    set file("manta_svimmer_graphTyper_chromosomes.vcf.gz"), file("manta_svimmer_graphTyper_chromosomes.vcf.gz") into allResult
    
    """
    bcftools concat *vcf.gz --output-type z --output manta_svimmer_graphTyper_unsorted.vcf.gz
    bcftools sort manta_svimmer_graphTyper_unsorted.vcf.gz -Oz -o manta_svimmer_graphTyper_chromosomes.vcf.gz   
    tabix manta_svimmer_graphTyper_chromosomes.vcf.gz
    """
}
// =============================================== MANTA ============================================== //