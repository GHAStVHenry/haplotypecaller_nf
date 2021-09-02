nextflow.enable.dsl = 2

params.sampleID = ""
params.fastqR1 = ""
params.fastqR2 = ""
params.ref = ""
params.fastaIn = ""
params.knownVariants = ""
params.dbsnp = ""
params.outDir = ""

sampleID = params.sampleID
fastqR1 = Channel.fromPath(params.fastqR1)
fastqR2 = Channel.fromPath(params.fastqR2)
ref = Channel.fromPath(params.ref)
fastaIn = Channel.fromPath(params.fastaIn)
knownVariants = Channel.fromPath(params.knownVariants)
dbsnp = Channel.fromPath(params.dbsnp)

workflow {
    main:
        mapping( sampleID, fastqR1, fastqR2, ref )
        markdup( sampleID, mapping.out.bam )
        fastaIndex( fastaIn )
        recal( sampleID, markdup.out.bam_md, markdup.out.bai_md, fastaIndex.out.fasta, fastaIndex.out.fai, knownVariants )
        bamIndex( recal.out.bam_recal )
        haplotypecaller( sampleID, recal.out.bam_recal, bamIndex.out.bai_recal, fastaIndex.out.fasta, fastaIndex.out.fai, recal.out.dict, dbsnp )
        genotype( sampleID, fastaIndex.out.fasta, fastaIndex.out.fai, recal.out.dict, dbsnp, haplotypecaller.out.gvcf )
    emit:
        genotype.out.vcf
        genotype.out.tbi
}

process mapping {
    machineType "mem2_ssd1_v2_x32"
    container "quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1"
    publishDir "${params.outDir}"
    input:
        val sampleID
        path fastqR1
        path fastqR2
        path ref
    output:
        path "${sampleID}.bam", emit: bam
        path "${sampleID}.bam.bai", emit: bai
    script:
        """
        mkdir -p genome
        zcat "${ref}" | tar xvf - -C genome
        genome_file=`ls ./genome/*.bwt`
        genome_file="\${genome_file%.bwt}"
        bwa mem -K 100000000 -t 16 -M \
            \${genome_file} \
            ${fastqR1} ${fastqR2} \
            -R "@RG\\tID:None\\tPL:None\\tPU:None\\tLB:None\\tSM:${sampleID}" | \
            samtools sort --threads 16 -m 2G - > ${sampleID}.bam
        samtools index ${sampleID}.bam
        """
}

process markdup {
    machineType "mem3_ssd1_v2_x8"
    container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    publishDir "${params.outDir}"
    input:
        val sampleID
        path bam
    output:
        path "./${sampleID}.md.bam", emit: bam_md
        path "./${sampleID}.md.bam.bai", emit: bai_md
    script:
        """
        gatk --java-options "-Xmx40g -Xms6000m" \
            MarkDuplicates \
                --INPUT ${bam} \
                --METRICS_FILE ${sampleID}.md.bam.metrics \
                --TMP_DIR . \
                --ASSUME_SORT_ORDER coordinate \
                --CREATE_INDEX true \
                --OUTPUT ${sampleID}.md.bam
        mv ${sampleID}.md.bai ${sampleID}.md.bam.bai
        """
}

process fastaIndex {
    machineType "mem1_ssd1_v2_x2"
    container "quay.io/biocontainers/mulled-v2-0560a8046fc82aa4338588eca29ff18edab2c5aa:c17ce694dd57ab0ac1a2b86bb214e65fedef760e-0"
    publishDir "${params.outDir}"
    input:
        path fastaIn
    output:
        path "fasta/*.fa", emit: fasta
        path "fasta/*.fai", emit: fai
    script:
        """
        mkdir -p fasta
        cp ${fastaIn} ./fasta/
        fasta=`ls ./fasta/*.fa.gz`
        zcat \${fasta} > \${fasta%.gz}
        fasta="\${fasta%.gz}"
        samtools faidx \${fasta}
        """
}

process recal {
    machineType "mem3_ssd1_v2_x8"
    container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    publishDir "${params.outDir}"
    input:
        val sampleID
        path bam
        path bai
        path fasta
        path fai
        path knownVariants
    output:
        path "*.recal.bam", emit: bam_recal
        path "fasta/*.dict", emit: dict
    script:
        """
        mkdir -p fasta
        cp ${fasta} ./fasta/
        fasta=`ls ./fasta/*.fa`
        cp ${fai} ./fasta/
        fai=`ls ./fasta/*.fai`
        gatk --java-options -Xmx40g \
            CreateSequenceDictionary \
                -R "\${fasta}" \
                -O "\${fasta%.fa}.dict"
        known=""
        for knownVariant in `echo ${knownVariants} | tr "," "\n"`; do
            knownVariantBase=`basename "\${knownVariant}"`
            zcat "\${knownVariant}" > ./\${knownVariantBase%.gz}
            gatk --java-options -Xmx40g \
                IndexFeatureFile \
                -I \${knownVariantBase%.gz}
            known="\${known} --known-sites ./\${knownVariantBase%.gz}"
        done
        gatk --java-options -Xmx40g \
            BaseRecalibrator \
                \${known} \
                -I ${bam} \
                -O ${sampleID}.recal.table \
                --tmp-dir . \
                -R \${fasta} \
                --use-original-qualities \
                --verbosity INFO
        gatk --java-options -Xmx40g \
            ApplyBQSR \
                -R \${fasta} \
                --input ${bam} \
                --output ${sampleID}.recal.bam \
                --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --add-output-sam-program-record --use-original-qualities \
                --bqsr-recal-file ${sampleID}.recal.table
        """
}

process bamIndex {
    machineType "mem1_ssd1_v2_x2"
    container "quay.io/biocontainers/mulled-v2-0560a8046fc82aa4338588eca29ff18edab2c5aa:c17ce694dd57ab0ac1a2b86bb214e65fedef760e-0"
    publishDir "${params.outDir}"
    input:
        path bam
    output:
        path "bam/*.bai", emit: bai_recal
    script:
        """
        mkdir -p bam
        cp ${bam} ./bam/
        bam=`ls ./bam/*.bam`
        samtools index \${bam}
        """
}

process haplotypecaller {
    machineType "mem3_ssd1_v2_x8"
    container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    publishDir "${params.outDir}"
    input:
        val sampleID
        path bam
        path bai
        path fasta
        path fai
        path dict
        path dbsnp
    output:
        path "${sampleID}.g.vcf", emit: gvcf
    script:
        """
        mkdir -p fasta
        cp ${fasta} ./fasta/
        fasta=`ls ./fasta/*.fa`
        cp ${fai} ./fasta/
        fai=`ls ./fasta/*.fai`
        cp ${dict} ./fasta/
        dict=`ls ./fasta/*.dict`
        mkdir -p dbsnp
        cp ${dbsnp} ./dbsnp/
        dbsnp=`ls ./dbsnp/*.vcf.gz`
        zcat "\${dbsnp}" > ./\${dbsnp%.gz}
        gatk --java-options -Xmx40g \
            IndexFeatureFile \
            -I \${dbsnp%.gz}
        gatk --java-options "-Xmx40g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
                -R \${fasta} \
                -I ${bam} \
                -O ${sampleID}.g.vcf \
                -ERC GVCF
        """
}

process genotype {
    machineType "mem3_ssd1_v2_x8"
    container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    publishDir "${params.outDir}"
    input:
        val sampleID
        path fasta
        path fai
        path dict
        path dbsnp
        path gvcf
    output:
        path "${sampleID}.vcf.gz", emit: vcf
        path "${sampleID}.vcf.gz.tbi", emit: tbi
    script:
        """
        mkdir -p fasta
        cp ${fasta} ./fasta/
        fasta=`ls ./fasta/*.fa`
        cp ${fai} ./fasta/
        fai=`ls ./fasta/*.fai`
        cp ${dict} ./fasta/
        dict=`ls ./fasta/*.dict`
        mkdir -p dbsnp
        cp ${dbsnp} ./dbsnp/
        dbsnp=`ls ./dbsnp/*.vcf.gz`
        zcat "\${dbsnp}" > ./\${dbsnp%.gz}
        dbsnp=\${dbsnp%.gz}
        gatk --java-options -Xmx40g \
            IndexFeatureFile \
            -I \${dbsnp}
        gatk --java-options -Xmx40g \
            GenotypeGVCFs \
                -R "\${fasta}" \
                --dbsnp "\${dbsnp}" \
                -V "${gvcf}" \
                -new-qual -G StandardAnnotation \
                -O ${sampleID}.vcf.gz \
                -isr INTERSECTION
        """
}