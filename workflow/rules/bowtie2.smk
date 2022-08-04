rule bowtie2_index:
    input:
        config['GENOME']
    output:
        touch(config['GENOME']+'.bowtie2.index.done')
    log:
        config['GENOME']+'.bowtie2.index.log'
    benchmark:
        config['GENOME']+'.bowtie2.index.benchmark'
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000
    threads:
        8
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build --version &> {log}
        bowtie2-build --threads {threads} {input} {input} &>> {log}
        """

rule bowtie2:
    input:
        genome=config['GENOME'],
        index_flag=config['GENOME']+'.bowtie2.index.done',
        reads=["results/trimmed_reads/{sample}.R1.fastq.gz", "results/trimmed_reads/{sample}.R2.fastq.gz"]
    output:
        temp("results/mapped_reads/{sample}.bam")
    params:
        reads="-1 fastq/{sample}.R1.fastq.gz -2 fastq/{sample}.R2.fastq.gz",
        options=config['BOWTIE2_PARAMS']
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1600 # human need 18G
    threads:
        12
    log:
        "results/mapped_reads/{sample}.bowtie2.log"
    benchmark:
        "results/mapped_reads/{sample}.bowtie2.benchmark"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2 --version &> {log}
        bowtie2 -x {input.genome} -p {threads} {params.reads} {params.options} 2>> {log}  | \
        samtools view -Sb -1 -@ 2 - -o {output} 2>> {log}
        """


rule sam_sort_index:
    input:
        "results/mapped_reads/{sample}.bam"
    output:
        bam="results/sorted_reads/{sample}.bam",
        bai="results/sorted_reads/{sample}.bam.bai"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2500
    threads:
        4
    log:
        "log/samtools_sort/{sample}.sort.log"
    benchmark:
        "log/samtools_sort/{sample}.sort.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools --version &> {log}
        samtools sort -@ {threads} -m 2G {input} -o {output.bam} &>> {log}
        samtools index {output.bam} {output.bai} &>> {log}
        """
