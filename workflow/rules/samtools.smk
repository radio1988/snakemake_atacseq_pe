rule samtools_location_sort:
    input:
        "results/mapped_reads/{s}.bam",
    output:
        bam="results/sorted_reads/{s}.bam",
    log:
        "results/sorted_reads/{s}.bam.log",
    benchmark:
        "results/sorted_reads/{s}.bam.benchmark"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2500,
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools --version &> {log}
        samtools sort -@ {threads} -m 2G {input} -o {output.bam} &>> {log}
        """


rule samtools_index:
    input:
        "{file}.bam",
    output:
        "{file}.bam.bai",
    log:
        "{file}.bam.bai.log",
    benchmark:
        "{file}.bam.bai.benchmark"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2500,
    threads: 1
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools --version &> {log}
        samtools index {input} {output} &>> {log}
        """
