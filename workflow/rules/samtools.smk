rule samtools_sort_index:
    input:
        "results/mapped_reads/{s}.bam",
    output:
        bam="results/sorted_reads/{s}.bam",
        bai="results/sorted_reads/{s}.bam.bai",
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
        samtools index {output.bam} {output.bai} &>> {log}
        """
