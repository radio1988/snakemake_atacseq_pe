rule trimmomatic:
    input:
        r1="fastq/{sample}.R1.fastq.gz",
        r2="fastq/{sample}.R2.fastq.gz"
    output:
        r1="results/trimmed_reads/{sample}.R1.fastq.gz",
        r2="results/trimmed_reads/{sample}.R2.fastq.gz",
        u1="results/trimmed_reads/{sample}.R1.unpaired.fastq.gz",
        u2="results/trimmed_reads/{sample}.R2.unpaired.fastq.gz"
    params:
        config['TRIMMOMATIC']
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 5000  # human need 18G
    threads:
        4
    log:
        "results/trimmed_reads/{sample}.log",
    benchmark:
        "results/trimmed_reads/{sample}.benchmark",
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        """
        trimmomatic \
        PE -phred33 \
        {input.r1} {input.r2} \
        {output.r1} {output.u1} \
        {output.r2} {output.u2} \
        {params}
        """
