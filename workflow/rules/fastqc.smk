rule fastqc_raw:
    input:
        expand("fastq/{sample}.{r}.fastq.gz", sample=config['SAMPLES'], r=["R1", "R2"])
    output:
        ["results/fastqc_raw/details/{sample}.R1_fastqc.html", 
        "results/fastqc_raw/details/{sample}.R2_fastqc.html"], 
        ["results/fastqc_raw/details/{sample}.R1_fastqc.zip", 
        "results/fastqc_raw/details/{sample}.R2_fastqc.zip"]
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        1
    log:
        'results/fastqc_raw/details/{sample}.log'
    benchmark:
        'results/fastqc_raw/details/{sample}.benchmark'
    shell:
        "which fastqc &> {log};"
        "fastqc -t {threads} {input} -o results/fastqc_raw/details &>> {log};"


rule fastqc_raw_multiqc:
    input:
        expand("results/fastqc_raw/details/{sample}.{r}_fastqc.zip", \
                sample=config['SAMPLES'], r=["R1", "R2"])
    output:
        "results/fastqc_raw/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        1
    log:
        "results/fastqc_raw/multiqc_report.log"
    benchmark:
        "results/fastqc_raw/multiqc_report.benchmark"
    shell:
        '''
        which multiqc &> {log};
        multiqc {input} -f -o results/fastqc_raw &>> {log};
        '''


rule fastqc_trimmed:
    input:
        expand("results/trimmed_reads/{sample}.{r}.fastq.gz", sample=config['SAMPLES'], r=["R1", "R2"])
    output:
        ["results/fastqc_trimmed/details/{sample}.R1_fastqc.html", 
        "results/fastqc_trimmed/details/{sample}.R2_fastqc.html"], 
        ["results/fastqc_trimmed/details/{sample}.R1_fastqc.zip", 
        "results/fastqc_trimmed/details/{sample}.R2_fastqc.zip"]
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        1
    log:
        'results/fastqc_trimmed/details/{sample}.log'
    benchmark:
        'results/fastqc_trimmed/details/{sample}.benchmark'
    shell:
        "which fastqc &> {log};"
        "fastqc -t {threads} {input} -o results/fastqc_trimmed/details &>> {log};"


rule fastqc_multiqc_trimmed:
    input:
        expand("results/fastqc_trimmed/details/{sample}.{r}_fastqc.zip", \
                sample=config['SAMPLES'], r=["R1", "R2"])
    output:
        "results/fastqc_trimmed/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        1
    log:
        "results/fastqc/multiqc_report.log"
    benchmark:
        "results/fastqc/multiqc_report.benchmark"
    shell:
        '''
        which multiqc &> {log};
        multiqc {input} -f -o results/fastqc_trimmed &>> {log};
        '''
