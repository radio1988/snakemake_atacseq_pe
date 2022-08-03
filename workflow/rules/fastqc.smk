rule fastqc:
    input:
        expand("fastq/{sample}.{r}.fastq.gz", sample=config['SAMPLES'], r=["R1", "R2"])
    output:
        ["results/fastqc/details/{sample}.R1_fastqc.html", 
        "results/fastqc/details/{sample}.R2_fastqc.html"], 
        ["results/fastqc/details/{sample}.R1_fastqc.zip", 
        "results/fastqc/details/{sample}.R2_fastqc.zip"]
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        1
    log:
        'results/fastqc/details/{sample}.log'
    benchmark:
        'results/fastqc/details/{sample}.benchmark'
    shell:
        "which fastqc &> {log};"
        "fastqc -t {threads} {input} -o results/fastqc/details &>> {log};"


rule fastqc_multiqc:
    input:
        expand("results/fastqc/details/{sample}.{r}_fastqc.zip", \
                sample=config['SAMPLES'], r=["R1", "R2"])
    output:
        "results/fastqc/multiqc_report.html"
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
        multiqc {input} -f -o results/fastqc &>> {log};
        '''
