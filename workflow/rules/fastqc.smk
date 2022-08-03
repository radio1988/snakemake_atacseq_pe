#rule fastqc:
#    input:
#    output:
#        "results/fastqc/multiqc_report.html"
#    log:
#        "results/fastqc/fastqc.log"
#    params:
#        odir="results/fastqc"
#    threads:
#        4
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt * 2000
#    conda:
#        "../envs/fastqc.yaml"
#    shell:
#        """
#        mkdir -p results/fastqc
#        mkdir -p results/fastqc/details
#        fastqc -t {threads} {input} -o {params.odir}/details &> {log}
#        multiqc {params.odir}/details -o {params.odir} &>> {log}
#        """

rule fastqc:
    input:
        expand("fastq/{sample}.{r}.fastq.gz", sample=config['SAMPLES'], r=["R1", "R2"])
    output:
        ["results/fastqc/details/{sample}.R1_results/fastqc.html", 
        "results/fastqc/details/{sample}.R2_results/fastqc.html"], 
        ["results/fastqc/details/{sample}.R1_results/fastqc.zip", 
        "results/fastqc/details/{sample}.R2_results/fastqc.zip"]
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
        expand("results/fastqc/details/{sample}.{r}_results/fastqc.zip", \
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
        """
        which multiqc &> {log};
        multiqc {input} -f -o fastqc &>> {log};
        """
