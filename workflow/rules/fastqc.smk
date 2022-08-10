o2i_fastqc = {
    # from odir to idir
    # idirs = ['fastq', 'results/trimmed_reads'] # locate fastq files
    # odirs = ['results/fastqc_raw', 'results/fastqc_trimmed'] # fastqc outputs
    "results/fastqc_raw": "fastq",
    "results/fastqc_trimmed": "results/trimmed_reads",
}

rs = ["R1", "R2"]  # suffix


rule fastqc:
    input:
        lambda wildcards: "{id}/{s}.{r}.fastq.gz".format(
            id=o2i_fastqc[wildcards.d], s=wildcards.s, r=wildcards.r
        ),
    output:
        "{d}/details/{s}.{r}_fastqc.zip",
    params:
        d="{d}/details",
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    threads: 1
    log:
        "{d}/fastqc/details/{s}.{r}_fastqc.log",
    benchmark:
        "{d}/fastqc/details/{s}.{r}_fastqc.benchmark"
    shell:
        "which fastqc &> {log};"
        "fastqc -t {threads} {input} -o {params.d} &>> {log};"


rule fastqc_multiqc:
    input:
        expand("{{d}}/details/{s}.{r}_fastqc.zip", s=samples, r=rs),
    output:
        "{d}/fastqc_report.html",
    params:
        d="{d}",
    conda:
        "../envs/multiqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    threads: 1
    log:
        "{d}/fastqc_report.log",
    benchmark:
        "{d}/fastqc_report.benchmark"
    shell:
        """
        which multiqc &> {log};
        multiqc {input} -f -o {params.d} --title {wildcards.d} -n fastqc_report &>> {log};
        """
