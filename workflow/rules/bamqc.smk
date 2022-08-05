SAMPLES = config["SAMPLES"]
dirs = ["sorted_reads", "clean_reads", "markDup"]


rule bamqc:
    input:
        bam="results/{d}/{sample}.bam",
        bai="results/{d}/{sample}.bam.bai",
    output:
        idxstats="results/{d}_qc/idxstats/{sample}.idxstats.txt",
        flagstat="results/{d}_qc/flagstat/{sample}.flagstat.txt",
        stats="results/{d}_qc/stats/{sample}.stats.txt",
    log:
        idxstats="results/{d}_qc/idxstats/{sample}.idxstats.log",
        flagstat="results/{d}_qc/flagstat/{sample}.flagstat.log",
        stats="results/{d}_qc/stats/{sample}.stats.log",
    benchmark:
        "results/{d}_qc/idxstats/{sample}.idxstats.benchmark"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3000,
    threads: 1
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools idxstats {input.bam} > {output.idxstats} 2> {log.idxstats} 
        samtools flagstat {input.bam} > {output.flagstat} 2> {log.flagstat} 
        samtools stats {input.bam} > {output.stats} 2> {log.stats} 
        """


rule bamqc_multiqc:
    input:
        stats=expand("results/{{d}}_qc/stats/{sample}.stats.txt", sample=SAMPLES),
        idxstats=expand(
            "results/{{d}}_qc/idxstats/{sample}.idxstats.txt", sample=SAMPLES
        ),
        flagstat=expand(
            "results/{{d}}_qc/flagstat/{sample}.flagstat.txt", sample=SAMPLES
        ),
    output:
        "results/{d}_qc/bamqc_report.html",
    log:
        "results/{d}_qc/bamqc_report.log",
    benchmark:
        "results/{d}_qc/bamqc_report.benchmark"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3000,
    threads: 1
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        multiqc -f {input.stats} {input.idxstats} {input.flagstat} -o results/{wildcards.d}_qc --title {wildcards.d} -n bamqc_report &> {log}
        """
