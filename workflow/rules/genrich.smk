rule genrich_sort_bam:
    input:
        "results/clean_reads/{sample}.bam",
    output:
        temp("results/genrich/cleanBamByName/{sample}.bam"),
    log:
        "results/genrich/cleanBamByName/{sample}.bam.log",
    benchmark:
        "results/genrich/cleanBamByName/{sample}.bam.benchmark"
    conda:
        "../envs/samtools.yaml"
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2500
    shell:
        "samtools sort -n -@ 2 -m 2G {input} -o {output} &> {log}"


rule genrich_sample:
    input:
        "results/genrich/cleanBamByName/{sample}.bam",
    output:
        peak="results/genrich/{sample}.narrowPeak",
        bedGraph=temp("results/genrich/{sample}.bedgraph_ish")
    log:
        "results/genrich/{sample}.narrowPeak.log",
    benchmark:
        "results/genrich/{sample}.narrowPeak.benchmark"
    params:
        options=config["GENRICH"], # default:  '-j -y -d 100 -q 0.05 -a 20.0'
        MQ_MIN=config["MQ_MIN"],
        chrM=config["chrM"],
        BLACKLIST=config['BLACKLIST']
    conda:
        "../envs/genrich.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000,
    shell:
        "pwd > {log}; Genrich  -v \
        -t {input} -o {output.peak} -k {output.bedGraph} \
        {params.options} \
        -m {params.MQ_MIN} -e {params.chrM} -E {params.BLACKLIST} &>>{log}"

# rule genrich_group:
# rule genrich_group_with_control:
