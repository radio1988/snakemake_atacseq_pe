rule markDup:
    input:
        bam="results/sorted_reads/{sample}.bam",
        bai="results/sorted_reads/{sample}.bam.bai",
    output:
        bam=temp("results/markDup/{sample}.bam"),
        metrics="results/markDup/{sample}.markDup_metrics.txt",
    log:
        "results/markDup/{sample}.bam.log",
    benchmark:
        "results/markDup/{sample}.bam.benchmark"
    conda:
        "../envs/picard.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 12000,
    shell:
        """
        p=`which picard` && echo $p &> {log};
        PICARD=`echo $p|sed 's/bin\/picard/share\/picard-2.27.4-0\/picard.jar/'` && echo $PICARD &>>{log}; 
        java -Xmx8g -jar $PICARD \
        MarkDuplicates \
        I={input.bam} \
        O={output.bam} \
        M={output.metrics} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true \
        &>> {log}
        """


rule mark_duplicates:
    input:
        bams="results/sorted_reads/{sample}.bam",
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam="results/dedup/{sample}.bam",
        metrics="results/dedup/{sample}.metrics.txt",
    log:
        "results/dedup/{sample}.bam.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000,
    wrapper:
        "v1.7.1/bio/picard/markduplicates"


rule atacseq_bam_filter:
    # remove mito reads
    # keep paired reads with MQ>20 and 38-2000nt fragment size only
    input:
        "results/markDup/{sample}.bam",
    output:
        bam="results/clean_reads/{sample}.bam",
    params:
        chrM=config["chrM"],
    log:
        "results/clean_reads/{sample}.bam.log",
    benchmark:
        "results/clean_reads/{sample}.bam.benchmark"
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3000,
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        which samtools &> {log}
        samtools view -h {input} 2>>{log} | \
        perl -lane 'print unless ($F[2] eq {params.chrM} and $_ != /\@/)' 2>>{log} | \
        awk \'{config[FILTER]}\' 2>>{log} | \
        samtools sort -@ {threads} -m 2G -o {output.bam}  2>> {log}
        """
