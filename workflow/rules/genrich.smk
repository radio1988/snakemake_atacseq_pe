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

              
rule genrich_group:
    input:
        lambda wildcards: ["results/genrich/cleanBamByName/{}.bam".format(s) for s in g2s[wildcards.group]],
    output:
        peak="results/genrich_group/{group}.narrowPeak",
        bedGraph=temp("results/genrich_group/{group}.bedgraph_ish")
    log:
        "results/genrich_group/{group}.narrowPeak.log",
    benchmark:
        "results/genrich_group/{group}.narrowPeak.benchmark"
    params:
        samples=lambda wildcards: ','.join(
            ["results/genrich/cleanBamByName/{}.bam".format(s) for s in g2s[wildcards.group]]
            ),
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
        -t {params.samples} -o {output.peak} -k {output.bedGraph} \
        {params.options} \
        -m {params.MQ_MIN} -e {params.chrM} -E {params.BLACKLIST} &>>{log}"

def genrich_contrast_input(wildcards):
    '''
    return: 
    [results/genrich/cleanBamByName/Lrig1.bam, 
    results/genrich/cleanBamByName/Lrig1b.bam, 
    results/genrich/cleanBamByName/WT_RSmad7.bam, 
    results/genrich/cleanBamByName/WT_RSmad7b.bam]
    '''
    files = []
    for group in c2g[wildcards['contrast']]:
        for sample in g2s[group]:
            files.append("results/genrich/cleanBamByName/{}.bam".format(sample))
    return files

def genrich_contrast_treatments(wildcards):
    '''
    return: "results/genrich/cleanBamByName/Lrig1.bam,results/genrich/cleanBamByName/Lrig1b.bam"
    '''
    files = []
    group = c2g[wildcards['contrast']][0]
    for sample in g2s[group]:
        files.append("results/genrich/cleanBamByName/{}.bam".format(sample))    
    param_str = ",".join(files)
    return param_str

def genrich_contrast_controls(wildcards):
    '''
    return: "results/genrich/cleanBamByName/WT_RSmad7.bam,results/genrich/cleanBamByName/WT_RSmad7b.bam"
    '''
    files = []
    group = c2g[wildcards['contrast']][1]
    for sample in g2s[group]:
        files.append("results/genrich/cleanBamByName/{}.bam".format(sample))    
    param_str = ",".join(files)
    return param_str

        
rule genrich_contrast:
    input:
        genrich_contrast_input
    output:
        peak="results/genrich_contrast/{contrast}.narrowPeak",
        bedGraph=temp("results/genrich_contrast/{contrast}.bedgraph_ish")
    log:
        "results/genrich_contrast/{contrast}.narrowPeak.log",
    benchmark:
        "results/genrich_contrast/{contrast}.narrowPeak.benchmark"
    params:
        treatments=genrich_contrast_treatments,
        controls=genrich_contrast_controls,
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
        -t {params.treatments} -c {params.controls} -o {output.peak} -k {output.bedGraph} \
        {params.options} \
        -m {params.MQ_MIN} -e {params.chrM} -E {params.BLACKLIST} &>>{log}"
