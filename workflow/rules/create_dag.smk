rule create_dag:
    output:
        dag='results/workflow.dag'
    log:
        'results/workflow.dag.log'
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    shell:
        'snakemake --dag targets > {output.dag} 2> {log}'


rule dag2svg:
    input:
        'results/workflow.dag'
    output:
        svg='results/workflow.svg'
    log:
        'results/workflow.svg.log'
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    conda:
        '../envs/graphviz.yaml'
    shell:
        'which dot &> {log} ;'
        'cat {input} | dot -Tsvg > {output.svg} 2>> {log}'
