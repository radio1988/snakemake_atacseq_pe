rule create_dag:
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        1
    output:
        svg='results/workflow.svg'
    log:
        'results/workflow.svg.log'
    conda:
        '../envs/create_dag.yaml'
    shell:
        'which dot &> {log};'
        'snakemake --dag targets | dot -Tsvg > {output.svg} 2>> {log}'
