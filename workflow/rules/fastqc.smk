dirs = ['fastq', 'results/trimmed_reads'] # locate fastq files
rs = ['R1', 'R2']  # suffix

rule fastqc:
    input:
        '{d}/{s}.{r}.fastq.gz' 
    output:
        '{d}/fastqc/details/{s}.{r}_fastqc.zip'
    params:
        d='{d}'
    conda:
        '../envs/fastqc.yaml'
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        1
    log:
        '{d}/fastqc/details/{s}.{r}_fastqc.log'
    benchmark:
        '{d}/fastqc/details/{s}.{r}_fastqc.benchmark'
    shell:
        'which fastqc &> {log};'
        'fastqc -t {threads} {input} -o results/fastqc_raw/details &>> {log};'


rule fastqc_multiqc:
    input:
        expand('{{d}}/fastqc/details/{s}.{r}_fastqc.zip', 
                s=config['SAMPLES'], 
                r=rs)
    output:
        '{d}/fastqc/multiqc_report.html'
    params:
        d='{d}'
    conda:
        '../envs/multiqc.yaml'
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        1
    log:
        '{d}/multiqc_report.log'
    benchmark:
        '{d}/multiqc_report.benchmark'
    shell:
        '''
        which multiqc &> {log};
        multiqc {input} -f -o {params.d} &>> {log};
        '''
