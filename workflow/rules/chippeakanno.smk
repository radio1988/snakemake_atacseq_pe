rule chippeakanno:
    input:
        peak="results/genrich_contrast/{contrast}.narrowPeak",
        gtf=config['GTF'],
    output:
        "results/genrich_contrast/{contrast}.narrowPeak.full_anno.xlsx",
    params:
        CHIPPEAKANNO_MODE=config['CHIPPEAKANNO_MODE'],
        BIDING_LEFT=config['BIDING_LEFT'],
        BIDING_RIGHT=config['BIDING_RIGHT'],
        PeakLocForDistance=config['PeakLocForDistance'],
        FeatureLocForDistance=config['FeatureLocForDistance'],
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    log:
        "results/genrich_contrast/{contrast}.narrowPeak.full_anno.xlsx.log",

    benchmark:
        "results/genrich_contrast/{contrast}.narrowPeak.full_anno.xlsx.benchmark",
    conda:
        "../envs/chippeakanno.yaml"
    shell:
        """
        Rscript workflow/scripts/chippeakanno.R {input.peak} {input.gtf} \
        {params.CHIPPEAKANNO_MODE} {params.BIDING_LEFT} {params.BIDING_RIGHT} \
        {params.PeakLocForDistance} {params.FeatureLocForDistance} \
        {wildcards.contrast} \
        &> {log}
        """
       
    
    
rule chippeakanno_group:
    input:
        peak="results/genrich_group/{group}.narrowPeak",
        gtf=config['GTF'],
    output:
        "results/genrich_group/{group}.narrowPeak.full_anno.xlsx"
    params:
        CHIPPEAKANNO_MODE=config['CHIPPEAKANNO_MODE'],
        BIDING_LEFT=config['BIDING_LEFT'],
        BIDING_RIGHT=config['BIDING_RIGHT'],
        PeakLocForDistance=config['PeakLocForDistance'],
        FeatureLocForDistance=config['FeatureLocForDistance'],
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    log:
        "results/genrich_group/{group}.narrowPeak.full_anno.xlsx.log"
    benchmark:
        "results/genrich_group/{group}.narrowPeak.full_anno.xlsx.benchmark"
    conda:
        "../envs/chippeakanno.yaml"
    shell:
        """
        Rscript workflow/scripts/chippeakanno.R {input.peak} {input.gtf} \
        {params.CHIPPEAKANNO_MODE} {params.BIDING_LEFT} {params.BIDING_RIGHT} \
        {params.PeakLocForDistance} {params.FeatureLocForDistance} \
        {wildcards.group} \
        &> {log}
        """