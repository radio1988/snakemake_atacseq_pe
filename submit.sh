# bsub -q long -W 144:00 -R rusage[mem=4000]  -R select[rh=8] 'bash submit.sh'
rm -f lsf.log

source activate snake7

snakemake -pk --jobs 99 \
--notemp \
--use-conda --conda-prefix ~/anaconda/envs/snakemake_atacseq_pe \
--latency-wait 20 --ri --restart-times 0 \
--cluster 'bsub -q long -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -R select[rh=8] -W 140:00' \
&> workflow.log
