pipeline=$1
# add -F will be the forcerun
snakemake -n all --jobname "${pipeline}.{jobid}" --jobs 1000 \
                --keep-going \
                --rerun-incomplete \
                --snakefile ../smk/${pipeline}.smk \
                --configfile ../smk/${pipeline}.yaml \
                --use-conda \
                --conda-prefix ~/miniconda3/envs \
                --use-singularity \
                --singularity-prefix "~/singularity-cache" \
                --singularity-args "--bind /scratch" \
                --use-envmodules \
                --printshellcmds \
                --cluster-config ../smk/${pipeline}_cluster.yaml \
                --cluster "sbatch --account {cluster.account} --output {cluster.output} --time {cluster.time} --mem {cluster.mem} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --parsable -p largemem,standard" \
                --cluster-status ../run_snakemake/slurm_status.py \
                --latency-wait 60 \
                --rerun-triggers mtime \
                --touch \
                > logs/${pipeline}_dryrun_snakemake.log 2>&1 &