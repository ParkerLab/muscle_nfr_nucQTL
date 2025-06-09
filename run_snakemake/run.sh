pipeline=$1
# add -F will be the forcerun
snakemake all  --jobname "${pipeline}.{jobid}" --jobs 1000 \
                --keep-going \
                --rerun-incomplete \
                --snakefile ../smk/${pipeline}.smk \
                --configfile ../smk/${pipeline}.yaml \
                --use-conda \
                --conda-prefix ~/miniconda3 \
                --use-singularity \
                --singularity-prefix "~/singularity-cache" \
                --singularity-args "--bind /scratch" \
                --use-envmodules \
                --printshellcmds \
                --cluster-config ../smk/${pipeline}_cluster.yaml \
                --cluster "sbatch --account {cluster.account} --output {cluster.output} --time {cluster.time} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mem-per-cpu {cluster.mem} --parsable -p largemem,standard" \
                --cluster-status ../run_snakemake/slurm_status.py \
                --latency-wait 480 \
                --rerun-triggers mtime \
                > logs/${pipeline}_snakemake.log 2>&1 &