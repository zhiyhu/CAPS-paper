module load python3-cbrg

snakemake -p  --rerun-incomplete --cluster "qsub -o std -e std" -j 32
