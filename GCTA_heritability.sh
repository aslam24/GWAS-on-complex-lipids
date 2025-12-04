MAX_JOBS=50

launch()
{
    while [ $(jobs | wc -l) -ge $MAX_JOBS ]
    do
        # at least $MAX_JOBS are still running.
        sleep 2 
    done
    "$@" &
}

for x in {1..1105}

 do
launch ~/software/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static --grm GCTA_GWAS_lipids_SAMPLES --pheno GWAS_LIPIDS_M3.txt --mpheno ${x} --reml --out GWAS_LIPID_HER_${x} --thread-num 5
done 
