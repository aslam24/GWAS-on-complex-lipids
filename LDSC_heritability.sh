MAX_JOBS=200

launch()
{
    while [ $(jobs | wc -l) -ge $MAX_JOBS ]
    do
        # at least $MAX_JOBS are still running.
        sleep 2 
    done
    "$@" &
}
for x in `awk '{print$1}' OFS='\t' RS1.EA.Nightingale.linkfile.20250226.Imtiaz.txt`

do
source /docker_environment 
source /opt/miniconda/bin/activate ldsc
launch  /opt/ldsc-master/ldsc.py --h2 ${x}.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ${x}_h2

done 
wait
