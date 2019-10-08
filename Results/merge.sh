#cat results_BEAST.txt | awk '{print $1,$2,$4,$5,$3,$6,$7,$8;}' > results_all.txt
#cat results_LSD_LF_LgD.txt >> results_all.txt

grep "rmse" results_BEAST.txt | awk '{print $1,$2,$4,$5,$3,$6,$8;}' > results_all_2019Oct06.txt
cat results_2019Oct06.txt >> results_all_2019Oct06.txt
