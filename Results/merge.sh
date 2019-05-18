cat results_BEAST.txt | awk '{print $1,$2,$4,$5,$3,$6,$7,$8;}' > results_all.txt
cat results_LSD_LF_LgD.txt >> results_all.txt
