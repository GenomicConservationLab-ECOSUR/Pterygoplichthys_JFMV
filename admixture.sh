#!/bin/bash

# create folders

#mkdir -p ./data/admixture
#mkdir -p ./data/Resultado_plecos

# we run admixture to k 1-5

#for K in 1 2 3 4 5 6 7 8 9 10; do 
#admixture --cv ./data/nm.plink.bed $K | tee ./data/admixture/103Inds_log${K}.out;
#done

# we move the Q and P outputs to an output (folder) 
#mv {*.P,*.Q} ./Resultado_plecos/

## we save the testability results for each k in a file

# grep -h CV ./admixture/Filtro_LD_log*.out > ./admixture/Pleco_Reducido_Kerror.tx

