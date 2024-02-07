#!/bin/bash
set -e; shopt -s expand_aliases

script=checkProtein.sh

PATH_TO_SCRIPTS=$(pwd)
cd $(lsvol | cut -d " " -f 4)/BAPLargeScale

# final value of sequence is 5311
for i in $(seq 1 10 5311); do
cp ${PATH_TO_SCRIPTS}/${script} ${script}_${i}
sed -i ${script}_${i} -e "s/#tbr"${i}"PDBids/PDBids/g"
sbatch ${script}_${i}
done
