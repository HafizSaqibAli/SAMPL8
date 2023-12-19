#!/bin/bash --login
#$ -cwd


system=complex

cp submit_analysis.sh submit_${system}_analysis.sh

sed -i "s/XXX/${system}/g" "submit_${system}_analysis.sh"



qsub submit_${system}_analysis.sh
