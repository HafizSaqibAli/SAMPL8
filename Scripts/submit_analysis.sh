#!/bin/bash --login
#$ -cwd
#$ -N run_XXX_analysis
#$ -pe smp.pe 4   # Each task will use X cores
#$ -l short
#$ -hold_jid run_XXX_frames



module load apps/anaconda3/5.2.0

### change path here
EXE=$HOME/bin/POSEIDON_beta/objectAnalysis/objectIterations.py
	

#arrays="1 2 3 4 5 6 7 8 9 10"
arrays="1"
folders=( `seq 1 50 500` )

mkdir frame_analysis
cd frame_analysis


last_paths=""

for a in ${arrays}; do
	dir=""
	#paths1=""
	path1="../${a}job_array/"


	for index in ${!folders[*]}; do

		dir+="${folders[$index]}",
		objs=`echo ${dir} | sed 's/.$//'` #remove last comma
		#echo ${objs}
		last=`echo ${last_paths} | tr '?' ' '` #replace ? with space

	done

	objs2=`echo ${dir} | sed 's/.$//'` #remove last comma
	last_paths+="${path1}{${objs2}}?"
done

EXE2="${EXE} -p ${last} ${path1}{${objs}} -SE > file.log"
echo ${EXE2}
$(eval ${EXE2} )


cd ..

