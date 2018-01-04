#/bin/bash
#$ -pe smp 56-144
#$ -cwd -V
#$ -o optimization.log
#$ -j y
#$ -N Optimization
#$ -q all.q



source "/scicomp/home/nft9/tools/python/python3-6.1/bin/activate"
NSLOTS=${NSLOTS:=45}

python start.py --d $NSLOTS --g ecoli --m 0.90 --f "./allele_calls/ecoli_allele_calls_all.csv" --t "./thresholds/ecoli_thresholds_multi_mod.txt" --v "./allele_calls/ecoli_views.json" --o "./ecoli_results" --s "core"


if [ $? -gt 0 ]; then
	echo -e "Error with optimzation\n";
	exit 1;
fi
