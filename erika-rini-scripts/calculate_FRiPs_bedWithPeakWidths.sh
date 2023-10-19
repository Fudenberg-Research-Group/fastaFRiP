#First argument is the path to the bam file
#Second argument is the label
#Third argument is the bed file with the peaks.


bed=$3

source /wynton/home/nora/erikaa/python_env/bin/activate

python3 /wynton/home/nora/erikaa/scripts/calculate_FRiPs.py $1 $2 $bed
