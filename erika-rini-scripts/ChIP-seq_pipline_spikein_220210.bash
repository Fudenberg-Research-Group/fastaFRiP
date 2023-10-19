# This script takes four inputs. 
# First, provide the path to the fastq file for the ChIP sample. 
# Second, provide the prefix for the ChIP output files.
# Third provide the path to the fastq file for the input sample.
# Fourth provide the prefix for the input files.
# The script must be run from the parent directory where the ChIP and input folders are housed.
# The index used is "/wynton/home/nora/erikaa/genome_sequences/mm10_hg38/mm10_hg38."
# Check the .stats file for information about duplicates, the number of reads, and the scaling factor as well as UCSC track definitions.


index=/wynton/home/nora/erikaa/genome_sequences/mm10_hg38/mm10_hg38
	
#Load the already-installed modules that I'll need.
module load CBI bowtie2 samtools fastqc

ChIP_mapping () {
#Define the commandline inputs as variables within the script.
local input=$1
local sample=$2

mkdir $sample
local sampleID=$sample/$sample

#Run fastqc.
fastqc --outdir $sample $input

#Print "aligning reads" (not necessary, but can sometimes be useful to see how far the script got before giving an error)
echo "aligning reads"

bowtie2 -p 6 -x $index -U $input > $sampleID.sam

echo "filtering, sorting, and indexing bam"
#Keep reads with mapq > 30. For a mapping with mapq 30, there's a 0.1% chance the read truly maps somewhere else.
# -h includes header.
samtools view -h -q 30 $sampleID.sam > $sampleID.q30.bam

#delete unfiltered sam file
rm $sampleID.sam

samtools sort $sampleID.q30.bam -o $sampleID.q30.sort.bam
samtools index $sampleID.q30.sort.bam

rm $sampleID.q30.bam

samtools markdup -f $sampleID.stats -r -d 100 $sampleID.q30.sort.bam $sampleID.q30.dedup.sam

rm $sampleID.q30.sort.bam

#Count the number of reads that map to each genome and print the ratio
local mmreads=$(( `grep -c '.*mm.*' $sampleID.q30.dedup.sam` - 22 ))
local hgreads=$(( `grep -c '.*hg.*' $sampleID.q30.dedup.sam` - 26 ))
echo -e "\nmm reads" $mmreads >> $sampleID.stats
echo "hg38 reads" $hgreads >> $sampleID.stats
echo "ratio of mouse to human reads is" >> $sampleID.stats
echo "scale=2; $mmreads/$hgreads" | bc >> $sampleID.stats

#Use grep to create a sam file with only the mm reads.
echo "separating reads by species"
grep -v '.*hg.*' $sampleID.q30.dedup.sam > $sampleID.q30.mm.sam

#In the mm file, remove mm from chromosome names.
#This line needs to be updated if the index is not mm10.
sed 's/mm10_chr/chr/' $sampleID.q30.mm.sam > $sampleID.q30.mm.chr.sam
#Delete big intermediate file
rm $sampleID.q30.mm.sam

samtools sort $sampleID.q30.mm.chr.sam -o $sampleID.q30.sort.bam
samtools index $sampleID.q30.sort.bam
}

scaling () {
	local ChIP_sampleID=$1/$1
	local input_sampleID=$2/$2

#Calculate the scaling factor based on the number of mouse and human reads from input and human reads from the ChIP sample.
#The formula comes from Fursova...Klose 2019. The 15000000 is equivalent to alpha in their formula -- it makes the scaling factors be roughly 1 for convenience.
	local ChIP_mmreads=$(( `grep -c '.*mm.*' $ChIP_sampleID.q30.dedup.sam` - 22 ))
	local ChIP_hgreads=$(( `grep -c '.*hg.*' $ChIP_sampleID.q30.dedup.sam` - 26 ))
	local input_mmreads=$(( `grep -c '.*mm.*' $input_sampleID.q30.dedup.sam` - 22 ))
	local input_hgreads=$(( `grep -c '.*hg.*' $input_sampleID.q30.dedup.sam` - 26 ))
	local factor=`echo "scale=20; $input_hgreads / $input_mmreads / $ChIP_hgreads * 15000000" | bc`
	echo "Scaling factor is" >> $ChIP_sampleID.stats
	echo $factor >> $ChIP_sampleID.stats

#Activate my python virtual environment, which has deepTools installed.
	source /wynton/home/nora/erikaa/python_env/bin/activate
	bamCoverage -b $ChIP_sampleID.q30.sort.bam -o $ChIP_sampleID.scale.bw -of bigwig --binSize 20 --scaleFactor $factor
	
	if [ $USER = erikaa ];
		then 
		name=Erika
		elif [ $USER = rini ];
		then 
		name=Rini
		elif [ $USER = epnora ];
		then 
		name=Elphege
		elif [ $USER = khansen1 ];
		then 
		name=Karissa
	fi
	echo 'track type=bigWig name="'$1'.scale" visibility=full description="'$1'" windowingFunction=maximum maxHeightPixels=75,75,75 bigDataUrl=https://storage.googleapis.com/ucscbrowsertracks/'$name'/'$1'.scale.bw' >> $ChIP_sampleID.stats
}


#If the input sample has not already been mapped, first process it, then process the ChIP sample, then create the scaled bigwig.
#Otherwise, process the ChIP sample and calculate the scaling with the already existing input sample.
if [ ! -e $4/$4.q30.mm.chr.sam ]
then
	ChIP_mapping $3 $4
	ChIP_mapping $1 $2
	scaling $2 $4
	else
	ChIP_mapping $1 $2
	scaling $2 $4
fi
