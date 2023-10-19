#has to be run from python virtual environment
#First argument is the path to the bam file
#Second argument is the label
#Third argument is the bed file with the peaks

import deeptools.countReadsPerBin as crpb
import pysam
import sys

bam = str(sys.argv[1])
label = str(sys.argv[2])
peaks = str(sys.argv[3])


cr = crpb.CountReadsPerBin([bam], bedFile=peaks, numberOfProcessors=10)
reads_at_peaks = cr.run()
total = reads_at_peaks.sum(axis=0)

a = pysam.AlignmentFile(bam)
frip = float(total[0])/a.mapped

file_name = label + "_frip.txt"

file = open(file_name,"a")
file.write(label + "_frip\n")
file.write(str(frip)+ "\n")
file.close()
