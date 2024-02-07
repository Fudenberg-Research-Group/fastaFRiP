from optparse import OptionParser
import pysam
import deeptools.countReadsPerBin as crpb

def main():
    usage = 'usage: %prog [options] <bam_file1> <bed_file>'
    parser = OptionParser(usage)
    parser.add_option('--nproc', dest='nproc',
        default=1,
        help='How many processes to use for calculation. Default: nproc=1')
    parser.add_option('--output-prefix', dest='prefix',
        default='',
        help='What is the prefix of output file')
    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide BAM and BED files.')
    else:
        bam = args[0]
        bed = args[1]
        
    bam_files = [bam]
    reads_counter = crpb.CountReadsPerBin(bam_files, bedFile=bed, numberOfProcessors=options.nproc)
    reads_at_peaks = reads_counter.run()
    total_reads_at_peaks = reads_at_peaks.sum(axis=0)

    alignment = pysam.AlignmentFile(bam)
    frip = float(total_reads_at_peaks[0])/alignment.mapped

    file_name = options.prefix + "_frip.txt"

    file = open(file_name,"a")
    file.write(options.prefix + "_frip\n")
    file.write(str(frip)+ "\n")
    file.close()
if __name__ == '__main__':
    main()
