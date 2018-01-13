#!/Users/Adrian/miniconda3/bin/python3
# Adrian Bubie
# 1/12/18

## PCR_Deduper: 
## Argparse arguments - sam file (abs path), paired-end flag (boolean), UMI file (abs path), and help
import textwrap
import random
import fileinput
import os
import argparse as ap

def get_arguments():
    parser = ap.ArgumentParser(prog='./PCR_Deduper.py', formatter_class=ap.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    PCR Duplicate Remover
    ---------------------
    Removes duplicate reads from given position sorted SAM file 
    that are product of PCR duplication.
    
    Note that the SAM file input is assumed to be pre-sorted by position, and all reads mapped.
    '''))
    parser.add_argument("-f", help="Sorted SAM file to be processed. Must include absolute path the file. <str>", required=True, type=str)
    parser.add_argument("-p", help="If passed as 'True', SAM data is considered paired-ended. <boolean> (def=False)", required=False, type=bool, default=False)
    parser.add_argument("-u", help="Optional file to define UMI sequences (one UMI per line). UMIs are assumed to be in the last field of the QNAME for each read. Reads with invalid UMIs are filtered from final output. If file not specified, SAM is assumed to use randomers. Must include absolute path to file <str>", required=False, type=str, default='')
    parser.add_argument("-dup_keep", help="Optional argument to designate which read to keep of a duplicate set. 'Q' will keep read with highest quality score, 'R' will choose a random read to keep. No specification will keep the first read encounted. <str>", required=False, type=str, default='')
    parser.add_argument("-s", help="Set to 'True' for additional summary file output <boolean> (def='False')", required=False, type=bool, default=False)
    return parser.parse_args()


class Read:
    """ Read object: object called to get read parts from SAM file lines, or change attributes """
    def __init__(self, line):
        self.chr = line.split('\t')[2] # Get chromosome
        self.flag = int(line.split('\t')[1]) # Get bitwise flag
        self.tag = line.split('\t')[0].split(':')[-1] # Get UMI/randomer tag from Read Name
        self.cig = line.split('\t')[5] # Get the cigar string

        if 'S' in self.cig and self.cig.find('S')+1 < len(self.cig): # If the read is softclipped on the 5'-end (given by CIGAR) we need to adjust the position
            adj = int(self.cig[:self.cig.find('S')])
            self.pos = int(line.split('\t')[3])-adj # Set the position with the adjusted value
        else:
            self.pos = int(line.split('\t')[3]) # No soft-clipping on 5'-end, keep the value of the POS column.

        self.qual = int(line.split('\t')[4]) # Get the quality score of the read
        self.strand = flag_stat(self.flag) # Call the bitwise flag checker and get the strand of the read

        self.chars = [str(self.tag),str(self.strand),str(self.chr),str(self.pos)] # Create read true characteristics list

        self.raw = line # Keep copy of the input line


def umi_list(umi_file):
    '''Returns list of UMI tags from UMI file argument'''

    with open(umi_file, 'r') as ulist:
        return [line.strip('\n') for line in ulist.readlines()]


def flag_stat(flag):
    '''Bitwise flag checker: checks read flag and returns the strand of the read. 
    If the paired-end argument passed as True, which read in the pair the read is 
    checked and returned.'''

    strand = '+'
    if (flag & 16) == 16: # Check strand:
        strand = '-'  # If SEQ is reverse complemented by flag definition, than "-" is returned. Otherwise, "+" is returned.

    return strand


def true_char(file_in):
    '''Takes in the input file and appends a new column with the "true characteristics"; the chr, umi/randomer, strand, and position. 
    Position may be adjusted from what is given in the read based on soft-clipping at the 5'-end (see Read object class) This modifies 
    the input file directly, but will be fixed in the final steps to keep the SAM format of the file correct.'''

    for line in fileinput.input(file_in, inplace=True): # inplace here meaning each line is written back to it's original place in the file
        if line.startswith('@'):  # If a header line, just write back to the file, don't modify
            print(line.strip('\n'))

        else:
            read = Read(line) # take the line and create a Read object from it
            print(line.strip('\n')+'\t'+':'.join(read.chars)) # Add the true characteristics to a new column at the end of each read


def sort_by_true(file_in):
    '''Takes in imput file and resorts it by the "true position" of each read, meaning the positions post-softclipping adjustment
    The SAM file is sorted 100000 lines at a time and written out to an 'Output.sam' temporary file, which will be deduped. The output
    of this function is the newly sorted Output file.'''

    with open(file_in, 'r') as SAM:
        line = SAM.readline() # Read first line of the file

        while line.startswith('@'): # Keep the header lines as they are (write the header to Output). Stop when you hit the first read line
            with open('Output.sam','a') as out:
                out.write(line)
            line = SAM.readline()

        reads_proc = 0 # Counter for number of reads (not including headers) processed

        while line:
            to_sort = [] # List to store read chunk to sort
            i = 0 # indexer

            while i < 100001: # Grab 100000 reads at a time from the file
                to_sort.append(line.split(':')) # Split each read into a list and add to the sort list
                line = SAM.readline() 
                if line == '': # Stop if you get to the end of the file (before you reach the 100000 size limit)
                    break

                i += 1
                reads_proc += 1

            sortee = sorted(sorted(to_sort, key=lambda read: int(read[-1].strip('\n'))), key=lambda read: int(read[-2].strip('\t'))) # Sort the reads based on the true position (last 2 items in each list)

            for i in range(0,len(sortee)):
                sortee[i] = ':'.join(sortee[i]) # Mend each read back together

            with open('Output.sam','a') as out: 
                out.writelines(sortee) # Write the "true" sorted reads to the temp output file Output.sam
    
    return reads_proc


def dup_remover(file_in, file_out, umis, dup_keep):
    ''' Takes in the temp output file from sort_by_true and removes all duplicates by checking the final column (characteristics)
    of adjacent reads, sequentially. If the umi flag was specified and a set of UMIs passed in, reads with invalid umis are removed from the
    final file. The read to keep among duplicate reads is specifed by the dup_keep flag.'''
    
    #As reads are now sorted by soft-clipping adjusted positions, PCR duplicates with given POS that
    #were previously sorted apart should now be grouped together, such that stepping through the file 
    #comparing reads next to each other will capture all PCR duplicates

    with open(file_in,'r') as SAM:

        removed_d = 0 # Counter for the number of reads removed as duplicates
        removed_u = 0 # Counter for the number of reads removed for invalid UMI tags
        line = SAM.readline()

        while line.startswith('@'): # Read lines until the header lines are all read in (write these to final output sam)
            with open(file_out,'a') as out_sam:
                out_sam.write(line)
            line = SAM.readline()

        while line: # While the EOF has not been read in:
    
            line = Read(line) # Create read object from current line
            next_line = SAM.readline() # Read in the next read line

            if next_line == '': # If the next line read in is the EOF, break
                break

            next_line = Read(next_line) # Create a read object from line (this is a separate step as to not raise an error passing '' to Read())

            if umis != [] and line.tag not in umis: # If specific UMIs have been specified, and the read does not have a valid ID, skip that read and do not write it to the output!
                line = next_line.raw
                removed_u += 1 # Increment invalid umi count

            else: # If UMIs have not been specified, or line.tag is in UMIs, proceed:

                if dup_keep == '': # If we want to keep the first of any duplicates:
                    while line.chars == next_line.chars: # While the reads match on duplicate characteristics, read in lines to next_line until they don't match
                        next_line = SAM.readline() 
                        if next_line == '': # Again, if the next line read in is the EOF, break
                            break
                        next_line = Read(next_line)
                        removed_d += 1 #Increment the duplicate count

                    with open(file_out,'a') as out_sam: # At the point the two reads are no longer duplicates, write the first read out
                        out_sam.write(line.raw)

                    line = next_line.raw # Set the last line read in as the next_read as the new comparitor line

                elif dup_keep == 'Q':
                    dups = [] # Create a list for the duplicates
                    while line.chars == next_line.chars: # While the reads match on duplicate characteristics, read in lines to next_line until they don't match
                        dups.append(next_line) # Add the duplicate read to the list
                        next_line = SAM.readline() 
                        if next_line == '': # Again, if the next line read in is the EOF, break
                            break
                        next_line = Read(next_line)
                        removed_d += 1 #Increment the duplicate count

                    if len(dups) > 0: # If there are any duplicates found:
                        dups.append(line) # Also need to append the first read, since it is a duplicate as well!
                        max_q = 0 # Start the max quality at 0
                        for read in dups:
                            if read.qual >= max_q: # For each read in the set of duplicates, see if it has a higher quality than the current max
                                max_q = read.qual # If so, set that read's quality as the nex max
                                keep = read # set that read to keep

                        line = keep # Whatever read ends up being the one to keep, set line to be that (as that will be written out)

                    with open(file_out,'a') as out_sam:
                        out_sam.write(line.raw)

                    line = next_line.raw # Set the last line read in as the next_read as the new comparitor line

                elif dup_keep == 'R':
                    dups = [] # Create a list for the duplicates
                    while line.chars == next_line.chars: # While the reads match on duplicate characteristics, read in lines to next_line until they don't match
                        dups.append(next_line) # Add the duplicate read to the list
                        next_line = SAM.readline() 
                        if next_line == '': # Again, if the next line read in is the EOF, break
                            break
                        next_line = Read(next_line)
                        removed_d += 1 #Increment the duplicate count

                    if len(dups) > 0: # If there are any duplicates found:
                        dups.append(line) # Also need to append the first read, since it is a duplicate as well!
                        line = random.choice(dups) # Pick a random choice among the duplicates to keep

                    with open(file_out,'a') as out_sam:
                        out_sam.write(line.raw)

                    line = next_line.raw # Set the last line read in as the next_read as the new comparitor line


        return [removed_d, removed_u]


def summary_getter(sum_flag, file_in, counts, reads):
    if sum_flag == True: # If specified that the user wants the summary output file:
        with open('Dedup_summary.out','w') as sum_out:
            sum_out.write("Summary of PCR duplicate removal for "+file_in.split('/')[-1]+":\n")
            sum_out.write(str("Number of Reads Processed: "+str(reads)+"\n"))
            sum_out.write(str("Number of Reads Removed as Dups: "+str(counts[0])+"\n"))
            sum_out.write(str("Number of Reads Removed for Incorrect UMIs: "+str(counts[1])+"\n"))


####################
## Main function: ##
####################
args = get_arguments() # Get passed in arguments

if args.p == True:
    raise ValueError('Cannot Currently Process Paired reads. Please do not specify flag at this time') # Can't handle paired end reads yet

if args.dup_keep not in ['Q','R','']:
    raise ValueError('Duplicate Keep argument is not valid, must use "Q", "R", or ""') # Throw error if improper keep flag given

umis = []
if args.u != '':
    umis = umi_list(args.u)

true_char(args.f)
reads_processed = sort_by_true(args.f)
those_removed = dup_remover('Output.sam',args.f+'_deduped',umis,args.dup_keep)

if args.s == True:
    summary_getter(args.s, args.f, those_removed, reads_processed)

## Clean up:
os.remove('Output.sam') # Remove Output temp file

for line in fileinput.input(args.f, inplace=True): # Open the input file and remove the extra column from each of the reads
    if line.startswith('@'):
        print(line.strip('\n'))

    else:
        line = '\t'.join(line.split('\t')[:-1])
        print(line)


for line in fileinput.input(args.f+'_deduped', inplace=True): # Open the output file and do the same
    if line.startswith('@'):
        print(line.strip('\n'))

    else:
        line = '\t'.join(line.split('\t')[:-1])
        print(line)


