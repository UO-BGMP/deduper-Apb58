#!/Users/Adrian/miniconda3/bin/python3
# Adrian Bubie
#12/14/17

## PCR_Deduper: 
## Argparse arguments - sam file (abs path), paired-end flag (boolean), UMI file (abs path), and help
import textwrap
import random
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
    parser.add_argument("-umi", help="Optional file to define UMI sequences (one UMI per line). UMIs are assumed to be in the last field of the QNAME for each read. Reads with invalid UMIs are filtered from final output. If file not specified, SAM is assumed to use randomers. Must include absolute path to file <str>", required=False, type=str, default='')
    parser.add_argument("-dup_keep", help="Optional argument to designate which read to keep of a duplicate set. 'Q' will keep read with highest quality score, 'R' will choose a random read to keep. No specification will keep the first read encounted. <str>", required=False, type=str, default='')
    parser.add_argument("-s", help="Set to 'True' for additional summary file output <boolean> (def='False')", required=False, type=bool, default=False)
    return parser.parse_args()

args = get_arguments()


class Read:
    """ Read object: object called to get read parts from SAM file lines, or change attributes """
    def __init__(self, line):
        self.chr = line.split('\t')[2]
        self.flag = int(line.split('\t')[1])
        self.tag = line.split('\t')[0].split(':')[-1]
        self.seq = line.split('\t')[9]
        self.cig = line.split('\t')[5]

        if 'S' in self.cig and self.cig.find('S')+1 < len(self.cig): # If the read is softclipped on the 5'-end (given by CIGAR) we need to adjust the position
            adj = int(self.cig[:self.cig.find('S')])
            self.pos = int(line.split('\t')[3])-adj # Set the position with the adjusted value
        else:
            self.pos = int(line.split('\t')[3]) # No soft-clipping on 5'-end, keep the value of the POS column.

        self.qual = int(line.split('\t')[4])
        self.raw = line


def umi_list(umi_file):
    '''Returns list of UMI tags from UMI file argument'''

    with open(umi_file, 'r') as ulist:
        return [line.strip('\n') for line in ulist.readlines()]


def flag_stat(flag, paired):
    '''Bitwise flag checker: checks read flag and returns the strand of the read. 
    If the paired-end argument passed as True, which read in the pair the read is 
    checked and returned.'''

    strand = '+'
    if (flag & 16) == 16: # Check strand:
        strand = '-'  # If SEQ is reverse complemented by flag definition, than "-" is returned. Otherwise, "+" is returned.

    return strand

    # If the file is specified as paired-ended:
    #if paired == True:
    #    if (flag & 40) == 40 and (flag & 80) != 80:
    #        pair = 1
    #    elif (flag & 40) != 40 and (flag & 80) == 80:
    #        pair = 2

    #    return pair


def dup_remover(reads, dups, dup_keep):
    '''Takes list of reads and list of duplicates, and removes the correct duplicates from the reads
    based on the dup_keep flag. Returns a subset of the total reads list.'''

    to_keep = [] # Initialize list of reads to keep; note that this will only ever be one item long, but is easy to use to remove the rest of the duplicates

    if dup_keep == 'Q': # If Quality flag has been specified, keep read only with the largest MAPQ 
        max_q = 0
        for read in dups:
            if Read(read).qual >= max_q:
                max_q = Read(read).qual
                keep = read
        to_keep.append(keep)

    elif dup_keep == 'R': # If Random flag has been specified, keep a random duplicate
        to_keep.append((random.choice(dups)))
        
    elif dup_keep == '': # Otherwise, keep the first duplicate
        to_keep.append(dups[-1]) # This is equivalent to the first read used to compare, since it is added to the duplicates last
        

    dups = [y for y in dups if y not in to_keep]
    sub_reads = [x for x in reads if x not in dups]

    return sub_reads


def dedup_chunk(chunk, dup_keep, paired):
    '''Takes in chunk of reads (list), and uses the first read in the chunk as the comparitor.
    Uses this read to determine any PCR duplicates of the read in the chunk. Each read's POS is adjusted
    if needed based on the CIGAR field (for soft clipping), and strand orientation is determined from the FLAG before
    reads are compared. Duplicates are then removed from list of reads, and deduped chunk is returned.'''

    set_of_reads = chunk
    removed = 0 # Initialize the number of reads removed as dups
    r_std = Read(chunk[0]) # Create comparitor read from first read in chunk
    dups = []

    for i in range(1,len(chunk)):
        pos_dup = Read(chunk[i])
        if flag_stat(r_std.flag, paired) == flag_stat(pos_dup.flag, paired): #First check if the reads have the same strand
            if pos_dup.pos == r_std.pos and pos_dup.tag == r_std.tag: #If tag (UMI/randomer) and positions match for reads:
                dups.append(pos_dup.raw)
    
    if len(dups) > 0:
        dups.append(r_std.raw)
        set_of_reads = dup_remover(chunk, dups, dup_keep)
        removed = int(len(chunk)-len(set_of_reads)) # Count of reads removed from chunk, for summary file

    return [set_of_reads, removed]


def read_chunker(line, umis, dup_keep, paired, summary, out_file):
    '''Takes line from SAM file and list of UMIs, if available, and checks that  
    line has valid UMI. If empty UMI list is provided, then reads assumed to use
    randomers. Once validated, lines are read from file and stored in list 
    until read positions exceed maximum read length is found.'''

    read1 = Read(line) # create Read object from SAM line
    lines_proc = 0 # Counter to track number reads read in (only used in summary file)
    removed_cd = 0 # Counter to track number of duplicates removed (only used in summary file)
    removed_cu = 0 # Counter to track number of reads removed for incorrect UMIs (only used in summary file)
    lines_added = 0 # Total number of non-header lines added to the file (only used in summary file)

    if umis != []: # If UMIs are used:
        # Check the umi tag on read1
        utag = read1.tag
        while utag not in umis: # if utag not in umis, read file lines until a read is found with a valid umi.
            read1 = Read(SAM.readline())
            utag = read1.tag
            lines_proc = lines_proc + 1
            removed_cu = removed_cu + 1 # Increment the number of lines removed for incorrect UMI tag by 1
        else:
            pass

    # Get length of the read, then add it to read1 pos to find cut-off for position chunk
    r_len = len(read1.seq)

    # Create list to store chunk of reads
    chunk = []
    chunk.append(read1.raw) #Add in first read
    read2 = ''

    while read1: # Until the end of the SAM file

        read2 = SAM.readline() # Read in the next line of the SAM
        lines_proc = lines_proc + 1 # Increment the number of reads read in

        if read2 == '':
            break
        
        read2 = Read(read2) # create Read object from read2

        while read2.pos <= read1.pos+r_len and read2.chr == read1.chr: # Check to see if the position of the next read is within the first read's position + read length, then add it to the chunk
            # If umis passed in, only add reads if they have a valid umi:
            if umis != []:
                if read2.tag in umis:
                    chunk.append(read2.raw)
                    read2 = SAM.readline()
                    lines_proc = lines_proc + 1
                    if read2 == '':
                        break

                    read2 = Read(read2) # read in the next line

                else:
                    read2 = SAM.readline()
                    lines_proc = lines_proc + 1
                    removed_cu = removed_cu + 1 # Increment the number of lines removed for incorrect UMI tag by 1
                    if read2 == '':
                        break

                    read2 = Read(read2) # if tag is wrong, just read in next line, and don't add read

            else: # If umis not passed in, don't do additional check
                chunk.append(read2.raw)
                read2 = SAM.readline()
                lines_proc = lines_proc + 1
                if read2 == '':
                    break

                read2 = Read(read2) # read in the next line
            

        while len(chunk) > 1 and read1.pos+r_len < read2.pos: # Once the chunk is determined, start with the first read and compare it to all the reads in the chunk as potential dups.
            #print([c.split('\t')[3] for c in chunk])
            process = dedup_chunk(chunk, dup_keep, paired)
            chunk = process[0]
            #print([c.split('\t')[3] for c in chunk])
            removed_cd = removed_cd + process[1]

            if len(chunk) > 1: # Once the comparison is done, and there is more than one read left in the chunk after deduping, remove the first read from the chunk and start the comparison again with the next read in the chunk. Continue until only one read remains in the chunk:
                with open(out_file,'a') as out:
                    if umis != [] and Read(chunk[0]).tag in umis:
                        out.write(chunk[0])
                        lines_added = lines_added + 1
                    elif umis == []:
                        out.write(chunk[0])

                read1 = Read(chunk[1])
                chunk = chunk[1:]

        if len(chunk) > 1:
            chunk.append(read2.raw)

        if len(chunk) == 1: # reset the chunk when there is only one read left, by removing the last read (nothing left to compare for duplicates) and starting a new chunk with the last read read in but NOT included (last read2):
            #print(chunk[0].split('\t')[3])
            with open(out_file,'a') as out:
                if umis != [] and Read(chunk[0]).tag in umis:
                    out.write(chunk[0])
                    lines_added = lines_added + 1
                elif umis == []:
                    out.write(chunk[0])

            chunk = []
            chunk.append(read2.raw)
            read1 = Read(chunk[0])


    while len(chunk) > 1: # This additional loop is required to catch the final, non-empty reads of the file if loop is broken because of a break from EOF:
        #print([c.split('\t')[3] for c in chunk])
        process = dedup_chunk(chunk, dup_keep, paired)
        chunk = process[0]
        #print([c.split('\t')[3] for c in chunk])
        removed_cd = removed_cd + process[1]

        if len(chunk) > 1:
            with open(out_file,'a') as out:
                if umis != [] and Read(chunk[0]).tag in umis:
                    out.write(chunk[0])
                    lines_added = lines_added + 1
                elif umis == []:
                    out.write(chunk[0])
                
            read1 = Read(chunk[1])
            chunk = chunk[1:]

    if len(chunk) == 1: 
        #print(chunk[0].split('\t')[3])
        with open(out_file,'a') as out:
            if umis != [] and Read(chunk[0]).tag in umis:
                out.write(chunk[0])
                lines_added = lines_added + 1
            elif umis == []:
                out.write(chunk[0])


    if summary == True: # If specified that the user wants the summary output file:
        with open('Dedup_summary.out','w') as sum_out:
            sum_out.write("Summary of PCR duplicate removal for "+out_file.split('/')[-1]+":\n")
            sum_out.write(str("Number of Reads Processed: "+str(lines_proc)+"\n"))
            sum_out.write(str("Number of Reads Removed as Dups: "+str(removed_cd)+"\n"))
            sum_out.write(str("Number of Reads Removed for Incorrect UMIs: "+str(removed_cu)+"\n"))
            sum_out.write(str("Number of Non-duplicate, remaining reads: "+str(lines_added)+"\n"))


### Main Function: ###

# Set an error for paired-end argument at this time:
if args.p == True:
    raise ValueError('Cannot Currently Process Paired reads. Please do not specify flag at this time')

if args.dup_keep not in ['Q','R','']:
    raise ValueError('Duplicate Keep argument is not valid, must use "Q", "R", or ""')

# If UMI file given, create list of tags; otherwise, leave it as an empty list:
umis = []
if args.umi != '':
    umis = umi_list(args.umi)

with open(args.f, 'r') as SAM:

    line = SAM.readline() # Read in first line of the SAM file

    # Read through all header lines and store for output file writing:
    header = []
    while line.startswith('@') == True:
        header.append(line)
        line = SAM.readline()

    if header != []: # If header is read in, write it to the deduped output file:
        with open((args.f.strip('.sam'))+'_deduped.sam','a') as out:
            out.write(''.join(header))

    # Read through Reads, remove duplicates, and write results out:
    output_file = args.f.strip('.sam')+'_deduped.sam'
    read_chunker(line, umis, args.dup_keep, args.p, args.s, output_file)



