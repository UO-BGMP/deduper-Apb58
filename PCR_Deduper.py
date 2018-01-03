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
    parser.add_argument("-dup_keep", help="Optional argument to designate which read to keep of a duplicate set. 'Q' will keep read with highest quality score, 'F' will keep the first read encountered. No specification will choose a random read. <str>", required=False, type=str, default='')
    #parser.add_argument("-s", help="Set to 'True' for additional summary file output <boolean> (def='False')")
    #parser.add_argument("-h","--help", help="Show this message and exit", required=False, type=str)
    return parser.parse_args()

args = get_arguments()


class Read:
    """ Read object: object called to get read parts from SAM file lines, or change attributes """
    def __init__(self, line):
        self.pos = int(line.split('\t')[3])
        self.chr = line.split('\t')[2]
        self.flag = int(line.split('\t')[1])
        self.tag = line.split('\t')[0].split(':')[-1]
        self.seq = line.split('\t')[9]
        self.cig = line.split('\t')[5]
        self.qual = line.split('\t')[4]
        self.raw = line

    def adj_pos(self, adj):
        self.pos = self.pos-adj


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
    if paired == True:
        if (flag & 40) == 40 and (flag & 80) != 80:
            pair = 1
        elif (flag & 40) != 40 and (flag & 80) == 80:
            pair = 2

        return pair


def dup_remover(reads, dups, dup_keep):
    '''Takes list of reads and list of duplicates, and removes the correct duplicates from the reads
    based on the dup_keep flag. Returns a subset of the total reads list.'''

    to_keep = []

    if dup_keep == 'Q': # If Quality flag has been specified, keep read only with the MAPQ 
        max_q = 0
        for read in dups:
            if Read(read).qual > max_q:
                max_q = Read(read).qual
                to_keep.append(read)

    elif dup_keep == 'F': # If First flag has been specified, keep the first duplicate
        to_keep.append(dups[-1]) # This is equivalent to the first read used to compare, since it is added to the duplicates last

    elif dup_keep == '': # Otherwise, keep a random duplicate
        to_keep.append((random.choice(dups)))

    dups = [y for y in dups if y not in to_keep]
    sub_reads = [x for x in reads if x not in dups]

    return sub_reads


def dedup_chunk(chunk, dup_keep, paired):
    '''Takes in chunk of reads (list), and uses the first read in the chunk as the comparitor.
    Uses this read to determine any PCR duplicates of the read in the chunk. Each read's POS is adjusted
    if needed based on the CIGAR field (for soft clipping), and strand orientation is determined from the FLAG before
    reads are compared. Duplicates are then removed from list of reads, and deduped chunk is returned.'''

    set_of_reads = chunk
    r_std = Read(chunk[0]) # Create comparitor read from first read in chunk
    dups = []

    for i in range(1,len(chunk)):
        pos_dup = Read(chunk[i])
        if flag_stat(r_std.flag, paired) == flag_stat(pos_dup.flag, paired): #First check if the reads have the same strand
            if 'S' in pos_dup.cig and pos_dup.cig.find('S')+1 < len(pos_dup.cig): #determing if there is soft-clipping at the 5'-end such that the pos has to be adjusted
                adj = int(pos_dup.cig[:pos_dup.cig.find('S')]) # get value to adjust by
                pos_dup.adj_pos(adj)

            if pos_dup.pos == r_std.pos and pos_dup.tag == r_std.tag: #If tag (UMI/randomer) and positions match for reads:
                dups.append(pos_dup.raw)
    
    if len(dups) > 0:
        dups.append(r_std.raw)
        set_of_reads = dup_remover(chunk, dups, dup_keep)

    return set_of_reads


def read_chunker(line, umis, dup_keep, paired, out_file):
    '''Takes line from SAM file and list of UMIs, if available, and checks that  
    line has valid UMI. If empty UMI list is provided, then reads assumed to use
    randomers. Once validated, lines are read from file and stored in list 
    until read positions exceed maximum read length is found.'''

    read1 = Read(line) # create Read object from SAM line
    if umis != []: # If UMIs are used:
        # Check the umi tag on read1
        utag = read1.tag
        while utag not in umis: # if utag not in umis, read file lines until a read is found with a valid umi.
            read1 = Read(SAM.readline())
            utag = read1.tag
        else:
            pass

    # Get length of the read, then add it to read1 pos to find cut-off for position chunk
    r_len = len(read1.seq)

    # Create list to store chunk of reads
    chunk = []
    chunk.append(read1.raw) #Add in first read

    while read1: # Until the end of the file

        read2 = SAM.readline() # Read in the next line of the SAM
        if read2 == '':
            break
        
        read2 = Read(read2) # create Read object from read2

        while read2.pos <= read1.pos+r_len and read2.chr == read1.chr: # Check to see if the position of the next read is within the first read's position + read length, then add it to the chunk
            # If umis passed in, only add reads if they have a valid umi:
            if umis != []:
                if read2.tag in umis:
                    chunk.append(read2.raw)
                    read2 = SAM.readline()
                    if read2 == '':
                        break

                    read2 = Read(read2) # read in the next line

                else:
                    read2 = SAM.readline()
                    if read2 == '':
                        break

                    read2 = Read(read2) # just read in next line, and don't add read

            else: # If umis not passed in, don't do any additional check
                chunk.append(read2.raw)
                read2 = SAM.readline()
                if read2 == '':
                    break

                read2 = Read(read2) # read in the next line
            

        while len(chunk) > 1 and read1.pos+r_len < read2.pos:
            chunk = dedup_chunk(chunk, dup_keep, paired)
            with open(out_file,'a') as out:
                if umis != [] and Read(chunk[0]).tag in umis:
                    out.write(chunk[0])
                elif umis == []:
                    out.write(chunk[0])

            read1 = Read(chunk[1])
            chunk = chunk[1:]

        if len(chunk) == 1: # reset the chunk when there is only one read in the list left
            with open(out_file,'a') as out:
                if umis != [] and Read(chunk[0]).tag in umis:
                    out.write(chunk[0])
                elif umis == []:
                    out.write(chunk[0])

            chunk = []
            chunk.append(read2.raw)
            read1 = Read(chunk[0])

    with open(out_file,'a') as out: 
        if read1 != '':
            if umis != [] and read1.tag in umis:
                out.write(read1.raw)   # This additional write is required to catch the final, non-empty read of the file if it is not otherwise captured in a chunk.
            elif umis == []:
                out.write(read1.raw)



### Main Function: ###

# If UMI file given, create list of tags:
if args.p == True:
    raise ValueError('Cannot Currently Process Paired reads. Please do not specify flag at this time')

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

    if header != []:
        with open((args.f.strip('.sam'))+'_deduped.sam','a') as out:
            out.write(''.join(header))

    # Read through Reads, remove duplicates, and write results out:
    output_file = args.f.strip('.sam')+'_deduped.sam'
    read_chunker(line, umis, args.dup_keep, args.p, output_file)



