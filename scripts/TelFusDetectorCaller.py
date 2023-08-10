import argparse
import pysam
import regex
import timeit
import subprocess
import re
import multiprocessing as mp
import distance
from collections import Counter
import os
import glob

'''
Create unmapped bam file
This bam file will be used for looking unmapped mate reads later
'''
def CreateUnmappedBam(bam,bam2,Cores,tmp,sample):
    if bam2 == None:
        bam2 = f"{tmp}/{sample}.unmapped.bam"
        command = f"samtools view -h -b -f 4 -@ {Cores} -o {bam2} {bam} && samtools index {bam2}"

        # Submit linux command
        try:
            subprocess.run(command, shell=True)
        except subprocess.CalledProcessError as error:
            print(error)

        return(bam2)
    else:
        return(bam2)

'''
Collect results from multiprocessing
'''
def collect_result(result):
    pass


'''
Function to know the chromosomes (or contigs) to check
'''
def extract_contigs(bam):
    # Open bam file
    samfile = pysam.AlignmentFile(bam, "rb")

    # Get bam statistics
    # Same as samtools idxstats
    A = samfile.get_index_statistics()

    # Get a list of contigs with at least one read (mapped or unmapped)
    contigs = [A[x].contig for x in range(len(A)) if A[x].total > 0]

    # Close bam file
    samfile.close()

    return(contigs)


'''
Function to extract fusions from stdin
'''
# Check if read contain telomere fusion events
def check_fusions(seq):

    # Telomere repeat sequence
    forward="TTAGGG"
    reverse="CCCTAA"
    
    # Twice the fwd/rev pattern
    forwardtwice="TTAGGGTTAGGG"
    reversetwice="CCCTAACCCTAA"

    # First screen of TTAGGG
    if (forward in seq and reverse in seq):
        # If at both telomere repeats are present
        # We check for two consecutive repeats (but allowing errors)
        # matches_forward = regex.findall("("+forwardtwice+"){s<=2}", seq)
        # matches_reverse = regex.findall("("+reversetwice+"){s<=2}", seq)
        matches_forward = regex.findall("("+forwardtwice+")", seq)
        matches_reverse = regex.findall("("+reversetwice+")", seq)
        if len(matches_forward) > 0 and len(matches_reverse) > 0:
            return('Yes')
        else:
            return(None)
    else:
        return(None)

# Function to extract the reads with fusions and their corresponding read mates
def ExtractReadsWithFusionsAndMates(chrom,bam,bam2,tmp_folder,sample):
    # Output files
    if (chrom == '*'):
        # Create temporary file paths
        outfile1 = f"{tmp_folder}/{sample}__None__fusions_with_mapped_reads.bam"
        outfile2 = f"{tmp_folder}/{sample}__None__fusions_without_mapped_reads.bam"
        mates1 = f"{tmp_folder}/{sample}__None__fusions_with_mapped_reads.mates.bam"
        mates2 = f"{tmp_folder}/{sample}__None__fusions_without_mapped_reads.mates.bam"
        out_coverage = f"{tmp_folder}/{sample}__None__coverage.txt"
    else:        
        # Create temporary file paths
        outfile1 = f"{tmp_folder}/{sample}__{chrom}__fusions_with_mapped_reads.bam"
        outfile2 = f"{tmp_folder}/{sample}__{chrom}__fusions_without_mapped_reads.bam"
        mates1 = f"{tmp_folder}/{sample}__{chrom}__fusions_with_mapped_reads.mates.bam"
        mates2 = f"{tmp_folder}/{sample}__{chrom}__fusions_without_mapped_reads.mates.bam"
        out_coverage = f"{tmp_folder}/{sample}__{chrom}__coverage.txt"

    #-------
    # Step 1
    # Extract reads with fusions
    #-------

    # Open bam file
    samfile = pysam.AlignmentFile(bam, "rb", require_index = True)

    bamfile_1 = pysam.AlignmentFile(outfile1, "wb",template=samfile)
    bamfile_2 = pysam.AlignmentFile(outfile2, "wb",template=samfile)

    mates_1 = pysam.AlignmentFile(mates1, "wb",template=samfile)

    # Check read length
    Read_length = 0
    # Check total number of reads analysed
    N_reads = 0
    # Potential fusion count
    N_fusions = 0
    # Iterate across reads in the bam
    for read in samfile.fetch(chrom, until_eof = True, multiple_iterators = True):
        
        # Remove supplementary and secondary alignment
        # Supplementary Alignment: A chimeric reads but not a representative reads
        # Primary Alignment and Secondary Alignment: A read may map ambiguously to multiple locations, e.g. due to repeats. Only one of the multiple read alignments is considered primary, and this decision may be arbitrary
        if read.is_secondary:
            continue
        
        # Update max read length
        Read_length_temp = read.infer_read_length()
        if (Read_length_temp != None):
            Read_length = max(Read_length,Read_length_temp)
        
        # Sum one more read
        N_reads += 1

        # Read sequence including soft clipped bases
        Seq = read.query_sequence
        
        # Fusion check
        Fusion = check_fusions(Seq)
        if (Fusion == "Yes"):
            N_fusions = N_fusions + 1
            try:
                mate = samfile.mate(read)
                mates_1.write(mate)
                bamfile_1.write(read)
            except:
                # When the mate is not mapped, the function samfile.mate is not possible
                # So we will extract these mates later
                bamfile_2.write(read)

    # Close bam file
    samfile.close()
    bamfile_1.close()
    bamfile_2.close()


    # If no reads in the contig, we remove the temp files and pass to the next chromosome
    if (N_reads == 0):
        os.remove(outfile1)
        os.remove(outfile2)
        os.remove(mates1)
        return ([Read_length,N_reads])

    # Save coverage information
    with open(out_coverage,'w') as out_cov:
        # Header
        Header_cov = ['Read_length','Total_reads']
        Header_cov = '\t'.join(Header_cov)+'\n'
        out_cov.write(Header_cov)
        # Values
        Info_cov = [str(Read_length),str(N_reads)]
        Info_cov = '\t'.join(Info_cov)+'\n'
        out_cov.write(Info_cov)

    # If not potential fusions, we remove the temp files and pass to the next chromosome
    if (N_fusions == 0):
        os.remove(outfile1)
        os.remove(outfile2)
        os.remove(mates1)
        return ([Read_length,N_reads])        

    # Index bam files
    pysam.index(outfile1)
    pysam.index(outfile2)

    #-------
    # Step 2
    # Extract mate reads for those that in the previous step was not possible (due to impossibility of pysam to look for unmapped mates)
    #-------

    # Get read ids for the fusion reads which mates were not mapped
    bamfile_2 = pysam.AlignmentFile(outfile2, "rb")
    READ_IDS = []
    for read in bamfile_2.fetch(until_eof = True):
        READ_IDS.append(read.qname)
    bamfile_2.close()
    READ_IDS = set(READ_IDS)

    # Extract unmapped mates for the reads with fusions
    # Open bam file
    samfile2 = pysam.AlignmentFile(bam2, "rb")

    # Mates 2 file
    mates_2 = pysam.AlignmentFile(mates2, "wb",template=samfile2)

    for read in samfile2.fetch(until_eof = True,multiple_iterators = True):
        # Don't want secondary reads
        if read.is_secondary:
            continue
        # Only looking for unmapped reads
        if read.is_unmapped:
            if (read.qname in READ_IDS):
                mates_2.write(read)

    # Close bam file
    samfile2.close()
    mates_2.close()

    # Return number of reads and read length
    return([Read_length,N_reads])

# Function to extract read count and length info from each conting and compute coverage
def ComputeCoverage(coverage_out,tmp,sample):
    # Find read count info in the temp folder
    cov_files = glob.glob(f"{tmp}/{sample}__*__coverage.txt")

    # Check if there are coverage files
    if len(cov_files) > 0:
        Multiprocessing_list = {'Read_length':[],'Total_reads':[]}
        # Gets info from each contig (chromosome)
        for cov_file in cov_files:
            # Get Read length, Total reads and Coverage
            with open(cov_file, 'r') as coverage_in:
                for line in coverage_in:
                    if not line.startswith('Read_length'):
                        line = line.rstrip('\n')
                        read_length,total_reads = line.split('\t')
                        Multiprocessing_list['Read_length'].append(int(read_length))
                        Multiprocessing_list['Total_reads'].append(int(total_reads))


        # Get max read length
        Read_length = max(Multiprocessing_list['Read_length'])
        # Get total reads
        Total_reads = sum(Multiprocessing_list['Total_reads'])
        # Compute coverage
        Coverage = round(int(Read_length)*int(Total_reads)/3200000000,2)

        # Saving results
        with open(coverage_out,'w') as out_cov:
            # Header
            Header_cov = ['Read_length','Total_reads','Coverage']
            Header_cov = '\t'.join(Header_cov)+'\n'
            out_cov.write(Header_cov)
            # Values
            Info_cov = [str(Read_length),str(Total_reads),str(Coverage)]
            Info_cov = '\t'.join(Info_cov)+'\n'
            out_cov.write(Info_cov)

        return('Yes')
    else:
        with open(coverage_out,'w') as out_cov:
            Info_cov = 'Warning: No useful reads found in the input bam file'+'\n'
            out_cov.write(Info_cov)

        return('No')

# Check if there are potential fusions files
def FusionsFlag(outfolder,tmp,sample):     
    # Fusion files
    fusion_files = glob.glob(f"{tmp}/{sample}__*__fusion*mapped_reads.bam")

    if len(fusion_files) > 0 :
        return('Yes')
    else:
        # Pass fusions out file
        pass_outf = f"{outfolder}/{sample}.summary_fusions.pass.tsv"
        # Filtered mutations out file
        filtered_outf = f"{outfolder}/{sample}.summary_fusions.filtered.tsv"

        with open(pass_outf,'w') as f:
            Info_f = 'Warning: No reads supporting potential telomere fusions found.'+'\n'
            f.write(Info_f)

        with open(filtered_outf,'w') as f:
            Info_f = 'Warning: No reads supporting potential telomere fusions found.'+'\n'
            f.write(Info_f)

        return('No')

'''
Collapse reads from all chromosomes
'''
def CollapseReads(bam,tmp_folder,sample):

    # Create temporary file paths
    outfile0 = f"{tmp_folder}/{sample}__header.sam"
    outfile1 = f"{tmp_folder}/{sample}__all_fusions.bam"
    outfile2 = f"{tmp_folder}/{sample}__all_mates.bam"
    outfile3 = f"{tmp_folder}/{sample}__all_fusion_mates.bam"

    #--------
    # 1. Merge fusions from all chromosomes in a single bam file
    #--------
    command0 = f"samtools view -H --no-PG {bam} > {outfile0}"

    # Submit linux command
    try:
        subprocess.run(command0, shell=True)
    except subprocess.CalledProcessError as error:
        print(error)

    #--------
    # 1. Merge fusions from all chromosomes in a single bam file
    #--------
    command1 = f"samtools merge -fpc {outfile1}.temp {tmp_folder}/{sample}__*__fusion*mapped_reads.bam && \
    samtools view {outfile1}.temp | sort | uniq | cat {outfile0} - | samtools view -bS | samtools sort -n > {outfile1} && rm -rf {outfile1}.temp"

    # Submit linux command
    try:
        subprocess.run(command1, shell=True)
    except subprocess.CalledProcessError as error:
        print(error)

    #--------
    # 2. Merge all mates from all chromosomes in a single bam file
    #--------
    command2 = f"samtools merge -fpc {outfile2}.temp {tmp_folder}/{sample}__*__fusion*mapped_reads.mates.bam && \
    samtools view {outfile2}.temp | sort | uniq | cat {outfile0} - | samtools view -bS | samtools sort -n > {outfile2} && rm -rf {outfile2}.temp"

    # Submit linux command
    try:
        subprocess.run(command2, shell=True)
    except subprocess.CalledProcessError as error:
        print(error)

    #--------
    # 3. Merge fusion and mate files and remove duplicated entries
    #--------
    command3 = f"samtools merge -fpc {outfile3}.temp {outfile1} {outfile2} && \
    samtools view {outfile3}.temp | sort | uniq | cat {outfile0} - | samtools view -bS | samtools sort -n > {outfile3} && rm -rf {outfile3}.temp {outfile0}"

    # Submit linux command
    try:
        subprocess.run(command3, shell=True)
    except subprocess.CalledProcessError as error:
        print(error)

    return(outfile1,outfile2,outfile3)

'''
Create a summary file with the reads
'''
# Build dictionary for pure sequences
# Create a dictionary with all possible pure breakpoint sequences
def PureBreakpointSequences():
    forward = 'TTAGGG'
    reverse = 'CCCTAA'

    Pure_breakpoint_sequences = {'Inward': [], 'Outward': []}

    # For inward fusions
    for i in range(len(forward)):
        donor_seq = forward[0:i]
        for j in range(len(reverse)):
            rev = reverse[::-1]
            receptor_seq = rev[0:j]
            receptor_seq = receptor_seq[::-1]

            # New breakpoint sequence
            BreakpointSequence = donor_seq+receptor_seq
            BreakpointSequence = regex.sub("^("+forward+")+",'',BreakpointSequence) # Remove telomere repeats (TTAGGG) at the beginning of the sequence
            BreakpointSequence = regex.sub("("+reverse+")+$",'',BreakpointSequence) # Remove telomere repeats (CCCTAA) at the end of the sequence

            if (BreakpointSequence not in Pure_breakpoint_sequences['Inward']):
                Pure_breakpoint_sequences['Inward'].append(BreakpointSequence)

    # For outward fusions
    for i in range(len(reverse)):
        donor_seq = reverse[0:i]
        for j in range(len(forward)):
            rev = forward[::-1]
            receptor_seq = rev[0:j]
            receptor_seq = receptor_seq[::-1]

            # New breakpoint sequence
            BreakpointSequence = donor_seq+receptor_seq
            BreakpointSequence = regex.sub("^("+reverse+")+",'',BreakpointSequence) # Remove telomere repeats (CCCTAA) at the beginning of the sequence
            BreakpointSequence = regex.sub("("+forward+")+$",'',BreakpointSequence) # Remove telomere repeats (TTAGGG) at the end of the sequence

            if (BreakpointSequence not in Pure_breakpoint_sequences['Outward']):
                Pure_breakpoint_sequences['Outward'].append(BreakpointSequence)

    return(Pure_breakpoint_sequences)

# Function to summarise each read
def overlaps(a, b):
    return min(a[1], b[1]) - max(a[0], b[0])

# Reverse complement
def ReverseComplement(SEQ):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    SEQ = regex.sub("__fw__",'TTAGGG',SEQ)
    SEQ = regex.sub("__rv__",'CCCTAA',SEQ)
    reverse_complement = "".join(complement.get(base.upper(), base.upper()) for base in reversed(SEQ))
    reverse_complement = regex.sub("TTAGGG",'__fw__',reverse_complement)
    reverse_complement = regex.sub("CCCTAA",'__rv__',reverse_complement)   
    return(reverse_complement)

# Function to simplify breakpoint sequence
# Getting sequences between expected telomere variant repeats based on the fusion orientation
def SimplifySequence(SEQ,Orientation):
    forward = 'TTAGGG'
    reverse = 'CCCTAA'

    if (Orientation == 'Inward'):
        TVR1 = forward
        tvr1 = "__fw__"
        TVR2 = reverse
        tvr2 = "__rv__"
    elif (Orientation == 'Outward'):
        TVR1 = reverse
        tvr1 = "__rv__"
        TVR2 = forward
        tvr2 = "__fw__"
    else:
        return()

    SEQ = TVR1 + SEQ + TVR2
    SEQ0 = regex.sub(forward,'__fw__',SEQ)
    SEQ0 = regex.sub(reverse,'__rv__',SEQ0)

    Breaks = regex.findall(tvr1 + "[A-Z]*" + tvr2 , SEQ0)
    Breaks = [regex.sub("^" + tvr1, '', x) for x in Breaks]
    Breaks = [regex.sub(tvr2 + "$", '', x) for x in Breaks]

    return(Breaks)

# Left-starting correction
def BreakpointSequenceLeftCorrection(SEQ1, max_dist = 1):
    forward = 'TTAGGG'
    reverse = 'CCCTAA'

    # New sequence
    SEQ_left = ''
    n = 0
    n_max = len(SEQ1)
    symbol_hold_n = 0
    symbol_hold = ''
    while n < len(SEQ1):
        # Extract 6-mer
        n_start = min(n,n_max) 
        n_end = min(n+6, n_max)
        kmer = SEQ1[n_start:n_end]

        if ('_' in kmer): # Ignore non-nucleotide bases edited during middle correction
            if kmer[-1] in '_fw_rv':
                SEQ_left = SEQ_left + kmer
                n = n + len(kmer)
            else:
                for BASE in kmer:
                    if BASE in '_fw_rv':
                        SEQ_left = SEQ_left + BASE
                        n = n + 1
        else:
            if (len(kmer) == 6):
                # Hamming distance for each 6-mer and the telomere repeats
                fwd_dist = distance.hamming(kmer, forward)
                rev_dist = distance.hamming(kmer, reverse)

                # Get new seq
                if (fwd_dist <= max_dist):
                    kmer = '__fw__'
                    n = n + len(kmer)
                    SEQ_left = SEQ_left + kmer
                elif (rev_dist <= max_dist):
                    kmer = '__rv__'
                    n = n + len(kmer)
                    SEQ_left = SEQ_left + kmer
                else:
                    n = n + 1
                    SEQ_left = SEQ_left + kmer[0]
            else:
                SEQ_left = SEQ_left + kmer
                n = n + len(kmer)
    return(SEQ_left)

# Right-starting correction
def BreakpointSequenceRightCorrection(SEQ1):
    SEQ1_rev = ReverseComplement(SEQ1)
    SEQ1_right = BreakpointSequenceLeftCorrection(SEQ1_rev)
    SEQ1_right = ReverseComplement(SEQ1_right)
    return(SEQ1_right)

# Clean sequence of isolated bases between telomere repeats
def CleanBreakpointSequence(SEQ):
    SEQ = regex.sub("^([A|C|T|G|N]?[(__fw__)|(__rv__)]+)+",'',SEQ)
    SEQ = regex.sub("([__fw__|__rv__]+[A|C|T|G|N]?)+$",'',SEQ)
    return(SEQ)

# Look for pure breakpoint sequences
def LookForPureBreakpointSequences(SEQ1,Orientation,Pure_breakpoint_sequences):
    forward = 'TTAGGG'
    reverse = 'CCCTAA'

    if (Orientation == 'Inward'):
        TVR1 = forward
        tvr1 = "__fw__"
        TVR2 = reverse
        tvr2 = "__rv__"
    elif (Orientation == 'Outward'):
        TVR1 = reverse
        tvr1 = "__rv__"
        TVR2 = forward
        tvr2 = "__fw__"

    Pures = []
    for pure in Pure_breakpoint_sequences[Orientation]:
        if pure in SEQ1:
            Pures.append(pure)

    # Return if empty
    if (len(Pures) < 1):
        return([SEQ1,'Alternative'])

    # Get the longest potential pure sequence
    Pure_long = max(Pures, key=len) 
    if (len(Pure_long) < 4):
        return([SEQ1,'Alternative'])

    # Get the location of the pure sequence in the original sequence
    Location_start = SEQ1.find(Pure_long)
    Location_end =  Location_start + len(Pure_long)

    # Extract left and right sequences around the pure candidate sequence
    Seq_left = SEQ1[0:Location_start]
    Seq_right = SEQ1[Location_end:len(SEQ1)]

    # Correct sequences
    Seq_left_corrected = BreakpointSequenceLeftCorrection(Seq_left,2)
    Seq_right_corrected = BreakpointSequenceLeftCorrection(Seq_right,2)

    # New sequence
    New_SEQ = tvr1 + Seq_left_corrected + Pure_long + Seq_right_corrected + tvr2

    # Simplify sequence
    New_SEQ2 = SimplifySequence(New_SEQ,Orientation)
    if (len(New_SEQ2) == 1 and New_SEQ2[0] == Pure_long):
        New_SEQ3 = [Pure_long,'Pure']
    else:
        New_SEQ3 = [SEQ1,'Alternative']

    # Return final sequence
    return(New_SEQ3)

# Breakpoint Sequence Correction (step 3)
def BreakpointSequenceCorrection(SEQ2):
    # Previous step sequence clean
    SEQ2_clean = CleanBreakpointSequence(SEQ2)

    # Left-starting corrected clean
    SEQ3_left = BreakpointSequenceLeftCorrection(SEQ2)
    SEQ3_left_clean = CleanBreakpointSequence(SEQ3_left)

    # Right-starting corrected clean
    SEQ3_right = BreakpointSequenceRightCorrection(SEQ2)
    SEQ3_right_clean = CleanBreakpointSequence(SEQ3_right)

    # Final Breakpoint sequence (the shortest of the previous)
    SEQ3 = min([SEQ2_clean,SEQ3_left_clean,SEQ3_right_clean], key=len)

    return(SEQ3)

def MainBreakpointSequenceCorrection(SEQ,Orientation):
    # Generate dicctionary for telomeres
    Pure_breakpoint_sequences = PureBreakpointSequences()

    # Expected telomere variant repeats
    forward = 'TTAGGG'
    reverse = 'CCCTAA'

    # Type of the sequence without correction
    if (SEQ in Pure_breakpoint_sequences[Orientation]):
        Type0 = 'Pure'
    elif (len(SEQ) <= 12):
        Type0 = 'Pure2'
    else:
        Type0 = 'Alternative'

    # Step 1. Simplify middle
    SEQ1 = SimplifySequence(SEQ,Orientation)

    if (len(SEQ1) > 1): # Different mini-breakpoints
        return ([SEQ1,Type0,'MultipleMiniBreakpoints'])
    else:
        SEQ1 = SEQ1[0]

    ## Step 2. Check the presence of pure sequences
    SEQ2, Type = LookForPureBreakpointSequences(SEQ1,Orientation,Pure_breakpoint_sequences)
    #print ('Done 2')

    if (Type == 'Pure'): # If pure, return the results
        return ([SEQ2,Type0,Type])

    # Step 3. Correcting breakpoint sequence allowing one mismatch
    # https://www.codespeedy.com/sequencematcher-in-python/
    # https://stackoverflow.com/questions/35517353/how-does-pythons-sequencematcher-work
    SEQ3= BreakpointSequenceCorrection(SEQ2)
    #print ('Done 3')

    # Step 4. Simplify sequence again
    SEQ4 = SimplifySequence(SEQ3,Orientation)[0]
    if (SEQ4 in Pure_breakpoint_sequences[Orientation]):
        Type2 = 'Pure'
    elif (len(SEQ4) <= 12):
        Type2 = 'Pure2'
    else:
        Type2 = 'Alternative'

    return([SEQ4,Type0,Type2])

def summary_read(read,Genome):

    # Read name
    Read_name = read.qname

    # Chrom
    Chrom = read.reference_name

    # Start
    Start = read.reference_start

    # End
    End = read.reference_end

    # MAPQ
    MQ = read.mapping_quality

    # Flag
    Flag = read.flag

    # Read 1/2 label
    if (read.is_read1):
        Read_pair = 'R1'
    else:
        Read_pair = 'R2'

    # Cigar string
    Cigar = read.cigarstring

    # Reverse/forward strand
    if (read.is_reverse):
        Strand = '-'
    else:
        Strand = '+'

    # Secondary
    if (read.is_secondary):
        Secondary = 'Secondary'
    else:
        Secondary = 'Primary'

    # Supplementary
    if (read.is_supplementary):
        Supplementary = 'Supplementary'
    else:
        Supplementary = 'Non-supplementary'

    # Duplicate
    if (read.is_duplicate):
        Duplicate = 'Duplicate'
    else:
        Duplicate = 'Non-duplicate'

    # SEQ
    SEQ = read.query_sequence

    # Read length
    Seq_length = len(SEQ)

    # Aligned length
    Aligned_length = read.query_alignment_length

    #--
    # Search for telomere repeats in the sequence
    #--
    # Telomere repeat sequence
    forward="TTAGGG"
    reverse="CCCTAA"

    # Twice the fwd/rev pattern
    forwardtwice="TTAGGGTTAGGG"
    reversetwice="CCCTAACCCTAA"

    # Forward 0 mismatches
    Forward0 = regex.findall("("+forward+")", SEQ)

    # Forward 1 mismatches
    Forward1 = regex.findall("("+forward+"){s<=1}", SEQ)

    # Forward 2 mismatches
    Forward2 = regex.findall("("+forward+"){s<=2}", SEQ)

    # Reverse 0 mismatches
    Reverse0 = regex.findall("("+reverse+")", SEQ)

    # Reverse 1 mismatches
    Reverse1 = regex.findall("("+reverse+"){s<=1}", SEQ)

    # Reverse 2 mismatches
    Reverse2 = regex.findall("("+reverse+"){s<=2}", SEQ)

    # Forward_twice 0 mismatches
    ForwardTwice0 = regex.findall("("+forwardtwice+")", SEQ)

    # Reverse_twice 0 mismatches
    ReverseTwice0 = regex.findall("("+reversetwice+")", SEQ)

    if (len(ForwardTwice0) > 0 and len(ReverseTwice0) > 0):
        # Get locations of telomere repeats
        ForwardTwice0_min_i = SEQ.find(forwardtwice) # Forward repeat more to the left
        ForwardTwice0_max_i = SEQ.rfind(forwardtwice) # Forward repeat more to the right
        ReverseTwice0_min_i = SEQ.find(reversetwice) # Reverse repeat more to the left
        ReverseTwice0_max_i = SEQ.rfind(reversetwice) # Reverse repeat more to the right

        # Decide orientation
        if (ForwardTwice0_min_i < min(ReverseTwice0_min_i,ReverseTwice0_max_i) and ForwardTwice0_max_i < min(ReverseTwice0_min_i,ReverseTwice0_max_i)):
            Orientation = 'Inward'
            BreakpointSequence0 = SEQ[ForwardTwice0_max_i:(ReverseTwice0_min_i+12)]
            BreakpointSequence = regex.sub("^("+forward+")+",'',BreakpointSequence0) # Remove telomere repeats (TTAGGG) at the beginning of the sequence
            BreakpointSequence = regex.sub("("+reverse+")+$",'',BreakpointSequence) # Remove telomere repeats (CCCTAA) at the end of the sequence

        elif (ReverseTwice0_min_i < min(ForwardTwice0_min_i,ForwardTwice0_max_i) and ReverseTwice0_max_i < min(ForwardTwice0_min_i,ForwardTwice0_max_i)):
            Orientation = 'Outward'
            BreakpointSequence0 = SEQ[ReverseTwice0_max_i:(ForwardTwice0_min_i+12)]
            BreakpointSequence = regex.sub("^("+reverse+")+",'',BreakpointSequence0) # Remove telomere repeats (CCCTAA) at the beginning of the sequence
            BreakpointSequence = regex.sub("("+forward+")+$",'',BreakpointSequence) # Remove telomere repeats (TTAGGG) at the end of the sequence
        else:
            Orientation = 'Other'
            BreakpointSequence0 = 'Other'
            BreakpointSequence = 'Other'
    else:
        Orientation = 'None'
        BreakpointSequence0 = 'None'
        BreakpointSequence = 'None'

    # Looking for endogenous fusions
    Endogenous_fusions = {'Hg19':{'chr2':[114355250,114365750],'chr4':[0,0],'chr9':[130911250,130922000]},
                    'Hg38':{'chr2':[113597750,113608250],'chr4':[190177250,190188000],'chr9':[128149000,128159750]}}

    if (Chrom != None):
        Chrom2 = 'chr'+regex.sub("^chr",'',Chrom) # Ensuring that the prefix is included in the chromosome name
    else:
        Chrom2 = 'none'

    #if (Chrom2 in Endogenous_fusions[Genome].keys() and int(MQ) >= 8 and Orientation == 'Inward'):
    if (Chrom2 in Endogenous_fusions[Genome].keys() and int(MQ) >= 8):
        Overalap = overlaps([Start,End],Endogenous_fusions[Genome][Chrom2]) # Check the bases of overlap
        if (Overalap >= 0):
            Chrom = Chrom2 + '_endogenous'

    # Check the mapping region
    Telomere_ends = {"Hg19":{"chr1": [249140621, float('inf')], "chr10": [135424747, float('inf')], "chr11": [134896516, float('inf')], "chr12": [133741895, float('inf')],
            "chr13": [115059878, float('inf')], "chr14": [107239540, float('inf')], "chr15": [102421392, float('inf')], "chr16": [90244753, float('inf')],
            "chr17": [81095210, float('inf')], "chr18": [77967248, float('inf')], "chr19": [59018983, float('inf')], "chr2": [243089373, float('inf')], 
            "chr20": [62915520, float('inf')],"chr21": [48019895, float('inf')], "chr22": [51194566, float('inf')], "chr3": [197912430, float('inf')], 
            "chr4": [191044276, float('inf')],"chr5": [180805260, float('inf')], "chr6": [171005067, float('inf')], "chr7": [159028663, float('inf')],
            "chr8": [146254022, float('inf')],"chr9": [141103431, float('inf')], "chrX": [155160560, float('inf')], "chrY": [59263566, float('inf')]},
        "Hg38":{"chr1": [248846422, float('inf')], "chr10": [133687422, float('inf')], "chr11": [134976622, float('inf')], "chr12": [133165309, float('inf')],
            "chr13": [114254328, float('inf')], "chr14": [106933718, float('inf')], "chr15": [101881189, float('inf')], "chr16": [90228345, float('inf')],
            "chr17": [83147441, float('inf')], "chr18": [80263285, float('inf')], "chr19": [58507616, float('inf')], "chr2": [242083529, float('inf')], 
            "chr20": [64334167, float('inf')], "chr21": [46599983, float('inf')], "chr22": [50708468, float('inf')], "chr3": [198185559, float('inf')],
            "chr4": [190104555, float('inf')], "chr5": [181428259, float('inf')], "chr6": [170695979, float('inf')], "chr7": [159235973, float('inf')],
            "chr8": [145028636, float('inf')], "chr9": [138284717, float('inf')], "chrX": [155930895, float('inf')], "chrY": [57117415, float('inf')]}} 

    if (int(MQ) >= 8 and Seq_length >= 50):
        if (Chrom2 in Telomere_ends[Genome].keys()):
            Region_end_chromosome = overlaps([Start,End],Telomere_ends[Genome][Chrom2])
            Region_start_chromosome = overlaps([Start,End],[0,100000]) # We consider as telomere the first 100,000 bps of each chromosome

            if (Region_end_chromosome >= 0 or Region_start_chromosome >= 0):
                Chromosomal_region = 'Subtelomeric'
            else:
                Chromosomal_region = 'Intra-chromosomal'

        else:
            Chromosomal_region = 'Other_contig'
    else:
        Chromosomal_region = 'Unknown'

    # Correct and clean breakpoint sequence
    if Orientation in ['Inward','Outward']:
        BreakpointSequenceCorrected,Type,TypeCorrected = MainBreakpointSequenceCorrection(BreakpointSequence,Orientation)
        if Chrom in ['chr9_endogenous','chr4_endogenous','chr9_endogenous']:
            TypeCorrected == Chrom
    else:
        BreakpointSequenceCorrected,Type,TypeCorrected = ['None','None','None']

    # We edit the None chromosome by *
    if (Chrom == None):
        Chrom = '*'

    # Info that we want to collect
    INFO = [Read_name,Chrom,Start,End,Chromosomal_region,Seq_length,Aligned_length,MQ,Flag,Read_pair,Cigar,Strand,Duplicate,Supplementary,Secondary,SEQ,len(Forward0),len(Forward1),len(Forward2),len(Reverse0),len(Reverse1),len(Reverse2),Orientation,BreakpointSequence,BreakpointSequenceCorrected,Type,TypeCorrected]
    return([Orientation,Read_pair,INFO])


'''
Collapse pair read info
To create a summary row for each
'''
# Function to sort and prioritize repeated reads (p.e supplementary alignments)
def SortReadsPriorization(READS): 
    Reads_sorted = [] # List to put the sorted reads
    Reads_priorization = {} # Dictionary used for read priorization
    Fusion_category = 0 # Variable to know if there are fusions in this set of reads
    for x in  READS:
        Priorization = 0
        Read_name,Chrom,Start,End,Chromosomal_region,Seq_length,Aligned_length,MQ,Flag,Read_pair,Cigar,Strand,Duplicate,Supplementary,Secondary,SEQ,Forward0,Forward1,Forward2,Reverse0,Reverse1,Reverse2,Orientation,BreakpointSequence,BreakpointSequenceCorrected,Type,TypeCorrected = x

        # Give values for priorization
        # Top importance for reads having a fusion
        # Then, based on if they are supplementary or not alignment
        if (Orientation != 'None'):
            Priorization = Priorization + 2
            Fusion_category = Fusion_category + 1
        if (Supplementary != 'Supplementary'):
            Priorization = Priorization + 1

        if (Priorization not in Reads_priorization.keys()):
            Reads_priorization[Priorization] = [x]
        else:
            Reads_priorization[Priorization].append(x)

    # Proper sorting of reads
    for i in sorted(Reads_priorization.keys(), reverse = True):
        Reads_sorted.extend(Reads_priorization[i])

    return ([Fusion_category,Reads_sorted])

# Function to summarise each pair in a single row
def SummarisePair(reads,Genome):
    # Collect read info
    # To check if the same read is represented more than once (p.e. to supplementary or secondary alignments)
    Read_ids = {'R1':[], 'R2' : []} 
    for read_i in range(len(reads)):
        # Extract read information
        Orientation, Read_pair, Read_info = summary_read(reads[read_i], Genome)
        
        # Cound the number of times are particular read in the pair is considered (R1 or R2)
        Read_ids[Read_pair].append(Read_info)

    # Append reads to the group of reads with and without fusions in the pair
    # Select only one read if the same R1 or R2 is represented more than once
    All_reads_info = { 'Fusion' : [], 'Non_fusion' : []}
    Main_read = None
    Mate_read = None
    # Sort and save results for R1
    R1f, R1s = SortReadsPriorization(Read_ids['R1'])
    if (R1f > 0):
        Main_read = R1s[0]
    elif len(R1s) > 0:
        Mate_read = R1s[0]
    else:
        Mate_read = []

    # Sort and save results for R2
    R2f, R2s = SortReadsPriorization(Read_ids['R2'])
    if (R2f > 0):
        if (Main_read == None): # Check if R1 has telomere fusions
            Main_read = R2s[0]
        else: # If R1 had fusion, save R2 as mate
            Mate_read = R2s[0]
    elif len(R2s) > 0:
        if (Main_read == None): # Check if R1 has fusion and if none (R1/R2) has a fusion, return a flag
            return ('No_fusion')
        else:
            Mate_read = R2s[0]
    else: 
        if (Main_read == None): # Check if R1 has fusion and if none (R1/R2) has a fusion, return a flag
            return ('No_fusion')
        else: # If R1 had fusion, save R2 as mate
            Mate_read = [] # Empty list

    # Collapse and return all info
    if (len(Main_read) == 27):
        Chrom = Main_read[1]
        if len(Mate_read) == 27: # Check if the mate is properly created
            Chrom_mate = Mate_read[1]
        else:
            Mate_read = ['NA'] * 27
            Chrom_mate = Mate_read[1]

        # Create a variable to see if the resulting pair is a telomere fusion or maps in a endogenous region
        Chrom_list = [Chrom,Chrom_mate]
        Fusion = [x for x in Chrom_list if 'endogenous' in x]
        if len(Fusion) > 0:
            Fusion = Fusion[0]
        else:
            Fusion = 'Telomere_fusion'

        # Concensus orientation
        # Orientation considering the sequences in both paired reads
        Orientation1 = Main_read[22]
        Orientation2 = Mate_read[22]

        Orientation_full_pair = list(set([x for x in [Orientation1,Orientation2] if x not in ['NA','None']]))

        if (len(Orientation_full_pair) == 0):
            Orientation_full_pair = 'None'
        elif(len(Orientation_full_pair) == 1):
            Orientation_full_pair = ''.join(Orientation_full_pair)
        else:
            if 'Other' in Orientation_full_pair:
                Orientation_full_pair = 'Other'
            else:
                Orientation_full_pair = 'In-out'

        # Save results
        Fusion_info = Main_read + Mate_read + [str(Orientation_full_pair)] + [str(Fusion)]
        Fusion_info2 = [str(x) for x in Fusion_info]
        Fusion_info2 = '\t'.join(Fusion_info2)+'\n'
        return(Fusion_info2)
    else:
        return('Error_in_paired_reads') 

def CreateSummary(outfile3,summary_out,Genome):
    # Open filtered bam file with reads
    samfile = pysam.AlignmentFile(outfile3, "rb")

    # Open output file
    out = open(summary_out,'w')

    # Header of the output file
    Header_main=["Read_name","Chrom","Start","End","Chromosomal_region","Seq_length","Aligned_length","MQ","Flag","Read_pair","Cigar","Strand","Duplicate","Supplementary","Secondary","SEQ","Forward0","Forward1","Forward2","Reverse0","Reverse1","Reverse2","Orientation","BreakpointSequence","BreakpointSequenceCorrected","Type","TypeCorrected"]
    Header_mate=["Read_name_mate","Chrom_mate","Start_mate","End_mate","Chromosomal_region_mate","Seq_length_mate","Aligned_length_mate","MQ_mate","Flag_mate","Read_pair_mate","Cigar_mate","Strand_mate","Duplicate_mate","Supplementary_mate","Secondary_mate","SEQ_mate","Forward0_mate","Forward1_mate","Forward2_mate","Reverse0_mate","Reverse1_mate","Reverse2_mate","Orientation_mate","BreakpointSequence_mate","BreakpointSequenceCorrected_mate","Type_mate","TypeCorrected_mate"]
    Header = Header_main + Header_mate + ['Orientation_full_pair','Fusion']
    out.write('\t'.join(Header)+'\n')

    # Read all files
    pair = ""
    reads = []
    for read in samfile.fetch(until_eof = True):

        if (read.qname == pair): # Group the reads in pairs
            reads.append(read)
        else: 
            if (pair == ""): # If the first read
                pair = read.qname
                reads.append(read)
            else: 
                # When pair is completed, get the summary
                # Create the summary for the current pair
                Summary_line =  SummarisePair(reads,Genome)
                if Summary_line not in ['No_fusion','Error_in_paired_reads']:   
                    out.write(Summary_line)
                else:
                    pass

                # Start saving the next pair
                reads = [read]
                pair = read.qname

    # Run the last pair of reads in the file
    Summary_line =  SummarisePair(reads,Genome)
    if Summary_line not in ['No_fusion','Error_in_paired_reads']:   
        out.write(Summary_line)
    else:
        pass

    # Close bam file
    samfile.close()
    out.close()


'''
Filtering fusion calls
'''
# Looking for duplicated sequences
def LookDuplicatedSequences(summary_out):
    Sequences = []
    with open(summary_out,'r') as f:
        for line in f:
            if line.startswith('Read_name'):
                pass
            else:
                line = line.rstrip('\n')

            # Split line in different colummns
            Read_name,Chrom,Start,End,Chromosomal_region,Seq_length,Aligned_length,MQ,Flag,Read_pair,Cigar,Strand,Duplicate,Supplementary,Secondary,SEQ,Forward0,Forward1,Forward2,Reverse0,Reverse1,Reverse2,Orientation,BreakpointSequence,BreakpointSequenceCorrected,Type,TypeCorrected,Read_name_mate,Chrom_mate,Start_mate,End_mate,Chromosomal_region_mate,Seq_length_mate,Aligned_length_mate,MQ_mate,Flag_mate,Read_pair_mate,Cigar_mate,Strand_mate,Duplicate_mate,Supplementary_mate,Secondary_mate,SEQ_mate,Forward0_mate,Forward1_mate,Forward2_mate,Reverse0_mate,Reverse1_mate,Reverse2_mate,Orientation_mate,BreakpointSequence_mate,BreakpointSequenceCorrected_mate,Type_mate,TypeCorrected_mate,Orientation_full_pair,Fusion_type = line.split('\t')

            # Screening sequences
            SEQ2 = SEQ + '>' + SEQ_mate
            Sequences.append(SEQ2)

    Sequences_c = Counter(Sequences)
    Sequences_c2 = {k: v for k, v in Sequences_c.items() if v >= 2}

    return(Sequences_c2)

# Filtering (reasons) read
def GettingFilteringVariable(line,Duplicated_sequences,Duplicated_sequence_already_included):
    # Split line in different colummns
    Read_name,Chrom,Start,End,Chromosomal_region,Seq_length,Aligned_length,MQ,Flag,Read_pair,Cigar,Strand,Duplicate,Supplementary,Secondary,SEQ,Forward0,Forward1,Forward2,Reverse0,Reverse1,Reverse2,Orientation,BreakpointSequence,BreakpointSequenceCorrected,Type,TypeCorrected,Read_name_mate,Chrom_mate,Start_mate,End_mate,Chromosomal_region_mate,Seq_length_mate,Aligned_length_mate,MQ_mate,Flag_mate,Read_pair_mate,Cigar_mate,Strand_mate,Duplicate_mate,Supplementary_mate,Secondary_mate,SEQ_mate,Forward0_mate,Forward1_mate,Forward2_mate,Reverse0_mate,Reverse1_mate,Reverse2_mate,Orientation_mate,BreakpointSequence_mate,BreakpointSequenceCorrected_mate,Type_mate,TypeCorrected_mate,Orientation_full_pair,Fusion_type = line.split('\t')

    #---
    # Evaluation of the read
    # Getting filter variable
    #---
    Filter = []
    No_mate_flag = 0
    # Filter 0
    if Read_name_mate == 'NA' or SEQ_mate == 'NA':
        Filter.append('MateNotFound')
        No_mate_flag = 1

    # Filter 1: Duplicates
    SEQ2 = SEQ + '>' + SEQ_mate
    if SEQ2 in Duplicated_sequences:
        if SEQ2 not in Duplicated_sequence_already_included:
            Duplicated_sequence_already_included.append(SEQ2)
        else:
            Filter.append('Duplicate')

    # Filter 2: Number of breakpoints
    if Orientation == 'Other': # Main read
        Filter.append('MainMultipleBreakpoints')
    if Orientation_mate == 'Other': # Mate read
        Filter.append('MateMultipleBreakpoints')

    # Filter 3: Chromosomal region
    if Chromosomal_region == 'Intra-chromosomal':
        Filter.append('MainIntraChromosomal')
    if Chromosomal_region_mate == 'Intra-chromosomal':
        Filter.append('MateIntraChromosomal')

    # Filter 4: No TRVs in mate
    if (Forward0_mate != 'NA' and int(Forward0_mate) == 0 and int(Reverse0_mate) == 0):
        Filter.append('NoTelomereRepeatMate')

    # Filter 6: Endogenous fusions
    if Chrom.endswith('_endogenous'):
        a = 'Main_' + Chrom
        Filter.append(a)
    if Chrom_mate.endswith('_endogenous'):
        b = 'Mate_' + Chrom
        Filter.append(b)
    
    # Filter 7: Wrong endo 9 orientation
    if (Fusion_type == 'chr9_endogenous' and Orientation_full_pair != 'Inward'):
        Filter.append('chr9_endogenous_WrongOrientation')

    # FOR THE FUTURE ? --> Filter 7: Check if fusion is found in both reads

    # Final decision filter
    Filter = ','.join(Filter)
    if Filter == '':
        Filter = 'Pass'

    return([Filter,Duplicated_sequence_already_included])

# Function to filter the final summaries
def FilteringFinalSummaries(pass_outf,filtered_outf,summary_out,coverage_out):
    # Get Read length, Total reads and Coverage
    with open(coverage_out, 'r') as coverage_out:
        for line in coverage_out:
            if not line.startswith('Read_length'):
                line = line.rstrip('\n')
                Read_length,Total_reads,Coverage = line.split('\t')

    # Get duplicated sequences
    Duplicated_sequences = LookDuplicatedSequences(summary_out)

    # Create files with new columns
    pass_out = open(pass_outf,'w')
    filtered_out = open(filtered_outf,'w')

    Duplicated_sequence_already_included = []
    with open(summary_out,'r') as f:
        for line in f:
            if line.startswith('Read_name'):
                line = line.rstrip('\n')

                Header = line.split('\t')
                Header.extend(['Sample_read_length','Sample_total_reads','Sample_coverage','Filter'])
                Header = '\t'.join(Header)+'\n'

                # Print header
                pass_out.write(Header)
                filtered_out.write(Header)
            else:
                line = line.rstrip('\n')

                # Getting filtering criteria
                Filter, Duplicated_sequence_already_included = GettingFilteringVariable(line,Duplicated_sequences,Duplicated_sequence_already_included)

                # Getting new line
                new_line = line.split('\t')
                new_line.extend([str(Read_length),str(Total_reads),str(Coverage),str(Filter)])
                new_line = '\t'.join(new_line)+'\n'

                # Split reads in pass and filtered
                Not_wanted_for_endo = ['MateNotFound','Duplicate','MultipleBreakpoints','chr9_endogenous_WrongOrientation']
                if (Filter == 'Pass' or ('chr9_endogenous' in Filter and not any([x in Filter for x in Not_wanted_for_endo]))):
                    pass_out.write(new_line)
                else:
                    filtered_out.write(new_line)


    # Close out files
    pass_out.close()
    filtered_out.close()

'''
Arguments required for the computation of TelFusDetector
'''
def initialize_parser():
    parser = argparse.ArgumentParser(description='TelFusDetector: This tool serves to detect telomere fusions in DNA sequencing data sets (paired-end sequencing data)')
    parser.add_argument('--bam', type=str, help='BAM file to be analysed (Sorted by coordinate and indexed)', required = True)
    parser.add_argument('--genome', type=str, choices =  ["Hg38","Hg19"], default = 'Hg38', help='Reference genome used for the alignment (Choose Hg38 or Hg19)', required = False)
    parser.add_argument('--bam2', type=str, default=None, help='BAM file with unmapped reads (Sorted by coordinate and indexed). If provided, it will speed up the computation', required = False)
    parser.add_argument('--outfolder', default = '.', help='Out directory', required = False)
    parser.add_argument('--sample', type=str, default = 'Sample', help='Sample ID. All output files will start with this string', required = False)
    parser.add_argument('--tmpfolder', default = None, help='Directory to save all temporary files. If it exits, please empty it before running TelFusDetector', required = False)
    parser.add_argument('--threads',default = 1, help='Number of threads to use [Default: 1]',required=False,type = int)
    return (parser)


'''
TelFusDetector main function
'''

def main():

    #------------
    # Get arguments
    #------------

    parser = initialize_parser()
    args = parser.parse_args()

    bam = args.bam
    bam2 = args.bam2
    sample = args.sample
    outfolder = args.outfolder
    tmp = args.tmpfolder
    Cores = args.threads
    Genome = args.genome

    # Create main output folder if it does not exist
    if (outfolder != '.'):
        try:
            # Create target Directory
            os.mkdir(outfolder)
        except FileExistsError:
            pass
    else:
        pass

    # Create temp folder if it does not exist
    if (tmp != None):
        try:
            # Create target Directory
            os.mkdir(tmp)
        except FileExistsError:
            pass
    else:
        tmp = f"{outfolder}/tmp_dir"
        try:
            # Create target Directory
            os.mkdir(tmp)
        except FileExistsError:
            pass

    # Printing arguments
    print ("\n----------------------")
    print ("Running TelFusDetectorCaller")
    print ("Detecting telomere fusions")
    print ("----------------------\n")
    print ("Input arguments:")
    print (f"   Input bam file: {bam}")
    print (f"   Reference genome: {Genome}")
    print (f"   Input unmapped bam file: {bam2}")
    print (f"   Out folder: {outfolder}")
    print (f"   Sample ID: {sample}")
    print (f"   Temporary folder: {tmp}")
    print (f"   Number of threads: {Cores}\n")

    #------------
    # 1. Get the list of chromosomes to be analysed
    #------------

    print ("1. Getting the list of contigs (chromosomes) to be analysed\n")

    # Get a list of contigs with at least one read (mapped or unmapped)
    Contigs = extract_contigs(bam)

    # Quality measure
    if (len(Contigs) < 1):
        print ('    ERROR: No contigs (chromosomes) with mapped or unmapped reads.')
        exit()

    # Add reads with NO chromosome attributed 
    Contigs = ['*'] + Contigs
    
    #------------
    # 2. Create the bam file with unmapped reads
    # If this is provided, do not re-generate it
    #------------

    print ("2. Checking (or creating) unmapped bam file")
    bam2 = CreateUnmappedBam(bam,bam2,Cores,tmp,sample)
    print (f"   - Unmapped bam file used: {bam2}\n")

    # ------------
    # 3. Extract reads with fusions
    # ------------

    print ("3. Extracting reads with fusions (and the mates)")
    print (f"   - Using {Cores} thread(s)\n")

    if (Cores > 1):
        pool = mp.Pool(Cores)

        # Step 3.1: Use loop to parallelize
        for chrom in Contigs:
            # Extract fusion reads plus mates for each chromosome (in parallel)
            pool.apply_async(ExtractReadsWithFusionsAndMates, args=(chrom,bam,bam2,tmp,sample), callback=collect_result)
                   
        # Step 3.2: Close Pool and let all the processes complete    
        pool.close()
        pool.join()
    else:
        for chrom in Contigs:
            # Extract fusion reads plus mates for each chromosome
            collect_result(ExtractReadsWithFusionsAndMates(chrom,bam,bam2,tmp,sample))

    # Get coverage file
    coverage_out = f"{outfolder}/{sample}.all_chromosomes.coverage.tsv"
    Cov_flag = ComputeCoverage(coverage_out,tmp,sample)
 
    # Fusion flag
    Fusion_flag = FusionsFlag(outfolder,tmp,sample)

    # Check if we stop computation or we continue
    if Cov_flag == 'No' and Fusion_flag == 'No':
        print (f"TelFusDetector finished:") 
        print (f"   - No useful reads found in the input bam file\n") 
        print (f"   - No reads supporting potential telomere fusions found\n")
        return()
    elif Cov_flag == 'No':
        print (f"TelFusDetector finished:") 
        print (f"   - No useful reads found in the input bam file\n")
        return()
    elif Fusion_flag == 'No':
        print (f"TelFusDetector finished:") 
        print (f"   - No reads supporting potential telomere fusions found\n")
        return()

    #------------
    # 4. Collapse all reads with fusions and all mates
    #------------
    
    print ("4. Collapsing all reads with fusions and all mates in a single file\n")
       
    outfile1,outfile2,outfile3 = CollapseReads(bam,tmp,sample)

    #-----------
    # 5. Create initial summaries
    #-----------

    print ("5. Creating initial fusion summary file\n")

    # Raw fusion summary output
    summary_out = f"{tmp}/{sample}.summary_fusions.temp.tsv"
    CreateSummary(outfile3,summary_out,Genome)

    #-----------
    # Filter fusion calls
    #-----------

    print ("6. Creating final summaries")

    # Pass fusions out file
    pass_outf = f"{outfolder}/{sample}.summary_fusions.pass.tsv"
    # Filtered mutations out file
    filtered_outf = f"{outfolder}/{sample}.summary_fusions.filtered.tsv"

    FilteringFinalSummaries(pass_outf,filtered_outf,summary_out,coverage_out)

    print (f"   - Passed fusions: {pass_outf}")
    print (f"   - Filtered fusion: {filtered_outf}\n")

#---------------
# Running TelFusDetector
#---------------

if __name__ == '__main__':
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    Seconds = round(stop - start)
    print(f"TelFusDetectorCaller computation time: {Seconds} seconds\n") 

