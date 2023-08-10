import pysam
import numpy as np
import random
import regex
import timeit
import argparse


'''
Useful functions for computation
'''
def Reverse(Seq):
	Dict_rev = {'A':'T','C':'G','T':'A','G':'C','N':'N'}
	Seq2 = [Dict_rev[x] for x in Seq]

	Seq3 = ''.join(Seq2)
	return(Seq3)

def ComplementaryReverse(Seq):
	Seq2 = Seq[::-1]
	Seq3 = Reverse(Seq2)
	return(Seq3)

def overlaps(a, b):
	return min(a[1], b[1]) - max(a[0], b[0])

def ExtractTelomereSequence(Donor_seq0,Min_no_nas,window):
	# Get the index of the starting sequence with at least 10 continous non NA bases
	L = len(Donor_seq0)
	No_NA = 0
	idx = 0
	while No_NA < Min_no_nas and idx < L:
		base = Donor_seq0[idx]
		if base == 'N':
			idx += 1
			No_NA = 0
		else:
			idx += 1
			No_NA += 1

	if idx == L:
		return(['No_clean_sequence',None])

	# Get the new sequence
	New_start = idx-Min_no_nas
	New_end = New_start + window

	Donor_seq = Donor_seq0[New_start:New_end]

	return([idx,Donor_seq])

# Function to generate sequences in the telomeres
def GenerateTelomereSequence(chrom,region,window,Min_no_nas,Random,inFasta,contigs):
	if (Random == 'Yes'):
		Telomere_random_start_end = random.randint(0,2500) # To force the script to not generate random telomere regions only at the end (or beginning) of the telomere
	else:
		Telomere_random_start_end = 0

	if region == 'start':
		# Initial indexes
		Start = 0 + Telomere_random_start_end
		End = 100000 

		# Extract reference sequence
		Donor_seq0 = inFasta.fetch(chrom, Start, End).upper()

		# Extract new sequence to be fused
		Donor_seq_new1 = ExtractTelomereSequence(Donor_seq0,Min_no_nas,window)

		# Get real coordinates in the chromosomes fasta file
		Real_coordenate_start = Start + Donor_seq_new1[0] - Min_no_nas 
		Real_coordenate_end = Start + Donor_seq_new1[0] + window - Min_no_nas

	else:
		# Get chromosome end coordinates
		End_coordinate = contigs[chrom] - Telomere_random_start_end
		Start = End_coordinate - 100000
		End = End_coordinate
		
		# Extract reference sequence
		Donor_seq0 = inFasta.fetch(chrom, Start, End).upper()

		# Reverse the sequence to move across the sequence
		Donor_seq0 = ComplementaryReverse(Donor_seq0)

		# Extract sequence
		Donor_seq_new1 = ExtractTelomereSequence(Donor_seq0,Min_no_nas,window)

		# Reverse back the sequence
		Donor_seq_new1[1] = ComplementaryReverse(Donor_seq_new1[1])

		# Get real coordinates in the chromosomes fasta file
		Real_coordenate_end = End_coordinate - Donor_seq_new1[0] + Min_no_nas
		Real_coordenate_start = Real_coordenate_end - window

	# Return real coordinates and sequences
	return(region,chrom,Real_coordenate_start,Real_coordenate_end,Donor_seq_new1[1])

# Extract breakpoint sequence
def GetBreakpointSequence(SEQ,orientation):
	forward = 'TTAGGG'
	reverse = 'CCCTAA'

	Forward_min_i = SEQ.find(forward) # Forward repeat more to the left
	Forward_max_i = SEQ.rfind(forward) # Forward repeat more to the right
	Reverse_min_i = SEQ.find(reverse) # Reverse repeat more to the left
	Reverse_max_i = SEQ.rfind(reverse) # Reverse repeat more to the right

	# Breakpoint sequence impossible to be found
	if (Forward_min_i == -1 or Forward_max_i == -1 or Reverse_min_i == -1 or Reverse_max_i == -1):
		return ('BreakpointSequence_not_found')

	# Check if there are crossing tsv repeats
	if (overlaps([Forward_min_i,Forward_max_i],[Reverse_min_i,Reverse_max_i]) > 0):
		return ('Overlapping')

	# Check for single repeats
	# Decide orientation
	if orientation == 'Inward':
		BreakpointSequence0 = SEQ[Forward_max_i:(Reverse_min_i+6)]
		BreakpointSequence = regex.sub("^("+forward+")+",'',BreakpointSequence0) # Remove telomere repeats (TTAGGG) at the beginning of the sequence
		BreakpointSequence = regex.sub("("+reverse+")+$",'',BreakpointSequence) # Remove telomere repeats (CCCTAA) at the end of the sequence
	else: # orientation == 'Outward'
		BreakpointSequence0 = SEQ[Reverse_max_i:(Forward_min_i+6)]
		BreakpointSequence = regex.sub("^("+reverse+")+",'',BreakpointSequence0) # Remove telomere repeats (CCCTAA) at the beginning of the sequence
		BreakpointSequence = regex.sub("("+forward+")+$",'',BreakpointSequence) # Remove telomere repeats (TTAGGG) at the end of the sequence

	return (BreakpointSequence)

# Function to generate telomere fusions
def GenerateTelomereFusion(Donor,Receptor,orientation,window,Min_no_nas,Count,Random,inFasta,contigs):
	# Get donor sequence
	chrom,region = Donor.split('_')
	Donor_info = GenerateTelomereSequence(chrom,region,window,Min_no_nas,Random,inFasta,contigs)
	region,chrom, Real_coordenate_start, Real_coordenate_end, Sequence = Donor_info

	# Get receptor sequence
	chrom2,region2 = Receptor.split('_')
	Receptor_info = GenerateTelomereSequence(chrom2,region2,window,Min_no_nas,Random,inFasta,contigs)
	region2,chrom2, Real_coordenate_start2, Real_coordenate_end2, Sequence2 = Receptor_info

	# Create telomere fusion and fasta entry
	# based on orientation
	if orientation == 'Inward':
		if (region == 'start' and region2 == 'start'):
			Seq1 = ComplementaryReverse(Sequence)
			Seq2 = Sequence2
			strand1 = '-'
			strand2 = '+'
		elif (region == 'start' and region2 == 'end'):
			Seq1 = ComplementaryReverse(Sequence)
			Seq2 = ComplementaryReverse(Sequence2)
			strand1 = '-'
			strand2 = '-'
		elif (region == 'end' and region2 == 'start'):
			Seq1 = Sequence
			Seq2 = Sequence2
			strand1 = '+'
			strand2 = '+'
		elif (region == 'end' and region2 == 'end'):
			Seq1 = Sequence
			Seq2 = ComplementaryReverse(Sequence2)
			strand1 = '+'
			strand2 = '-'
	else: # orientation == 'Outward'
		if (region == 'start' and region2 == 'start'):
			Seq1 = Sequence
			Seq2 = ComplementaryReverse(Sequence2)
			strand1 = '+'
			strand2 = '-'
		elif (region == 'start' and region2 == 'end'):
			Seq1 = Sequence
			Seq2 = Sequence2
			strand1 = '+'
			strand2 = '+'
		elif (region == 'end' and region2 == 'start'):
			Seq1 = ComplementaryReverse(Sequence)
			Seq2 = ComplementaryReverse(Sequence2)
			strand1 = '-'
			strand2 = '-'
		elif (region == 'end' and region2 == 'end'):
			Seq1 = ComplementaryReverse(Sequence)
			Seq2 = Sequence2
			strand1 = '-'
			strand2 = '+'

	# Telomere fusion sequence
	Telomere_fusion_seq = Seq1 + Seq2

	# Junction sequences
	if (len(Seq1) < 30):
		Seq1_junction = Seq1
	else:
		Seq1_junction = Seq1[-30:]
	if (len(Seq2) < 30):
		Seq2_junction = Seq2
	else:
		Seq2_junction = Seq2[0:30]

	Raw_junction = f"...{Seq1_junction}>{Seq2_junction}..."

	# Breakpoint_sequence
	Breakpoint_sequence = GetBreakpointSequence(Telomere_fusion_seq,orientation)

	# Text label (entry in the output fasta file)
	Label = f">TelomereFusion_{orientation}_{Count}\t{orientation};{chrom}:{Real_coordenate_start}-{Real_coordenate_end},{strand1},{region};{chrom2}:{Real_coordenate_start2}-{Real_coordenate_end2},{strand2},{region2};{Breakpoint_sequence};{Raw_junction}"

	return([Label,Telomere_fusion_seq,Breakpoint_sequence])


'''
Arguments required for the random generation of Telomere Fusions
'''
def initialize_parser():
	parser = argparse.ArgumentParser(description='Tool for generating random telomere fusions using as template a reference genome.')
	parser.add_argument('--genome', type=str, help='Reference genome used as template for the generation of telomere fusions. Be sure that the indexed file (.fai) exists', required = True)
	parser.add_argument('--outp', default = 'Simulated_telomere_fusions', help='Out prefix file', required = False)
	parser.add_argument('--n_tfs', default = 10, help='Number of telomere fusions to simulate (for both Inward and Outward orientations) [Default = 10]', required = False,type = int)
	parser.add_argument('--size', default = 500, help='Number of bases taken from each part of the telomeres to be fused [Default = 500]',required=False,type = int)
	parser.add_argument('--seed', default = 1111, help='Random seed used to run the tool [Default = 1111]',required=False,type = int)
	parser.add_argument('--unique', choices = ['Yes','No'], default = 'Yes', help='Yes: Each BreakpointSequence is presented only once; No: The same BreakpointSequence can be presented more than once. [Default = Yes]', required=False)	
	parser.add_argument('--random', choices = ['Yes','No'], default = 'Yes', help='Yes: Random region inside the telomere; No: Take initial bases of the telomere. [Default = Yes]', required=False)

	return (parser)


'''
Generating telomere fusions
using all functions described above
'''

def main():

	#------------
	# Get arguments
	#------------

	parser = initialize_parser()
	args = parser.parse_args()

	# To be passed as parameters
	FASTA = args.genome
	N_fusions = args.n_tfs
	window = args.size
	seed = args.seed
	Random = args.random
	out_prefix = args.outp
	Min_no_nas = 20
	Unique = args.unique

	# Printing arguments
	print ("\n----------------------")
	print ("Simulating telomere fusions")
	print ("----------------------\n")
	print ("Input arguments:")
	print (f"   Reference genome: {FASTA}")
	print (f"   Number of telomere fusions to be simulated: {N_fusions}")
	print (f"   Size of each telomere fragment fused: {window}")
	print (f"   Out prefix: {out_prefix}")
	print (f"   Random seed: {seed}")
	print (f"   Unique: {Unique}")
	print (f"   Allow random selection of regions inside the telomeres: {Random}\n")

	# Random seed
	random.seed(seed)

	# List of expected chromosomes in human genomce
	# We ignore chrY
	Expected_chroms = ['chr' + str(i) for i in range(1,23)]
	Expected_chroms.extend(['chrX'])

	# Load raw fasta file
	inFasta = pysam.FastaFile(FASTA)

	# List of contigs
	# Filter to the ones found in the Expected_chroms list
	contigs_names = inFasta.references
	contigs_lengths = inFasta.lengths

	contigs = {contigs_names[i] : contigs_lengths[i] for i in range(len(contigs_names)) if contigs_names[i] in Expected_chroms}

	# Generate random combinatios of start and end chromosomes
	contigs0 = [x+'_start' for x in contigs.keys()] # Starting telomere
	contigs1 = [x+'_end' for x in contigs.keys()] # End telomere
	contigs_all = contigs0 + contigs1

	#-----
	# 1. Generate X number of chormosome end pairs
	#-----
	outfile1 = f"{out_prefix}.fasta"
	Outfile1 = open(outfile1, 'w')

	Telomere_fusions_ids = []
	for orientation in ['Inward','Outward']:
		BreakpointSequences_used = []
		Count = 1
		Real_count = 0
		while Count <= N_fusions and Real_count < (N_fusions + 100):

			Real_count = Real_count + 1

			# Select randomly a telomere region to be used as donor
			Donor_i = random.randrange(len(contigs_all))
			Donor = contigs_all[Donor_i]

			# Select randomly a telomere region to be used as receptor
			Receptor_i = random.randrange(len(contigs_all))
			Receptor = contigs_all[Receptor_i]

			# Generate telomere fusion
			try:
				Tel_id,Telomere_fusion_sequence,Breakpoint_sequence = GenerateTelomereFusion(Donor,Receptor,orientation,window,Min_no_nas,Count,Random,inFasta,contigs)
				if Unique == 'Yes':
					if (Breakpoint_sequence not in BreakpointSequences_used and Breakpoint_sequence != 'Overlapping'):
						BreakpointSequences_used.append(Breakpoint_sequence)
			
						# Save the results in the fasta file
						Outfile1.write(Tel_id + '\n')
						Outfile1.write(Telomere_fusion_sequence + '\n')

						# Save telomere fusion info
						Telomere_fusions_ids.append(Tel_id)
						Count = Count + 1
				else:
					if (Breakpoint_sequence != 'Overlapping'):
						BreakpointSequences_used.append(Breakpoint_sequence)
			
						# Save the results in the fasta file
						Outfile1.write(Tel_id + '\n')
						Outfile1.write(Telomere_fusion_sequence + '\n')

						# Save telomere fusion info
						Telomere_fusions_ids.append(Tel_id)
						Count = Count + 1
			except: 
				continue

	inFasta.close()
	Outfile1.close()

	#-----
	# 2. Get report with telomere fusions
	#-----
	outfile2 = f"{out_prefix}.summary.txt"
	Outfile2 = open(outfile2, 'w')

	# Header of the summary file
	Header = ["TF_id","Orientation","Breakpoint_sequence"," Junction_sequence","Donor_coordinates","Donor_RevComp","Donor_chromosome_part","Receptor_coordinates","Receptor_RevComp","Receptor_chromosome_part"]
	Header = '\t'.join(Header)
	Outfile2.write(Header + '\n')   

	for tf in Telomere_fusions_ids:
		# Extract information from telomere fusion id
		TF_id, info = tf.split('\t')

		# Clean TF id
		TF_id = regex.sub("^>",'',TF_id)

		# Get other information
		Orientation, donor_info, receptor_info, Breakpoint_sequence, Junction_sequence = info.split(';')
		Donor_coordinates,Donor_RevComp,Donor_chromosome_part = donor_info.split(',')
		Receptor_coordinates,Receptor_RevComp,Receptor_chromosome_part = receptor_info.split(',')

		# Line to be printed in the summary file
		Line = [TF_id,Orientation,Breakpoint_sequence, Junction_sequence,Donor_coordinates,Donor_RevComp,Donor_chromosome_part,Receptor_coordinates,Receptor_RevComp,Receptor_chromosome_part]
		Line = '\t'.join(Line)
		Outfile2.write(Line + '\n')

	Outfile2.close()

	# Printing the output files
	print (f"\nNew generated files:")
	print (f"   - Fasta file with telomere fusions: {outfile1}")
	print (f"   - Summary file with telomere fusions information: {outfile2}")


'''
Running telomere fusion functions
'''

if __name__ == '__main__':
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    Seconds = round(stop - start)
    print(f"Telomere fusion simulator computation time: {Seconds} seconds\n") 

