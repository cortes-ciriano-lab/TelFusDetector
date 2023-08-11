import argparse
import pandas as pd
import timeit

'''
Check if file is correct
'''
def CheckInputFile(fusion_file):
    with open(fusion_file, 'r') as f:
        for line in f:
            if line.startswith('Warning'):
                return('No_fusions')
            elif line.startswith('Read_name'):
                return('Pass')
            else:
                return('Weird_format')

'''
Arguments required for the computation of TelFusDetector
'''
def initialize_parser():
    parser = argparse.ArgumentParser(description='TelFusDetectorRate: Tool that takes the output generated by TelFusDetector to calculate telomere fusion rates')
    parser.add_argument('--fusion_file', type=str, help='Telomere fusion summary file generated by TelFusDetector (*summary_fusions.pass.tsv)', required = True)
    parser.add_argument('--variables', nargs='+', help='Input column names to be used for grouping the telomere fusions (follow this format --variables VAR1 VAR2 VAR3 ...)', required = False)
    parser.add_argument('--purity', type=float, default = 1, help='Tumour purity [Default = 1]',required=False)
    parser.add_argument('--outfile', default = 'Telomere_fusion_rates.tsv', help='Output file with the calculated telomere fusion rates', required = False)
    return (parser)


'''
TelFusDetectorRates main function
'''

def main():

    #------------
    # Get arguments
    #------------

    parser = initialize_parser()
    args = parser.parse_args()

    fusion_file = args.fusion_file
    Variables = args.variables
    outfile = args.outfile
    Purity = args.purity

    if (Variables != None):
        Variables_str = ','.join(Variables)
        Variables2 = [x for x in Variables if x not in ['Fusion','Sample_coverage','Purity']]
        Variables2 = ['Fusion','Sample_coverage'] + Variables2
    else:
        Variables_str = 'None'
        Variables2 = ['Fusion','Sample_coverage']

    # Printing arguments
    print ("\n----------------------")
    print ("Running TelFusDetectorRates")
    print ("Telomere fusion rate estimation")
    print ("----------------------\n")
    print ("Input arguments:")
    print (f"   Input fusion file: {fusion_file}")
    print (f"   Variables used for grouping: {Variables_str}")
    print (f"   Purity: {Purity}")
    print (f"   Output file with telomere fusion rates: {outfile}\n")

    # 0. Check if the file is empty
    File_flag = CheckInputFile(fusion_file)

    # 1. Open fusion file
    fusions = pd.read_csv(fusion_file, sep='\t')

    Fix_variables = Variables

    N_rows = len(fusions.index)
    if N_rows > 0 and File_flag == 'Pass':
        # 2. Check if the grouping variable provided by the user are valid
        # Expected variables (or columns)
        Header_main=["Read_name","Chrom","Start","End","Chromosomal_region","Seq_length","Aligned_length","MQ","Flag","Read_pair","Cigar","Strand","Duplicate","Supplementary","Secondary","SEQ","Forward0","Forward1","Forward2","Reverse0","Reverse1","Reverse2","Orientation","BreakpointSequence","BreakpointSequenceCorrected","Type","TypeCorrected"]
        Header_mate=["Read_name_mate","Chrom_mate","Start_mate","End_mate","Chromosomal_region_mate","Seq_length_mate","Aligned_length_mate","MQ_mate","Flag_mate","Read_pair_mate","Cigar_mate","Strand_mate","Duplicate_mate","Supplementary_mate","Secondary_mate","SEQ_mate","Forward0_mate","Forward1_mate","Forward2_mate","Reverse0_mate","Reverse1_mate","Reverse2_mate","Orientation_mate","BreakpointSequence_mate","BreakpointSequenceCorrected_mate","Type","TypeCorrected_mate"]
        Expected_variables = Header_main + Header_mate + ['Fusion','Sample_read_length','Sample_total_reads','Sample_coverage', 'Filter'] 

        Variables_not_valid = [x for x in Variables2 if x not in fusions.columns]
        if (len(Variables_not_valid) > 0):
            Not_valid = ', '.join(Variables_not_valid)
            print ("- These set of variables are not valid:\n")
            print (f"{Variables_not_valid}\n")
            print ("- Use variables found in this list:\n")
            print (f"{Expected_variables}\n")
            exit()

        # 3. Add purity column to the table
        def AddPurity(row,Purity):
            if 'endogenous' in row['Fusion']:
                val = 1
            else:
                val = Purity
            return val

        fusions['Purity'] = fusions.apply(AddPurity, args = [Purity], axis=1)
        Variables2 = Variables2 + ['Purity']

        # 4. Count telomere fusions and compute rates
        Counting = fusions.groupby(Variables2).size().reset_index(name='TelomereFusion_counts')
        # Compute telomere fusion rate
        Counting['TelomereFusion_rate'] = Counting['TelomereFusion_counts'] / Counting['Sample_coverage']  / Counting['Purity']
        Counting.TelomereFusion_rate  = Counting.TelomereFusion_rate.round(4)
    elif File_flag == 'No_fusions':
        # Expected variables (or columns)
        Header_main=["Read_name","Chrom","Start","End","Chromosomal_region","Seq_length","Aligned_length","MQ","Flag","Read_pair","Cigar","Strand","Duplicate","Supplementary","Secondary","SEQ","Forward0","Forward1","Forward2","Reverse0","Reverse1","Reverse2","Orientation","BreakpointSequence","BreakpointSequenceCorrected","Type","TypeCorrected"]
        Header_mate=["Read_name_mate","Chrom_mate","Start_mate","End_mate","Chromosomal_region_mate","Seq_length_mate","Aligned_length_mate","MQ_mate","Flag_mate","Read_pair_mate","Cigar_mate","Strand_mate","Duplicate_mate","Supplementary_mate","Secondary_mate","SEQ_mate","Forward0_mate","Forward1_mate","Forward2_mate","Reverse0_mate","Reverse1_mate","Reverse2_mate","Orientation_mate","BreakpointSequence_mate","BreakpointSequenceCorrected_mate","Type","TypeCorrected_mate"]
        Expected_variables = Header_main + Header_mate + ['Fusion','Sample_read_length','Sample_total_reads','Sample_coverage', 'Filter'] 
        
        Variables_not_valid = [x for x in Variables2 if x not in Expected_variables]
        if (len(Variables_not_valid) > 0):
            Not_valid = ', '.join(Variables_not_valid)
            print ("- These set of variables are not valid:\n")
            print (f"{Variables_not_valid}\n")
            print ("- Use variables found in this list:\n")
            print (f"{Expected_variables}\n")
            return()

        # Add common variables
        Variables2 = Variables2 + ['Purity']
        if Variables == None:
            Variables = []
        Counting = pd.DataFrame(columns = Variables2 + ['TelomereFusion_counts','TelomereFusion_rate'])
        Counting.loc[0] = ['chr9_endogenous','NA',1] + ['NA'] * len(Variables) + [0] * 2
    else:
        print (f"Corrupted file: {fusion_file}\n")
        return()

    # Save results
    Counting.to_csv(outfile, sep = '\t', index=False)  

#---------------
# Running TelFusDetectorRates
#---------------

if __name__ == '__main__':
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    Seconds = round(stop - start)
    print(f"TelFusDetectorRates computation time: {Seconds} seconds\n") 





