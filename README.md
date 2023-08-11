# TelFusDetector
TelFusDetector provides functionalities for the detection of telomere fusions using whole-genome sequencing data.

# Usage

### TelFusDetectorCaller
This script detects telomere fusions in WGS DNA sequencing data sets (paired-end sequencing data). It can be run with multiple processes to speed up the computation. 

```
python scripts/TelFusDetectorCaller.py --help
usage: TelFusDetectorCaller.py [-h] --bam BAM [--genome {Hg38,Hg19}]
                               [--bam2 BAM2] [--outfolder OUTFOLDER]
                               [--sample SAMPLE] [--tmpfolder TMPFOLDER]
                               [--threads THREADS]

TelFusDetector: This tool serves to detect telomere fusions in DNA sequencing
data sets (paired-end sequencing data)

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM             BAM file to be analysed (Sorted by coordinate and
                        indexed)
  --genome {Hg38,Hg19}  Reference genome used for the alignment (Choose Hg38
                        or Hg19)
  --bam2 BAM2           BAM file with unmapped reads (Sorted by coordinate and
                        indexed). If provided, it will speed up the
                        computation
  --outfolder OUTFOLDER
                        Out directory
  --sample SAMPLE       Sample ID. All output files will start with this
                        string
  --tmpfolder TMPFOLDER
                        Directory to save all temporary files. If it exits,
                        please empty it before running TelFusDetector
  --threads THREADS     Number of threads to use [Default: 1]
```

<br>



**Out files generated**

*TelFusDetectorCaller.py* will generate the next files:
- *Sample.all_chromosomes.coverage.tsv* : File listing the read length, total number of reads and mean coverage in the sample.
- *Sample.summary_fusions.pass.tsv* : File listing the *PASS* telomere fusion calls and all characteristics. It includes the read-pairs supporting the somatic telomere fusions and the read-pairs supporting the chromosome 9 endogenous fusion.
- *Sample.summary_fusions.filtered.tsv* : File listing the *Filtered* telomere fusion calls and the reason for being filtered. 

# License
**TelFusDetector is free for academic use only.** If you are not a member of a public funded academic and/or education and/or research institution you must obtain a commercial license from EMBL Enterprise Management GmbH (EMBLEM); please email EMBLEM (info@embl-em.de).

# Contact
If you have any comments or suggestions please raise an issue or contact us:
Francesc Muyas: fmuyas@ebi.ac.uk
Isidro Cortes-Ciriano: icortes@ebi.ac.uk
