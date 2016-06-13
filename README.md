# pseudo-it
An approach that iteratively generates pseudoreferences, incorporating sample-specific variation and reducing mapping biases. Coded for Python 3.5+ but ought to work under other versions. This script manages the (potentially large) number of files and handles arguments appropriately.

This is being converted into a true (and easier to follow/extend) Python package and will eventually be replaced.

***I am currently tracking down I/O issues and/or bugs in various applications used for the nocall masking, possibly with bedtools sort. Be aware that using --nocall is memory intensive for full genomes after reverting back to a non-streaming approach. I am converting masking to its own module for ease in logging and exception handling. You could also pull the ```devel``` branch which uses the old stream-based processing.***

Requires: 

1. [Biopython][1]
2. [Samtools][2]
3. [The Genome Analysis Toolkit][3]
4. [Picard][4]
5. [BWA][5]
6. [bedtools][6]

The GATK and Picard jarfiles are expected to be in /usr/local/bin; the rest of the executables are expected to be in your $PATH. The system calls might require some modification depending on how your system is configured. I will be incorporating the ability to easily override default paths (if your administrator has things set up differently, for example) in a later version.

Please report issues via 'Issues' above, or send me an email.

***
```
usage: pseudo-it.py [-h] [--PE1 PE1] [--PE2 PE2] [--SE SE] [--proc PROC]
                    [--bed BED] [--haplotype] [--nocall] [--nocall-filter NCF]
                    [--soft-masking] [--iupac] [--keep-haploid-reference]
                    [--filter FIL] [--nct NCT] [--nt NT]
                    iterations reference prefix

Iterative pseudoreference generation with BWA and the GATK

positional arguments:
  iterations            number of iterations. one iteration will not inject
                        IUPAC ambiguities
  reference             path to the reference/contigs/scaffolds used for the
                        first iteration
  prefix                prefix to use on output files; this is also added to
                        the SM field with AddOrReplaceReadGroups

optional arguments:
  -h, --help            show this help message and exit
  --proc PROC, -np PROC
                        number of cores to use for multithreaded applications
                        (default: 1)
  --bed BED, -b BED     a BED file of regions to call genotypes via GATK's -L
                        (default: None)
  --haplotype           invoke to use HaplotypeCaller instead of
                        UnifiedGenotyper. runtime will increase dramatically.
                        indels are still ignored. HaplotypeCaller cannot be
                        threaded (default: False)
  --nocall              identify nocall and low-quality sites and mask these
                        in the final reference. has the effect of changing
                        bases that cannot be called to N, by default. requires
                        more than a single iteration at currently; this
                        functionality may be introduced in subsequent versions
                        (default: False)
  --nocall-filter NCF, -ncf NCF
                        additional filtering threshold for low-quality bases
                        to be used for the masking step (default:
                        --filterExpression "MQ < 30.0 || DP < 10 || DP > 60")
  --soft-masking        soft mask (i.e., replace with lowercase) instead of
                        hard mask (i.e., replace with N). requires `--nocall`
                        (default: False)
  --iupac               invoke to inject IUPAC ambiguity codes for
                        heterozygotes into the final reference (default:
                        False)
  --keep-haploid-reference
                        if using '--iupac', this argument also keeps a haploid
                        reference this reference is not masked (default:
                        False)
  --filter FIL, -f FIL  overwrite the default filter used to select variants.
                        you MUST specify --filterName and might want to
                        consider selecting something meaningful if you plan to
                        use the VCFs again. you can also specify multiple
                        filters by passing multiple --filterExpression and
                        --filterName arguments (will need a --filterExpression
                        for each additional filter) (default: "MQ < 30.0 || DP
                        < 5 || DP > 60" --filterName "mq30-5dp60")
  --nct NCT             number of compute threads for the GATK's
                        UnifiedGenotyper. total CPU usage is nct*nt (default:
                        1)
  --nt NT               number of data threads for the GATK's
                        UnifiedGenotyper. total CPU usage is nct*nt (default:
                        1)

required arguments:
  --PE1 PE1, -1 PE1     Data: PE1. PE, SE, OR PE+SE DATA IS REQUIRED (default:
                        None)
  --PE2 PE2, -2 PE2     Data: PE2. PE, SE, OR PE+SE DATA IS REQUIRED (default:
                        None)
  --SE SE, -s SE        Data: SE. PE, SE, OR PE+SE DATA IS REQUIRED (default:
                        None)
```
***
#FAQs:


**How many iterations should I perform?**

This depends on the sequence divergence of your sample relative to your reference and the type of data you have. For exome data, I found that three iterations performs well for samples with ~7.5 million years of divergence. If you expect more or have quickly evolving loci (noncoding, etc.), you might need more.


**No-call masking takes a long time! Can this be sped up?**

I'm looking into this, but it will require a multiprocessing approach and possibly the use of more sophisticated data structures.


**You don't account for indels? Why not?**

They shouldn't be included you want to keep the same coordinate system as your reference! Easy enough to do afterwards using standard variant-calling aproaches.

**I can't remember if I allowed ambiguities in my final reference!**

Besides looking at the name of the file, you can
```grep "[insert IUPAC ambiguity code of choice here]" [your FINAL fasta reference]```

[1]:http://biopython.org/wiki/Biopython
[2]:http://www.htslib.org
[3]:http://www.broadinstitute.org/gatk/
[4]:http://broadinstitute.github.io/picard/
[5]:http://bio-bwa.sourceforge.net
[6]:http://bedtools.readthedocs.io/en/latest/
