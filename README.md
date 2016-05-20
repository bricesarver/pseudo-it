# pseudo-it
An approach that iteratively generates pseudoreferences, incorporating sample-specific variation and reducing mapping biases. Coded for Python 3.5+ but ought to work under other versions. This script manages the (potentially large) number of files and handles arguments appropriately.

This is being converted into a true (and easier to follow/extend) Python package and will eventually be replaced.

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
                    [--bed BED] [--haplotype HAPLO] [--nocall] [--iupac]
                    [--keep-haploid-reference] [--filter FIL] [--nct NCT]
                    [--nt NT]
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
  --PE1 PE1, -1 PE1     Data: PE1. PE, SE, OR PE+SE DATA IS REQUIRED (default:
                        None)
  --PE2 PE2, -2 PE2     Data: PE2. PE, SE, OR PE+SE DATA IS REQUIRED (default:
                        None)
  --SE SE, -s SE        Data: SE. PE, SE, OR PE+SE DATA IS REQUIRED (default:
                        None)
  --proc PROC, -np PROC
                        number of cores to use for multithreaded applications
                        (default: 1)
  --bed BED, -b BED     a BED file of regions to call genotypes via GATK's -L
                        (default: None)
  --haplotype HAPLO, -haplo HAPLO
                        set to 1 to use HaplotypeCaller instead of
                        UnifiedGenotyper. runtime will increase dramatically.
                        indels are still ignored. HaplotypeCaller cannot be
                        threaded (default: 0)
  --nocall              identify no-call sites and inject these into the final
                        reference. has the effect of changing bases that
                        cannot be called to Ns. WARNING: if disabled, all
                        bases that cannot be called will default to the
                        reference allele in the final iteration. This will not
                        work if you are just performing a single iteration;
                        consider running the commands sequentially (one
                        EMIT_ALL_SITES, identify nocalls and subset VCF, and a
                        FastaAltenateReferenceMaker); this functionality may
                        be introduced in subsequent versions. this DOES NOT
                        happen by default and needs to be invoked (default:
                        False)
  --iupac               invoke to inject IUPAC ambiguity codes for
                        heterozygotes into the final reference (default:
                        False)
  --keep-haploid-reference
                        if using '--iupac', this argument also keeps a haploid
                        reference. this reference is not masked (default:
                        False)
  --filter FIL, -f FIL  overwrite the default filter used to select variants.
                        you MUST specify --filterName and might want to
                        consider selecting something meaningful if you plan to
                        use the VCFs again. you can also specify multiple
                        filters by passing multiple --filterExpression and
                        --filterName arguments (will need a --filterExpression
                        for each additional filter) (default: "MQ < 30 || DP <
                        5 || DP > 60" --filterName "mq30-5dp60")
  --nct NCT             number of compute threads for the GATK's
                        UnifiedGenotyper. total CPU usage is nct*nt (default:
                        1)
  --nt NT               number of data threads for the GATK's
                        UnifiedGenotyper. total CPU usage is nct*nt (default:
                        1)
```
***
#FAQs:


**How many iterations should I perform?**

This depends on the sequence divergence of your sample relative to your reference and the type of data you have. For exome data, I found that three iterations performs well for samples with ~7.5 million years of divergence. If you expect more or have quickly evolving loci (noncoding, etc.), you might need more.


**No-call masking takes a long time! Can this be sped up?**

I'm looking into this, but it will require a multiprocessing approach and possibly the use of more sophisticated data structures.


**You don't account for indels? Why not?**

They shouldn't be included you want to keep the same coordinate system as your reference! Easy enough to do afterwards using standard variant-calling aproaches.

**I can't remember if I allowed multiallelic calls in my final reference!**


```grep "[insert IUPAC ambiguity code of choice here]" [your FINAL fasta reference]```

[1]:http://biopython.org/wiki/Biopython
[2]:http://www.htslib.org
[3]:http://www.broadinstitute.org/gatk/
[4]:http://broadinstitute.github.io/picard/
[5]:http://bio-bwa.sourceforge.net
[6]:http://bedtools.readthedocs.io/en/latest/
