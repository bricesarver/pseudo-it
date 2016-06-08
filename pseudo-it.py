#pseudo-it.py v1.0.4
#iterative pseudoreference generation using BWA and the GATK
#Brice A. J. Sarver
#v0.1.0 completed 3 June 2015
#v0.8.0 release 18 Jan 2016; bug fixes, minor functionality improvements
#v1.0.1 release 8 May 2016; major functionality improvements, rewrite with Python 3.5+, bug fixes

import sys
import os
import argparse
import itertools
import subprocess
import fnmatch
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Iterative pseudoreference generation with BWA and the GATK", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('iterations', help='number of iterations. one iteration will not inject IUPAC ambiguities', type=int)
parser.add_argument('reference', help='path to the reference/contigs/scaffolds used for the first iteration')
parser.add_argument('prefix', help="prefix to use on output files; this is also added to the SM field with AddOrReplaceReadGroups")
required = parser.add_argument_group('required arguments')
required.add_argument('--PE1', '-1', dest='pe1', help="Data: PE1. PE, SE, OR PE+SE DATA IS REQUIRED", default=None)
required.add_argument('--PE2', '-2', dest='pe2', help="Data: PE2. PE, SE, OR PE+SE DATA IS REQUIRED", default=None)
required.add_argument('--SE', '-s', dest='se', help="Data: SE. PE, SE, OR PE+SE DATA IS REQUIRED", default=None)
parser.add_argument('--proc', '-np', dest='proc', type=int, help='number of cores to use for multithreaded applications', default=1)
parser.add_argument('--bed', '-b', dest='bed', help="a BED file of regions to call genotypes via GATK's -L", default=None)
parser.add_argument('--haplotype', dest='haplo', help="invoke to use HaplotypeCaller instead of UnifiedGenotyper. runtime will increase dramatically. indels are still ignored. HaplotypeCaller cannot be threaded", action='store_true')
parser.add_argument('--nocall', dest='nocall', help="identify no-call sites and inject these into the final reference. has the effect of changing bases that cannot be called to Ns. WARNING: if disabled, all bases that cannot be called will default to the reference allele in the final iteration. This will not work if you are just performing a single iteration; consider running the commands sequentially (one EMIT_ALL_SITES, identify nocalls and subset VCF, and a FastaAltenateReferenceMaker); this functionality may be introduced in subsequent versions. this DOES NOT happen by default and needs to be invoked", action='store_true')
parser.add_argument('--nocall-filter', '-ncf', dest='ncf', help='additional filtering to be used for the masking step', default='--filterExpression "MQ < 30.0 || DP < 10 || DP > 60"')
parser.add_argument('--iupac', dest='iupac', help='invoke to inject IUPAC ambiguity codes for heterozygotes into the final reference', action='store_true')
parser.add_argument('--keep-haploid-reference', dest='haploid', help="if using '--iupac', this argument also keeps a haploid reference this reference is not masked", action='store_true')
parser.add_argument('--filter', '-f', dest='fil', help='overwrite the default filter used to select variants. you MUST specify --filterName and might want to consider selecting something meaningful if you plan to use the VCFs again. you can also specify multiple filters by passing multiple --filterExpression and --filterName arguments (will need a --filterExpression for each additional filter)', default='"MQ < 30.0 || DP < 5 || DP > 60" --filterName "mq30-5dp60"')
parser.add_argument('--nct', dest='nct', help="number of compute threads for the GATK's UnifiedGenotyper. total CPU usage is nct*nt", default=1)
parser.add_argument('--nt', dest='nt', help="number of data threads for the GATK's UnifiedGenotyper. total CPU usage is nct*nt", default=1)
#parser.add_argument('--resume', dest="resumeRun", help="resume a previously failed run [WORK IN PROGRESS]")
#parser.add_argument('--clean', dest="cleanRun", help="blow away all intermediate files")

args = parser.parse_args()

#argument sanity checks
if args.iterations == 1 and args.nocall:
    sys.exit("One iteration and injection of nocalls is not currently supported. If iterations == 1, please do not invoke --nocall")
elif args.haploid and not args.iupac:
    sys.exit("You cannot specify an additional haploid reference without --iupac (the process normally generates a haploid reference)")
elif not args.pe1 and not args.pe2 and not args.se:
    sys.exit("You need to specify data with the --PE1, --PE2, and --SE options")
elif args.pe1 and not args.pe2:
    sys.exit("You specified PE1 but not its mate (PE2)")
elif not args.pe1 and args.pe2:
    sys.exit("You specified PE2 but not its mate (PE1)")

#print(args)

##############################
#FUTURE TO-DO: have python attempt to locate the executables itself and pass these to subprocess.
#              as-is, this script is not really portable because executable paths are hard-coded (but
#              probably housed in /usr/local/bin on many systems. easy to change here with sed or equivalent)
#              perhaps pass as a file (config.txt)?

#picardpath = ""
#gatkpath = ""
#bwapath = ""
#samtoolspath = ""
#gatkpath = ""
#picard = ""
#samtools = ""
#bwa = ""
#GATK = ""
# 
#               additionally, all of this was coded modularly on purpose so I can easily extend it to a true package
#               especially, the contig renaming function needs to be kicked out as its own thing
#               perhaps still get IDs and store that list as a global since the footprint can be large for eukaryotes?
#
#               the haploid in addition to diploid reference output was an afterthought and can probably be
#               implemented better...
##############################

def make_indices(reference):
    print("Generate dictionaries and indices (if needed)...")
    #for ease of use, this picks everything before the first period
    #consider symlinking to 'species.fa' or something if this is an issue (e.g., lots of '.')
    dictprefix = os.path.splitext(reference)[0]
    if not os.path.isfile('{}.dict'.format(dictprefix)):
        subprocess.check_call('java -jar /usr/local/bin/picard.jar CreateSequenceDictionary R={} O={}.dict'.format(reference, dictprefix), shell=True)
    if not os.path.isfile('{}.fai'.format(reference)):
        print("faidx...")
        subprocess.check_call('samtools faidx {}'.format(reference), shell=True)
    if not os.path.isfile('{}.bwt'.format(reference)):
        print("BWA index...")
        subprocess.check_call('bwa index {}'.format(reference), shell=True)

def first_iteration(iterations, reference, prefix, proc, bed, haplo, fil, pe1, pe2, se, nct, nt):
    finalseqs = []
    #if bed file exists, define it as a variable now
    if bed:
        bedoption = "-L " + bed
    else:
        bedoption = ""

    if nct:
        nct = '-nct {}'.format(nct)
    else:
        nct = ""

    if nt:
        nt = '-nt {}'.format(nt)
    else:
        nt = ""

    #create dictionaries, etc.
    make_indices(reference)
    #map PE and SE reads | convert to bam
    if pe1 and pe2:
        print("PE BWA map | convert to BAM...")
        subprocess.check_call('bwa mem -M -t {} {} {} {} | samtools view -Sb - > {}.iteration1.pe.bam 2> {}.iteration1.pe.bam.stderr'.format(proc, reference, pe1, pe2, prefix, prefix), shell=True)
    if se:
        print("SE BWA map | conver to BAM...")
        subprocess.check_call('bwa mem -M -t {} {} {} | samtools view -Sb - > {}.iteration1.se.bam 2> {}.iteration1.se.bam.stderr'.format(proc, reference, se, prefix, prefix), shell=True)
    
    if pe1 and pe2 and se:
        print("Merge BAMs...")
        subprocess.check_call('java -jar /usr/local/bin/picard.jar MergeSamFiles I={}.iteration1.pe.bam I={}.iteration1.se.bam O={}.iteration1.merged.bam USE_THREADING=TRUE VALIDATION_STRINGENCY=LENIENT'.format(prefix, prefix, prefix), shell=True)
    elif pe1 and pe2 and not se:
        os.rename('{}.iteration1.pe.bam'.format(prefix), '{}.iteration1.merged.bam'.format(prefix))
    elif se and not pe1 and not pe2:
        os.rename('{}.iteration1.se.bam'.format(prefix), '{}.iteration1.merged.bam'.format(prefix))
    
    #add readgroups and mark dups
    print("AddOrReplaceReadGroups...")
    subprocess.check_call('java -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups I={}.iteration1.merged.bam O={}.iteration1.merged.RG.bam SO=coordinate LB=spret_exome PL=illumina PU=misc SM={} VALIDATION_STRINGENCY=LENIENT'.format(prefix, prefix, prefix), shell=True)
    print("Mark duplicates...")
    subprocess.check_call('java -jar /usr/local/bin/picard.jar MarkDuplicates I={}.iteration1.merged.RG.bam O={}.iteration1.merged.RG_dedup.bam VALIDATION_STRINGENCY=LENIENT M=iteration1.dup_metrics'.format(prefix, prefix), shell=True)
    
    #indel realignment
    print("Index BAMs...")
    subprocess.check_call('samtools index {}.iteration1.merged.RG_dedup.bam'.format(prefix), shell=True)
    print("Identify targets to realign...")
    subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {} -I {}.iteration1.merged.RG_dedup.bam -o iteration1.indel_intervals.list {} -nt {}'.format(reference, prefix, bedoption, proc), shell=True)
    print("Realign Indels...")
    subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T IndelRealigner -R {} -I {}.iteration1.merged.RG_dedup.bam -targetIntervals iteration1.indel_intervals.list -o {}.iteration1.realigned.bam --filter_bases_not_stored'.format(reference, prefix, prefix), shell=True)
    
    #variant calling
    if haplo:
        print("HaplotypeCaller...")
        subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T HaplotypeCaller -R {} -I {}.iteration1.realigned.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o {}.iteration1.raw.vcf {}'.format(reference, prefix, prefix, bedoption), shell=True)
    else:
        print("UnifiedGenotyper...")
        subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T UnifiedGenotyper -R {} -I {}.iteration1.realigned.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o {}.iteration1.raw.vcf {} {} {}'.format(reference, prefix, prefix, bedoption, nct, nt), shell=True)
    
    #selecting SNPs and filtering
    print("Select SNPs from VCF...")
    subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T SelectVariants -R {} -V {}.iteration1.raw.vcf -o {}.iteration1.snps.vcf --selectTypeToInclude SNP'.format(reference, prefix, prefix), shell=True)
    print("Filter variants...")
    subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T VariantFiltration -R {} -V {}.iteration1.snps.vcf --filterExpression {} -o {}.iteration1.filtered.vcf'.format(reference, prefix, fil, prefix), shell=True)
    
    #consensus calling and filtering
    print("Generate consensus from reference and variants...")
    subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R {} -o {}.gatk.iteration1.consensus.fa -V {}.iteration1.filtered.vcf'.format(reference, prefix, prefix), shell=True)

    with open("{}.gatk.iteration1.consensus.fa".format(prefix), "rU") as consensus:
        finalseqs = list(SeqIO.parse(consensus, "fasta"))
        for i in range(0, len(ids), 1):
            finalseqs[i].id = finalseqs[i].name = finalseqs[i].description = ids[i]
    with open("{}.gatk.iteration1.consensus.fa".format(prefix), "w") as outfile:
        SeqIO.write(finalseqs, outfile, "fasta")

 
def other_iterations(iterations, prefix, proc, totalIterations, bed, haplo, ncall, iupac, fil, pe1, pe2, se, nct, nt, haploid, ncf):
    finalseqs = []
    #again, define BED argument
    if bed:
        bedoption = "-L " + bed
    else:
        bedoption = ""

    if nct:
        nct = '-nct {}'.format(nct)
    else:
        nct = ""

    if nt:
        nt = '-nt {}'.format(nt)
    else:
        nt = ""

    previousIteration = iterations - 1
    reference = '{}.gatk.iteration{}.consensus.fa'.format(prefix, previousIteration)

    make_indices(reference)
    
    if pe1 and pe2:
        print("PE BWA map | covert to BAM...")
        subprocess.check_call('bwa mem -M -t {} {} {} {} | samtools view -Sb - > {}.iteration{}.pe.bam 2> {}.iteration{}.pe.bam.stderr'.format(proc, reference, pe1, pe2, prefix, iterations, prefix, iterations), shell=True)

    if se:
        print("SE BWA map | convert to BAM...")
        subprocess.check_call('bwa mem -M -t {} {} {} | samtools view -Sb - > {}.iteration{}.se.bam 2> {}.iteration{}.se.bam.stderr'.format(proc, reference, se, prefix, iterations, prefix, iterations), shell=True)
    
    if pe1 and pe2 and se:
        print("Merge BAMs...")
        subprocess.check_call('java -jar /usr/local/bin/picard.jar MergeSamFiles I={}.iteration{}.pe.bam I={}.iteration{}.se.bam O={}.iteration{}.merged.bam USE_THREADING=TRUE VALIDATION_STRINGENCY=LENIENT'.format(prefix, iterations, prefix, iterations, prefix, iterations), shell=True)
    elif pe1 and pe2 and not se:
        os.rename('{}.iteration1.pe.bam'.format(prefix), '{}.iteration{}.merged.bam'.format(prefix, iterations))
    elif se and not pe1 and not pe2:
        os.rename('{}.iteration1.se.bam'.format(prefix), '{}.iteration{}.merged.bam'.format(prefix, iterations))

    #sort bam
    print("Sorting BAMs...")
    subprocess.check_call('samtools sort -o {}.iteration{}.merged.sorted.bam -T hold.sorting -@ {} {}.iteration{}.merged.bam'.format(prefix, iterations, proc, prefix, iterations), shell=True)
    
    #add readgroups and mark dups
    print("AddOrReplaceReadGroups...")
    subprocess.check_call('java -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups I={}.iteration{}.merged.bam O={}.iteration{}.merged.RG.bam SO=coordinate LB=spret_exome PL=illumina PU=misc SM={} VALIDATION_STRINGENCY=LENIENT'.format(prefix, iterations, prefix, iterations, prefix), shell=True)
    print("Mark duplicates...")
    subprocess.check_call('java -jar /usr/local/bin/picard.jar MarkDuplicates I={}.iteration{}.merged.RG.bam O={}.iteration{}.merged.RG_dedup.bam VALIDATION_STRINGENCY=LENIENT M=iteration1.dup_metrics'.format(prefix, iterations, prefix, iterations), shell=True)
    
    #indel realignment
    print("Index BAMs...")
    subprocess.check_call('samtools index {}.iteration{}.merged.RG_dedup.bam'.format(prefix, iterations), shell=True)
    print("Identify targets to realign...")
    subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {} -I {}.iteration{}.merged.RG_dedup.bam -o iteration{}.indel_intervals.list {} -nt {}'.format(reference, prefix, iterations, iterations, bedoption, proc), shell=True)
    print("Realign Indels...")
    subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T IndelRealigner -R {} -I {}.iteration{}.merged.RG_dedup.bam -targetIntervals iteration{}.indel_intervals.list -o {}.iteration{}.realigned.bam --filter_bases_not_stored'.format(reference, prefix, iterations, iterations, prefix, iterations), shell=True)
    
    #variant calling
    if haplo:
        print("HaplotypeCaller...")
        subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T HaplotypeCaller -R {} -I {}.iteration{}.realigned.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o {}.iteration{}.raw.vcf {}'.format(reference, prefix, iterations, prefix, iterations, bedoption), shell=True)
    else:
        print("UnifiedGenotyper...")
        subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T UnifiedGenotyper -R {} -I {}.iteration{}.realigned.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o {}.iteration{}.raw.vcf {} {} {}'.format(reference, prefix, iterations, prefix, iterations, bedoption, nct, nt), shell=True)
    
    #selecting SNPs and filtering
    print("Select SNPs from VCF...")
    subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T SelectVariants -R {} -V {}.iteration{}.raw.vcf -o {}.iteration{}.snps.vcf --selectTypeToInclude SNP'.format(reference, prefix, iterations, prefix, iterations), shell=True)
    print("Filter variants...")
    subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T VariantFiltration -R {} -V {}.iteration{}.snps.vcf --filterExpression {} -o {}.iteration{}.filtered.vcf'.format(reference, prefix, iterations, fil, prefix, iterations), shell=True)
    
    #consensus calling and filtering; write out IUPAC ambiguities on last iteration if indicated
    if totalIterations == iterations:
        if ncall:
            if iupac:
                print("Generate consensus from reference and variants...")
                subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R {} -o {}.gatk.iteration{}.consensus.FINAL.fa -V {}.iteration{}.filtered.vcf -IUPAC {}'.format(reference, prefix, iterations, prefix, iterations, prefix), shell=True)
                if haploid:
                    subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R {} -o {}.gatk.iteration{}.consensus.FINAL.haploid.fa -V {}.iteration{}.filtered.vcf'.format(reference, prefix, iterations, prefix, iterations), shell=True)
            else:
                print("Generate consensus from reference and variants...")
                subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R {} -o {}.gatk.iteration{}.consensus.FINAL.fa -V {}.iteration{}.filtered.vcf'.format(reference, prefix, iterations, prefix, iterations), shell=True)

            with open("{}.gatk.iteration{}.consensus.FINAL.fa".format(prefix, iterations), "rU") as consensus:
                finalseqs = list(SeqIO.parse(consensus, "fasta"))
                for i in range(0, len(ids), 1):
                    finalseqs[i].id = finalseqs[i].name = finalseqs[i].description = ids[i]
            with open("{}.gatk.iteration{}.consensus.FINAL.fa".format(prefix, iterations), "w") as outfile:
                SeqIO.write(finalseqs, outfile, "fasta")

            print("Emitting all sites...")
            subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T UnifiedGenotyper -R {} -I {}.iteration{}.realigned.bam --genotyping_mode DISCOVERY --output_mode EMIT_ALL_SITES -stand_emit_conf 10 -stand_call_conf 30 -o {}.allcalls.vcf {} {}'.format(reference, prefix, iterations, prefix, nct, nt), shell=True)
            subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T VariantFiltration -R {} -V {}.allcalls.vcf {} --filterName "allcallfilter" -o {}.allcalls.filtered.vcf'.format(reference, prefix, ncf, prefix)) 
            
            ##FOR LATER START HERE
            print("filtering of nocalls...")
            #whip up a quick BED from the VCF using awk; this appears fastest even though it's a system call
            subprocess.check_call('''grep "\./\." {}.allcalls.filtered.vcf | awk '{{OFS="\t"; if ($0 !~ /\#/); print $1, $2-1, $2}}' > {}.nocalls.positions.bed'''.format(prefix, prefix), shell=True)

            subprocess.check_call("split -d -l 100000000 {}.nocalls.positions.bed {}.nocalls.split".format(prefix, prefix), shell=True)

            files = os.listdir('.')
            matches = fnmatch.filter(files, "{}.nocalls.split*".format(prefix))
            maskiteration = 0
            if len(matches) == 1:
                subprocess.check_call("bedtools maskfasta -fi {}.gatk.iteration{}.consensus.FINAL.fa -fo {}.masked.fa -bed {}.nocalls.split00".format(prefix, iterations, prefix, prefix), shell=True)
            else:
                for match in matches:
                    if maskiteration == 0:
                        subprocess.check_call("bedtools maskfasta -fi {}.gatk.iteration{}.consensus.FINAL.fa -fo {}.masked.fa -bed {}".format(prefix, iterations, prefix, match), shell=True)
                        maskiteration = 1
                    else:
                        subprocess.check_call("bedtools maskfasta -fi {}.masked.fa -fo tmp.fa -bed {}".format(prefix, match), shell=True)
                        os.rename("tmp.fa", "{}.masked.fa".format(prefix))
        elif iupac:
            subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R {} -o {}.gatk.iteration{}.consensus.FINAL.fa -V {}.iteration{}.filtered.vcf -IUPAC {}'.format(reference, prefix, iterations, prefix, iterations, prefix), shell=True)
            if haploid:
                print("Generate consensus from reference and variants...")
                subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R {} -o {}.gatk.iteration{}.consensus.FINAL.haploid.fa -V {}.iteration{}.filtered.vcf'.format(reference, prefix, iterations, prefix, iterations), shell=True)

            with open("{}.gatk.iteration{}.consensus.FINAL.fa".format(prefix, iterations), "rU") as consensus:
                finalseqs = list(SeqIO.parse(consensus, "fasta"))
                for i in range(0, len(ids), 1):
                    finalseqs[i].id = finalseqs[i].name = finalseqs[i].description = ids[i]
            with open("{}.gatk.iteration{}.consensus.FINAL.fa".format(prefix, iterations), "w") as outfile:
                SeqIO.write(finalseqs, outfile, "fasta")
        else:
            subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R {} -o {}.gatk.iteration{}.consensus.FINAL.fa -V {}.iteration{}.filtered.vcf'.format(reference, prefix, iterations, prefix, iterations), shell=True)

            with open("{}.gatk.iteration{}.consensus.FINAL.fa".format(prefix, iterations), "rU") as consensus:
                finalseqs = list(SeqIO.parse(consensus, "fasta"))
                for i in range(0, len(ids), 1):
                    finalseqs[i].id = finalseqs[i].name = finalseqs[i].description = ids[i]
            with open("{}.gatk.iteration{}.consensus.FINAL.fa".format(prefix, iterations), "w") as outfile:
                SeqIO.write(finalseqs, outfile, "fasta")

    else:
        print("Generate consensus from reference and variants...")
        subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R {} -o {}.gatk.iteration{}.consensus.fa -V {}.iteration{}.filtered.vcf'.format(reference, prefix, iterations, prefix, iterations), shell=True)

        with open("{}.gatk.iteration{}.consensus.fa".format(prefix, iterations), "rU") as consensus:
            finalseqs = list(SeqIO.parse(consensus, "fasta"))
            for i in range(0, len(ids), 1):
                finalseqs[i].id = finalseqs[i].name = finalseqs[i].description = ids[i]
        with open("{}.gatk.iteration{}.consensus.fa".format(prefix, iterations), "w") as outfile:
            SeqIO.write(finalseqs, outfile, "fasta")

    #if present, rename the haploid references, too        
    if os.path.isfile('{}.gatk.iteration{}.consensus.FINAL.haploid.fa'.format(prefix, iterations)):
        with open("{}.gatk.iteration{}.consensus.FINAL.haploid.fa".format(prefix, iterations), "rU") as consensus:
            hapfinalseqs = list(SeqIO.parse(consensus, "fasta"))
            for i in range(0, len(ids), 1):
                hapfinalseqs[i].id = hapfinalseqs[i].name = hapfinalseqs[i].description = ids[i]
        with open("{}.gatk.iteration{}.consensus.FINAL.haploid.fa".format(prefix, iterations), "w") as outfile:
            SeqIO.write(hapfinalseqs, outfile, "fasta")

       
#iteration 1 function call
#FastaAlternateReferenceMaker in the GATK renames the contigs. 'ids' faciliates renaming after injection
ids = []
with open(args.reference) as ref:
    for record in SeqIO.parse(ref, "fasta"):
        ids.append(record.id)

first_iteration(iterations=args.iterations, reference=args.reference, prefix=args.prefix, proc=args.proc, bed=args.bed, haplo=args.haplo, fil=args.fil, pe1=args.pe1, pe2=args.pe2, se=args.se, nct=args.nct, nt=args.nt)

#and the other iterations
if args.iterations > 1:
    for i in range(2, args.iterations + 1, 1):
       other_iterations(iterations=i, prefix=args.prefix, proc=args.proc, totalIterations=args.iterations, bed=args.bed, haplo=args.haplo, ncall=args.nocall, iupac=args.iupac, fil=args.fil, pe1=args.pe1, pe2=args.pe2, se=args.se, nct=args.nct, nt=args.nt, haploid=args.haploid, ncf=args.ncf)

