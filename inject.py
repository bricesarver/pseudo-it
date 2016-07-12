#usage python[3+] inject.py [reference for the vcf] [prefix] [iterations] [iupac; 1:True, 0:False]
#splitting out the reference function for analysis of incomplete pseudoreferences
#output called [prefix].iteration[iteration].pseudo.fa and contigs renamed to ref

import subprocess
import sys
from Bio import SeqIO


def main():
    ids = []
    with open(reference) as ref:
        for record in SeqIO.parse(ref, "fasta"):
            ids.append(record.id)

    if iupac == 1:
        subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R {} -o {}.gatk.iteration{}.pseudo.fa -V {}.iteration{}.filtered.vcf -IUPAC {}'.format(reference, prefix, iterations, prefix, iterations, prefix), shell=True)
    else:
        subprocess.check_call('java -jar /usr/local/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R {} -o {}.gatk.iteration{}.pseudo.fa -V {}.iteration{}.filtered.vcf'.format(reference, prefix, iterations, prefix, iterations), shell=True)

    with open("{}.gatk.iteration{}.pseudo.fa".format(prefix, iterations), "rU") as consensus:
        finalseqs = list(SeqIO.parse(consensus, "fasta"))
        for i in range(0, len(ids), 1):
            finalseqs[i].id = finalseqs[i].name = finalseqs[i].description = ids[i]
    with open("{}.gatk.iteration{}.pseudo.fa".format(prefix, iterations), "w") as outfile:
        SeqIO.write(finalseqs, outfile, "fasta")


if __name__ == "__main__":
	reference = sys.argv[1]
    prefix = sys.argv[2]
    iterations = sys.argv[3]
    iupac = sys.argv[4]
    print(reference, prefix, iterations, iupac)
	main()