import sys
import subprocess

#as a stand-alone: [python3+] mask.py prefix iterations soft; soft is a 1 (true) or 0 (false)

def main():
    subprocess.check_call('''grep "\./\." {}.allcalls.filtered.vcf | awk '{{OFS="\t"; if ($0 !~ /\#/); print $1, $2-1, $2}}' | bedtools merge -i - > nocalls.combined.bed'''.format(prefix), shell=True)
    subprocess.check_call('''grep "allcallfilter" {}.allcalls.filtered.vcf | awk '{{OFS="\t"; if ($0 !~ /\#/); print $1, $2-1, $2}}' | bedtools merge -i - > filtered.combined.bed'''.format(prefix), shell=True)
    subprocess.check_call("cat nocalls.combined.bed filtered.combined.bed > both.bed", shell=True)
    subprocess.check_call("bedtools sort -i both.bed | bedtools merge -i - > all_positions_to_mask.bed", shell=True)
    if soft == 1:
        subprocess.check_call("bedtools maskfasta -fi {}.gatk.iteration{}.consensus.FINAL.fa -fo {}.masked.fa -bed all_positions_to_mask.bed -soft".format(prefix, iterations, prefix), shell=True)
    else:
        subprocess.check_call("bedtools maskfasta -fi {}.gatk.iteration{}.consensus.FINAL.fa -fo {}.masked.fa -bed all_positions_to_mask.bed".format(prefix, iterations, prefix), shell=True)

if __name__=="__main__":
    prefix = sys.argv[1]
    iterations = sys.argv[2]
    soft = sys.argv[3]
    print(prefix, iterations, soft)
    main()