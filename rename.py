import sys
from Bio import SeqIO

#as a stand-alone: [python3+] rename.py [original_reference] [reference_to_rename]

def main():
    #get original IDs
    ids = []
    with open(original, "rU") as ref:
        for record in SeqIO.parse(ref, "fasta"):
            ids.append(record.id)

    #get rename IDs
    with open(toRename, "rU") as new:
        newseqs = list(SeqIO.parse(new, "fasta"))
        for i in range(0, len(ids), 1):
            newseqs[i].id = newseqs[i].name = newseqs[i].description = ids[i]

    #write new names and seqs to file
    with open(toRename, "w") as outfile:
        SeqIO.write(newseqs, outfile, "fasta")


if __name__=="__main__":
    original = sys.argv[1]
    toRename = sys.argv[2]
    print("Renaming: ", toRename, " based on: ", original)
    main()