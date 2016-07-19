library(Biostrings)
library(gtools)
library(parallel)

#usage: Rscript --nosave [genome_file.fa] [gtf_file.gtf] [transcript_ids.txt] [number of processors] [sample ID (for fasta output)]

args <- commandArgs(TRUE)
#a genome file to subset
genomefile <- args[1]
#a gtf associated with the genome file
gtffile <- args[2]
#a list of canonical transcript IDs (or any transcript IDs) in the GTF. header="Ensembl.Transcript.ID"
canonfile <- args[3]
#number of processors for multithreaded applications
proc <- args[4]
#the name of the individuals
libraryid <- args[5]

print(args)

print("reading genome...")
genome <- readDNAStringSet(genomefile, format="fasta")
genome <- genome[mixedsort(names(genome))]
print("reading gtf...")
gtf <- read.table(gtffile, sep="\t", stringsAsFactors=FALSE)
print("reading IDs...")
canon <- read.table("canonical_transcript_ids.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

ids <- unique(canon$Ensembl.Transcript.ID)

cds <- gtf[which(gtf$V3 == "CDS"), ]

obtainCanonicalSequence <- function(transcriptID){
	short <- cds[grep(pattern=transcriptID, cds$V9), ]
	if (nrow(short) == 0) return(1)
	len <- nrow(short)
	hld <- list()
	for (i in 1:len){
		wrk <- short[i, ]
		hld[[i]] <- subseq(genome[paste("chr", wrk$V1, sep="")], start=as.numeric(wrk$V4), end=as.numeric(wrk$V5))
		if (wrk$V7 == "-") hld[[i]] <- reverseComplement(hld[[i]])
	}

	combined <- unlist(as(sapply(hld, "[[", 1L), "DNAStringSet"))
	combined <- DNAStringSet(combined)
	if ((width(combined) %% 3) != 0) return(2) 
	chr <- short[1, 1]
	start <- short[1, 4]
	names(combined) <- libraryid
	#return(combined)
	writeXStringSet(combined, filepath=paste(transcriptID, chr, start, "fa", sep="."), append=FALSE, format="fasta")
	return(0)
}

print("principal execution...")
mclapply(ids, obtainCanonicalSequence, mc.cores=proc)