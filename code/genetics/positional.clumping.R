#!/usr/bin/env Rscript
# ===========================
# === positional clumping ===
# ===========================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=10) {
  stop(paste0('expected 10 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
inFile = args[1] # inFile="results/gap_gm/replicate/mrmega.all/mrmega.weights.gz"
outFile = args[2] # outFile="results/gap_gm/replicate/mrmega.all/mrmega.weights.clumped.gz"
idCol = args[3] # idCol="ID"
chrCol = args[4] # chrCol="CHR"
bpCol = args[5] # bpCol="BP"
pCol = args[6] # pCol="P.value_association"
maxP = as.numeric(args[7]) # maxP=5e-8
decrease = args[8] # decrease=FALSE
windowSize = as.numeric(args[9]) # windowSize=500000
printLeadOnly = args[10] # printLeadOnly=TRUE

message(paste0('\n--- positional clumping ---',
               '\ninFile: ', inFile,
               '\noutFile: ', outFile,
               '\nidCol: ', idCol,
               '\nchrCol: ', chrCol,
               '\nbpCol: ', bpCol,
               '\npCol: ', pCol,
               '\nmaxP: ', maxP,
               '\ndecrease ', decrease,
               '\nwindowSize: ', windowSize,
               '\nprintLeadOnly: ', printLeadOnly,'\n'))

# load file and sort
message(paste0('Loading ', inFile))
if (stringr::str_sub(inFile,-3,-1)=='.gz') {
  df = data.frame(data.table::fread(cmd=paste0("gzip -dc ", inFile), tmpdir = getwd(), header=T, stringsAsFactors=FALSE))
} else {
  df = read.delim(targetFile, header=T, stringsAsFactors=FALSE)
}
df = df[df[,pCol] < maxP,]
df = df[order(df[,pCol],decreasing = decrease),]

# clump by position
message('Starting clumping procedure')
df$LEAD = 1
df$LOCUS = 0
k = 0
for (i in 1:nrow(df)) {
	if (df$LEAD[i] == 1) {
		k = k + 1
		df$LOCUS[i] = k
		if (i != nrow(df)) {
			chr = df[i,chrCol]
			bp = df[i,bpCol]
			for (j in (i+1):nrow(df)) {
				if (df$LEAD[j] == 1 & df[j,chrCol] == chr & abs(bp-df[j,bpCol]) < windowSize) {
					df$LEAD[j] = 0
					df$LOCUS[j] = k
				}
			}
		}
	}	
}

# print lead only
if (printLeadOnly == TRUE) {
	df = df[df$LEAD==1,]
}

# write output
message(paste0('Saving ', outFile))
if (stringr::str_sub(outFile,-3,-1)=='.gz') {
	data.table::fwrite(df, outFile, sep = '\t', compress = 'gzip')
} else {
	data.table::fwrite(df, outFile, sep = '\t', compress = 'none')
}
system(sprintf('chmod 770 %s',outFile))
message('--- Completed: positional clumping ---')
