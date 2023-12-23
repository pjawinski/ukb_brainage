#!/usr/bin/env Rscript

# =======================
# === ApoE-e4 scoring ===
# =======================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop(paste0('expected 2 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# attach required packages
for (pkg in c('dplyr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# set arguments
input = args[1] # input="results/gap_gm/apoe/apoe.vcf"
output = args[2] # output="results/gap_gm/apoe/apoe.e4.vcf"
message(paste0('\n--- Settings for ApoE-e4 scoring ---',
               '\ninput: ', input,
               '\noutput: ', output,'\n'))

# load data
df = read.table(input, skip = 8, comment.char = "", header = T) 
df.t = data.frame(t(df[,10:ncol(df)]))
names(df.t) = paste0(df$ID,'_',df$ALT)
for (i in 1:ncol(df.t)) {
	df.t[,i] = factor(df.t[,i], levels = c('0/0','0/1','1/0','1/1'), labels = c(paste0(df$REF[i],'/',df$REF[i]),paste0(df$REF[i],'/',df$ALT[i]),paste0(df$ALT[i],'/',df$REF[i]),paste0(df$ALT[i],'/',df$ALT[i])))
}

# convert to e2,e3,e4
df.t$hap = paste0(df.t[,1],' ',df.t[,2])
df.t$hap = as.character(factor(df.t$hap, 
	levels = c('C/C C/C','C/C C/T','C/C T/C','C/C T/T','C/T C/C','C/T C/T','C/T T/C','C/T T/T',
			   'T/C C/C','T/C C/T','T/C T/C','T/C T/T','T/T C/C','T/T C/T','T/T T/C','T/T T/T'),
	labels = c('e4/e4','e4/e1','e1/e4','e1/e1','e4/e3','e4/e2','e1/e3','e1/e2',
			   'e3/e4','e3/e1','e2/e4','e2/e1','e3/e3','e3/e2','e2/e3','e2/e2')))

# replace e4 by 1 and other alleles as 0
df.t$e4 = df.t$hap %>% gsub(pattern = 'e4', replacement = '1') %>%
	gsub(pattern = 'e1', replacement = '0') %>%
	gsub(pattern = 'e2', replacement = '0') %>%
	gsub(pattern = 'e3', replacement = '0')

# replace e2 by 1 and other alleles as 0
df.t$e2 = df.t$hap %>% gsub(pattern = 'e4', replacement = '0') %>%
	gsub(pattern = 'e1', replacement = '0') %>%
	gsub(pattern = 'e2', replacement = '1') %>%
	gsub(pattern = 'e3', replacement = '0')

# add snp-line
df = rbind(df,c(19,0,'apoe_e4','other','e4','.','.','PR','GT',df.t$e4))
df = rbind(df,c(19,0,'apoe_e2','other','e2','.','.','PR','GT',df.t$e2))
df = df[-c(1,2),]

# get header
header = read.delim(input, nrow = 9, comment.char = "", sep = '%', header = F, quote = "") 
header = gsub("_.*?\t","\t", header$V1) 
header = data.frame(V1 = header)
header = gsub("[_].*","", header$V1)
header = data.frame(V1 = header)

# write data
write.table(header, output, col.names = F, row.names = F, quote = F, sep = '\t')
write.table(df, output, col.names = F, row.names = F, quote = F, sep = '\t', append = T)

# ==== code snippets
# df = read.table('/slow/projects/ukb_brainage/results/gap_gm/apoe/apoe.raw', header = T) 
# df$hap=paste0(df$rs429358_C,df$rs7412_T)
# table(df$hap)
# df$apoe = as.numeric(as.character(factor(df$hap, levels = c('00','01','02','10','11','20'), labels = c(0,0,0,1,NA,2))))
# table(df$apoe)
# df.raw = df