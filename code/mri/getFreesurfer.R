#!/usr/bin/env Rscript

# ================================================================================
# === get FreeSurfer information from UKB T1 surface model files (field 20263) ===
# ================================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop(paste0('expected 2 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
surfacePath=args[1] # surfacePath="data/surface/"
outFile=args[2] # outFile="results/mri/freesurfer.tsv"

# attach required packages
for (pkg in c('doParallel','dplyr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# get list of freesurfer files
files = list.files(path = "surfacePath", pattern = ".zip", recursive = TRUE)

# define cortical structures of interest
cortical = data.frame(StructName = c('bankssts','caudalanteriorcingulate','caudalmiddlefrontal','cuneus','entorhinal','fusiform','inferiorparietal','inferiortemporal','isthmuscingulate','lateraloccipital',
	'lateralorbitofrontal','lingual','medialorbitofrontal','middletemporal','parahippocampal','paracentral','parsopercularis','parsorbitalis','parstriangularis','pericalcarine','postcentral',
	'posteriorcingulate','precentral','precuneus','rostralanteriorcingulate','rostralmiddlefrontal','superiorfrontal','superiorparietal','superiortemporal','supramarginal','frontalpole',
	'temporalpole','transversetemporal','insula'))

# define subcortical structures of interest
subcortical = data.frame(StructName = c('Left-Accumbens-area','Left-Amygdala','Left-Caudate','Left-Hippocampus','Left-Pallidum','Left-Putamen','Left-Thalamus-Proper','Left-Lateral-Ventricle',
	'Right-Accumbens-area','Right-Amygdala','Right-Caudate','Right-Hippocampus','Right-Pallidum','Right-Putamen','Right-Thalamus-Proper','Right-Lateral-Ventricle'))

# create empty data frame
header = cortical$StructName
header = c(paste0('SurfArea_',header),paste0('GrayVol_',header),paste0('ThickAvg_',header))
header = c(paste0('LH_',header),paste0('RH_',header))
header = c(header,subcortical$StructName)
header = c('IID', 'Filename', 'eTIV', header)

# set up parallel pool
n.cores = 100 # parallel::detectCores() - 1
my.cluster = parallel::makeCluster(n.cores, type = "FORK")
doParallel::registerDoParallel(cl = my.cluster) # foreach::getDoParRegistered()
    
# get data
iterationsPerWorker = floor(length(files)/n.cores) # must be an integer, check by typing: iterationsPerWorker == as.integer(iterationsPerWorker)
start_time = Sys.time()
freesurfer = foreach(i = 1:n.cores, .combine='rbind') %dopar% {

	# define block for worker
	j.start = i*iterationsPerWorker-iterationsPerWorker+1
	j.end = i*iterationsPerWorker
	if (i == n.cores) { j.end = length(files) }
	    
 	# get data
    tmp = data.frame(matrix(NA, nrow = j.end-j.start+1, ncol = length(header)))
    names(tmp) = header

    # start iterations
    k = 0
    for (j in j.start:j.end) {

        # message
        skip_to_next = FALSE
        k = k+1
        filename = gsub('.*/','',files[j])
        iid = gsub('_.*','',filename)
        message(paste0('Starting with ',iid, ' (', k, '/', j.end-j.start+1, ')'))
        tmp[k,1:2] = c(iid,filename)

		# open files                        
        tryCatch(lh <- read.table(unz(paste0("surfacePath",files[j]),"FreeSurfer/stats/lh.aparc.stats")), error = function(e) { skip_to_next <<- TRUE})
        tryCatch(rh <- read.table(unz(paste0("surfacePath",files[j]),"FreeSurfer/stats/rh.aparc.stats")), error = function(e) { skip_to_next <<- TRUE})
        tryCatch(aseg <- read.table(unz(paste0("surfacePath",files[j]),"FreeSurfer/stats/aseg.stats")), error = function(e) { skip_to_next <<- TRUE})
        if(skip_to_next) { next } 

        # name columns
        names(lh) = names(rh) = c("StructName","NumVert","SurfArea","GrayVol","ThickAvg","ThickStd","MeanCurv","GausCurv","FoldInd","CurvInd")
        names(aseg) = c('Index','SegId','NVoxels','Volume_mm3','StructName','normMean','normStdDev','normMin','normMax','normRange')

        # extract information
        lh = left_join(cortical, lh, by = 'StructName')
        rh = left_join(cortical, rh, by = 'StructName')
        aseg = left_join(subcortical, aseg, by = 'StructName')
        etiv = read.delim(unz(paste0("surfacePath",files[j]),"FreeSurfer/stats/aseg.stats"), sep = '\n')
        etiv = as.numeric(gsub(', .*','',gsub('.*Volume, ','',etiv[grep('Total Intracranial Volume',etiv[,1]),1])))/1000

        # concatenate
        tmp[k,3:ncol(tmp)] = c(etiv,lh$SurfArea,lh$GrayVol,lh$ThickAvg,rh$SurfArea,rh$GrayVol,rh$ThickAvg,aseg$Volume_mm3)
    }
    tmp
}
Sys.time() - start_time

# get freesurfer variables for ENIGMA
stopCluster(my.cluster)

# count NA
# - 666 in a matrix of 43075 subjects x 222 variables
# - 4 subjects
sum(is.na(freesurfer))
sum(!complete.cases(freesurfer))

# write results
write.table(freesurfer, outFile, sep = '\t', quote = F, col.names = T, row.names = F)
system(sprintf('gzip -f %s; chmod 770 %s', outFile, outFile)

