#!/usr/bin/env Rscript
# ====================================
# === prepare for machine learning ===
# ====================================

# load discovery data
message('Loading discovery data.')
gm = data.frame(data.table::fread('results/mri/prepML.discovery.gm.txt', sep='\t', header=F, quote="\"", stringsAsFactors=FALSE))
wm = data.frame(data.table::fread('results/mri/prepML.discovery.wm.txt', sep='\t', header=F, quote="\"", stringsAsFactors=FALSE))
u = read.delim('results/mri/prepML.discovery.u.txt', sep='\t', header=F, quote="\"", stringsAsFactors=FALSE)
v = read.delim('results/mri/prepML.discovery.v.txt', sep='\t', header=F, quote="\"", stringsAsFactors=FALSE)
covs = read.delim('results/mri/prepML.discovery.covs.txt', sep='\t', header=T, quote="\"", stringsAsFactors=FALSE)
meta = read.delim('results/mri/prepML.discovery.meta.txt', sep='\t', header=T, quote="\"", stringsAsFactors=FALSE)

# save discovery data and clean up
message('Saving discovery data in .RData file and cleaning up.')
save(gm, wm, u, v, covs, meta, file = 'results/mri/prepML.discovery.RData')
system('rm -f results/mri/prepML.discovery.*.txt')

# load replication data
message('Loading replication data.')
gm = data.frame(data.table::fread('results/mri/prepML.replication.gm.txt', sep='\t', header=F, quote="\"", stringsAsFactors=FALSE))
wm = data.frame(data.table::fread('results/mri/prepML.replication.wm.txt', sep='\t', header=F, quote="\"", stringsAsFactors=FALSE))
covs = read.delim('results/mri/prepML.replication.covs.txt', sep='\t', header=T, quote="\"", stringsAsFactors=FALSE)
meta = read.delim('results/mri/prepML.replication.meta.txt', sep='\t', header=F, quote="\"", stringsAsFactors=FALSE)

# save replication data and clean up
message('Saving replication data and cleaning up.')
save(gm, wm, covs, meta, file = 'results/mri/prepML.replication.RData')
system('rm -f results/mri/prepML.replication.*.txt')
message('Completed.')

# load retest data
message('Loading retest data.')
gm = data.frame(data.table::fread('results/mri/prepML.retest.gm.txt', sep='\t', header=F, quote="\"", stringsAsFactors=FALSE))
wm = data.frame(data.table::fread('results/mri/prepML.retest.wm.txt', sep='\t', header=F, quote="\"", stringsAsFactors=FALSE))
covs = read.delim('results/mri/prepML.retest.covs.txt', sep='\t', header=T, quote="\"", stringsAsFactors=FALSE)
meta = read.delim('results/mri/prepML.retest.meta.txt', sep='\t', header=F, quote="\"", stringsAsFactors=FALSE)

# save replication data and clean up
message('Saving retest data and cleaning up.')
save(gm, wm, covs, meta, file = 'results/mri/prepML.retest.RData')
system('rm -f results/mri/prepML.retest.*.txt')
message('Completed.')
