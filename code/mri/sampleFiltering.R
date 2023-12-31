#!/usr/bin/env Rscript

# ========================
# === Sample Filtering === 
# ========================
message('\n--- Sample Filtering ---')

# -------------------------------------------
# --- Filter data released until Feb 2020 ---
# -------------------------------------------
message('\nSelecting MRI subjects released in Feb 2020.')

# load required packages
for (pkg in c('dplyr','igraph','data.table')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# load basket data and CAT12 statistics of preprocessed T1 files
message(' - loading r2020 basket data.')
load('data/basket/20200217_2005558/data/ukb40487.RData')
message(sprintf(' - %d columns and %d rows.', dim(bd)[1],dim(bd)[2]))

# get list of usable T1 files
message(' - loading r2020 cat12 data.')
cat2020 = read.table("results/mri/cat12.r2020.txt", header = TRUE, sep = "\t")
cat2020 = data.frame(cat2020$IID[cat2020$IQR_poor==0])
names(cat2020) = "f.eid"
message(sprintf(' - %d cases', dim(cat2020)[1]))

# crop basket dataset to individuals with usable T1 files
message(' - intersect basket and cat12 data.')
bd2020 = left_join(x = cat2020, y = bd,  by = "f.eid")
message(sprintf(' - %d cases remaining.', dim(bd2020)[1]))

# apply filter criteria
# - reported and genetic sex mismatch
# - sex aneuploidy
# - outliers in heterozygosity and missingness rates
# - kinship information available
message(' - applying filter criteria.')
sex = as.numeric(bd2020$f.31.0.0)==as.numeric(bd2020$f.22001.0.0)
sex[is.na(sex)] = FALSE
noAneuploidy = is.na(as.numeric(bd2020$f.22019.0.0))
noHetMiss = is.na(as.numeric(bd2020$f.22027.0.0))
kinship = rep(TRUE, nrow(bd2020)); kinship[is.na(bd2020$f.22021)] = FALSE; kinship[as.numeric(bd2020$f.22021)==1] = FALSE
bd2020 = bd2020[sex & noAneuploidy & noHetMiss & kinship,]
message(sprintf(' - %d cases remaining.', dim(bd2020)[1]))

# define function to select unrelated individuals
select.unrelated <- function(kinshipTable, iid, seed, keep = NULL) {

  # set random number seed
  set.seed(seed)

  # shrink kinshipTable to pairs of individuals included in sample of interest
  kinshipTable <- kinshipTable[(kinshipTable$ID1 %in% iid) & (kinshipTable$ID2 %in% iid),]

  # identify dyadic pairs
  cat <- c(kinshipTable$ID1,kinshipTable$ID2)
  nodeSubjects <- cat[duplicated(cat)]
  dyadic = kinshipTable[!(kinshipTable$ID1 %in% nodeSubjects) & !(kinshipTable$ID2 %in% nodeSubjects), ]

  # remove one out of two subjects in dyadic pair
  exclude = sample(1:2, nrow(dyadic), replace=T)
  dyadicExclude = c()
  for (i in 1:nrow(dyadic)) {
    if (!is.null(keep) & sum(dyadic[i,c('ID1','ID2')] %in% keep) > 0) { dyadicExclude = c(dyadicExclude, dyadic[i,which(!(dyadic[i,c('ID1','ID2')] %in% keep))]) }
    else { dyadicExclude = c(dyadicExclude, dyadic[i,exclude[i]]) }
  }

  # identify subjects with multiple relationships
  multi = kinshipTable[kinshipTable$ID1 %in% nodeSubjects | kinshipTable$ID2 %in% nodeSubjects, ]
  multi$ID1 = factor(multi$ID1)
  multi$ID2 = factor(multi$ID2)
  adj = get.adjacency(graph.edgelist(as.matrix(multi[,c("ID1","ID2")]), directed = FALSE))

  # Identify trios
  adj = as.array(adj[,]) # create complete adjacency matrix (upper and lower triangle identical)
  for (i in 1:length(adj[,1])) {
    for (j in 1:length(adj[,1])) {
      if (adj[i,j] == 1) {
        adj[j,i] = 1
      }
    }
  }

  trios = c(NULL) # identify trios
  degrees = c(NULL)
  for (i in 1:nrow(adj)) {
    trio = c(row.names(adj)[i])
    degree = c(sum(adj[,i]))
    
    for (j in 1:nrow(adj)) {
      if (adj[j,i] == 1 && !(row.names(adj)[j] %in% trio) && !(row.names(adj)[j] %in% trios)) {
        trio = c(trio, row.names(adj)[j])
        degree = c(degree, sum(adj[,j]))
        
        for (k in 1:nrow(adj)) {
          if (adj[k,j] == 1 && !(row.names(adj)[k] %in% trio) && !(row.names(adj)[k] %in% trios)) {
            trio = c(trio, row.names(adj)[k])
            degree = c(degree, sum(adj[,k]))
            
            for (l in 1:nrow(adj)) {
              if (adj[l,k] == 1 && !(row.names(adj)[l] %in% trio) && !(row.names(adj)[l] %in% trios)) {
                trio = c(trio, row.names(adj)[l])
                degree = c(degree, sum(adj[,l]))
              }
            }
          }
        }
      }
    }
    if (length(trio) == 3) {
      trios = rbind(trios, trio)
      degrees = rbind(degrees, degree)
    }
  }

  # remove "hub" subject from type 1 trios (A related with B, B related with C) 
  triosType1 = matrix(trios[rowSums(degrees) == 4,], ncol = 3)
  degreesType1 = matrix(degrees[rowSums(degrees) == 4,], ncol = 3)
  triosType1Exclude = c()
  for (i in 1:nrow(triosType1)) {
    if (!is.null(keep) & sum(triosType1[i,] %in% keep) > 0 & identical(as.numeric(degreesType1[i,which(triosType1[i,] %in% keep)]),2)) { 
        triosType1Exclude = c(triosType1Exclude, triosType1[i,which(!(triosType1[i,] %in% keep))]) }
    else {
        triosType1Exclude = c(triosType1Exclude, triosType1[i, degreesType1[i,]==2]) }
  }

  # randomly remove 2 subjects from type 2 trios (a, b and c related to each other)
  triosType2 = matrix(trios[rowSums(degrees) == 6,], ncol = 3)
  select = sample(1:3, nrow(triosType2), replace=T)
  triosType2Exclude = c()
  for (i in 1:nrow(triosType2)) {
    if (!is.null(keep) & sum(triosType2[i,] %in% keep) > 0) { triosType2Exclude = c(triosType2Exclude, triosType2[i,which(!(triosType2[i] %in% keep))]) }
    else { triosType2Exclude = c(triosType2Exclude, triosType2[i, -select[i]]) }
  }

  # get list of remaining subjects
  kinshipExclude = c(dyadicExclude, triosType1Exclude, triosType2Exclude)
  iid = iid[!(iid %in% kinshipExclude)]

  # identify remaining kinships
  kinshipTable <- kinshipTable[(kinshipTable$ID1 %in% iid) & (kinshipTable$ID2 %in% iid),]
  cat = c(kinshipTable$ID1,kinshipTable$ID2)
  multi = kinshipTable[kinshipTable$ID1 %in% cat | kinshipTable$ID2 %in% cat, ]
  multi$ID1 = factor(multi$ID1)
  multi$ID2 = factor(multi$ID2)
  adj = get.adjacency(graph.edgelist(as.matrix(multi[,c("ID1","ID2")]), directed = FALSE))

  # iteratively remove subject with highest degree until no kinship remains
    # a) remove subjects related to 'keep' subjects
    if (!is.null(keep) & sum(row.names(adj) %in% keep) > 0) {
      idx = c()
      for (i in which(row.names(adj) %in% keep)) {
        idx = c(idx, as.numeric(which(adj[i,] == 1)))
      }
      idx = unique(idx)
      kinshipExclude = c(kinshipExclude,row.names(adj)[idx])
      adj = adj[-idx,-idx]
      adj = adj[which(row.names(adj) %in% keep),which(row.names(adj) %in% keep)]
    }
    # b) remove remaining relationships
    while (sum(adj) > 0) {
      idx = sample.int(nrow(adj),nrow(adj)) # randomly shuffle participants in case of ties
      adj = adj[idx,idx]
      subExclude = row.names(adj)[rowSums(as.array(adj)) %>% order(decreasing = T)][1]
      kinshipExclude = c(kinshipExclude,subExclude)
      idx = which(row.names(adj) %in% subExclude)
      adj = adj[-idx,-idx]
    }

  # return list of individuals to keep
  iid = iid[!(iid %in% kinshipExclude)]
  return(iid)
}

# apply function and get list of unrelated individuals
message(' - removing related individuals.')
kinshipTable = read.table('data/genetics/meta/ukb42032_rel_s488264.dat', head=TRUE)
unrelated = select.unrelated(kinshipTable,bd2020$f.eid,86609)
bd = bd[bd$f.eid %in% unrelated,]

# Sanity check: Test for two subjects in a row in relatedness table
if (sum(kinshipTable$ID1 %in% bd$f.eid & kinshipTable$ID2 %in% bd$f.eid) == 0) {
  message(sprintf(' - relatedness successfully removed\n - %d cases remaining.', dim(bd)[1]))
} else {
  message(' - relatedness has not been removed, exiting script.')
  quit()
}

# determine ancestry
caucasian = as.numeric(bd$f.22006)
caucasian[is.na(caucasian)] = 0

# creating output files
message(' - writing r2020 output files.')
save(bd, file = "data/basket/20200217_2005558/data/ukb40487_MRI.RData")
system('mkdir -p results/sample')
write.table(data.frame(FID = bd$f.eid, IID = bd$f.eid), file = 'results/sample/iid.r2020.unrelated.txt', quote = F, sep = '\t', row.names = F)
write.table(data.frame(FID = bd$f.eid[caucasian==1], IID = bd$f.eid[caucasian==1]), file = 'results/sample/iid.discovery.txt', quote = F, sep = '\t', row.names = F)


# ---------------------------------------------
# === add MRI subjects released in Feb 2021 ---
# ---------------------------------------------
message('\nAdding MRI subjects released in Feb 2021.')

# load feb 2020 and 2021 basket data
# load("data/basket/20200409_2007685/data/ukb41573_MRI.RData")
bd2020 = bd
message(' - loading r2021 basket data.')
load("data/basket/20210205_2007685/data/ukb45233.RData")
bd2021 = bd

# get list of usable T1 files
message(' - loading cat12 data.')
cat2020 = read.table("results/mri/cat12.r2020.txt", header = TRUE, sep = "\t")
cat2021 = read.table("results/mri/cat12.r2021.txt", header = TRUE, sep = "\t")
cat2021 = cat2021[!(cat2021$IID %in% cat2020$IID),]
cat2021 = data.frame(cat2021$IID[cat2021$IQR_poor==0])
names(cat2021) = "f.eid"

# bd2021: only keep individuals with usable T1 files
message(' - intersect basket and cat12 dataset.')
bd2021 = inner_join(x = cat2021, y = bd2021,  by = "f.eid")
message(sprintf(' - %d additional individuals with usable T1 files.', dim(bd2021)[1]))

# apply filter criteria
# - reported and genetic sex mismatch
# - sex aneuploidy
# - outliers in heterozygosity and missingness rates
# - kinship information available
message(' - applying filter criteria.')
sex = as.numeric(bd2021$f.31.0.0)==as.numeric(bd2021$f.22001.0.0)
sex[is.na(sex)] = FALSE
noAneuploidy = is.na(as.numeric(bd2021$f.22019.0.0))
noHetMiss = is.na(as.numeric(bd2021$f.22027.0.0))
kinship = rep(TRUE, nrow(bd2021)); kinship[is.na(bd2021$f.22021)] = FALSE; kinship[as.numeric(bd2021$f.22021)==1] = FALSE
bd2021 = bd2021[sex & noAneuploidy & noHetMiss & kinship,]
message(sprintf(' - %d cases remaining.', dim(bd2021)[1]))

# get list of unrelated individuals (keep individuals selected from 2020 release)
message(' - removing relatedness (taking into account relatedness with r2020 subjects).')
kinshipTable = read.table('data/genetics/meta/ukb42032_rel_s488264.dat', head=TRUE)
unrelated = select.unrelated(kinshipTable, iid = c(bd2020$f.eid,bd2021$f.eid), seed = 81262, keep = bd2020$f.eid)

# sanity check: test for two subjects in a row in relatedness table
if (sum(kinshipTable$ID1 %in% unrelated & kinshipTable$ID2 %in% unrelated) == 0) {
  message(sprintf(' - relatedness successfully removed.\n - %d additional cases remaining (%d in total).', sum(bd2021$f.eid %in% unrelated),length(unrelated)))
} else {
  message(' - relatedness has not been removed, exiting script.')
  quit()
}

# sanity check: make sure that all bd2020 individuals have been kept
if (sum(bd2020$f.eid %in% unrelated) == nrow(bd2020)) { 
  message(' - successfully kept all r2020 individuals in overall list of unrelated individuals.')
} else { 
  message(' - failed to keep all r2020 individualds in overall list of unrelated individuals.')
  quit()
}

# create datasets with subjects selected in 2020 and 2021
message(' - creating combined dataset.')
bd2021 = bd[bd$f.eid %in% unrelated,]

  # identify and add subjects missing in feb 2021 release
  # - 1 white-British ancestry individual 
  # bd2020[!(bd2020$f.eid %in% bd2021$f.eid), c("f.eid","f.22006.0.0")]
  bd2020crop = bd2020[!(bd2020$f.eid %in% bd2021$f.eid), names(bd2020) %in% names(bd2021)]
  bdTemp = suppressWarnings(data.table::rbindlist(list(bd2021, bd2020crop), fill = TRUE)) # Column 13685 of item 2 is an ordered factor but level 5 ['IC Scottish deaths (2017 onwards)'] is missing ... same for 13686 | numeric values did not change

    # sanity check
    # names(bd2020crop)[13685] # f.40018.0.0
    # names(bd2020crop)[13686] # f.40018.1.0
    # table(as.numeric(bd2021$f.40018.0.0)) 
    # table(as.numeric(bdTemp$f.40018.0.0)) # numeric values did not change
    # table(as.numeric(bd2021$f.40018.1.0)) 
    # table(as.numeric(bdTemp$f.40018.1.0)) # no values at all

  bd2021 = bdTemp[order(bdTemp$f.eid),]
  message(sprintf(' - %d rows and %d columns in combined dataset.', dim(bd2021)[1],dim(bd2021)[2]))

# save dataset
message(' - writing output files.')
bd = bd2021
save(bd, file = "data/basket/20210205_2007685/data/ukb45233_MRI.RData")
system('mkdir -p results/mri')
write.table(data.frame(FID = bd$f.eid, IID = bd$f.eid), file = 'results/mri/iid.r2021.unrelated.txt', quote = F, sep = '\t', row.names = F)
iid.discovery = read.table('results/mri/iid.discovery.txt', header = T)$IID
iid.replication = bd$f.eid[!(bd$f.eid %in% iid.discovery)]
write.table(data.frame(FID = iid.replication, IID = iid.replication), file = 'results/mri/iid.replication.txt', quote = F, sep = '\t', row.names = F)

# --------------------------------------------
# === add MRI retest data released in 2022 ---
# --------------------------------------------
message('\nAdding retest data released in 2022.')

# load r2021 and r2022 basket data
# load("data/basket/20210205_2007685/data/ukb45233_MRI.RData")
bd2021 = data.frame(bd)
message(' - loading r2022 basket data.')
load("data/basket/20220914_2016290/data/ukb669463.RData")
bd2022 = bd

# remove individualds not in MRI cohort from r2022
bd2022 = bd2022[bd2022$f.eid %in% bd2021$f.eid,]

  # identify and add subjects missing in r2022
  # - 11 individuals
  # bd2021[!(bd2021$f.eid %in% bd2022$f.eid), c("f.eid","f.22006.0.0")]
  bd2021crop = bd2021[!(bd2021$f.eid %in% bd2022$f.eid), names(bd2021) %in% names(bd2022)]
  bdTemp = suppressWarnings(data.table::rbindlist(list(bd2022, bd2021crop), fill = TRUE)) 

    # sanity check
    # names(bd2020crop)[13685] # f.40018.0.0
    # names(bd2020crop)[13686] # f.40018.1.0
    # table(as.numeric(bd2021$f.40018.0.0)) 
    # table(as.numeric(bdTemp$f.40018.0.0)) # numeric values did not change
    # table(as.numeric(bd2021$f.40018.1.0)) 
    # table(as.numeric(bdTemp$f.40018.1.0)) # no values at all

  bd2022 = bdTemp[order(bdTemp$f.eid),]
  message(sprintf(' - %d rows and %d columns in combined dataset.', dim(bd2022)[1],dim(bd2022)[2]))

# save dataset
message(' - writing output files.')
bd = data.frame(bd2022)
save(bd, file = "data/basket/20220914_2016290/data/ukb669463_MRI.RData")

# =============================
# === add pan ancestry data ===
# =============================
message('\nAdding pan ancestry data to basket mri data.')

# load dataset
# load("data/basket/20210205_2007685/data/ukb45233_MRI.RData")
pan = read.delim('data/basket/20210327_2008261/data/Files for retman/all_pops_non_eur_pruned_within_pop_pc_covs.tsv', sep = '\t', header = TRUE)
bridge = read.delim('data/basket/20210327_2008261/data/ukb42032bridge31063.txt', sep = ' ', header = FALSE)

# prepare pan-ancestry data
# remove duplicates with related flag 'true'
dupids = pan$s[duplicated(pan$s)]
pan = pan[!(pan$s %in% dupids & pan$related == 'true'),]
names(bridge) = c('f.eid', 's')
pan = inner_join(pan, bridge, 's')

# add pan to MRI and save dataset
message(' - writing output files.')
bd = left_join(bd, pan, 'f.eid')
save(bd, file = "data/basket/20220914_2016290/data/ukb669463_MRI_pan.RData")

# ==============================================================
# === Create variable file and ancestry-stratified iid files ===
# ==============================================================
message('Creating variable file and ancestry-stratified iid files.')

# load dataset
# load("data/basket/20220914_2016290/data/ukb669463_MRI_pan.RData")

# discovery or replication?
iid.discovery = read.table('results/mri/iid.discovery.txt', header = T)$IID
iid.replication = read.table('results/mri/iid.replication.txt', header = T)$IID
bd$discovery = NA
bd$discovery[bd$f.eid %in% iid.discovery] = 1
bd$discovery[bd$f.eid %in% iid.replication] = 0
message(sprintf(' - discovery n = %d | replication n = %d',sum(bd$discovery==1),sum(bd$discovery==0)))

# get caucasian ancestry and pan-ancestry population data
message(' - preparing output variables.')
caucasian = as.numeric(bd$f.22006)
caucasian[is.na(caucasian)] = 0
pan = bd$pop

# calculate exact age
birth_year_month = paste(as.numeric(bd$f.34.0.0),as.numeric(bd$f.52.0.0),"01", sep="-")
t1.age = as.numeric((bd$f.53.2.0 - as.Date(birth_year_month,"%Y-%m-%d"))/365.24219052)
t2.age = as.numeric((bd$f.53.3.0 - as.Date(birth_year_month,"%Y-%m-%d"))/365.24219052)

# prepare assessment center
ac.f = factor(bd$f.54.2.0)
ac.dummies = model.matrix(~ac.f)
ac = matrix(NA, nrow = length(ac.f), ncol = ncol(ac.dummies)-1)
ac[!is.na(ac.f),] = ac.dummies[,c(2:ncol(ac.dummies))]
ac = data.frame(ac)
names(ac) = paste0('t1.ac',1:ncol(ac))
t1.ac.dummy = ac

ac.f = factor(bd$f.54.3.0)
ac.dummies = model.matrix(~ac.f)
ac = matrix(NA, nrow = length(ac.f), ncol = ncol(ac.dummies)-1)
ac[!is.na(ac.f),] = ac.dummies[,c(2:ncol(ac.dummies))]
ac = data.frame(ac)
names(ac) = paste0('t2.ac',1:ncol(ac))
t2.ac.dummy = ac

# genotyping array dichotom
array = bd$f.22000.0.0
array[array>0] = 1
array[array<0] = 0

# get cat12 data
t1.cat12 = read.table('results/mri/cat12.r2021.txt', sep = '\t', header = T)
t2.cat12 = read.table('results/mri/cat12.r2022.retest.txt', sep = '\t', header = T)
t1.cat12 = t1.cat12[,-which(names(t1.cat12) == 'IQR_poor')]
t2.cat12 = t2.cat12[,-which(names(t2.cat12) == 'IQR_poor')]
t1.cat12 = t1.cat12[t1.cat12$IQR < 3,]
t2.cat12 = t2.cat12[t2.cat12$IQR < 3,]
names(t1.cat12) = c('IID',paste0('t1.',names(t1.cat12)[-1]))
names(t2.cat12) = c('IID',paste0('t2.',names(t2.cat12)[-1]))

# Covs data.frame
df = data.frame(
          FID = bd$f.eid,
          IID = bd$f.eid,
          discovery = bd$discovery,
          wba = caucasian,
          pan = pan,
          sex = as.numeric(bd$f.31.0.0),
          t1.age = t1.age,
          t1.age2 = t1.age^2,
          t1.ac.dummy,
          t1.x = bd$f.25756.2.0,
          t1.y = bd$f.25757.2.0, 
          t1.z = bd$f.25758.2.0,
          t1.rfMRI = bd$f.25741.2.0, 
          t1.tfMRI = bd$f.25742.2.0) %>%
      left_join(t1.cat12, by = 'IID') %>%
      cbind(data.frame(
          t2.age = t2.age,
          t2.age2 = t2.age^2,
          t2.ac.dummy,
          t2.x = bd$f.25756.3.0, 
          t2.y = bd$f.25757.3.0, 
          t2.z = bd$f.25758.3.0,
          t2.rfMRI = bd$f.25741.3.0, 
          t2.tfMRI = bd$f.25742.3.0)) %>%
      left_join(t2.cat12, by = 'IID') %>%
      cbind(data.frame(
          array = array,
          PanC1 = bd$PC1, PanC2 = bd$PC2, PanC3 = bd$PC3, PanC4 = bd$PC4, PanC5 = bd$PC5, 
          PanC6 = bd$PC6, PanC7 = bd$PC7, PanC8 = bd$PC8, PanC9 = bd$PC9, PanC10 = bd$PC10,
          PanC11 = bd$PC11, PanC12 = bd$PC12, PanC13 = bd$PC13, PanC14 = bd$PC14, PanC15 = bd$PC15,
          PanC16 = bd$PC16, PanC17 = bd$PC17, PanC18 = bd$PC18, PanC19 = bd$PC19, PanC20 = bd$PC20))

# check whether pan sex is the same as covs sex
# - one mismatch
# sum(df$sex[!is.na(bd$sex)]-1 == bd$sex[!is.na(bd$sex)])
# sum(df$sex[!is.na(bd$sex)]-1 != bd$sex[!is.na(bd$sex)])

# write subject IIDs and vars
message(' - writing output files.')
system('mkdir -p results/mri')
write.table(df, file = 'results/mri/r2022.vars.txt', sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(df[df$discovery == 0,] , file = 'results/mri/r2022.vars.replication.txt', sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

# write subject IIDs for each ancestry (replication sample only)
repl = df[df$discovery == 0,]
for (anc in names(table(repl$pan))) {
  tmp = data.frame(FID = repl$FID[repl$pan==anc & !is.na(repl$pan)], IID = repl$IID[repl$pan==anc & !is.na(repl$pan)])
  write.table(tmp, file = sprintf('results/mri/iid.replication.%s.txt',anc), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
}
anc = 'NA'
tmp = data.frame(FID = repl$FID[is.na(repl$pan)], IID = repl$IID[is.na(repl$pan)])
write.table(tmp, file = sprintf('results/mri/iid.replication.%s.txt',anc), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
message('\n--- Completed: Sample Filtering ---')

