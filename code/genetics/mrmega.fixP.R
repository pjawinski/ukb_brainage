#Manhattan plot script for MR-MEGA
#Written by Joshua C Randall & Reedik Magi
for (e in commandArgs(trailingOnly=TRUE))
{
  ta = strsplit(e,"=",fixed=TRUE)
  if(!is.null(ta[[1]][2]))
  {
    assign(ta[[1]][1],ta[[1]][2])
  } else {
    assign(ta[[1]][1],TRUE)
  }
}
if(!exists("input"))
{
  input <- paste("mrmega.result")
}
if(!exists("out")) {
  out <- paste(input,".fixedP",sep="")
}
data<-data.frame(data.table::fread(input,stringsAsFactors=FALSE,header=TRUE,sep = "\t")) # PJ: added data.table::fread for faster reading
data$P.value_association <- pchisq(data$chisq_association, data$ndf_association, lower.tail=FALSE)
data$P.value_ancestry_het <- pchisq(data$chisq_ancestry_het, data$ndf_ancestry_het, lower.tail=FALSE)
data$P.value_residual_het <- pchisq(data$chisq_residual_het, data$ndf_residual_het, lower.tail=FALSE)
data.table::setnames(data, gsub("\\.", "-", names(data))) # PJ: Restore original column names by replacing dots with dashes
data.table::fwrite(data, file=out, row.names=F, quote=F, sep = "\t", na = "NA", compress='none'); # PJ: added data.table::fwrite for faster writing