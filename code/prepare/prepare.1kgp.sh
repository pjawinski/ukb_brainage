#!/bin/bash

# ===================================
# === get 1k genomes project data ===
# ===================================

# download phase3 v5a
mkdir -p /fast/software/1kgp/v5a && cd /fast/software/1kgp/v5a
for i in {1..22}; do (
	echo starting with chromosome ${i}
	wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi ) & 
done
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1a.20130502.genotypes.vcf.gz &
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1a.20130502.genotypes.vcf.gz.tbi &
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz &
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz.tbi
wait
chmod 550 ALL*
rm -f wget-log*

# keep meta data only (to speed up 1kgp allele + frequency matching)
rm -f *meta*
FILES=$(ls *.vcf.gz)
N=24
time (
for FILE in ${FILES}; do
	((i=i%N)); ((i++==0)) && wait
	(echo "Extracting meta data of $FILE"
	outfile=$(echo $FILE | sed 's/vcf.gz/meta.vcf/g')
	zcat $FILE | awk '{ print $1, $2, $3, $4, $5, $6, $7, $8, $9 }' > ${outfile}
	gzip -f ${outfile}) &
done
wait
)


wget $pgen
mv 'all_phase3.pgen.zst?dl=1' all_phase3.pgen.zst
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen

wget $pvar
mv 'all_phase3.pvar.zst?dl=1' all_phase3.pvar.zst
plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar

wget $sample
mv 'phase3_corrected.psam?dl=1' all_phase3.psam

# create list of duplicate snps - takes about 6 minutes with 12 cpus
plink2 --pfile $source/all_phase3 vzs \
	--write-snplist \
	--out $source/snps

gsort $source/snps.snplist --parallel=24 | uniq -d > $source/snps_duplicate.txt

# Convert 1000 Genomes phase 3 data to plink 1 binary format
# keep European subjects
# keep biallelic snps
# exclude duplicate rs ids

mkdir -p bed_EUR_nodups
plink2 --pfile $source/all_phase3 vzs \
	--max-alleles 2 \
	--exclude $source/snps_duplicate.txt \
	--keep ${parent}/04_EUR_individuals.txt \
	--make-bed \
	--out $source/bed_EUR_nodups/all_phase3_EUR_nodups

# Convert 1000 Genomes phase 3 data to plink 1 binary format
# keep biallelic snps
# exclude duplicate rs ids

mkdir -p bed_nodups
plink2 --pfile $source/all_phase3 vzs \
	--max-alleles 2 \
	--exclude $source/snps_duplicate.txt \
	--make-bed \
	--out $source/bed_nodups/all_phase3_nodups



# download phase3 v5b
mkdir -p /ifs/loni/faculty/dhibar/EEG/pwr/data/1kgp/v5b
cd /ifs/loni/faculty/dhibar/EEG/pwr/data/1kgp/v5b
qsub -cwd -N downl <<eof
	#!/bin/bash
	for i in {1..22}; do (
	echo starting with chromosome \${i}
	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr\${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr\${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi ) & 
	done
	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz &
	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz.tbi &
	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz &
	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz.tbi &
	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz &
	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz.tbi
	wait
eof
chmod 550 ALL*
rm -f downl*
