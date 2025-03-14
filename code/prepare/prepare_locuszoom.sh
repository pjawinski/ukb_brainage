#!/bin/bash

# ============================================
# === prepare locuszoom for regional plots ===
# ============================================

# set working directory
cd /slow/projects/ukb_brainage

# create conda environment for locuszoom
conda create -p envs/locuszoom -c conda-forge python=2.7 r-base=4.1.1 r-gridExtra r-lattice r-cowplot r-dplyr r-ggplot2 r-magick r-patchwork r-pdftools
conda activate envs/locuszoom

# save python environment in yml file
conda env export --no-builds -p envs/locuszoom > envs/locuszoom.yml

# add name and remove prefix
awk 'NR==1 { print "name: locuszoom"; next } $1=="prefix:" { $2="locuszoom"; print; next} { print }' envs/locuszoom.yml > envs/locuszoom.yml.tmp; \mv envs/locuszoom.yml.tmp envs/locuszoom.yml

# install tabix
cd /fast/software
wget -O https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
chmod 770 tabix-0.2.6.tar.bz2
tar jxvf tabix-0.2.6.tar.bz2
cd tabix-0.2.6
make
ln -s /fast/software/tabix-0.2.6/tabix /fast/software/bin/tabix
ln -s /fast/software/tabix-0.2.6/bgzip /fast/software/bin/bgzip

# create chromosome .vcf files
cd /slow/projects/ukb_brainage
geneticsDir="data/genetics/"

for i in {1..21} X XY Y MT; do
	mkdir -p "${geneticsDir}/chr${i}/imp_mri_qc_EURjoined/vcf"
	plink --bfile "${geneticsDir}/chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc" \
	--recode vcf bgz tabx \
	--out "${geneticsDir}/chr${i}/imp_mri_qc_EURjoined/vcf/chr${i}_mri_qc"
chmod 770 "${geneticsDir}/chr${i}/imp_mri_qc_EURjoined/vcf/chr${i}_mri_qc"*
done

# create tabix file
task () {
echo "Starting with chr$i ..."
cd "${geneticsDir}/chr${i}/imp_mri_qc_EURjoined/vcf"
tabix -f "chr${i}_mri_qc.vcf.gz"
echo "Finished chr$i."
}

N=26; (
for i in {1..22} X XY Y MT; do
   ((j=j%N)); ((j++==0)) && wait
   task "$i" &
done
)

# download locuszoom and create symbolic links in bin folder
cd /fast/software
wget https://statgen.sph.umich.edu/locuszoom/download/locuszoom_1.4.tgz
tar -xvzf locuszoom_1.4.tgz
ln -s /fast/software/locuszoom/bin/locuszoom /fast/software/bin/locuszoom
ln -s /fast/software/locuszoom/bin/lzupdate.py /fast/software/bin/lzupdate.py
ln -s /fast/software/locuszoom/bin/locuszoom.R /fast/software/bin/locuszoom.R

# update locuszoom
mkdir /fast/software/locuszoom/data/database/deprecated
mv /fast/software/locuszoom/data/database/* /fast/software/locuszoom/data/database/deprecated/
cd /fast/software/locuszoom/

	# change RS_MERGE_ARCH_URL = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/database/data/organism_data/RsMergeArch.bcp.gz" 
	# into   RS_MERGE_ARCH_URL = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/database/organism_data/RsMergeArch.bcp.gz"
	cat bin/lzupdate.py | sed s%data/organism_data/RsMergeArch.bcp.gz%organism_data/RsMergeArch.bcp.gz%g > bin/lzupdate.py.tmp; \mv bin/lzupdate.py.tmp bin/lzupdate.py
	chmod 770 bin/lzupdate.py

	# before updating: get HapMap recombination rate
	mkdir -p data/build/hg19/recomb_rate/GRCh37
	cd data/build/hg19/recomb_rate/GRCh37
	wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
	tar -zxvf genetic_map_HapMapII_GRCh37.tar.gz
	cd ..
	echo chr$'\t'pos$'\t'recomb$'\t'cm_pos > recomb_rate.tab
	awk -F'\t' 'FNR==1 { next } { gsub(/chr/,""); gsub(/X_par1/,"25", $1); gsub(/X_par2/,"25", $1); gsub(/X/,"23", $1); print }' OFS='\t' GRCh37/genetic_map_GRCh37_chr*.txt | sort -k1,1g -k2,2g >> recomb_rate.tab
	
	# update locuszoom
	bin/lzupdate.py --build hg19 # if it fails, downloaded files may be corrupt; if so, delete them / download them manually and call 'bin/lzupdate.pj.py --build hg19'

	# recombination rate may also be added after updating
	bin/dbmeister.py --db db.backup/locuszoom_hg19.db --recomb_rate data/build/hg19/recomb_rate/recomb_rate.tab

# Change positions of XY variations from X coordinates to Y coordinates
cd /fast/software/locuszoom/data/database
sqlite3 -separator $'\t' locuszoom_hg19.db "select * from snp_pos;" > locuszoom_hg19_snp_pos.txt
cd ${geneticsDir}/chrXY/imp_mri_qc_EURjoined/vcf
mkdir original
cp chrXY* original/
awk 'NR==FNR && $2==24 { snp[$1]=$3; next } NR==FNR { next } FNR <=7 { print; next} $3 in snp { sub(23,24,$1); $2=snp[$3]; print $0; next } { print $0; next }' OFS="\t" /fast/software/locuszoom/data/database/locuszoom_hg19_snp_pos.txt <(gzip -dc chrXY_mri_qc.vcf.gz) > chrXY_mri_qc.vcf 
awk 'NR <= 7 { print }' OFS="\t" chrXY_mri_qc.vcf > chrXY_mri_qc.tmp.vcf
awk 'NR > 7 { print }' OFS="\t" chrXY_mri_qc.vcf | sort -k1,1 -k2,2n >> chrXY_mri_qc.tmp.vcf
\mv chrXY_mri_qc.tmp.vcf chrXY_mri_qc.vcf
rm -f chrXY_mri_qc.vcf.gz; bgzip chrXY_mri_qc.vcf
rm -f *.tbi; tabix -f "chrXY_mri_qc.vcf.gz"
chmod -R 770 *

