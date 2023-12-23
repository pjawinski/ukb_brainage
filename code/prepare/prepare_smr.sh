#!/bin/bash

# =================================
# === prepare GTEx eQTL mapping ===
# =================================

# add smr software
parent="/home/groups/markett" # parent="/fast"
mkdir -p $parent/software/smr
cd $parent/software/smr
wget https://yanglab.westlake.edu.cn/software/smr/download/smr_Linux.zip
unzip *.zip
ln -s $parent/software/smr/smr_linux_x86_64 $parent/software/bin/smr

# create smr link in project folder
cd /home/groups/markett/ukb_brainage # cd /slow/projects/ukb_brainage
ln -s $parent/software/smr/ data/smr

# download westra
cd $parent/software/smr/
mkdir westra; cd westra
wget https://yanglab.westlake.edu.cn/data/SMR/westra_eqtl_hg19.zip
unzip westra_eqtl_hg19.zip
cd ..

# download CAGE
mkdir cage; cd cage
wget https://yanglab.westlake.edu.cn/data/SMR/cage_eqtl_data_lite_hg19.tar.gz
tar -xzvf cage_eqtl_data_lite_hg19.tar.gz
cd ..

# download BrainMeta v2 eqtl and sqtl
mkdir BrainMeta; cd BrainMeta
wget https://yanglab.westlake.edu.cn/data/SMR/BrainMeta_cis_eqtl_summary.tar.gz
tar -xzvf BrainMeta_cis_eqtl_summary.tar.gz

wget https://yanglab.westlake.edu.cn/data/SMR/BrainMeta_cis_sqtl_summary.tar.gz
tar -xzvf BrainMeta_cis_sqtl_summary.tar.gz
cd ..

# download GTEx
mkdir GTEx; cd GTEx
wget https://yanglab.westlake.edu.cn/data/SMR/GTEx_V8_cis_eqtl_summary_lite.tar
tar -xvf GTEx_V8_cis_eqtl_summary_lite.tar
awk 'NR==FNR {id=$1; next} $1==id { print id, $1, "md5 is equal"; next} { print id, $1, "md5 is NOT equal"; next} ' <(md5sum GTEx_V8_cis_eqtl_summary_lite/eQTL_besd_lite.zip) GTEx_V8_cis_eqtl_summary_lite/eQTL_besd_lite.zip.md5sum
unzip GTEx_V8_cis_eqtl_summary_lite/eQTL_besd_lite.zip
rm -rf GTEx_V8_cis_eqtl_summary_lite/

wget https://yanglab.westlake.edu.cn/data/SMR/GTEx_V8_cis_sqtl_summary_lite.tar
tar -xvf GTEx_V8_cis_sqtl_summary_lite.tar
awk 'NR==FNR {id=$1; next} $1==id { print id, $1, "md5 is equal"; next} { print id, $1, "md5 is NOT equal"; next} ' <(md5sum GTEx_V8_cis_sqtl_summary_lite/sQTL_besd_lite.zip) GTEx_V8_cis_sqtl_summary_lite/sQTL_besd_lite.zip.md5sum
unzip GTEx_V8_cis_sqtl_summary_lite/sQTL_besd_lite.zip
rm -rf GTEx_V8_cis_sqtl_summary_lite/

# meta-analyze gtex
cd eQTL_besd_lite
rm -f eqtl*
ls -l *.besd | awk '{ print $NF}' | sed 's/.besd//g' > eqtl.besd.list 
smr --besd-flist eqtl.besd.list --mecs --thread-num 100 --out eqtl.mecs

smr eQTL_besd_lite/


# change access rights
chmod -R 770 *

# create symbolic link in project folder
cd /home/groups/markett/ukb_brainage
ln -s /home/groups/markett/software/gtex data/gtex


