# This script downloads a few BED files for the experiment


# Given a bed file hash, download the file from bedbase and unzip it into the
# /data folder
download_from_bedbase() {
  fileid=$1
  wget http://bedbase.org/api/bed/$fileid/file/bedfile -O data/$fileid.bed.gz
  gunzip data/$fileid.bed.gz > temp
  cut -f1-3 data/$fileid.bed > data/$fileid_cut.bed
  mv data/$fileid_cut.bed data/$fileid.bed
}

download_from_bedbase "713f58a6497a9168a326123919672ebe"
download_from_bedbase "0f84fea95b736ec99914bc66e74ab6e0"
download_from_bedbase "c75ea5133f825d779a02be41a529342e"

# Download screen universes:
wget https://api.wenglab.org/screen_v13/fdownloads/GRCh38-ccREs.bed -O data/GRCh38-ccREs.bed
wget https://api.wenglab.org/screen_v13/fdownloads/GRCh38-ccREs.CTCF-only.bed  -O data/GRCh38-ccREs.CTCF-only.bed
wget https://api.wenglab.org/screen_v13/fdownloads/GRCh38-ccREs.PLS.bed  -O data/GRCh38-ccREs.PLS.bed
wget https://api.wenglab.org/screen_v13/fdownloads/GRCh38-ccREs.DNase-H3K4me3.bed  -O data/GRCh38-ccREs.DNase-H3K4me3.bed
