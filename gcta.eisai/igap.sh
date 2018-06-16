#!/bin/sh

# CD33: 19:51728320-51747115
# CD33 variants
# flanking: +/-50K

igap="/home/xulong/HBDS/IGAP/IGAP_stage_1.txt"

awk -F '\t' '$1 == 19 { print $0 }' $igap | \
awk -F '\t' '$2 > 51678320 { print $0 }' | \
awk -F '\t' '$2 < 51797115 { print $0 }' > ./IGAP_stage_1.txt

# awk -F '\t' '$11 < 1e-5 { print $0 }' ./Suhre2017/3166-92_1_one.out

plink="/home/xulong/HBDS/xulong/1000G/plink"
./gcta64 --bfile $plink --cojo-file input.ma --out test0
./gcta64 --bfile $plink --cojo-file input.ma --cojo-cond cond.snplist --out test1

