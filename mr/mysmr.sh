#!/bin/bash

mysmr="/Users/xwang/Git/merck/mr/smr_Mac"

mydata="/Users/xwang/Git/merck/gcta/test"
mygwas="/Users/xwang/Git/merck/mr/IGAP_stage_1.smr"
myeqtl="/Users/xwang/Git/merck/mr/westra_eqtl_hg19/westra_eqtl_hg19"

$mysmr --bfile $mydata --gwas-summary $mygwas --beqtl-summary $myeqtl --out mysmr --thread-num 3

