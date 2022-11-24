#!/bin/bash

inputFold=$1
outFold=$2
testdevgeno=0
devgeno=0.01
split_tot=$3


cd ${outFold}
mkdir -p devgeno${devgeno}_testdevgeno${testdevgeno}/

gene_info=${inputFold}resPrior_regEval_allchr.txt

# test dev geno
tmp="$(head -1 ${gene_info} | tr '\t' '\n' | cat -n | grep -w "test_dev_geno")"
echo ${tmp}
id_tdg=($tmp)
id_tdg=${id_tdg[0]}

# dev geno
tmp="$(head -1 ${gene_info} | tr '\t' '\n' | cat -n | grep -w "dev_geno")"
echo ${tmp}
id_dg=($tmp)
id_dg=${id_dg[0]}

# ensembl_gene_id
tmp="$(head -1 ${gene_info} | tr '\t' '\n' | cat -n | grep -w "ensembl_gene_id")"
echo ${tmp}
id_egi=($tmp)
id_egi=${id_egi[0]}

# total number of columns
n_col="$(awk '{print NF; exit}' ${gene_info})"
echo ${n_col}

awk -v id1="$id_egi" -v id2="$id_dg" -v id3="$id_tdg" '{print $id1, $id2, $id3}' ${gene_info} > geneList_tmp
# print gene names to remove
awk -v tdg=${testdevgeno} -v dg=${devgeno} '{ if(($2 == "NA") || ($3 == "NA") || ($3 < tdg) || ($2 < dg)) { print } }' geneList_tmp > rm_geneList

# first split, keep gene info
zcat split1_predictedExpression.txt.gz > tmp.txt
awk -v id1="$id_egi" 'NR==FNR{a[$1];next} !($id1 in a)' rm_geneList tmp.txt > devgeno${devgeno}_testdevgeno${testdevgeno}/split1_predictedExpression_filt.txt
rm tmp.txt
cp devgeno${devgeno}_testdevgeno${testdevgeno}/split1_predictedExpression_filt.txt devgeno${devgeno}_testdevgeno${testdevgeno}/predictedExpression_filt.txt

# from each split file remove the gene list
for i in $(seq 2 1 ${split_tot})
do
	
	zcat split${i}_predictedExpression.txt.gz > tmp.txt
	awk -v id1="$id_egi" 'NR==FNR{a[$1];next} !($id1 in a)' rm_geneList tmp.txt > devgeno${devgeno}_testdevgeno${testdevgeno}/split${i}_predictedExpression_filt.txt
	rm tmp.txt	
	cut --complement -f1-"$n_col" devgeno${devgeno}_testdevgeno${testdevgeno}/split${i}_predictedExpression_filt.txt > devgeno${devgeno}_testdevgeno${testdevgeno}/split${i}_predictedExpression_filt_nogene.txt
	paste -d "\t" devgeno${devgeno}_testdevgeno${testdevgeno}/predictedExpression_filt.txt  devgeno${devgeno}_testdevgeno${testdevgeno}/split${i}_predictedExpression_filt_nogene.txt > devgeno${devgeno}_testdevgeno${testdevgeno}/tmp_tot
	mv devgeno${devgeno}_testdevgeno${testdevgeno}/tmp_tot devgeno${devgeno}_testdevgeno${testdevgeno}/predictedExpression_filt.txt
	rm devgeno${devgeno}_testdevgeno${testdevgeno}/split${i}_predictedExpression_filt_nogene.txt

done

gzip devgeno${devgeno}_testdevgeno${testdevgeno}/predictedExpression_filt.txt


