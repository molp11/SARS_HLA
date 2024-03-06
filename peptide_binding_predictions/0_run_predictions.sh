#!/bin/bash

echo 'Preparing UM inputs for NetMHCpan and NetMHCIIpan...'
python 1_prepare_prediction_inputs.py
mkdir -p input/UM_peptides/MHCI
for I in {8..11}; do cp input/UM_peptides/*_${I}mer.pep input/UM_peptides/MHCI; done
mkdir -p input/UM_peptides/MHCII
for I in {9..25}; do cp input/UM_peptides/*_${I}mer.pep input/UM_peptides/MHCII; done
rm input/UM_peptides/*.pep
mkdir -p output

echo 'Calculating UM NetMHCpan predictions...'
mkdir -p output/netMHCpan
declare -i N=$(ls -l input/UM_peptides/MHCI | grep -c '.pep')
for pep in input/UM_peptides/MHCI/*.pep; do 
    while read -r allele; do 
        ~/netMHCpan-4.1/netMHCpan -p $pep -BA -xls -a $allele -xlsfile output/netMHCpan/$(basename $pep .pep)"_"${allele}".xls" > tmp
        rm tmp
    done<input/MHCI_alleles.txt
    printf .
done | pv -pt -i0.2 -s$N -w 100 > /dev/null
python 2_combine_MHCI_output.py

echo 'Calculating UM NetMHCIIpan predictions...'
mkdir -p output/netMHCIIpan
declare -i N=$(ls -l input/UM_peptides/MHCII | grep -c '.pep')
for pep in input/UM_peptides/MHCII/*.pep; do
    while read -r allele; do
        ~/netMHCIIpan-4.0/netMHCIIpan -BA -f $pep -inptype 1 -a $allele -xls -xlsfile output/netMHCIIpan/$(basename $pep .pep)"_"${allele}".xls" > tmp
        rm tmp
    done<input/MHCII_alleles.txt
    printf .
done | pv -pt -i0.2 -s$N -w 100 > /dev/null
python 3_combine_MHCII_output.py

echo 'Realigning and annotating predictions...'
mkdir -p output/netMHCpan/xls
mv output/netMHCpan/*.xls output/netMHCpan/xls
mkdir -p output/netMHCIIpan/xls
mv output/netMHCIIpan/NSP1*.xls output/netMHCIIpan/xls
mv output/netMHCIIpan/NSP*.xls output/netMHCIIpan/xls
mv output/netMHCIIpan/NS*.xls output/netMHCIIpan/xls
mv output/netMHCIIpan/*.xls output/netMHCIIpan/xls
mkdir -p output/netMHCpan/pkl
mkdir -p output/netMHCIIpan/pkl
python 4_realign_and_annotate.py

echo 'Filtering binding interactions...'
python 5_filter_interactions.py

echo 'Preparing genbank inputs for NetMHCpan...'
mkdir -p ./input/genbank_peptides
python 6_prepare_genbank_peptides.py

echo 'Calculating genbank NetMHCpan predictions...'
mkdir -p output/netMHCpan_genbank
declare -i N=$(ls -l input/genbank_peptides | grep -c '.pep')
for pep in input/genbank_peptides/*.pep; do 
    while read -r allele; do 
        ~/netMHCpan-4.1/netMHCpan -p $pep -BA -xls -a $allele -xlsfile output/netMHCpan_genbank/$(basename $pep .pep)"_"${allele}".xls" > tmp
        rm tmp
    done<input/MHCI_alleles.txt
    printf .
done | pv -pt -i0.2 -s$N -w 100 > /dev/null
python 7_combine_genbank_output.py
mkdir -p output/netMHCpan_genbank/xls
mv output/netMHCpan_genbank/*.xls output/netMHCpan_genbank/xls

echo 'Done!'