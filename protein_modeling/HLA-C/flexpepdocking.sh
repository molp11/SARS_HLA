#!/bin/bash
BIN="/mnt/f/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin"
DB="/mnt/f/rosetta_bin_linux_2021.16.61629_bundle/main/database"
mkdir -p ../top_models
for DIR in C*
do
  MODEL=$(echo $DIR | cut -f1-3 -d'_')
  PDB=${DIR}_unrelaxed_rank_1_model_1.pdb
  echo 'Starting refinement for '${MODEL}'...'
  cd $DIR
  $BIN/FlexPepDocking.static.linuxgccrelease -database $DB -s $PDB -flexpep_prepack -ex1 -ex2aro
  $BIN/FlexPepDocking.static.linuxgccrelease -database $DB -s $(basename $PDB .pdb)_0001.pdb -out:file:silent decoys_lowres.silent -out:file:silent_struct_type binary -pep_refine -ex1 -ex2aro -use_input_sc -nstruct 100 -lowres_preoptimize
  TOP_LOWRES=$(cat decoys_lowres.silent | grep SCORE | awk '{print $2"\t"$NF}' | grep -v description | sort -n -k1 | head -1 | awk '{print $2}')
  $BIN/extract_pdbs.static.linuxgccrelease -in:file:silent decoys_lowres.silent -in:file:tags $TOP_LOWRES
  $BIN/FlexPepDocking.static.linuxgccrelease -database $DB -s ${TOP_LOWRES}.pdb -out:file:silent decoys_highres.silent -out:file:silent_struct_type binary -pep_refine -ex1 -ex2aro -use_input_sc -nstruct 100
  TOP_HIGHRES=$(cat decoys_highres.silent | grep SCORE | awk '{print $2"\t"$NF}' | grep -v description | sort -n -k1 | head -1 | awk '{print $2}')
  $BIN/extract_pdbs.static.linuxgccrelease -in:file:silent decoys_highres.silent -in:file:tags $TOP_HIGHRES
  #cp -v ${TOP_HIGHRES}.pdb ../../top_models
  cd ../
  echo 'Refinement for '${MODEL}' finished!'
done
