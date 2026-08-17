#!/bin/bash

## 1J3K chain A crystal ligand redocking
## Ligand: WRA, chain A
## Other HETATM in chain A: HOH and NDP are removed from receptor

## Download 1J3K
wget -N https://files.rcsb.org/download/1J3K.pdb.gz
gunzip -kf 1J3K.pdb.gz

## Extract crystal ligand WRA from chain A
awk 'substr($0,1,6)=="HETATM" && substr($0,18,3)=="WRA" && substr($0,22,1)=="A" {print}' \
1J3K.pdb > crystal_ligand_1J3K_A.pdb

echo "Crystal ligand atom count:"
wc -l crystal_ligand_1J3K_A.pdb
head crystal_ligand_1J3K_A.pdb

## Prepare receptor: chain A protein only
## This removes HOH, NDP and WRA because they are HETATM records
awk 'substr($0,1,4)=="ATOM" && substr($0,22,1)=="A" {print}' \
1J3K.pdb > 1J3K_chainA_clean.pdb

echo "END" >> 1J3K_chainA_clean.pdb

echo "Receptor atom count:"
wc -l 1J3K_chainA_clean.pdb
head 1J3K_chainA_clean.pdb

## Convert receptor to PDBQT
/home/mehnaz/app/mgltools_x86_64Linux2_1.5.7/bin/pythonsh \
/home/mehnaz/app/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py \
-r 1J3K_chainA_clean.pdb \
-o 1J3K_chainA_clean.pdbqt

## Prepare ligand: add hydrogens using OpenBabel
obabel crystal_ligand_1J3K_A.pdb -O ligand_1J3K_A.pdb -h

## Convert ligand to PDBQT
/home/mehnaz/app/mgltools_x86_64Linux2_1.5.7/bin/pythonsh \
/home/mehnaz/app/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py \
-l ligand_1J3K_A.pdb \
-o ligand_1J3K_A.pdbqt

## Run docking
vina --config 1J3K_config.txt --out 1J3K_docked.pdbqt --cpu 8

## Convert docking result to PDB
obabel 1J3K_docked.pdbqt -O 1J3K_docked.pdb

## Show docking score
head -30 1J3K_docked.pdbqt
