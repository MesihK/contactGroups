#!/bin/bash
RED='\033[0;31m'
NC='\033[0m' 

echo 'This script compares given PDB file with MSA made DCA calculations'
echo ''

if [ $# -ne 4 ]; then
    echo 'usage: eval.sh pfam.file prot.pdb subs.mat clustThreshold'
    echo '   ex: eval.sh PF18-3.pfam 1shg.pdb protsub.mat 0.70'
    exit 0
fi

PFAM=$1
PROT=$2
MAT=$3
CLUST=$4
SIM=0.7

# remove gaps from MSA file
echo -e ${RED}sed remove gaps: ${NC}
sed 's/\.//g' $PFAM > $PFAM.ng

# get protein fasta file
echo -e ${RED}python utils_protein2.py writeseqfa $PROT: ${NC}
python utils_protein2.py writeseqfa $PROT

# paste the protein fasta to the end of the MSA file
cat $PROT.fa >> $PFAM.ng
echo -e ${RED}sed uppercase $PFAM.ng.up: ${NC}
sed -e 's/\(.*\)/\U\1/' $PFAM.ng > $PFAM.ng.up

# Get the name of protein from fasta file, and convert it to uppercase
NAME=$(cat $PROT.fa | head -n1 | cut -c 2- | sed -e 's/\(.*\)/\U\1/')


# Recalculate MSA with muscle
if [ ! -f $PFAM.$MAT.msa ]; then
    echo -e ${RED}muscle -in $PFAM.ng.up -out $PFAM.$MAT.msa -matrix protsub.mat -gapopen -12.0 -gapextend -1.0 -center 0.0: ${NC}
    muscle -in $PFAM.ng.up -out $PFAM.$MAT.msa -matrix $MAT -gapopen -6.0 -gapextend -0.5 -center 0.0
    #muscle -in $PFAM.ng.up -out $PFAM.ng.up.bl -matrix blosum62.mat -gapopen -12.0 -gapextend -1.0 -center 0.0
else
    echo Previous Muscle MSA found, use it.
fi

# Get our protein alignment from MSA
echo -e ${RED}python utils_pfammsa.py getsinglemsa $PFAM.$MAT.msa $NAME $PFAM.$PROT: ${NC}
python utils_pfammsa.py getsinglemsa $PFAM.$MAT.msa $NAME $PFAM.$PROT.pb

# Calculate map file between the PDB and aligned sequence 
echo -e ${RED}python utils_resimap.py posmapvec4 $PROT.fa $PFAM.${PROT}.pb_MSA.fa $PROT $PROT.pb.map: ${NC}
python utils_resimap.py posmapvec4 $PROT.fa $PFAM.${PROT}.pb_MSA.fa $PROT $PROT.pb.map

# Cluster the result of MSA with given threshold
echo -e ${RED}python utils_pfammsa.py getsinglemsacluster $PFAM.$MAT.msa $PROT.pb.map $NAME aa 0.7 $PFAM.$PROT: ${NC}
python utils_pfammsa.py getsinglemsacluster $PFAM.$MAT.msa $PROT.pb.map $NAME aa $CLUST $PFAM.$PROT

#echo -e ${RED}python utils_pfammsa.py msareduce_withmap $PFAM.$MAT.msa $PROT.pb.map aa $PFAM.$PROT.aa: ${NC}
#python utils_pfammsa.py msareduce_withmap $PFAM.$MAT.msa $PROT.pb.map aa $PFAM.$PROT.aa

# Calculate score by reducing msa
echo -e ${RED}python utils_pfammsa.py msareduce_withmap $PFAM.${PROT}_${CLUST}_cluster.fa $PROT.pb.map aa $PFAM.$PROT.aa: ${NC}
python utils_pfammsa.py msareduce_withmap $PFAM.${PROT}_${CLUST}_cluster.fa $PROT.pb.map aa $PFAM.$PROT.aa

# Calculate hammingweight
echo -e ${RED}python proc_hammingweight.py $PFAM.$PROT.aa.scoremat $SIM '>' $PROT.weight: ${NC}
python proc_hammingweight.py $PFAM.$PROT.aa.scoremat $SIM > $PROT.weight

# Convert the weight from [1, ... 3] to 1 ... 3
sed 's/\[//g' $PROT.weight > $PROT.weight1
sed 's/\]//g' $PROT.weight1 > $PROT.weight2
sed 's/,//g' $PROT.weight2 > $PROT.weight
rm $PROT.weight1 $PROT.weight2

# Calculate DCA
echo -e ${RED}python utils_dca.py foo $PFAM.$PROT.aa.scoremat $PROT.weight $PFAM.pydca: ${NC}
python utils_dca.py dca $PFAM.$PROT.aa.scoremat $PFAM.$PROT.aa.rcol $PROT.weight $PFAM.pydca

# Calculate residue distance using PDB file
echo -e ${RED}python utils_protein2.py writeresdists $PROT $PROT.dist : ${NC}
python utils_protein2.py writeresdists $PROT $PROT.dist 

# Compare DCA and DIST and calculate RMSD 
echo -e ${RED}python utils_mesih.py evaldistdca $PROT.pb.map $PROT.dist $PFAM.pydca: ${NC}
python utils_mesih.py evaldistdca $PROT.pb.map $PROT.dist $PFAM.pydca 5

mv roc_curve.png $PROT.$MAT.roc.png

#rm $PFAM.ng $PFAM.ng.up $PROT.fa $PFAM.$MAT.msa $PFAM.${PROT}* $PROT.weight # $PROT.map 

