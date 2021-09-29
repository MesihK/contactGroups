#!/bin/bash
RED='\033[0;31m'
NC='\033[0m' 

echo 'This script compares given PDB file with MSA made DCA calculations'

if [ $# -lt 3 ]; then
    echo 'usage: eval.sh pfam.file prot.pdb subs.mat gapopen gapextention'
    echo '   ex: eval.sh PF18-3.pfam 1shg.pdb protsub.mat 10 1'
    exit 0
fi

PFAM=$1
PROT=$2
MAT=$3
echo -e ${RED}Matrix: $MAT ${NC}

if [ -z $4 ]; then
    GAPOPEN=10
else
    GAPOPEN=$4
    echo -e ${RED}Gap open: $GAPOPEN ${NC}
fi

if [ -z $5 ]; then
    GAPEXT=1
else
    GAPEXT=$5
    echo -e ${RED}Gap extend: $GAPEXT ${NC}
fi
SIM=0.7
DIR=$PFAM.$PROT.$MAT.$GAPOPEN.$GAPEXT

if [ -d $DIR ]; then
    echo -e ${RED}The result directory: $DIR exists.$NC
    exit 0
fi

# remove gaps from MSA file
if [ ! -f $PFAM.ng ]; then 
    echo -e ${RED}sed remove gaps: ${NC}
    sed 's/\.//g' $PFAM > $PFAM.ng
fi

# get protein fasta file
if [ ! -s $PROT.fa ]; then
    echo -e ${RED}python utils_protein2.py writeseqfa $PROT: ${NC}
    python utils_protein2.py writeseqfa $PROT
fi

# paste the protein fasta to the end of the MSA file
if [ ! -s $PFAM.ng.up ]; then
    echo -e '\n' >> $PFAM.ng
    cat $PROT.fa >> $PFAM.ng
    echo -e ${RED}sed uppercase $PFAM.ng.up: ${NC}
    sed -e 's/\(.*\)/\U\1/' $PFAM.ng > $PFAM.ng.up
fi

# Get the name of protein from fasta file, and convert it to uppercase
NAME=$(cat $PROT.fa | head -n1 | cut -c 2- | sed -e 's/\(.*\)/\U\1/')


# Recalculate MSA with muscle
if [ ! -s $PFAM.$MAT.msa ]; then
    
    echo -e ${RED}mafft --op $GAPOPEN --ep $GAPEXT --textmatrix $MAT --thread 8 --quiet $PFAM.ng.up '>' $PFAM.$MAT.msa ${NC}
    #mafft --op $GAPOPEN --ep $GAPEXT --textmatrix $MAT --thread 8 --quiet $PFAM.ng.up > $PFAM.$MAT.msa
    mafft --op $GAPOPEN --ep $GAPEXT --textmatrix $MAT --thread 8 $PFAM.ng.up > $PFAM.$MAT.msa
    if [ $? -ne 0 ]; then
        echo -e ${RED} MSA failed $NC
        exit 1
    fi
else
    echo -e ${RED}Previous MSA found, use it.$NC
fi

# Get our protein alignment from MSA
echo -e ${RED}python utils_pfammsa.py getsinglemsa $PFAM.$MAT.msa $NAME $PFAM.$PROT: ${NC}
python utils_pfammsa.py getsinglemsa $PFAM.$MAT.msa $NAME $PFAM.$PROT

# Calculate map file between the PDB and aligned sequence 
echo -e ${RED}python utils_resimap.py posmapvec4 $PROT.fa $PFAM.${PROT}_MSA.fa $PROT $PROT.map: ${NC}
python utils_resimap.py posmapvec4 $PROT.fa $PFAM.${PROT}_MSA.fa $PROT $PROT.map

# Calculate score by reducing msa
echo -e ${RED}python utils_pfammsa.py msareduce_withmap $PFAM.$MAT.msa $PROT.map aa $PFAM.$PROT.aa: ${NC}
python utils_pfammsa.py msareduce_withmap $PFAM.$MAT.msa $PROT.map aa $PFAM.$PROT.aa

# Calculate hammingweight
echo -e ${RED}python proc_hammingweight.py $PFAM.$PROT.aa.scoremat $SIM '>' $PROT.weight: ${NC}
python proc_hammingweight.py $PFAM.$PROT.aa.scoremat $SIM > $PROT.hamm.weight

# Convert the weight from [1, ... 3] to 1 ... 3
weight_convert () {
    sed 's/\[//g' $1 > ${1}.1
    sed 's/\]//g' ${1}.1 > ${1}.2
    sed 's/,//g' ${1}.2 > $1
    rm ${1}.1 ${1}.2
}

weight_convert $PROT.hamm.weight

# Calculate DCA
echo -e ${RED}python utils_dca.py foo $PFAM.$PROT.aa.scoremat $PROT.hamm.weight $PFAM.pydca: ${NC}
python utils_dca.py dca $PFAM.$PROT.aa.scoremat $PFAM.$PROT.aa.rcol $PROT.hamm.weight $PFAM.pydca

# Calculate residue distance using PDB file
echo -e ${RED}python utils_protein2.py writeresdists $PROT $PROT.dist : ${NC}
python utils_protein2.py writeresdists $PROT $PROT.dist 

# Compare DCA and DIST and calculate RMSD 
echo -e ${RED}python utils_mesih.py evaldistdca $PROT.map $PROT.dist $PFAM.pydca: ${NC}
python utils_mesih.py evaldistdca $PROT.map $PROT.dist $PFAM.pydca 5

mv roc_curve.png hamm.png

# Calculate dummyweight
echo -e ${RED}python proc_hammingweight.py $PFAM.$PROT.aa.scoremat $SIM '>' $PROT.weight: ${NC}
python utils_mesih.py dummyweight $PFAM.$PROT.aa.scoremat > $PROT.dummy.weight

weight_convert $PROT.dummy.weight

# Calculate DCA
echo -e ${RED}python utils_dca.py foo $PFAM.$PROT.aa.scoremat $PROT.dummy.weight $PFAM.pydca: ${NC}
python utils_dca.py dca $PFAM.$PROT.aa.scoremat $PFAM.$PROT.aa.rcol $PROT.dummy.weight $PFAM.pydca

# Compare DCA and DIST and calculate RMSD 
echo -e ${RED}python utils_mesih.py evaldistdca $PROT.map $PROT.dist $PFAM.pydca: ${NC}
python utils_mesih.py evaldistdca $PROT.map $PROT.dist $PFAM.pydca 5

mv roc_curve.png dummy.png

# Calculate clusterweight
echo -e ${RED}python proc_hammingweight.py $PFAM.$PROT.aa.scoremat $SIM '>' $PROT.weight: ${NC}
python utils_mesih.py clusterweight $PFAM.$PROT.aa.scoremat > $PROT.clust.weight

weight_convert $PROT.clust.weight

# Calculate DCA
echo -e ${RED}python utils_dca.py foo $PFAM.$PROT.aa.scoremat $PROT.clust.weight $PFAM.pydca: ${NC}
python utils_dca.py dca $PFAM.$PROT.aa.scoremat $PFAM.$PROT.aa.rcol $PROT.clust.weight $PFAM.pydca

# Compare DCA and DIST and calculate RMSD 
echo -e ${RED}python utils_mesih.py evaldistdca $PROT.map $PROT.dist $PFAM.pydca: ${NC}
python utils_mesih.py evaldistdca $PROT.map $PROT.dist $PFAM.pydca 5

mv roc_curve.png clust.png

mkdir $DIR

#mv roc_curve.png $DIR/$PROT.$MAT.$CLUST.roc.png
#mv $PROT.fa $PFAM.$MAT.msa $PFAM.${PROT}_* $PFAM.$PROT.aa* $PROT.weight $PROT.map $PROT.dist $PFAM.pydca $DIR
mv $PROT.fa $PFAM.$MAT.msa $PFAM.${PROT}_* $PFAM.$PROT.aa* *.weight $PROT.map $PFAM.pydca *.png $DIR
echo ./eval.sh $@ > $DIR/cmd

