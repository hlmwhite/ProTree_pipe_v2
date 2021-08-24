#!/bin/bash

# ProTree is a simple automoated pipeline for generating phylogenetic trees based on protein data sets.
#	- It should suit quick and simple phylogenetic inferences without much input, however if there are any modifications to the commands required,
#	- you may do so within this script at their respective command.

# the program is run within a conda environment, with most tools installed through conda. Some other scripts are installed through github:

# conda install orthofinder
# conda install -c bioconda gblocks
# conda install -c bioconda mafft
# faSomeRecords  -  https://github.com/santiagosnchez/faSomeRecords
# fasta_to_phylip.py  -  https://github.com/audy/bioinformatics-hacks/blob/master/bin/fasta-to-phylip
# prottest3 - https://github.com/ddarriba/prottest3
#	- prottest3 also requires java!
# conda install -c genomedk raxml-ng 
# conda install -c bioconda fasttree

# 0. activate the conda environment and install the necessary software

# e.g   source activate env

# 1. copy all protein sets to compare into a single directory 
#	- it is helpful to rename each with an organism/specific identifier to rename the proteins within the file. 
#	- The name must also be short, where Gblocks will complain of protein fasta headers being too long and disregarding the alilgnment.

#	- ALL FASTAS MUST END IN " .faa ".

# execute via : ' ProTree_pipe.sh <threads to use> <option - "ft" > ' 

#  using fasttree removes the use of protest3 and raxml, is useful for quick trees, or for focusing on 4 or less organisms
#  this requires the flag 'ft' at the end, but leave blank to use prottest and raxml as a default

if [ "$#" == 0 ]; then
	echo ''
	echo 'usage: ProTree_pipe.sh <threads to use> <raxml threads> <option - "ft" >'
	echo '	- run ProTree in current working directory with protein fastas to use.'
	echo '	- all protein fastas must end in ".faa"!'
	echo ''
	echo 'install:'
	echo '	conda install orthofinder'
        echo '  conda install -c bioconda gblocks'
        echo '  conda install -c bioconda mafft'
        echo '  faSomeRecords  -  https://github.com/santiagosnchez/faSomeRecords'
        echo '  fasta_to_phylip.py  -  https://github.com/audy/bioinformatics-hacks/blob/master/bin/fasta-to-phylip'
        echo '  prottest3 - https://github.com/ddarriba/prottest3'
        echo '        - prottest3 also requires java!'
        echo '  conda install -c genomedk raxml-ng '
        echo '  conda install -c bioconda fasttree'

	exit
fi

threads=$1
tree=$3
rax_threads=$2

echo ''
echo 'running ProTree_pipe...'
echo ''
echo '   WARNING'
echo ''
echo 'ensure command is the following format: (10 seconds to cancel)'
echo ''
echo 'ProTree_pipe.sh <threads to use> <rax_threads> <option - "ft">'

sleep 5
echo '(5 seconds to cancel)'
sleep 5

echo 'renaming fastas...'
sleep 2
#mkdir work_dir

#for file in *faa ; do
#	name=${file%.faa}
#	awk -v var="$name" '/^>/{print ">"var "_" ++i "_"; next}{print}' "$file" > work_dir/"$name"_mod.faa
#done

echo 'generating header files...'
sleep 2
cd work_dir

#for file in *faa ; do 
#org=${file%_mod.faa}
#grep $org $file > $org.headers
#done

# 2 run orthofinder to infer single copy orthogroups

echo 'running orthofinder...'
sleep 2
#orthofinder -f . -t $threads -og

cd OrthoFinder/Results_*/
#cd Results*

#grep -f SingleCopyOrthogroups.txt Orthogroups.csv > single_orthos.csv

#while IFS= read -r pattern ; do
#	printf "$pattern" > OG.txt
#	OG=$( awk '{print $1}' OG.txt )
#	cat OG.txt | sed 's/\t/\n/g' | grep -v 'OG' > $OG.list
#done < single_orthos.csv

mkdir single_OGs
cd single_OGs/
cp ../Single_Copy_Orthologue_Sequences/*fa . 

#cp ../OrthoFinder/Results_*/
#mv ../OG*.list .

#cat ../../*faa > all.faa

#for file in *.list ; do
#	name=${file%.list}
#	faSomeRecords.py --fasta all.faa --list $file --outfile $name.faa
#done

# 3. align sequences with mafft and trim erroneous alignments using Gblocks

echo 'running mafft...'
sleep 2
for file in OG*.fa ;
do
num=$(ls -hal OG*fa | wc -l)
  mafft --auto $file > $file.mafft.out 2> $file.mafft.err
	count=$(ls -hal *mafft.out | wc -l)
	echo $count'/'$num
done 

final_num=$(ls -hal OG*fa | wc -l)

#rename "s/faa.mafft.out/fasta/" *.faa.mafft.out
rename "s/fa.mafft.out/fasta/" *.fa.mafft.out

echo 'running Gblocks...'
sleep 2
for file in *.fasta ;
do
	Gblocks $file -t=p -p=y
done

rename "s/fasta-gb/fas/" *-gb

# 4. concatenating alignments

echo 'concatenating alignments...'
sleep 2
mkdir gblocks_fas
cd gblocks_fas

cat ../*fas > all.fas

#mv ../../../*headers .
mv ../../../../*headers .

for file in *headers ; do 
	name=${file%.headers}
	faSomeRecords.py --fasta all.fas --list $file --outfile $name.faa
done 

for file in *.faa ; do
	name=${file%.faa}
	grep -v '>' $file > 1.file
	printf ">$name\n" > fas_head
	cat fas_head 1.file > "$name"_all.fa
done

cat *_all.fa > ../../../alignment.fasta

cd ../../../

# 5. converting alignment fasta to phylip format

echo 'running fasta_to_phylip...'

fasta_to_phylip.py --input-fasta alignment.fasta --output-phy alignment.phy

# 6. generate tree with fasttree

sleep 2

if [ $tree == "ft" ]
then
        echo 'running fasttree...'

	mkdir fasttree_out
	cd fasttree_out

	model="N/A"
	fasttree ../alignment.fasta > fasttree.nwk

	tree_file="fasttree.nwk"
else 
# 6. choosing the best model for tree inference with prottest3 - this is the default!

echo 'running prottest3...'

sleep 2

java -Xmx32g -jar /pub28/markw/bin/prottest3/prottest-3.4.2/prottest-3.4.2.jar -i alignment.phy -o prottest3.out -all-distributions -F -BIC -tc 0.5 -threads $threads

model=$(grep Best prottest3.out | awk '{print $6}' | head -n 1 )


# 7. phylogenetic tree inference with raxml-ng

echo 'running raxml-ng...'

mkdir raxml-ng-out
cd raxml-ng-out 

raxml-ng --all --msa ../alignment.phy --model $model --prefix alignment_raxml --threads $rax_threads --tree pars{10} --bs-trees 1000

tree_file="alignment_raxml.raxml.support"

fi
echo '######################'

echo '	'

echo 'ProTree complete, the final tree with support from bootstraps (alignment_raxml.raxml.support, or fasttree.nwk) is in newick format and can now be viewed with your favourite tree viewer'

echo ' ' >> ../../../model.file

echo 'alignment based on: ' $final_num '  Single copy orthologues' >> ../../../model.file
echo '' >> ../../../model.file
echo 'BIC model used = '$model >> ../../../model.file

echo ' ' >> ../../../model.file

cat ../../../model.file

echo 'Final tree with support:' 

cp $tree_file ../../../

cat $tree_file

cd ../../../

tar -zcf 'ProTree_out.tar.gz' work_dir

echo ''

echo 'use " rm -r work_dir " to delete the working directory'

echo ' '

