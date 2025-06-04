# Download containers

MCScanX from seqera (https://seqera.io/containers/)

```{r}
singularity pull mcscanx_1.0.0.sif oras://community.wave.seqera.io/library/mcscanx:1.0.0--a4560d3bd6065987
```

# Download/copy java plotting scripts

dual_synteny_plotter.java can be found at:
https://github.com/wyp1125/MCScanX/blob/master/downstream_analyses/dual_synteny_plotter.java 

circle_plotter.java can be found at:
https://github.com/wyp1125/MCScanX/blob/master/downstream_analyses/circle_plotter.java 

Cubic.java can be found at:
https://github.com/wyp1125/MCScanX/blob/master/downstream_analyses/Cubic.java 

IMPORTANT!!!! – put circle_plotter.java and Cubic.java in the same file

# Make blast databases

```{r}
# make one for each assembly you want to blast
# example command (local installation) – whole genome protein sequence
makeblastdb -in synteny/protein/NC.fasta -out synteny/db/NC -dbtype prot
# output is three files - .phr .pin .psq
```
```{r}
# loop
for f in *.fasta; 
do 
	makeblastdb -in $f -out ../db/${f//.fasta//} -dbtype prot; 
done
```

# Run blastp

```{r}
# example command (local installation)
blastp -db synteny/db/NC -query synteny/protein/FC_p.fasta -evalue 1e-10 -num_alignments 5 -outfmt 6 -out synteny/blast_out/NC_FC.blast
```
run for all the combinations of database and query for your desired comparison
in this case: FCxFC, FCxNC, NCxNC, NCxFC

blast takes time to run – use a screen session
especially if you want to run many comparisons
```{r}
screen -S blastp_runs
```

write a bash script with all the commands and execute it
```{r}
bash blast.sh
```

then you can detach screen and do something else
```{r}
screen -d blastp_runs
```

simply reattach when you want to check
```{r}
screen -r blastp_runs
```

# Create “gff” with right format for MCScanX

```(r}
# this uses augustus output as the starting point
# remove #lines, lines with start_codon, intron, CDS, stop_codon
# remove columns 2,6,7 and 8
sed -e '/^#/d' -e '/start_codon/d' -e '/intron/d' -e '/CDS/d' -e '/stop_codon/d' synteny/augustus/NC.fasta_augustus | cut --complement -f2,6,7,8 > synteny/NCtest
# get gene ids from transcript lines, everything else from gene lines
cut -f5 synteny/NCtest | sed '/t/!d' > synteny/NCtest1
cut --complement -f5 synteny/NCtest | sed '/transcript/d' > synteny/NCtest2
# put them together, get only the needed columns and reorder them
paste synteny/NCtest1 synteny/NCtest2 | awk -F '\t' 'BEGIN { OFS=FS } { print $2, $1, $4, $5 }' > synteny/NCtest3
# remove temporary files
rm NCtest NCtest1 NCtest2
```

or full script with loops
```{r}
for f in *_augustus; 
do 
	sed -e '/^#/d' -e '/start_codon/d' -e '/intron/d' -e '/CDS/d' -e '/stop_codon/d' $f | cut --complement -f2,6,7,8 > ../mcscanx_gff/${f}_temp; 
done
cd ../mcscanx_gff/
for f in *_temp; 
do 
	cut -f5 $f | sed '/t/!d' > ${f}1; 
	cut --complement -f5 $f | sed '/transcript/d' > ${f}2; 
	paste ${f}1 ${f}2 | awk -F '\t' 'BEGIN { OFS=FS } { print $2, $1, $4, $5 }' > ${f//.fasta_augustus_temp/.gff}; 
rm $f ${f}1 ${f}2; 
done
```

# Gene and contig ids

make sure gene and contig ids are different for the different assemblies – in all files
make sure the ids in .gff and .blast match for a given assembly

suggest adding a distinct prefix or suffix for each assembly
note: in blast output 1st column = query, 2nd column = database
add assembly id to gene ids for FCxFC .blast
also works for .gff and .fasta – change file path and substitution
```{r}
sed -i -e 's/t1/t1FC/g' synteny/misc/FC_FCtest.blast
```

add assembly id to gene ids for FC=database
```{r}
awk -F'\t' 'sub(/t1/,"t1FC",$2)' OFS='\t' synteny/blast_out/FC_NC.blast > synteny/misc/FC_NCtest.blast
```

add assembly id to gene ids for FC=query
```{r}
awk -F'\t' 'sub(/t1/,"t1FC",$1)' OFS='\t' synteny/blast_out/NC_FC.blast > synteny/misc/NC_FCtest.blast
```

change contig names
```{r}
sed -i -e 's/utg/FC/g' FC.gff
```

# Concatenate blast output and gffs

put in same directory, with same prefix

concatenate blast output
```{r}
cat synteny/blast_out/FC_FC.blast synteny/blast_out/FC_NC.blast synteny/blast_out/NC_FC.blast synteny/blast_out/NC_NC.blast > synteny/data/FC_NC_all.blast
```

concatenate gffs
```{r}
cat synteny/NCtest3 synteny/FCtestall > synteny/data/FC_NCtest.gff
```

# Run MCScanX

```{r}
apptainer exec containers/mcscanx_1.0.0.sif MCScanX synteny/data/FC_NC_all
```
output is two files - .collinearity and .tandem, and a directory - .html

# Make .ctl

for dual_synteny_plotter
```{r}
echo 600 > synteny/misc/FC_NC.ctl # dimensions – may need to change
echo 2500 >> synteny/misc/FC_NC.ctl # dimensions – may need to change
cut -f1 synteny/FCtestall | sort | uniq | awk 1 ORS=, >> synteny/misc/FC_NC.ctl
echo "" >> synteny/misc/FC_NC.ctl
cut -f1 synteny/NCtest3 | sort | uniq | awk 1 ORS=, >> synteny/misc/FC_NC.ctl
```

# Plot

dual_synteny_plotter
```{r}
java containers/dual_synteny_plotter.java -g synteny/data/FC_NCtest.gff -s synteny/data/FC_NCtest.collinearity -c synteny/misc/FC_NC.ctl -o synteny/plots/FC_NCtest_synteny.png
```

# Filter out small contigs

filter genome sequence with BBmap
```{r}
for f in *.fasta; 
do 
apptainer exec ../../containers/bbmap_39.15.sif reformat.sh in=$f out=${f/.fasta/_100kb_filtered.fasta} minlength=100000; 
apptainer exec ../../containers/bbmap_39.15.sif reformat.sh in=$f out=${f/.fasta/_50kb_filtered.fasta} minlength=50000; 
done
```

get contig ids
```{r}
for f in *50*; 
do 
    grep ^\> $f | cut -d' ' -f1 | sed -e 's/>//g' >> ${f/_filtered.fasta/_ids.txt}; 
done
for f in *100*; 
do 
    grep ^\> $f | cut -d' ' -f1 | sed -e 's/>//g' >> ${f/_filtered.fasta/_ids.txt}; 
done
```

rename contig ids
```{r}
sed -i -e 's/utg/FC/g' synteny/genomes/FC_50kb_ids.txt
```

use list of contig ids to filter .gff
```{r}
grep -f synteny/genomes/FC_50kb_ids.txt synteny/misc/FC_testnames > synteny/misc/FC50kbtestnames
grep -f synteny/genomes/NC_50kb_ids.txt synteny/misc/NCtest3 > synteny/misc/NC50kbtest3
cat synteny/misc/FC50kbtestnames synteny/misc/NC50kbtest3 > synteny/data/FC_NC_testnames_50kb.gff
```

use filtered .gff to remake .ctl
```{r}
echo 600 > synteny/misc/FC_NC_100.ctl # dimensions – may need to change
echo 1500 >> synteny/misc/FC_NC_100.ctl # -II-
cut -f1 synteny/FC100kbtestnames | sort | uniq | awk 1 ORS=, >> synteny/misc/FC_NC_100.ctl
echo "" >> synteny/misc/FC_NC_100.ctl
cut -f1 synteny/NCtest3 | sort | uniq | awk 1 ORS=, >> synteny/misc/FC_NC_100.ctl
```

copy .blast from previous – rename
```{r}
cp synteny/data/FC_NC_all.blast synteny/data/FC_NC_100kb.blast
```
loop
```{r}
for f in *.blast; do cp $f ../${f/.blast/_50kb.blast}; done
```

rerun MCScanX
```{r}
apptainer exec containers/mcscanx_1.0.0.sif MCScanX synteny/data/FC_NC_100kb
```
loop
```{r}
for f in *.blast; 
do 
	apptainer exec ../../containers/mcscanx_1.0.0.sif MCScanX ${f/.blast/}; 
done
```

plot
```{r}
java containers/dual_synteny_plotter.java -g synteny/data/FC_NC_100kb.gff -s synteny/data/FC_NC_100kb.collinearity -c synteny/misc/FC_NC_100kb.ctl -o synteny/plots/FC_NC_100kb.png
```

# Test on another set

concatenate blast output
```{r}
cat synteny/blast_mcscanx/FC_FC.blast synteny/blast_mcscanx/FC_PSWKBGH_1.blast synteny/blast_mcscanx/PSWKBGH_1_FC.blast synteny/blast_mcscanx/PSWKBGH_1_PSWKBGH_1.blast > synteny/data/FC_PSWKBGH_1_all.blast
```
concatenate gffs
```{r}
cat synteny/mcscanx_gff/PSWKBGH_1.gff synteny/mcscanx_gff/FC_p.gff > synteny/data/FC_PSWKBGH_1_all.gff
```
run MCScanX
```{r}
apptainer exec containers/mcscanx_1.0.0.sif MCScanX synteny/data/FC_PSWKBGH_1_all
```
make .ctl for dual synteny
```{r}
echo 600 > synteny/misc/FC_PSWKBGH_1.ctl
echo 2500 >> synteny/misc/FC_PSWKBGH_1.ctl
cut -f1 synteny/mcscanx_gff/FC_p.gff | sort | uniq | awk 1 ORS=, >> synteny/misc/FC_PSWKBGH_1.ctl
echo "" >> synteny/misc/FC_PSWKBGH_1.ctl
cut -f1 synteny/mcscanx_gff/PSWKBGH_1.gff | sort | uniq | awk 1 ORS=, >> synteny/misc/FC_PSWKBGH_1.ctl
```
plot before filter – WORKED!!
```{r}
java containers/dual_synteny_plotter.java -g synteny/data/FC_PSWKBGH_1_all.gff -s synteny/data/FC_PSWKBGH_1_all.collinearity -c synteny/misc/FC_PSWKBGH_1.ctl -o synteny/plots/FC_PSWKBGH_1.png
```
rename contig ids in .txt
```{r}
sed -e 's/utg/FC/g' synteny/genomes/FC_50kb_ids.txt > synteny/genomes/FC_50kb_ids_renamed.txt
```
filter gff
```{r}
grep -f synteny/genomes/FC_50kb_ids.txt synteny/mcscanx_gff/FC_p.gff > synteny/mcscanx_gff/FC_50kb.gff
grep -f synteny/genomes/PSWKBGH_1_50kb_ids.txt synteny/mcscanx_gff/PSWKBGH_1.gff > synteny/mcscanx_gff/PSWKBGH_1_50kb.gff
cat synteny/mcscanx_gff/FC_50kb.gff synteny/mcscanx_gff/PSWKBGH_1_50kb.gff > synteny/data/FC_PSWKBGH_1_50kb.gff
```
copy _all.blast from previous 
```{r}
cp synteny/data/FC_PSWKBGH_1_all.blast synteny/data/FC_PSWKBGH_1_50kb.blast
```
new .ctl
```{r}
echo 600 > synteny/misc/FC_PSWKBGH_1_50kb.ctl
echo 2000 >> synteny/misc/FC_PSWKBGH_1_50kb.ctl
cut -f1 synteny/mcscanx_gff/FC_50kb.gff | sort | uniq | awk 1 ORS=, >> synteny/misc/FC_PSWKBGH_1_50kb.ctl
echo "" >> synteny/misc/FC_PSWKBGH_1_50kb.ctl
cut -f1 synteny/mcscanx_gff/PSWKBGH_1_50kb.gff | sort | uniq | awk 1 ORS=, >> synteny/misc/FC_PSWKBGH_1_50kb.ctl
```
plotting with copied .collinearity doesn’t work
rerun MCScanX
```{r}
apptainer exec containers/mcscanx_1.0.0.sif MCScanX synteny/data/FC_PSWKBGH_1_50kb
```
plot – WORKED!!!
```{r}
java containers/dual_synteny_plotter.java -g synteny/data/FC_PSWKBGH_1_50kb.gff -s synteny/data/FC_PSWKBGH_1_50kb.collinearity -c synteny/misc/FC_PSWKBGH_1_50kb.ctl -o synteny/plots/FC_PSWKBGH_1_50kb.png
```
FC has 36,62% collinearity with PSWKBGH_1 (outgroup), 84,46% with NC

