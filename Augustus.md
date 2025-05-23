Annotates (Tilletia) genomes using Augustus (container) - to use on other species, check for appropriate --species=

Run from inside directory containing the genome assemblies

Creates output files in output directory called augustus

```
for f in *.fasta; 
do apptainer exec ../containers/augustus_3.5.0.sif augustus --species=ustilago_maydis ${f} > ../augustus/${f}_augustus; 
done
```

Use getAnnoFasta.pl to extract amino acid sequences etc
```{r}
perl containers/getAnnoFasta.pl augustus/NC_augustus --seqfile=genomes/NC.fasta
```

Move amino acid seq into protein dir, then

Extract total predicted genes
```{r}
for f in protein/*.fasta*; 
do 
echo $f >> gene_counts.txt; 
grep ^\> $f | cut -d' ' -f1 | sort | wc -l >> gene_counts.txt; 
done
# fix file – remove newlines, etc
sed -i '$!N;s/\n/ /g' gene_counts.txt
sed -i -e 's/protein\///g' -e 's/.fasta//g' gene_counts.txt
```
