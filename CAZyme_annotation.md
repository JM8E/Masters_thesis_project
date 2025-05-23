# CAZyme annotation

Run all whole genome protein sequences through the dbCAN3 webserver (http://dbcan.unl.edu/dbcan/)

Save all the output in directories for each genome

## Count cazymes

Remove newlines from fasta
```{r}
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' caz_annotation/FC.fa >> caz_annotation/FC.fasta
```
Filter fasta to remove seq only id by 1 tool
```{r}
awk '/^>/ {P=index($0,":1")==0} {if(P) print} ' FC.fasta > FC_filtered.fasta
```
Count total CAZymes
```{r}
for dir in ls *; 
do 
  echo $dir >> caz_count.txt; 
  grep ^\> ${dir}/${dir}_filtered.fasta | wc -l >> caz_count.txt; 
done
```
## Get CAZyme substrates
Get lists of gene ids
```{r}
for dir in *; 
do 
  grep ^\> ${dir}/*_filtered.fasta | cut -d' ' -f1 | sort | uniq | sed -i -e 's/>//g' >> ${dir}_gene_ids.txt; 
done
```

Cross ref with HMMERdbcansub

For each gene id get field containing substrate
```{r}
mkdir subs
for f in *_gene_ids.txt; 
do
  fgrep -f $f ${f//_gene_ids.txt//}/*sub.txt | cut -f4 >> subs/${f//_gene_ids.txt/_subs.txt}; 
done
```

Concatenate results
```{r}
cd subs/
for f in *; 
do
  echo $f >> pred_caz_sub.txt; 
  cat $f >> pred_caz_sub.txt; 
done
```

Tidy up file
```{r}
sed -i -e 's/_subs.txt//g' pred_caz_sub.txt
```
## Cross ref cazymes and signal peptides
Copy code into .sh file and run
```{r}
# make new dir
mkdir caz_sigp
# copy signalP output (processed_entries.fasta) to caz_sigp
cd effector_annotation/
for dir in *.aa_; 
do 
    cp ${dir}/signalp/processed_entries.fasta ${dir}/1processed_entries.fasta;
    mv ${dir}/1processed_entries.fasta ../caz_sigp/${dir//.aa_\//}_processed_entries.fasta;
done 
cd ..
# copy all caz output (_filtered.fasta) into caz_sigp
cd caz_annotation
for dir in *; 
do 
    cp ${dir}/*_filtered.fasta ../caz_sigp/; 
done 
cd ..
# cat signalp and caz output fastas 
for f1 in caz_sigp/*_filtered.fasta;
do
    f2=${f1//_filtered.fasta/.aa__processed_entries.fasta}
    cat ${f1} ${f2} >> ${f1//_filtered.fasta/_caz_sigp.fasta}; 
done 
# extract gene id that appears in both and count
for f in caz_sigp/*caz_sigp*; 
do 
    echo $f >> caz_sigp/count_caz_sigp.txt; 
    grep ^\> $f | cut -d' ' -f1 | sort | uniq -d | wc -l >> caz_sigp/count_caz_sigp.txt; 
done
# fix names, remove newlines
sed -i -e 's/caz_sigp\///g' -e 's/_caz_sigp.fasta//g' caz_sigp/count_caz_sigp.txt
sed -i -e '$!N;s/\n/ /g' caz_sigp/count_caz_sigp.txt
```
