# Effector annotation loop (fungal genome)

SignalP, EffectorP and Deeploc2 are locally installed, Phobius uses a container

Works from inside directory with protein fasta files

Creates a directory for each assembly (in the output directory), containing all the results for that assembly

```
for f in *aa*;
do signalp6 --fastafile ${f} --organism eukarya --output_dir ../effector_annotation/${f}_/signalp/ --model_dir /bioinfo/signalp6_fast/signalp-6-package/models;
  apptainer exec ../containers/phobius.sif phobius ../effector_annotation/${f}_/signalp/processed_entries.fasta > ../effector_annotation/${f}_/phobius_output_long;
  python /bioinfo/EffectorP-3.0/EffectorP.py -f -i ../effector_annotation/${f}_/signalp/processed_entries.fasta -E ../effector_annotation/${f}_/${f}_effectors.fasta > ../effector_annotation/${f}_/${f}_effectorp;
  deeploc2 --fasta ../effector_annotation/${f}_/${f}_effectors.fasta --output ../effector_annotation/${f}_/${f/.aa/_effectors_deeploc};
done
```
# Get count of effectors
gives you total effectors, apoplastic and cytoplasmic effectors
```{r}
for file in *effectorp; 
do 
echo $file >> count_effectors.txt; 
cat $file | grep 'Number of' >> effectorp/count_effectors.txt; 
done
```

# Count signal peptides
```{r}
for dir in ls *.aa_; 
do 
echo ${dir//.aa_/ } >> count_signalp.txt; 
grep ^\> ${dir}/signalp/processed_entries.fasta | wc -l >> count_signalp.txt; 
done
```
