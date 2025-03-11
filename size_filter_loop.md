Loop to filter contigs by size (5000 bp - can easily be changed), using BBMap (container)

Run from inside directory with the assemblies

Creates a filtered file for each assembly with _filtered.fasta amended to the end of the name, in a subdirectory called filtered5k

```
for f in *.fasta;
do apptainer exec ../containers/bbmap_39.15.sif reformat.sh in=$f out=./filtered5k/${f/.fasta/_filtered.fasta} minlength=5000;
done
```
