Annotates (Tilletia) genomes using Augustus (container) - to use on other species, check for appropriate --species=
Run from inside directory containing the genome assemblies
Creates output files in output directory called augustus

```
for f in *.fasta; 
do apptainer exec ../containers/augustus_3.5.0.sif augustus --species=ustilago_maydis ${f} > ../augustus/${f}_augustus; 
done
```
