# Phylogenetic network workflow

Convert the vcf from variant_calling_workflow.md to nexus format. Used vcf_to_nexus.py (https://github.com/bmichanderson/RAD_scripts/blob/master/vcf_to_nexus.py)

Make file executable
```{r}
chmod 700 vcf_to_nexus.py
```

Make sure your vcf doesn't have any pipes in the genotypes
```{r}
sed -i -e 's/0|0/0\/0/g' -e 's/0|1/0\/1/g' -e 's/1|0/1\/0/g' -e 's/1|1/1\/1/g' FC_filtered_biallelic.recode.vcf
```

Run file conversion - command specificly intended for SplitsTree
```{r}
python ./vcf_to_nexus.py -o splits FC_filtered_biallelic.recode.vcf
```

Load your file into SplitsTree6 (app), you can colour your samples by species for example

Also did one without controversa

Used VCFtools to remove controversa samples from the vcf
```{r}
apptainer run containers/vcftools_0.1.16.sif vcftools --vcf mapping/vcf/FC_filtered_biallelic.recode.vcf --remove-indv OA2_FC --remove-indv OL14_FC --remove-indv OR_FC --remove-indv OV_FC --remove-indv OW_FC --remove-indv DAOMC251939_FC --remove-indv CBS121952_FC --remove-indv DAOMC236426_FC --remove-indv DAOMC238052_FC --out mapping/vcf/FC_nocon -–recode
```
And converted that to nexus format too and loaded it into SpitsTree6
