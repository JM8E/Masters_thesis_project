# Population structure workflow

Containers used:

VCFtools v.0.1.16

BCFtools v.1.21

## Additional filtering of vcf file produced in variant_calling_workflow.md

Set thinning distance at 10000
```{r}
apptainer exec containers/vcftools_0.1.16.sif vcftools --vcf mapping/vcf/FC_filtered_biallelic.recode.vcf --out mapping/vcf/FC_thin_pop_struc --recode --thin 10000
```

Remove missing sites and those without at least one sample with an alternate allele
```{r}
apptainer exec containers/vcftools_0.1.16.sif vcftools --vcf mapping/vcf/FC_thin_pop_struc.recode.vcf --out mapping/vcf/FC_thin_struc_use –recode --max-missing 1 --non-ref-ac-any 1
```

## Convert vcf to geno format

Copy the following script to a .sh file and bash it
```{r}
# vcf2geno4LEA
# this script can be used to convert vcf to geno
# to be used in snmf() function from R package LEA

VCF=input.vcf # name of your input vcf

OUTFILE=vcf_4_snmf.geno
# file with only genotypes in 0, 1, 2 format
# can be used for snmf() function in R package LEA

# start OUTFILE with only genotypes
apptainer exec <path_to_bcftools_container> bcftools query -e 'STRLEN(REF)>1' -f '[\t%GT]\n' $VCF | sed -e 's/\.\/\./N/g' >> $OUTFILE
# change genotypes to 0, 1, 2
sed -i -e 's/0\/0/0/g' -e 's/0\/1/1/g' -e 's/1\/0/1/g' -e 's/1\/1/2/g' -e 's/0|0/0/g' -e 's/0|1/1/g' -e 's/1|0/1/g' -e 's/1|1/2/g' $OUTFILE
# remove tabs
sed -i -e "s|\t||g" $OUTFILE
```

## Population structure

snmf function in the R package LEA

Had to run on a server since it requires a lot of memory

Make a .Rmd file with the commands below (ex tilletia_snmf.Rmd)

Install / load the following packages: LEA, here, vcfR, adegenet, openxlsx, dplyr, tidyr, tibble, ggplot2

Example:
```{r}
library(LEA)
library(here)
library(vcfR)
library(adegenet)
library(openxlsx)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
```

Set location to your .Rmd with the commands:
```{r}
here::i_am("tilletia_snmf.Rmd")
```

Start a new snmf project, using your .geno file
```{r}
project.snmf = snmf("tilletia.geno",
K = 1:20,
entropy = TRUE,
ploidy = 2,
repetitions = 10,
project = "new")
# plot cross-entropy criterion of all runs of the project
pdf("ce_tilletia.pdf")
plot(project.snmf, cex = 1.2, col = "darkcyan", pch = 19)
dev.off()
```

This is what I used for my structure analysis including T. controversa, assess your cross entropy plot and choose the best K for your situation

You also need a file with your sample names as well as one with metadata, such as species and geographical origin
```{r}
# get the cross-entropy of the 10 runs for K = 5
ce = cross.entropy(project.snmf, K = 5)
# select the run with the lowest cross-entropy for K = 4
best = which.min(ce)
#add names
metadata = read.xlsx("sample_meta.xlsx")
indv_snmf = read.csv("sample_ids.csv", header = FALSE)
names(indv_snmf) = "Sample"
datalist = list()
for (i in c(2, 3, 4, 5, 6, 7, 8, 9)){
best = which.min(cross.entropy(project.snmf, K = i))
temp = as.data.frame(Q(project.snmf, i, best))
temp = cbind(indv_snmf, temp)
temp = temp %>%
gather("Cluster", "Admix_coef", -"Sample") %>%
mutate(K=i)
datalist[[i]] = as_tibble(temp)
}

snmf_results_per_K = bind_rows(datalist) %>%
inner_join(., metadata, by = c("Sample" = "Sample")) %>%
unite(Sample, Location, Region, Country, Morphology, col = "for_display", remove = F)

pdf("struc_country_withcon_s.pdf", width = 20, height = 12)
ggplot(snmf_results_per_K, aes(x = reorder(Sample, ID_sandL), y = Admix_coef, fill = Cluster,
text = Sample)) +
geom_bar(position = "stack", stat = "identity", show.legend = F) +
facet_grid(K~.) +
theme_bw() +
theme(axis.title = element_blank(),
axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
legend.title = element_blank())
dev.off()
```
