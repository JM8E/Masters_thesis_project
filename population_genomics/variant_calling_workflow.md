# Variant calling workflow
unless otherwise stated - containers are used
## Quality control
Fastqc
```{r}
# run all read files in directory
apptainer exec containers/fastqc_0.12.1.sif fastqc -o fastqc/ reads/*.fastq.gz
```
multiqc
```{r}
apptainer exec containers/multiqc_1_27_1.sif multiqc fastqc/ -o multiqc/
```
Look at the output - copy overrepresented sequences into a txt file

## Trim adapters, etc
Use custom trim file when running Trimmomatic

See Trimming_loop.md

## Redo quality control
```{r}
# combo
apptainer exec containers/fastqc_0.12.1.sif fastqc -o fastqc/trimmomatic_out trimmomatic/trimmomatic_out/*.fastq.gz
apptainer exec containers/multiqc_1_27_1.sif multiqc fastqc/trimmomatic_out/ -o multiqc/trimmomatic_out/
```
## Variant calling
Index the reference
```{r}
apptainer exec ../containers/bwa-mem2_2.2.1.sif bwa-mem2 index -p FC_p.fa FC_p.fa
```
Mapping:

loop creates script to map the reads
```{r}
for R1 in trimmomatic/trimmomatic_out/*R1_001_paired*;
do
   R2=${R1//R1_001_paired.fastq.gz/R2_001_paired.fastq.gz}
   echo ”apptainer exec containers/bwa-mem2_2.2.1.sif bwa-mem2 mem -o mapping/${R1}_FCp.sam mapping/index/FC_p.fa trimmomatic/trimmomatic_out/${R1} trimmomatic/trimmomatic_out/${R2}” >> scripts/bwa-mem2_mapping.sh;
done;
echo “echo complete” >> scripts/bwa-mem2_mapping.sh
```
Execute the script
```{r}
bash scripts/bwa-mem2_mapping.sh
```
Convert to bam and remove sam
```{r}
for f1 in mapping/*sam;
do 
apptainer exec containers/samtools_1.21.sif samtools view -o ${f1}.bam ${f1};
rm ${f1};
done
```
Sort mapped reads
```{r}
for f1 in mapping/*bam; 
do 
apptainer exec containers/samtools_1.21.sif samtools sort ${f1} > ${f1}.sorted; 
done
```
Remove duplicates
```{r}
for f1 in mapping/*sorted; 
do 
apptainer exec containers/picard.3.1.1.sif picard MarkDuplicates -O ${f1}_dups -M ${f1}_metrics -I ${f1} --REMOVE_DUPLICATES true; 
done
```
Add read groups
```{r}
for f1 in mapping/*dups; 
do 
apptainer exec containers/picard.3.1.1.sif picard AddOrReplaceReadGroups R=mapping/index/FC_p.fa I=${f1} O=${f1}_rg \
       RGID=${f1} \
       RGLB=${f1} \
       RGPL=ILLUMINA \
       RGPU=${f1} \
       RGSM=${f1} 
done
```
Index mapped reads
```{r}
for f1 in mapping/*rg; 
do 
apptainer exec containers/samtools_1.21.sif samtools index ${f1}; 
done
```
Index reference
```{r}
apptainer exec containers/samtools_1.21.sif samtools faidx mapping/index/FC_p.fa
apptainer exec containers/gatk4_4.5.0.0.sif gatk CreateSequenceDictionary R=mapping/index/FC_p.fa O=mapping/index/FC_p.dict
```
Make genomic vcf files
```{r}
for f1 in mapping/*FC*rg; 
do 
apptainer exec containers/gatk4_4.5.0.0.sif gatk HaplotypeCaller -R mapping/index/FC_p.fa -I ${f1} -O ${f1}.gvcf -ERC GVCF;  
done
```
Combine gvcf files:

loop writes script to combine all gvcf files
```{r}
echo 'apptainer exec containers/gatk4_4.5.0.0.sif gatk CombineGVCFs -R mapping/index/FC_p.fa -O mapping/combined/FC_combined \' >> combinegvcf.sh;
for file in *.gvcf;
do echo -V ${file} \\ >> combinegvcf.sh;
done;
for file in mapping/*.gvcf;
do echo -V ${file} \\ >> combinegvcf.sh;
done
```
Execute script
```{r}
bash combinegvcf.sh
```
Genotype gvcf file
```{r}
apptainer exec containers/gatk4_4.5.0.0.sif gatk GenotypeGVCFs -R mapping/index/FC_p.fa -V mapping/combined/FC_combined -O mapping/vcf/FC_raw.vcf
```
Filter by depth and quality
```{r}
apptainer exec containers/picard.3.1.1.sif picard FilterVcf -I mapping/vcf/FC_raw.vcf -O mapping/vcf/FC_filtered.vcf --MIN_DP 5 --MIN_GQ 10
```
Retain only biallelic sites
```{r}
apptainer exec containers/vcftools_0.1.16.sif vcftools --vcf mapping/vcf/FC_filtered.vcf --out mapping/vcf/FC_filtered_biallelic --recode --min-alleles 2 --max-alleles 2 --remove-indels
```
Filter again - thinning, missingness, minimum 1 with non-reference allele
```{r}
apptainer exec containers/vcftools_0.1.16.sif vcftools --vcf mapping/vcf/FC_filtered_biallelic.recode.vcf --out mapping/vcf/FC_thin_pop_struc --recode --thin 10000
apptainer exec containers/vcftools_0.1.16.sif vcftools --vcf mapping/vcf/FC_thin_pop_struc.recode.vcf --out mapping/vcf/FC_thin_struc_use –recode --max-missing 1 --non-ref-ac-any 1
```
