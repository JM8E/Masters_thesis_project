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

