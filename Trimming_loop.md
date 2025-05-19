how to run Trimmomatic on multiple sets of paired end read files without having to manually write and execute a command for each set.

create a file with a command for each set of files.

The following loop creates an executable file with a command for each set of paired read files:
```
for R1 in *R1*;
do
   R2=${R1//R1_001.fastq.gz/R2_001.fastq.gz}
   R1paired=${R1//.fastq.gz/_paired.fastq.gz}
   R1unpaired=${R1//.fastq.gz/_unpaired.fastq.gz}	
   R2paired=${R2//.fastq.gz/_paired.fastq.gz}
   R2unpaired=${R2//.fastq.gz/_unpaired.fastq.gz}	
   echo "apptainer exec ../containers/trimmomatic_0.39.sif trimmomatic PE -phred33 $R1 $R2 ../trimmomatic/trimmomatic_out/$R1paired ../trimmomatic/trimmomatic_out/$R1unpaired ../trimmomatic/trimmomatic_out/$R2paired ../trimmomatic/trimmomatic_out/$R2unpaired ILLUMINACLIP:../trimmomatic/overrep_seq.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" >> ../scripts/trimmomatic.sh;
done;
echo “echo “complete”” >> ../scripts/trimmomatic.sh
```
Make sure the pattern "R1_001.fastq.gz" matches all the files you want to run - the code wont change the names to R2 etc if it doesn't match.

settings for Trimmomatic are very general, but can easily be changed - just edit the loop.

The file with sequences to trim (overrep_seq.fa) is specific and will need to be replaced, either with a general option available with Trimmomatic (for example TruSeq3-PE.fa) or with your own file.

This uses a container with Trimmomatic (apptainer system) - installed versions have other commands.

copy the loop into a new .sh file and execute it with:
```
bash <file_name>.sh
```

Then execute the file with the commands:
```
bash ../scripts/trimmomatic.sh
```

It may take some time to run. The output will be four files per paired set.
run QC after to make sure the trimming worked.

```{r}
echo "apptainer exec containers/fastqc_0.12.1.sif fastqc -o fastqc/trimmomatic_out trimmomatic/trimmomatic_out/*.fastq.gz
apptainer exec containers/multiqc_1_27_1.sif multiqc fastqc/trimmomatic_out/ -o multiqc/trimmomatic_out/" >> scripts/qc.sh
```
```{r}
bash ../scripts/qc.sh
```
