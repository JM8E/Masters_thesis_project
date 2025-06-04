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
Figures in R Studio
```{r}
# repeat commands for all assemblies until otherwise stated
# divide into 1 dataframe / assembly - change numbers for each assembly (check your file for the right numbers)
library(tidyverse)
ATCC42080_c <- pred_caz_sub %>% 
    slice(2:215)

# rename col
ATCC42080_c <- ATCC42080_c %>% 
    dplyr::rename('Substrate' = 'X1')

# count nr of CAZymes that work on each substrate
ATCC42080_subc <- ATCC42080_c %>% count(Substrate)

# add columns for assembly and species
ATCC42080_subc <- ATCC42080_subc %>% add_column(Assembly = "ATCC42080")
ATCC42080_subc <- ATCC42080_subc %>% add_column(Species = "T.laevis")

# combine all dataframes - from now use commands only once
sub_counts <- rbind(ATCC42080_subc, QB_1_subc, WJ_1_subc, TX6_subc, TC_1_MSG_1_subc, DAOMC238038_subc, CBS122995_subc, DAOMC238055_subc, DAOMC238035_subc, DAOMC238041_subc, CBS122992_subc, DAOMC238030_subc, DAOMC238049_subc, FC_subc, GC_subc, IC_subc, OA2_subc, PSWKBGH_1_subc, PSWKBGH_2_subc, DAOMC238032_subc, DAOMC236426_subc, DAOMC236422_subc, Umaydis_subc, KD_subc, LC_subc, NC_subc, QC_subc)

# rename col n
sub_counts <- sub_counts %>%
    dplyr::rename('Frequency' = 'n')

# turn frequency to numeric
sub_counts$Frequency <- as.numeric(sub_counts$Frequency)

# rename "-"
sub_counts <- sub_counts %>% 
    mutate(Substrate = recode(Substrate, '-' = 'Unknown'))

# group all identified substrates
sub_counts_filtered_t <- mutate(sub_counts, Substrate = case_when(
  Substrate != "Unknown" ~ "Identified", 
  TRUE   ~ Substrate 
))

# try plot with unknown / identified
#pdf("cazyme_subs_inc_unknown.pdf")
sub_counts_filtered_t %>% mutate(Species = factor(Species, levels = c("T.controversa","T.caries","T.laevis","T.bromi","T.goloskokovii","T.fusca","T.walkeri","T.indica","T.rugispora","T.horrida","T.setariae","U.maydis"))) %>% 
 ggplot(aes(fill = Substrate, x = Assembly, y = Frequency)) + 
 geom_bar(position = "stack", stat = "identity") + scale_y_continuous(expand = c(0, 0), limits = c(0, 275)) + theme(legend.text = element_text(size = 10)) + scale_fill_manual(values = c("darkolivegreen4", "skyblue")) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size=10)) + facet_grid(~ Species, switch = "x", space = "free_x", scales = "free_x") + theme(strip.text.x = element_text(angle=90, hjust=1), strip.placement = "outside",
        strip.background = element_rect(fill = c("honeydew"), color = "grey"),
        panel.spacing = unit(-.01,"cm"))
#dev.off()

# remove "-" values
sub_counts_4 <- sub_counts %>%
  filter(!(Substrate %in% c("Unknown")))

# filter 
sub_counts_4 <- filter(sub_counts_4, Frequency >= 3)

# plot
#pdf("cazyme_subs_=>3.pdf")
sub_counts_4 %>% mutate(Species = factor(Species, levels = c("T.controversa","T.caries","T.laevis","T.bromi","T.goloskokovii","T.fusca","T.walkeri","T.indica","T.rugispora","T.horrida","T.setariae","U.maydis"))) %>% 
 ggplot(aes(fill = Substrate, x = Assembly, y = Frequency)) + 
 geom_bar(position = "stack", stat = "identity") + scale_y_continuous(expand = c(0, 0), limits = c(0, 75)) + labs(title = "CAZyme substrates", size = 14) + theme(plot.title = element_text(hjust = 0.5, vjust = 1)) + theme(legend.text = element_text(size = 10)) + scale_fill_manual(values = c("darkolivegreen", "darkmagenta", "darkcyan", "darkolivegreen3", "orchid2", "aquamarine2", "darkolivegreen1", "mediumpurple1", "grey","royalblue", "lightgreen", "lavender", "darkseagreen4", "deepskyblue", "palevioletred" )) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size=8)) + facet_grid(~ Species, switch = "x", space = "free_x", scales = "free_x") + theme(strip.text.x = element_text(angle=90, hjust=1), strip.placement = "outside",
        strip.background = element_rect(fill = c("honeydew"), color = "grey"),
        panel.spacing = unit(-.01,"cm"))
#dev.off()

# plots side by side 
library(patchwork)
pdf("caz_sub_plots_new_guess.pdf", width=12, height=10)
all_sub <- sub_counts_filtered_t %>% mutate(Species = factor(Species, levels = c("T.controversa","T.caries","T.laevis","T.bromi","T.goloskokovii","T.fusca","T.walkeri","T.indica","T.rugispora","T.horrida","T.setariae","U.maydis"))) %>% 
 ggplot(aes(fill = Substrate, x = Assembly, y = Frequency)) + 
 geom_bar(position = "stack", stat = "identity") + scale_y_continuous(expand = c(0, 0), limits = c(0, 275)) + theme(legend.text = element_text(size = 10)) + scale_fill_manual(values = c("darkolivegreen4", "skyblue")) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size=7)) + facet_grid(~ Species, switch = "x", space = "free_x", scales = "free_x") + theme(strip.text.x = element_text(angle=90, hjust=1, size = 7), strip.placement = "outside",
        strip.background = element_rect(fill = c("honeydew"), color = "grey"),
        panel.spacing = unit(-.01,"cm"))

zoomed_sub <- sub_counts_4 %>% mutate(Species = factor(Species, levels = c("T.controversa","T.caries","T.laevis","T.bromi","T.goloskokovii","T.fusca","T.walkeri","T.indica","T.rugispora","T.horrida","T.setariae","U.maydis"))) %>% 
 ggplot(aes(fill = Substrate, x = Assembly, y = Frequency)) + 
 geom_bar(position = "stack", stat = "identity") + scale_y_continuous(expand = c(0, 0), limits = c(0, 75)) + theme(plot.title = element_text(hjust = 0.5, vjust = 1)) + theme(legend.text = element_text(size = 10)) + scale_fill_manual(values = c("darkolivegreen", "darkmagenta", "darkcyan", "darkolivegreen3", "orchid2", "aquamarine2", "darkolivegreen1", "mediumpurple1", "grey","royalblue", "lightgreen", "lavender", "darkseagreen4", "deepskyblue", "palevioletred" )) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size=7)) + facet_grid(~ Species, switch = "x", space = "free_x", scales = "free_x") + theme(strip.text.x = element_text(angle=90, hjust=1, size = 7), strip.placement = "outside",
        strip.background = element_rect(fill = c("honeydew"), color = "grey"),
        panel.spacing = unit(-.01,"cm"))

sub_combo <- all_sub + zoomed_sub

sub_combo + plot_annotation(title = "CAZyme substrates", theme = theme(plot.title = element_text(hjust = 0.5, size = 16, vjust = 1))) + plot_layout(widths = c(80, 80), heights = c(100, 100))
dev.off()

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
## CAZyme families

In R studio

Load overview tables from dbCAN3
```{r}
library(tidyverse)
# repeat commands for all genomes until otherwise stated
# rename column #ofTools to nrofTools and EC# to EC 
FCtable <- FCtable %>% rename("nrofTools" = "#ofTools", "EC" = "EC#")

# filter out all genes that are only predicted by one tool
FCtable_filtered <- filter(FCtable, !nrofTools == 1)

# new dataframe with all rows including '+' in the DIAMOND col
FCctemp <- FCtable_filtered %>% filter(grepl('\\+', DIAMOND))

# split the two cazyme family names into separate cols (with + as delimiter)
FCctemp <- FCctemp %>% separate(DIAMOND, c('DIAMOND', 'Diamond'), sep = '\\+')

# merge DIAMOND and Diamond, creating new rows duplicating the other cols
df1 <- FCctemp %>% select(`Gene ID`, EC, HMMER, dbCAN_sub, DIAMOND, Signalp, nrofTools)
df2 <- FCctemp %>% select(`Gene ID`, EC, HMMER, dbCAN_sub, Diamond, Signalp, nrofTools) %>% rename(DIAMOND = Diamond)
FCctemp <- bind_rows(df1, df2)

# remove rows from FCtable_filtered with + in DIAMOND col
FCtable_filtered <- FCtable_filtered %>% filter(!grepl('\\+', DIAMOND))

# merge FCtable_filtered and FCtest
FCtable_filtered <- bind_rows(FCtable_filtered, FCctemp)

# count nr of genes in families
FC_cazf <- FCtable_filtered %>% count(DIAMOND)

# add cols: assembly, species
FC_cazf <- FC_cazf %>% add_column(Assembly = "FC")
FC_cazf <- FC_cazf %>% add_column(Species = "T.caries")

# merge caz count dataframes - from now on do commands only once
caz_family_merged <- bind_rows(ATCC42080_cazf, QB_1_cazf, WJ_1_cazf, TX6_cazf, TC_1_MSG_1_cazf, DAOMC238038_cazf, CBS122995_cazf, DAOMC238055_cazf, DAOMC238035_cazf, DAOMC238041_cazf, CBS122992_cazf, DAOMC238030_cazf, DAOMC238049_cazf, FC_cazf, GC_cazf, IC_cazf, OA2_cazf, PSWKBGH_1_cazf, PSWKBGH_2_cazf, DAOMC238032_cazf, DAOMC236426_cazf, DAOMC236422_cazf, Umaydis_cazf, KD_cazf, LC_cazf, NC_cazf, QC_cazf)

# change col names
caz_family_merged <- caz_family_merged %>% dplyr::rename('Frequency' = 'n', 'CAZyme_family' = 'DIAMOND')

# rename N
caz_family_merged <- caz_family_merged %>% 
    mutate(CAZyme_family = recode(CAZyme_family, 'N' = 'Unknown'))

# group identified cazyme families
CAZfam_all <- mutate(caz_family_merged, CAZyme_family = case_when(
  CAZyme_family != "Unknown" ~ "Identified", 
  TRUE   ~ CAZyme_family 
))

# plot with unknown and identified
#pdf("cazyme_fams_inc_unknown.pdf")
CAZfam_all %>% mutate(Species = factor(Species, levels = c("T.controversa","T.caries","T.laevis","T.bromi","T.goloskokovii","T.fusca","T.walkeri","T.indica","T.rugispora","T.horrida","T.setariae","U.maydis"))) %>% 
 ggplot(aes(fill = CAZyme_family, x = Assembly, y = Frequency)) + 
 geom_bar(position = "stack", stat = "identity") + scale_y_continuous(expand = c(0, 0), limits = c(0, 275)) + theme(legend.text = element_text(size = 10)) + scale_fill_manual(values = c("darkolivegreen4", "skyblue")) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size=10)) + facet_grid(~ Species, switch = "x", space = "free_x", scales = "free_x") + theme(strip.text.x = element_text(angle=90, hjust=1), strip.placement = "outside",
        strip.background = element_rect(fill = c("honeydew"), color = "grey"),
        panel.spacing = unit(-.01,"cm"))
#dev.off()

# remove rows with unknown in cazyme family col
CAZfam_filtered <- filter(caz_family_merged, !CAZyme_family == 'Unknown')

# top 5 or 10 per assembly... Freq >= 
caz_family_counts <- filter(CAZfam_filtered, Frequency >= 4)

# stacked bars 
#pdf("caz_fam_zoomed.pdf")
caz_family_counts %>% mutate(Species = factor(Species, levels = c("T.controversa","T.caries","T.laevis","T.bromi","T.goloskokovii","T.fusca","T.walkeri","T.indica","T.rugispora","T.horrida","T.setariae","U.maydis"))) %>%  
 ggplot(aes(fill = CAZyme_family, x = Assembly, y = Frequency)) + 
 geom_bar(position = "stack", stat = "identity") + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values = c("orchid3", "skyblue", "violetred3", "navy", "royalblue", "darkseagreen", "darkmagenta", "magenta3", "gray10", "darkolivegreen3", "darkcyan",  "turquoise3", "slategray3", "darkseagreen1", "darkolivegreen", "blueviolet", "coral1", "plum3", "lavender", "aquamarine2", "purple4", "coral", "olivedrab", "firebrick")) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size=8)) + facet_grid(~ Species, switch = "x", space = "free_x", scales = "free_x") + theme(strip.text.x = element_text(angle=90, hjust=1), strip.placement = "outside",
        strip.background = element_rect(fill = c("honeydew"), color = "grey"),
        panel.spacing = unit(-.01,"cm"))
#dev.off()

# plots side by side
library(patchwork)
#pdf("caz_fam_plots2.pdf", width=10, height=8)
cazall <- CAZfam_all %>% mutate(Species = factor(Species, levels = c("T.controversa","T.caries","T.laevis","T.bromi","T.goloskokovii","T.fusca","T.walkeri","T.indica","T.rugispora","T.horrida","T.setariae","U.maydis"))) %>% 
 ggplot(aes(fill = CAZyme_family, x = Assembly, y = Frequency)) + 
 geom_bar(position = "stack", stat = "identity") + scale_y_continuous(expand = c(0, 0), limits = c(0, 275)) + theme(legend.text = element_text(size = 10)) + scale_fill_manual(values = c("darkolivegreen4", "skyblue")) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size=10)) + theme(legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=10), legend.text = element_text(size=8)) + facet_grid(~ Species, switch = "x", space = "free_x", scales = "free_x") + theme(strip.text.x = element_text(angle=90, hjust=1), strip.placement = "outside",
        strip.background = element_rect(fill = c("honeydew"), color = "grey"),
        panel.spacing = unit(-.01,"cm"))

cazzoomed <- caz_family_counts %>% mutate(Species = factor(Species, levels = c("T.controversa","T.caries","T.laevis","T.bromi","T.goloskokovii","T.fusca","T.walkeri","T.indica","T.rugispora","T.horrida","T.setariae","U.maydis"))) %>%  
 ggplot(aes(fill = CAZyme_family, x = Assembly, y = Frequency)) + 
 geom_bar(position = "stack", stat = "identity") + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values = c("orchid3", "skyblue", "violetred3", "navy", "royalblue", "darkseagreen", "darkmagenta", "magenta3", "gray10", "darkolivegreen3", "darkcyan",  "turquoise3", "slategray3", "darkseagreen1", "darkolivegreen", "blueviolet", "coral1", "plum3", "lavender", "aquamarine2", "purple4", "coral", "olivedrab", "firebrick")) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size=8)) + theme(legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=10), legend.text = element_text(size=8)) + guides(color = guide_legend(override.aes = list(size = 0.5))) + guides(shape = guide_legend(override.aes = list(size = 0.5))) + theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) + facet_grid(~ Species, switch = "x", space = "free_x", scales = "free_x") + theme(strip.text.x = element_text(angle=90, hjust=1), strip.placement = "outside",
        strip.background = element_rect(fill = c("honeydew"), color = "grey"),
        panel.spacing = unit(-.01,"cm")) 

cazfam_combo <- cazall + cazzoomed

cazfam_combo + plot_annotation(title = "CAZyme families", theme = theme(plot.title = element_text(hjust = 0.5, size = 16, vjust = 1))) + plot_layout(widths = c(80, 80), heights = c(100, 100))
#dev.off()

```
