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
# first move all effectorp output to one directory
for file in *effectorp; 
do 
echo $file >> count_effectors.txt; 
cat $file | grep 'Number of' >> effectorp/count_effectors.txt; 
done
```
# Count those proteins predicted to be both effectors and CAZymes
Cat the effector and CAZyme fastats for each assembly, _eandc suffix

loop count effectors and cazymes by gene id
```{r}
for f in *eandc*; 
do 
echo ${f//_eandc.fasta/ } >> counts/count_effcaz_nr.txt; 
grep ^\> $f | cut -d' ' -f1 | sort | uniq -d | wc -l >> counts/count_effcaz_nr.txt; 
# next line adds seq ids - if you want them, otherwise remove
grep ^\> $f | cut -d' ' -f1 | sort | uniq -d >> /counts/count_effcaz.txt; 
done
```


# Figures
In R Studio

Put all effector (and CAZyme) prediction values in excel - columns with total effectors, cytoplasmic, apoplastic, both effector and CAZyme, CAZymes - gene_pred - load into R Studio

numbers of cols depend on how you wrote them into excel - change to fit your situation

Effectors - cytoplasmic and apoplastic
```{r}
library(tidyverse)

addnew_cyt <- cbind(gene_pred[3], gene_pred[5], gene_pred[6])
addnew_apo <- cbind(gene_pred[4], gene_pred[5], gene_pred[6])
addnew_cyt <- addnew_cyt %>% add_column(Type = "Cytoplasmic effectors")
addnew_apo <- addnew_apo %>% add_column(Type = "Apoplastic effectors")
addnew_apo <- rename(addnew_apo, Number_of_effectors = 'Apoplastic effectors')
addnew_cyt <- rename(addnew_cyt, Number_of_effectors = 'Cytoplasmic effectors')
addnew_eff <- rbind(addnew_apo,addnew_cyt)

# stacked bar plot with effectors divided into cytoplasmic and apoplastic
#pdf("effectors.pdf")
addnew_eff %>% mutate(addnew_eff, Type = factor(Type, levels=c("Apoplastic effectors", "Cytoplasmic effectors"))) %>% mutate(addnew_eff, Species=factor(Species,levels=c("T.controversa", "T.caries", "T.laevis", "T.bromi", "T.goloskokovii", "T.fusca", "T.walkeri", "T.indica", "T.rugispora", "T.horrida", "T.setariae", "U.maydis"))) %>% ggplot(aes(fill = Type, x = Assembly, y = Number_of_effectors)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 450)) + ggtitle("Predicted apoplastic and \n cytoplasmic effectors") + theme(plot.title = element_text(hjust = 0.5)) +
 geom_bar(position = "stack", stat = "identity") + scale_fill_manual(values = c("darkmagenta","darkcyan")) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) + facet_grid(~ Species, switch = "x", space = "free_x", scales = "free_x") + theme(strip.text.x = element_text(angle=90, hjust=1), strip.placement = "outside",
        strip.background = element_rect(fill = c("honeydew"), color = "grey"),
        panel.spacing = unit(-.01,"cm")) 
#dev.off()

# -II- but in percentages
#pdf("effplot_perc1.pdf")
addnew_eff %>% 
 ggplot(aes(fill = Type, x = Assembly, y = Number_of_effectors)) + 
 geom_bar(position = "fill", stat = "identity") + ylab("Percentage of effectors") + scale_fill_manual(values = c("olivedrab", "turquoise3")) + theme(plot.margin = unit(c(1,1,1,1), "cm")) + theme(axis.title.y=element_text(hjust=0.8)) + scale_y_continuous(labels = scales::percent) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) + facet_grid(~ Species, switch = "x", space = "free_x", scales = "free_x") + theme(strip.text.x = element_text(angle=90, hjust=1), strip.placement = "outside",
        strip.background = element_rect(fill = c("honeydew"), color = "grey"),
        panel.spacing = unit(-.01,"cm")) #scale_x_discrete(limits=c("DAOMC238035","DAOMC238055","TC-1-MSG-1","DAOMC238032","FC","IC","ATCC42080","DAOMC238038","OA2","DAOMC236426","GC","CBS122995","CBS122992","DAOMC238041","DAOMC238049","DAOMC236422","PSWKBGH_1_1","PSWKBGH_1_2","DAOMC238030","QB-1","TX6","WJ-1","Umaydis")) 
#dev.off()
```

Effectors and CAZymes
```{r}
library(tidyverse)

# new col (nr) with nr of predicted Effectors, CAZymes and Both effector and CAZyme stacked, new col with type of protein (Type)
effectors <- cbind(gene_pred[11],gene_pred[5],gene_pred[6])
effectors <- effectors %>% add_column(Type = "Effectors")
effectors <- effectors %>% dplyr::rename('nr' = 'Effectors')
cazymes <- cbind(gene_pred[10],gene_pred[5],gene_pred[6])
cazymes <- cazymes %>% add_column(Type = "CAZymes")
cazymes <- cazymes %>% dplyr::rename('nr' = 'CAZymes')
bothec <- cbind(gene_pred[9],gene_pred[5],gene_pred[6])
bothec <- bothec %>% add_column(Type = "Both effector and CAZyme")
bothec <- bothec %>% dplyr::rename('nr' = 'Both effector and CAZyme')
pred_temp <- rbind(effectors,bothec,cazymes)

# remove NA
pred_temp1 <- pred_temp[complete.cases(pred_temp),]

# stacked bar plot showing nr of predicted Effectors, CAZymes and Both effector and CAZyme
#pdf("pred_effandcaz_stacked_new_guess.pdf")
pred_temp1 %>% mutate(Type = factor(Type, levels = c("Effectors", "Both effector and CAZyme", "CAZymes"))) %>% mutate(Species = factor(Species, levels = c("T.controversa","T.caries","T.laevis","T.bromi","T.goloskokovii","T.fusca","T.walkeri","T.indica","T.rugispora","T.horrida","T.setariae","U.maydis"))) %>%
 ggplot(aes(fill = Type, x = Assembly, y = nr)) + ggtitle("Predicted effectors and CAZymes") + theme(plot.title = element_text(hjust = 0.5, size = 16)) +
 geom_bar(position = "stack", stat = "identity") + scale_y_continuous(expand = c(0, 0), limits = c(0, 700)) + theme(axis.title.y=element_text(size = 14)) + scale_fill_manual(values = c("olivedrab4", "aquamarine2", "darkmagenta")) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size = 10)) + facet_grid(~ Species, switch = "x", space = "free_x", scales = "free_x") + theme(strip.text.x = element_text(angle=90, hjust=1), strip.placement = "outside",
        strip.background = element_rect(fill = c("honeydew"), color = "grey"),
        panel.spacing = unit(-.01,"cm"))
#dev.off()

# stacked bar plot of percentages - predicted Effectors, CAZymes and Both effector and CAZyme
#pdf("predicted_effandcaz_perc.pdf")
pred_temp1 %>% mutate(Type = factor(Type, levels = c("Effectors", "Both effector and CAZyme", "CAZymes"))) %>% mutate(Species = factor(Species, levels = c("T.controversa","T.caries","T.laevis","T.bromi","T.goloskokovii","T.fusca","T.walkeri","T.indica","T.rugispora","T.horrida","T.setariae","U.maydis"))) %>%
ggplot(aes(fill = Type, x = Assembly, y = nr)) + 
 geom_bar(position = "fill", stat = "identity") + ylab("Percentage") + scale_fill_manual(values = c("darkred", "darkcyan", "darkorchid")) + theme(plot.margin = unit(c(1,1,1,1), "cm")) + theme(axis.title.y=element_text(hjust=0.8)) + scale_y_continuous(labels = scales::percent) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) + facet_grid(~ Species, switch = "x", space = "free_x", scales = "free_x") + theme(strip.text.x = element_text(angle=90, hjust=1), strip.placement = "outside",
        strip.background = element_rect(fill = c("honeydew"), color = "grey"),
        panel.spacing = unit(-.01,"cm"))
#dev.off()
```
Effector localizations

```{r}
# load Deeploc 2 results into R Studio

library(tidyverse)

# until otherwise stated - repeat commands for all assemblies
# count localizations
#ATCC42080_eloc <- ATCC42080_results_20250220_144123 %>% count(Localizations)

#add cols: assembly, species - matching the assembly in question
ATCC42080_eloc <- ATCC42080_eloc %>% add_column(Assembly = "ATCC42080")
ATCC42080_eloc <- ATCC42080_eloc %>% add_column(Species = "T.laevis")

# following command only once
# bind rows
Eff_loc <- rbind(ATCC42080_eloc,QB_1_eloc,WJ_1_eloc,TX6_eloc,TC_1_MSG_1_eloc,DAOMC238038_eloc,CBS122995_eloc,DAOMC238055_eloc,DAOMC238035_eloc,DAOMC238041_eloc,CBS122992_eloc,DAOMC238030_eloc,DAOMC238049_eloc,FC_eloc,GC_eloc,IC_eloc,OA2_eloc,PSWKBGH_1_eloc,PSWKBGH_2_eloc,DAOMC238032_eloc,DAOMC236426_eloc,DAOMC236422_eloc,Umaydis_eloc,KD_eloc,LC_eloc,NC_eloc,QC_eloc)

# rename Var1 and freq
Eff_loc <- Eff_loc %>% dplyr::rename('Localization' = 'Var1')
Eff_loc <- Eff_loc %>% dplyr::rename('Frequency' = 'n')

# filter 
Eff_loc_filtered <- filter(Eff_loc, Frequency >=6 )

# stacked bar plot
#pdf("effector_loc1.pdf")
Eff_loc_filtered %>% mutate(Species = factor(Species, levels = c("T.controversa","T.caries","T.laevis","T.bromi","T.goloskokovii","T.fusca","T.walkeri","T.indica","T.rugispora","T.horrida","T.setariae","U.maydis"))) %>% 
 ggplot(aes(fill = Localizations, x = Assembly, y = Frequency)) + 
 geom_bar(position = "stack", stat = "identity") + scale_y_continuous(expand = c(0, 0), limits = c(0, 425)) + labs(title = "Effector localizations", size = 14) + scale_fill_manual(values = c("darkolivegreen", "darkmagenta", "darkcyan",  "orchid2", "aquamarine2", "darkolivegreen3", "mediumpurple1", "darkolivegreen1", "grey", "skyblue")) + theme(plot.title = element_text(hjust = 0.5, vjust = 1)) + theme(legend.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size=8)) + facet_grid(~ Species, switch = "x", space = "free_x", scales = "free_x") + theme(strip.text.x = element_text(angle=90, hjust=1), strip.placement = "outside",
        strip.background = element_rect(fill = c("honeydew"), color = "grey"),
        panel.spacing = unit(-.01,"cm"))
#dev.off()

```
# Count signal peptides
```{r}
for dir in ls *.aa_; 
do 
echo ${dir//.aa_/ } >> count_signalp.txt; 
grep ^\> ${dir}/signalp/processed_entries.fasta | wc -l >> count_signalp.txt; 
done
```
