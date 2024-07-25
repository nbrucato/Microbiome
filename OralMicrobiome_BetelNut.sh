#!/bin/bash

##This script calculates microbial diversity indices based on Kraken results
##It can idenitify differences between groups if specified
#############################################
##Fill in information
#############################################
input_folder=  ##A folder including all Kraken files. One subfolder per individual
output_folder=  ##An output folder
output_filename=  ##Output filename
samplelist=  ##A 5-colums list of all studied individuals: ID Population Sex Age Group
microb_class=  ## a 2-colums list of all identified bacteria: Genus Phylum 
min_reads= ##Minimum number of reads to include the identification
LDA_thres=  ##threshold to detect results with the LDA. The lower the stricter. default 1000000

##############################################################################################################################################################################################################################################################################
######################################
##Data preparation
######################################
##Sort and calculate abundance of bacteria per ind
echo -e "ID\tgroup\tpop\tgender\tage\tpercent_read\treads1\treads2\ttype\ttype_nb\tgenus\tclass\tabundance" > $output_folder/$output_filename.summary.abundance.txt
while read Ind2; do
Ind=$(echo $Ind2 | awk '{print $1}')
Population=$(echo $Ind2 | awk '{print $2}')
Gender=$(echo $Ind2 | awk '{print $3}')
Group=$(echo $Ind2 | awk '{print $5}')
Age=$(echo $Ind2 | awk '{print $4}')
grep $'\tG\t' $input_folder/$Ind/"$Ind"_report.microbiome.txt | sed 's/ \+/ /g' | sed -e 's/^ //g' | awk 'FNR==NR{a[$1]=$0;next};{ if ($6 in a) {print $1,$2,$3,$4,$5,a[$6]}}' $microb_class - | awk -v min_reads="$min_reads" '{if($2>min_reads) print $1,$2,$3,$4,$5,$6,$7}' | sort -k2,2nr | sed -e 's/ /\t/g' -e 's/^\t//g'> $output_folder/$Ind.microbiome.txt
sum=$(cat $output_folder/$Ind.microbiome.txt | datamash sum 2)
awk -v sum="$sum" '{if(sum>0) print $0,$2/sum}' $output_folder/$Ind.microbiome.txt | sed -e 's/^/'$Ind'\t'$Group'\t'$Population'\t'$Gender'\t'$Age'\t/g' -e 's/ /\t/g' >> $output_folder/$output_filename.summary.abundance.txt
done < <(tail -n+2 $samplelist)

##Make the matrix of abundance
awk '{print $11}' $output_folder/$output_filename.summary.abundance.txt | tail -n+2 | sort | awk '!seen[$0]++' > $output_folder/temp.Genuslist.txt
awk '{print $1,$2,$3,$4,$5}' $output_folder/$output_filename.summary.abundance.txt | tail -n+2 | sort | awk '!seen[$0]++' > $output_folder/temp.Indlist.txt
cat $output_folder/temp.Indlist.txt > $output_folder/temp.Matrix.2.txt

while read K; do
grep "$K" $output_folder/$output_filename.summary.abundance.txt | awk 'FNR==NR{a[$1]=$13;next}{ if ($1 in a) {print $0,a[$1]*100} else {print $0,"0"}}' - $output_folder/temp.Matrix.2.txt > $output_folder/temp.Matrix.txt
cp $output_folder/temp.Matrix.txt $output_folder/temp.Matrix.2.txt
done < $output_folder/temp.Genuslist.txt

sed '1 i\ID Group Population Sex Age' $output_folder/temp.Genuslist.txt | tr "\n" " " | awk '{print $0}' | cat - $output_folder/temp.Matrix.txt > $output_folder/Matrix.microbiome.txt
##sed '1 i\ID Group Population Sex Age' $output_folder/temp.Genuslist.txt | tr "\n" " " | sed 's/ $/#/'| tr "#" "\n" | cat - $output_folder/temp.Matrix.txt > $output_folder/Matrix.microbiome.txt

##Calculate percentage of Class per Group
tail -n+2 $samplelist | awk '{print $5}' | sort | awk '!seen[$1]++' > $output_folder/temp.Grouplist.txt
echo -e "group\tclass\tpercent" > $output_folder/Group.summary.abundance.txt
while read Group; do
N_Group=$(grep "$Group" $output_folder/temp.Indlist.txt | wc -l)
grep $Group $output_folder/$output_filename.summary.abundance.txt | sort -k12,12 | awk -v N_Group="$N_Group" '{a[$12] += $13} END{for (i in a) print $2,i, a[i]/N_Group}'| sort -k3,3gr | sed -e 's/ /\t/g' -e 's/NA/Other/g' >> $output_folder/Group.summary.abundance.txt
done < $output_folder/temp.Grouplist.txt

##Calculate percentage of Class per Individual
echo -e "ID\tgroup\tclass\tpercent" > $output_folder/Individual.summary.classabundance.txt
while read Ind; do
grep $Ind $output_folder/$output_filename.summary.abundance.txt | sort -k12,12 | awk '{a[$12] += $13} END{for (i in a) print $1,$2,i, a[i]}'| sort -k4,4gr | sed -e 's/ /\t/g' -e 's/NA/Other/g' >> $output_folder/Individual.summary.classabundance.txt
done < <(tail -n+2 $samplelist | awk '{print $1}')

################################################
##Alpha and Beta Div
################################################
##PCoA of all samples
awk '!($2=$3=$4=$5="")' $output_folder/Matrix.microbiome.txt | sed 's/ \+/ /g' | sed 's/ID/ /g' > $output_folder/temp.Matrix.txt

echo -e "library(ape)
library(vegan)
##Compute PCoA
matrix <- read.table(\"$output_folder/temp.Matrix.txt\", header=TRUE)
matrixdist <- vegdist(matrix, \"bray\")
res <- pcoa(matrixdist)
sink(file=\"$output_folder/temp.pcoa.txt\")
res$vectors
sink()

##Compute MDS
mds <- metaMDS(matrixdist, maxit=60)
mds_data <- as.data.frame(mds\$points)
sink(file=\"$output_folder/temp.beta.txt\")
mds_data
sink()
" > $output_folder/microbiome.PCoA_Beta.R
Rscript $output_folder/microbiome.PCoA_Beta.R

##Compute Alpha diversity
echo -e "library(ape)
library(vegan)
matrix <- read.table(\"$output_folder/temp.Matrix.txt\", header=TRUE)
alphaInvS <- diversity(matrix,MARGIN = 1,index = \"invsimpson\")
write.csv(alphaInvS, file=\"$output_folder/temp.Matrix.alpha.InvS.txt\")
alphaShan <- diversity(matrix,MARGIN = 1,index = \"shannon\")
write.csv(alphaShan, file=\"$output_folder/temp.Matrix.alpha.Shan.txt\")
alphaSimp <- diversity(matrix,MARGIN = 1,index = \"simpson\")
write.csv(alphaSimp, file=\"$output_folder/temp.Matrix.alpha.Simp.txt\")
" > $output_folder/microbiome.Alpha.R

Rscript $output_folder/microbiome.Alpha.R

##PCoA and NMDS files
paste $output_folder/temp.Matrix.alpha.InvS.txt $output_folder/temp.Matrix.alpha.Simp.txt $output_folder/temp.Matrix.alpha.Shan.txt |  sed -e 's/\"//g' -e 's/,/ /g'  -e 's/\t/ /g' | tail -n+2 | awk '{print $1,$2,$4,$6}' > $output_folder/temp.Matrix.alpha.txt
head -1 $samplelist | sed 's/$/\tmicrobiome_Alpha_InvS\tmicrobiome_Alpha_Shan\tmicrobiome_Alpha_Simp/g' > $output_folder/microbiome.alpha.txt
awk 'FNR==NR{a[$1]=$0;next}{ if ($1 in a) {print a[$1],$2,$3,$4}}' $samplelist $output_folder/temp.Matrix.alpha.txt  >> $output_folder/microbiome.alpha.txt

N_tot=$(tail -n+2 $samplelist | awk '!seen[$1]++' | wc -l)
head -1 $samplelist | sed 's/$/\taxis1\taxis2\taxis3\taxis4/g' > $output_folder/temp.pcoa.2.txt
awk -v N_tot="$N_tot" '/vector/{x=NR+N_tot+1}(NR<=x){print}' $output_folder/temp.pcoa.txt | tail -n+2 | sed 's/ \+/\t/g' | awk 'FNR==NR{a[$1]=$0;next}{ if ($1 in a) {print a[$1],$2,$3,$4,$5}}' $samplelist - | tr " " "\t" >> $output_folder/temp.pcoa.2.txt

head -1 $samplelist | sed 's/$/\tNMDS1\tNMDS2/g' > $output_folder/temp.beta.2.txt
tail -n+2 $output_folder/temp.beta.txt | sed 's/ \+/\t/g' | awk 'FNR==NR{a[$1]=$0;next}{ if ($1 in a) {print a[$1],$2,$3}}' $samplelist - | tr " " "\t" >> $output_folder/temp.beta.2.txt

cp $output_folder/temp.pcoa.2.txt $output_folder/$output_filename.microbiome.PCoA.txt
cp $output_folder/temp.beta.2.txt $output_folder/$output_filename.microbiome.NMDS.txt

######################################
##Plots
######################################
##Pie charts
echo -e "library(ggplot2)
png( \"$output_folder/$output_filename.microbiome.abundance.piechart.png\", height=500, width=1000)
df <- read.table(\"$output_folder/Group.summary.abundance.txt\", header=TRUE)
ggplot(df, aes(x=\"\", y=percent, fill=class)) + geom_col(color = \"white\") + coord_polar(\"y\", start=0) + facet_wrap(vars(group), ncol=4) + theme_void() 
dev.off()
" > $output_folder/microbiome.piechart.R
Rscript $output_folder/microbiome.piechart.R

##Bar Plot
echo "library(ggplot2)
library(RColorBrewer)
png(\"$output_folder/$output_filename.microbiome.abundance.barplot.png\", width=1000, height=400)
admix<-read.table(\"$output_folder/Individual.summary.classabundance.txt\", header=T)
ggplot(admix, aes(x = ID, y = percent, fill = factor(class), width=1)) + geom_bar(stat = \"identity\", position = \"stack\") + facet_grid(~ group, scales=\"free\", space=\"free\") + theme(axis.line = element_blank(),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),axis.text.x=element_blank(), strip.background =element_rect(fill=\"white\"),axis.ticks.x=element_blank(),panel.spacing=unit(0.1,\"lines\"))
dev.off()
" > $output_folder/microbiome.barplot.R
Rscript $output_folder/microbiome.barplot.R

##PCoA plot
echo -e "library(ggplot2)
##Plot PCoA
png( \"$output_folder/$output_filename.microbiome.PCoA.1-2.png\", height=500, width=750)
pcoa <- read.table(\"$output_folder/temp.pcoa.2.txt\", header=TRUE)
ggplot(pcoa, aes(x=axis1, y=axis2)) + geom_point(aes(color = factor(Group),shape = factor(Sex), size=2)) + theme_classic()
dev.off()

png( \"$output_folder/$output_filename.microbiome.PCoA.1-3.png\", height=500, width=750)
pcoa <- read.table(\"$output_folder/temp.pcoa.2.txt\", header=TRUE)
ggplot(pcoa, aes(x=axis1, y=axis3)) + geom_point(aes(color = factor(Group),shape = factor(Sex), size=2)) + theme_classic()
dev.off()

png( \"$output_folder/$output_filename.microbiome.PCoA.1-4.png\", height=500, width=750)
pcoa <- read.table(\"$output_folder/temp.pcoa.2.txt\", header=TRUE)
ggplot(pcoa, aes(x=axis1, y=axis4)) + geom_point(aes(color = factor(Group),shape = factor(Sex), size=2)) + theme_classic()
dev.off()

##Plot Alpha Diversity
png( \"$output_folder/$output_filename.microbiome.Alpha_InvS.Boxplot.png\", height=500, width=750)
alpha_InvS <- read.table(\"$output_folder/microbiome.alpha.txt\", header=TRUE)
ggplot(alpha_InvS, aes(x = Group, y = microbiome_Alpha_InvS)) + geom_violin(aes(fill = Group), trim = FALSE) + theme(axis.text.x = element_text(angle = 90), panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background= element_blank(),axis.line = element_line(colour = \"black\")) + geom_boxplot(width = 0.1)
dev.off()

png( \"$output_folder/$output_filename.microbiome.Alpha_Simp.Boxplot.png\", height=500, width=750)
alpha_Simp <- read.table(\"$output_folder/microbiome.alpha.txt\", header=TRUE)
ggplot(alpha_Simp, aes(x = Group, y = microbiome_Alpha_Simp)) + geom_violin(aes(fill = Group), trim = FALSE) + theme(axis.text.x = element_text(angle = 90), panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background= element_blank(),axis.line = element_line(colour = \"black\")) + geom_boxplot(width = 0.1)
dev.off()

png( \"$output_folder/$output_filename.microbiome.Alpha_Shan.Boxplot.png\", height=500, width=750)
alpha_Shan <- read.table(\"$output_folder/microbiome.alpha.txt\", header=TRUE)
ggplot(alpha_Shan, aes(x = Group, y = microbiome_Alpha_Shan)) + geom_violin(aes(fill = Group), trim = FALSE) + theme(axis.text.x = element_text(angle = 90), panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background= element_blank(),axis.line = element_line(colour = \"black\")) + geom_boxplot(width = 0.1)
dev.off()

##MDS Beta Diversity
png( \"$output_folder/$output_filename.microbiome.Beta_MDS.png\", height=500, width=750)
mds <- read.table(\"$output_folder/temp.beta.2.txt\", header=TRUE)
ggplot(mds, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(color = factor(Group),shape = factor(Sex), size=2)) + theme_classic()
dev.off()
" > $output_folder/microbiome.PCoA_Alpha.R

Rscript $output_folder/microbiome.PCoA_Alpha.R

###############################
##PERMANOVA and T TEST
###############################

echo -e "library(stats)
matrix <- read.table(\"$output_folder/microbiome.alpha.txt\", header=TRUE)
sink(file=\"$output_folder/$output_filename.microbiome.alpha.Ttest.txt\")
t.test(microbiome_Alpha_InvS ~ Group, data = matrix)
t.test(microbiome_Alpha_Simp ~ Group, data = matrix)
t.test(microbiome_Alpha_Shan ~ Group, data = matrix)
sink()
" > $output_folder/betel.alpha.Ttest.R
Rscript $output_folder/betel.alpha.Ttest.R

sed '1 i\  Group Population Sex Age' $output_folder/temp.Indlist.txt > $output_folder/temp.ID
echo -e "library(vegan)
matrix <- read.table(\"$output_folder/temp.Matrix.txt\", header=TRUE)
covar <- read.table(\"$output_folder/temp.ID\", header=TRUE)
sink(file=\"$output_folder/$output_filename.microbiome.beta.permanova.txt\")
adonis2(matrix ~ Group,data = covar,method = \"bray\",by = \"margin\",permutations = 9999)
sink()
" > $output_folder/betel.beta.permanova.R
Rscript $output_folder/betel.beta.permanova.R

##################################
##LEFSE
##################################
##Integrate phylogeny to matrix
awk '{for (i=1; i<=NF; i++){ a[NR,i] = $i}} NF>p { p = NF } END { for(j=1; j<=p; j++) {str=a[1,j]; for(i=2; i<=NR; i++){str=str" "a[i,j];} print str}}' $output_folder/Matrix.microbiome.txt | head -2 > $output_folder/temp.1.$output_filename.Matrix.All.txt
cut -d" " -f6- $output_folder/Matrix.microbiome.txt  | awk '{for (i=1; i<=NF; i++){ a[NR,i] = $i}} NF>p { p = NF } END { for(j=1; j<=p; j++) {str=a[1,j]; for(i=2; i<=NR; i++){str=str" "a[i,j];} print str}}' | awk 'FNR==NR{a[$1]=$3;next}{ if ($1 in a) {print a[$1],$0} else  {print $0}}' $microb_class - | awk '!($2="")'  | sed 's/  / /g' | cat $output_folder/temp.1.$output_filename.Matrix.All.txt - | tr " " "\t" > $output_folder/$output_filename.Matrix.All.txt

##LefSE
format_input.py $output_folder/$output_filename.Matrix.All.txt $output_folder/temp.$output_filename.Matrix.All.in -c 2 -s -1 -u 1 -o $LDA_thres
run_lefse.py $output_folder/temp.$output_filename.Matrix.All.in $output_folder/temp.$output_filename.Matrix.All.res
plot_res.py $output_folder/temp.$output_filename.Matrix.All.res $output_folder/$output_filename.Matrix.All.png
plot_cladogram.py $output_folder/temp.$output_filename.Matrix.All.res $output_folder/$output_filename.Matrix.All.cladogram.pdf --format pdf
if [ ! -d "$output_folder/lefse_plot" ]; then
mkdir $output_folder/lefse_plot
fi
plot_features.py $output_folder/temp.$output_filename.Matrix.All.in $output_folder/temp.$output_filename.Matrix.All.res $output_folder/lefse_plot/ | tee | sed 's/^.*\.//g' | awk 'FNR==NR{a[$1]=$1;next}{ if ($1 in a) {print a[$1]}}' $microb_class - > $output_folder/$output_filename.lefse.Genus.list

##Make the matrix of abundance per specie from Genus significant with LDA
echo -e "ID\tgroup\tpop\tgender\tage\tpercent_read\treads1\treads2\ttype\ttype_nb\tgenus\tclass\tabundance" > $output_folder/$output_filename.lefse_specie.abundance.txt
while read Ind2; do
Ind=$(echo $Ind2 | awk '{print $1}')
Population=$(echo $Ind2 | awk '{print $2}')
Gender=$(echo $Ind2 | awk '{print $3}')
Group=$(echo $Ind2 | awk '{print $5}')
Age=$(echo $Ind2 | awk '{print $4}')
grep $'\tS\t' $input_folder/$Ind/"$Ind"_report.microbiome.txt | sed 's/ \+/ /g' | sed -e 's/^ //g' | awk 'FNR==NR{a[$1]=$0;next};{ if ($6 in a) {print $0}}' $output_folder/$output_filename.lefse.Genus.list - | awk -v min_reads="$min_reads" '{if($2>min_reads) print $0}' | sort -k2,2nr | sed -e 's/ /\t/g' -e 's/^\t//g' > $output_folder/temp.2.$Ind.specie_lefse.txt
cut -f7- $output_folder/temp.2.$Ind.specie_lefse.txt | tr "\t" "_" > $output_folder/temp.$Ind.specie_name.txt
cut -f1-5 $output_folder/temp.2.$Ind.specie_lefse.txt | paste - $output_folder/temp.$Ind.specie_name.txt > $output_folder/temp.$Ind.specie_lefse.txt
sum=$(sed -e 's/\t/ /g' -e 's/ \+/ /g' $input_folder/$Ind/"$Ind"_report.microbiome.txt | grep "D .*Bacteria" | cut -d" " -f3)
awk -v sum="$sum" '{if(sum>0) print $0,$2/sum}' $output_folder/temp.$Ind.specie_lefse.txt | sed -e 's/^/'$Ind'\t'$Group'\t'$Population'\t'$Gender'\t'$Age'\t/g' -e 's/ /\t/g' >> $output_folder/$output_filename.specie_lefse.abundance.txt
done < <(tail -n+2 $samplelist)

awk '{print $11}' $output_folder/$output_filename.specie_lefse.abundance.txt | tail -n+2 | sort | awk '!seen[$0]++' > $output_folder/temp.specie_lefselist.txt
awk '{print $1,$2,$3,$4,$5}' $output_folder/$output_filename.specie_lefse.abundance.txt | tail -n+2 | sort | awk '!seen[$0]++' > $output_folder/temp.specie_lefse.Indlist.txt
cat $output_folder/temp.specie_lefse.Indlist.txt > $output_folder/temp.specie_lefse.Matrix.2.txt

while read K; do
grep "$K" $output_folder/$output_filename.specie_lefse.abundance.txt | awk 'FNR==NR{a[$1]=$12;next}{ if ($1 in a) {print $0,a[$1]*100} else {print $0,"0"}}' - $output_folder/temp.specie_lefse.Matrix.2.txt > $output_folder/temp.specie_lefse.Matrix.txt
cp $output_folder/temp.specie_lefse.Matrix.txt $output_folder/temp.specie_lefse.Matrix.2.txt
done < $output_folder/temp.specie_lefselist.txt
sed '1 i\ID Group Population Sex Age' $output_folder/temp.specie_lefselist.txt | tr "\n" " " | awk '{print $0}' | cat - $output_folder/temp.specie_lefse.Matrix.txt > $output_folder/Matrix.specie_lefselist.microbiome.txt

##Find significant Species
> $output_folder/Specie_lefselist.microbiome.Ttest.txt
while read K; do
echo -e "library(stats)
matrix <- read.table(\"$output_folder/Matrix.specie_lefselist.microbiome.txt\", header=TRUE)
sink(file=\"$output_folder/temp.specie_lefselist.microbiome.txt\")
t.test($K ~ Group, data = matrix)
sink()
" > $output_folder/specie_lefselist.Ttest.R
Rscript $output_folder/specie_lefselist.Ttest.R
grep "p-value" $output_folder/temp.specie_lefselist.microbiome.txt | awk -v specie="$K" '{if ($9<0.05) print specie,$0}' >> $output_folder/Specie_lefselist.microbiome.Ttest.txt
done < $output_folder/temp.specie_lefselist.txt

##Violin plot of each significant species
if [ ! -d "$output_folder/Specie_ViolinPlot" ]; then
mkdir $output_folder/Specie_ViolinPlot
fi

while read K; do
echo -e "library(ggplot2)
png( \"$output_folder/Specie_ViolinPlot/$output_filename.microbiome.$K.LDA_Specie_ViolinPlot.png\", height=500, width=750)
LEFSE <- read.table(\"$output_folder/Matrix.specie_lefselist.microbiome.txt\", header=TRUE)
ggplot(LEFSE, aes(x = Group, y = $K)) + geom_violin(aes(fill = Group), trim = FALSE) + theme(axis.text.x = element_text(angle = 90), panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background= element_blank(),axis.line = element_line(colour = \"black\")) + geom_boxplot(width = 0.1)
dev.off()
" > $output_folder/microbiome.LDA_Specie_ViolinPlot.R
Rscript $output_folder/microbiome.LDA_Specie_ViolinPlot.R
done < <(awk '{print $1}' $output_folder/Specie_lefselist.microbiome.Ttest.txt)

