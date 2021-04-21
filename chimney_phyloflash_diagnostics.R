# Plan ----
# 1. 
# Run phyloflash_compare.pl on the 2020 and 2019 sequencing projects together
# data for these found in /export/dahlefs/work/Emily

# 2. 
# Get duplication, human decontam file sizes. add to metadata.txt
# found in the multiQC from the project, grab from cmd line

# 3. 
# Look at typical kit & extraction microbial contaminants from paper Glassing_etal_2016:
# https://link.springer.com/article/10.1186/s13099-016-0103-7
# supplementary table, saved as: 13099_2016_103_MOESM1_ESM.xlsx 

# 4. 
# Rarify and / or normalize phyloflash hits

# 5.
# Write code to fetch # reads + sample % of each contaminants / sample
# Likely want to do this at the genera level.

# 6. 
# Make histograms of suspected contaminants / sample in grid across each gradient

# 7.
# Make histograms of Eukaryotic types / sample in grid across each gradient
# Maybe divide this into unknown, Fungi, Metazoa, Other.

# Environment setup  ----
library("reshape")
library("vegan")
library("ape")
library("janitor")

setwd("Desktop/Metagenomes_AMOR_2020/")
# Data Import ----

# NTU table long form to wide
FR <- read.delim(file='Metagenomes_AMOR_2020/ALL_withLokintu_table3.csv', sep='\t', header=FALSE)
df2 <- FR

# metadata
env_char <- read.csv("Metagenomes_AMOR_2020/SampleMasterAchim.csv", sep="\t", header=T, row.names=1)
humans <- read.delim(file='Metagenomes_AMOR_2020/human_clean_ratios2.txt', sep='\t', header=T, row.names = 1)

# contaminant genera
x <- scan("Metagenomes_AMOR_2020/ContamEval/contam_genera.txt", what="", sep="\n")
z <- scan("Metagenomes_AMOR_2020/myCommon_Contams.csv", what="", sep="\n")
# Extract the first vector element and set it as the list element name
names(x) <- sapply(x, `[[`, 1)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
contams <- lapply(x, `[`, -1)
#y <- lapply(y, function(x) x[-1]) # same as above

# Pre-processing ----

# Give column headers
names(df2)<- c("Taxonomy", "SampleNames","Count")

# Convert to wide format
df1 <- reshape(df2, idvar = "Taxonomy", timevar = "SampleNames", direction = "wide")

# Remove prefix added from wide-format transform
# Give column headers
colnames(df1) <- sub("Emily/MetaGen_chimneys_2019/PHYLOFLASH2/", "", colnames(df1))
colnames(df1) <- sub("Count.", "", colnames(df1))

# Change NAs to zeros
df1[is.na(df1)] <- as.numeric(0)

# Transpose the table, remove the first column before doing so.
spe <- t(df1[,-1])
names(spe) <- df1$Taxonomy

# Reorder the abundance table by sample names
spe <- spe[order(rownames(spe)),]
# Reorder environmental table by sample names
env_char <- env_char[order(rownames(env_char)),]
humans <- humans[order(rownames(humans)),]
env_char$humanRatios <- humans$HumanRatio

# Test if you have samplename matches between meta and NTU tables
# should just return TRUE if it is.
test <-  all.equal(rownames(spe), rownames(env_char), rownames(humans))

colnames(spe) <- df1$Taxonomy

# Normalizing ----
spe.tot <- decostand(spe, "total")

# ANCOM normalize instead? logratio transform, addresses sparse dataA ----
  # Had to do some dependency installs and download file and then "source on save"
# It's not a package just a function ...*eye roll*

spe.logratio <- 
# Functions ----

# Get list of groups over 1% in a sample
return_over_1perc <- function(spe.tot, smp_forgrep) {
  r <- grep(smp_forgrep, row.names(spe.tot))
  sort(spe.tot[r,which(spe.tot[r,]>0.01)]) 
}

# Return suspected contaminants only 
find_contams <- function(spe.tot, contams) {
  found_contams <- as.data.frame(spe.tot[,FALSE])
for (contam_name in names(contams)) {
  c <- grep(contam_name,colnames(spe.tot))
  ifelse(length(c) > 1, found_contams[,contam_name] <- rowSums(spe.tot[,c]), found_contams[,contam_name] <- spe.tot[,c])
}
  found_contams <<- found_contams   # save dataframe to the workspace
  return(found_contams)}    # Also return the dataframe to the command line

# Return dataframe of contaminants ----
found_contams <- find_contams(spe.tot, contams)

# Charts ----
# Overall, each contaminant
library(Hmisc)
found_contams_without_zeros <- found_contams
found_contams_without_zeros[found_contams_without_zeros== 0] <- NA
hist.data.frame(found_contams_without_zeros, mtitl = "Genera Frequency")


# Transform to long form for plotting in ggplot2
df <- melt(t(found_contams) ,  id.vars = colnames(found_contams))
new_df <- na.omit(df, c("X1", "X2", "value")) 


# Plot all potential contaminants line graph
install.packages("viridis")
library(viridis)
col_vector <- viridis_pal(option = "D")(180) 

env_char$SampleNames <- rownames(env_char)
new_df2 <- merge(x = new_df, y = env_char, by.x = "X2", by.y= "SampleNames", all = TRUE)

library(plotly)
p <- ggplot(new_df2, aes(x = X2, y = value, color = X1, group = X1)) +
  geom_point() +
  geom_line() +
  facet_wrap(Site ~ SampleType, scales = "free_x")+
  scale_color_manual(values= col_vector) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),
        strip.background = element_rect(colour="white", fill="white"))

g<- ggplotly(p)
# Save it locally
library(htmlwidgets)
htmlwidgets::saveWidget(g, "suspect_contam_genera.html")
dev.off()

# Plot of Human / Total Reads bar charts

f <- ggplot(new_df2) +
  geom_col(aes(x=X2, y = humanRatios), size = 1, color = "grey") +
  facet_wrap(Site ~ SampleType, scales = "free_x")+
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),
        strip.background = element_rect(colour="white", fill="white"))

f <- f +         scale_x_discrete(name= "Sample Names") +
  scale_y_continuous(name="Human / Total QC Reads", limits=c(0, 1))

f <- ggplotly(f)

htmlwidgets::saveWidget(f, "human_contam.html")
dev.off()


# Junk ----
par(mar=c(10,4,4,4))
par(mfrow=c(1,1))
barplot(as.matrix(found_contams),beside=TRUE, las=2, cex.axis = 1, cex.names= 0.3)

dev.off()


f <- ggplot(new_df2) +
  geom_point(aes(x=X2, y =value, color = X1, group = X1)) +
  geom_line(aes(x=X2, y =value, color = X1, group = X1)) +
  geom_hline(yintercept=.01, linetype='dotted', col = 'red')+
  annotate("text", x = "7_17ROV19_HD11", y = .01, label = "Cutoff 1 %?", vjust = -0.5)+
  geom_col(aes(x=X2, y = humanRatios), size = 1, color = "grey", fill = "white") +
  scale_y_continuous(
    
    # Features of the first axis
    name = "% abundance contamination",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*7, name="human reads ratio")
  )+
  #  facet_wrap(Site ~ SampleType, scales = "free_x")+
  scale_color_manual(values= col_vector) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),
        strip.background = element_rect(colour="white", fill="white"))

f