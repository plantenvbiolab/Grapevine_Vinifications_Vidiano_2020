#Vidiano 2020 Anticancer#
MUST_Fungi_Vidiano_2020_Merged

MUST_Fungi_Vidiano_2020_Merged_Genus

#MERGED REPLICATE##
MUST_Fungi_Vidiano_2020_Merged_Vin <- merge_samples(MUST_Fungi_Vidiano_2020_Merged, "vinification")

MUST_Fungi_Vidiano_2020_Merged_Vin <- prune_taxa(taxa_sums(MUST_Fungi_Vidiano_2020_Merged_Vin)>0,MUST_Fungi_Vidiano_2020_Merged_Vin)

View(data.frame(sample_data(MUST_Fungi_Vidiano_2020_Merged_Vin)))

write.table(data.frame(sample_data(MUST_Fungi_Vidiano_2020_Merged_Vin)), file="MUST_Fungi_Vidiano_2020_Merged_Vin.txt", quote = F,col.names = NA, sep="\t")

SampleDataNew67 <- read.table("MUST_Fungi_Vidiano_2020_Merged_Vin.txt", header=T,sep = "\t",row.names = 1)

sample_data(MUST_Fungi_Vidiano_2020_Merged_Vin) <- SampleDataNew67

MUST_Fungi_Vidiano_2020_Merged_Vin

View(data.frame(sample_data(MUST_Fungi_Vidiano_2020_Merged_Vin)))

MUST_Fungi_Vidiano_2020_Merged_Vin <- prune_taxa(taxa_sums(MUST_Fungi_Vidiano_2020_Merged_Vin)>0,MUST_Fungi_Vidiano_2020_Merged_Vin)

MUST_Fungi_Vidiano_2020_Merged_Vin_100 <- transform_sample_counts(MUST_Fungi_Vidiano_2020_Merged_Vin, function(OTU) 100*OTU/sum(OTU)) 

ord.nmds.bray1000 <- ordinate(MUST_Fungi_Vidiano_2020_Merged_Vin_100, method="CCA", distance="bray")

plot_ordination(MUST_Fungi_Vidiano_2020_Merged_Vin_100, ord.nmds.bray1000, color="vinification", shape ="", label = "", title=paste("Vidiano 2020")) + geom_point(size = 4)

VolatileVidiano2020 <- read.table("VolatileVidiano2020.txt", header=T,sep = "\t",row.names = 1)

fitVidiano <- envfit(ord.nmds.bray1000, VolatileVidiano2020[,])
#####


###Merged Replicate Genus#
MUST_Fungi_Vidiano_2020_Merged_Genus

#MERGED REPLICATE##
MUST_Fungi_Vidiano_2020_Merged_Genus_Vin <- merge_samples(MUST_Fungi_Vidiano_2020_Merged_Genus, "vinification")

MUST_Fungi_Vidiano_2020_Merged_Genus_Vin <- prune_taxa(taxa_sums(MUST_Fungi_Vidiano_2020_Merged_Genus_Vin)>0,MUST_Fungi_Vidiano_2020_Merged_Genus_Vin)

View(data.frame(sample_data(MUST_Fungi_Vidiano_2020_Merged_Genus_Vin)))

write.table(data.frame(sample_data(MUST_Fungi_Vidiano_2020_Merged_Genus_Vin)), file="MUST_Fungi_Vidiano_2020_Merged_Genus_Vin.txt", quote = F,col.names = NA, sep="\t")

SampleDataNew67 <- read.table("MUST_Fungi_Vidiano_2020_Merged_Genus_Vin.txt", header=T,sep = "\t",row.names = 1)

sample_data(MUST_Fungi_Vidiano_2020_Merged_Genus_Vin) <- SampleDataNew67

MUST_Fungi_Vidiano_2020_Merged_Genus_Vin

View(data.frame(sample_data(MUST_Fungi_Vidiano_2020_Merged_Genus_Vin)))

MUST_Fungi_Vidiano_2020_Merged_Genus_Vin <- prune_taxa(taxa_sums(MUST_Fungi_Vidiano_2020_Merged_Genus_Vin)>0,MUST_Fungi_Vidiano_2020_Merged_Genus_Vin)

MUST_Fungi_Vidiano_2020_Merged_Genus_Vin_100 <- transform_sample_counts(MUST_Fungi_Vidiano_2020_Merged_Genus_Vin, function(OTU) 100*OTU/sum(OTU)) 

ord.nmds.bray1000 <- ordinate(MUST_Fungi_Vidiano_2020_Merged_Genus_Vin_100, method="CCA", distance="bray")

plot_ordination(MUST_Fungi_Vidiano_2020_Merged_Genus_Vin_100, ord.nmds.bray1000, color="vinification", shape ="", label = "", title=paste("Vidiano 2020")) + geom_point(size = 4)

VolatileVidiano2020 <- read.table("VolatileVidiano2020.txt", header=T,sep = "\t",row.names = 1)

fitVidiano <- envfit(ord.nmds.bray1000, VolatileVidiano2020[,])


####______________________TRIPLOT_______________________________###
# prepare the names for the to be plotted taxa and save them in a new column
mytax_tbl <- data.frame(tax_table(MUST_Fungi_Vidiano_2020_Merged_Genus_Vin))
# the following command pastes the Phylum column with the Genus one
# !!! make sure to get rid of the NAs with the for loop as we discussed before !!!

# For ITS - Remove letter from taxonomy

mytax_tbl$forplt <- for (i in c(1:nrow(mytax_tbl))) {
  for(j in c(1:ncol(mytax_tbl))) {
    mytax_tbl[i,j] <- gsub("[a-z]__","",mytax_tbl[i,j])
  }
}


# prepare the plotting colours
myplotcols <- RColorBrewer::brewer.pal(n = 4, name = 'RdBu')
# prep the plot points
mynmdssit <- scores(ord.nmds.bray1000, display = "site")

# prep the species points
mynmdsspe <- scores(ord.nmds.bray1000, display = "species")

# prep the variable of interest (e.g. vinification) has to be a factor
my_sel_var <- factor(MUST_Fungi_Vidiano_2020_Merged_Vin_100@sam_data$vinification)
# convert it also into numbers in orde to use it for colours etc
myterr_sel <- as.numeric(my_sel_var)
# get the arrow data that came out of the fitting of the environmental variables
arrowdata <- data.frame(fitVidiano$vectors$arrows)

# start the graphics device
cairo_pdf("Fungi Vidiano 2020 Volatile.pdf", height = 7, width = 7)

# draw and empty plot that you will start to populate with points, arrows, ellipses, names etc.
plot(mynmdssit, frame = F, cex = 0, pch = 21, xlim = c(min(1.5*mynmdssit[,1]),max(1.5*mynmdssit[,1])))
# first add the ellipses in order tto keep them in the background and prevent them to obstruct information
vegan::ordiellipse(mynmdssit, groups = my_sel_var, kind = "ehull", lty = 2, lwd=1)
# add the sample points and give them colors one for each treatment
points(mynmdssit, bg = myplotcols[myterr_sel], pch = 21, cex = 1.5
       ##gia keimeno sta simeia##
       , text(mynmdssit, labels = row.names(mynmdssit), font = 6, col = rgb(55,55,55, max = 255, alpha = 100)))

# write the stress value in the bottom left of the plot (you need to change the global "adj" parameter from 0.5 to 0 in order to add the text in the left instead of the middle (if you gave the value of 1 it would add the text in the rightmost side... I provide some text as example))
par(adj = 0)
title(sub = paste("stress ", round(ord.nmds.bray1000$stress,2), sep = ""), cex.sub = 1.2)

par(adj = 1)
title(sub = "Fungi Vidiano Volatile")
par(adj = .5)
# add the env parameter arrows (in the rgb colouring mode we have the red, green and blue values from 0 to 255 and also the alpha value or transparency which also ranges between 0 and 255)

arrows(0,0,0.3*arrowdata[,1] , 0.3*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)
# also the labeling text
text(0.5*arrowdata[,1], 0.5*arrowdata[,2], labels = row.names(arrowdata), cex = 1.3, font = 2, col = rgb(20,20,20, max = 255, alpha = 255)) 

##simantiko species arrows and labels##
-----------------------------------------------------------------
  ##Epilego ta ASVs pou thelo ite aritmitika ite onomastika
  
  ##Aritmitika##
  arrows(0,0,mynmdsspe[1:22,1] , mynmdsspe[1:22,2], angle = 25, length = 0.15, col = rgb(40,40,40, max = 255, alpha = 100))


plotrix::thigmophobe.labels(1.2*mynmdsspe[1:20,1], 1.2*mynmdsspe[1:20,2], labels = mytax_tbl[row.names(mynmdsspe)[1:20],"forplt"], cex = .6, font = 2, col = rgb(120,120,120, max = 255, alpha =200))


##Onomastika anti gia diaforika afthona##
arrows(0,0,0.5*mynmdsspe[c("ASV0001","ASV0005","ASV0012","ASV0209","ASV0230","ASV2579","ASV5102"),1] , 0.5*mynmdsspe[c("ASV0001","ASV0005","ASV0012","ASV0209","ASV0230","ASV2579","ASV5102"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))


("ASV0001","ASV0005","ASV0012","ASV0209","ASV0230","ASV2579","ASV5102")

("Saccharomyces_sp.*","Lachancea_sp.**","Torulaspora_sp.*","Candida_sp*","Metschnikowia_sp.*","Saccharomycetales_sp.**","Saccharomycetales_sedis*")

##LABELS Short-cut##
a <- c("Saccharomyces_sp.*","Lachancea_sp.**","Torulaspora_sp.*","Candida_sp*","Wickerhamomyces_sp.*","Saccharomycetales_sp.**","Saccharomycetales_sedis*")


install.packages("plotrix")

plotrix::thigmophobe.labels(1*mynmdsspe[c("ASV0001","ASV0005","ASV0012","ASV0209","ASV0230","ASV2579","ASV5102"),1], 1*mynmdsspe[c("ASV0001","ASV0005","ASV0012","ASV0209","ASV0230","ASV2579","ASV5102"),2], labels = a, cex = .6, font = 2, col = rgb(120,120,120,max = 255, alpha =200)) 

#simantiko# plotrix::thigmophobe.labels(1.2*mynmdsspe[1:15,1], 1.2*mynmdsspe[1:15,2], labels = mytax_tbl[row.names(mynmdsspe)[1:15],"forplt"], cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
# the following prints the legend
graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()