###Merged Replicates
MUST_Bacteria_Vidiano_2020_Merged

MUST_Bacteria_Vidiano_2020_Merged_Genus

MUST_Bacteria_Vidiano_2020_Merged_Vin <- merge_samples(MUST_Bacteria_Vidiano_2020_Merged_Genus, "vinification")

MUST_Bacteria_Vidiano_2020_Merged_Vin <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_Merged_Vin)>0,MUST_Bacteria_Vidiano_2020_Merged_Vin)

View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_Vin)))

write.table(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_Vin)), file="MUST_Bacteria_Vidiano_2020_Merged_Vin.txt", quote = F,col.names = NA, sep="\t")

SampleDataNew67 <- read.table("MUST_Bacteria_Vidiano_2020_Merged_Vin.txt", header=T,sep = "\t",row.names = 1)

sample_data(MUST_Bacteria_Vidiano_2020_Merged_Vin) <- SampleDataNew67

MUST_Bacteria_Vidiano_2020_Merged_Vin

View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_Vin)))

MUST_Bacteria_Vidiano_2020_Merged_Vin <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_Merged_Vin)>0,MUST_Bacteria_Vidiano_2020_Merged_Vin)

MUST_Bacteria_Vidiano_2020_Merged_Vin_100 <- transform_sample_counts(MUST_Bacteria_Vidiano_2020_Merged_Vin, function(OTU) 100*OTU/sum(OTU)) 

ord.nmds.bray1000 <- ordinate(MUST_Bacteria_Vidiano_2020_Merged_Vin, method="CCA", distance="bray")

plot_ordination(MUST_Bacteria_Vidiano_2020_Merged_Vin_100, ord.nmds.bray1000, color="vinification", shape ="", label = "", title=paste("VIDIANO Bacteria 2020")) + geom_point(size = 4) + scale_color_manual(values = mycols) + theme_classic()

plot_ordination(MUST_Bacteria_Vidiano_2020_Merged_Vin_100, ord.nmds.bray1000, color="stage", shape ="", label = "", title=paste("VIDIANO Bacteria 2020")) + geom_point(size = 4) + scale_color_manual(values = mycols) + theme_classic()


VolatileVidiano2020 <- read.table("VolatileVidiano2020.txt", header=T,sep = "\t",row.names = 1)

fitVidiano <- envfit(ord.nmds.bray1000, VolatileVidiano2020[,])


####______________________TRIPLOT_______________________________###
# prepare the names for the to be plotted taxa and save them in a new column
mytax_tbl <- data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_Vin_100))
# the following command pastes the Phylum column with the Genus one
# !!! make sure to get rid of the NAs with the for loop as we discussed before !!!

# For ITS - Remove letter from taxonomy

mytax_tbl$forplt <- for (i in c(1:nrow(mytax_tbl))) {
  for(j in c(1:ncol(mytax_tbl))) {
    mytax_tbl[i,j] <- gsub("[a-z]__","",mytax_tbl[i,j])
  }
}


# prepare the plotting colours
myplotcols <- RColorBrewer::brewer.pal(n = 6, name = 'RdBu')
# prep the plot points
mynmdssit <- scores(ord.nmds.bray1000, display = "site")

# prep the species points
mynmdsspe <- scores(ord.nmds.bray1000, display = "species")

# prep the variable of interest (e.g. vinification) has to be a factor
my_sel_var <- factor(MUST_Bacteria_Vidiano_2020_Merged_Vin_100@sam_data$vinification)
# convert it also into numbers in orde to use it for colours etc
myterr_sel <- as.numeric(my_sel_var)
# get the arrow data that came out of the fitting of the environmental variables
arrowdata <- data.frame(fitVidiano$vectors$arrows)

# start the graphics device
cairo_pdf("Bacteria Vidiano 2020 Volatile15.pdf", height = 7, width = 7)

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
title(sub = "Bacteria 22 Vidiano Volatile2")
par(adj = .5)
# add the env parameter arrows (in the rgb colouring mode we have the red, green and blue values from 0 to 255 and also the alpha value or transparency which also ranges between 0 and 255)

#big
arrows(0,0,1.5*arrowdata[,1] , 1.5*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)
# also the labeling text
text(0.5*arrowdata[,1], 0.5*arrowdata[,2], labels = row.names(arrowdata), cex = 1.3, font = 2, col = rgb(20,20,20, max = 255, alpha = 255)) 

##simantiko species arrows and labels##
-----------------------------------------------------------------
  ##Epilego ta ASVs pou thelo ite aritmitika ite onomastika
  
  ##Aritmitika##
  arrows(0,0,mynmdsspe[1:22,1] , mynmdsspe[1:22,2], angle = 25, length = 0.15, col = rgb(40,40,40, max = 255, alpha = 100))


plotrix::thigmophobe.labels(1.2*mynmdsspe[1:20,1], 1.2*mynmdsspe[1:20,2], labels = mytax_tbl[row.names(mynmdsspe)[1:20],"forplt"], cex = .6, font = 2, col = rgb(120,120,120, max = 255, alpha =200))


##Onomastika anti gia diaforika afthona##
arrows(0,0,3.5*mynmdsspe[c("ASV00051","ASV00061","ASV00062","ASV00136","ASV00202","ASV00210","ASV00216","ASV00225","ASV00237","ASV00286","ASV00347","ASV00355"),1] , 3.5*mynmdsspe[c("ASV00051","ASV00061","ASV00062","ASV00136","ASV00202","ASV00210","ASV00216","ASV00225","ASV00237","ASV00286","ASV00347","ASV00355"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))


("ASV00051","ASV00061","ASV00062","ASV00136","ASV00202","ASV00210","ASV00216","ASV00225","ASV00237","ASV00286","ASV00347","ASV00355")

a <- c("Gluconobacter_sp.","Komagataeibacter_sp.","Burkholderia_sp.","Bacillales_sp.","Comamonadaceae_sp.","Asinibacterium_sp.","Sphingomonas_sp.","Xanthobacteraceae_sp.","Acetobacter_sp.","Sphingomonadaceae_sp.","Magnetospirillaceae_sp.","Obscuribacteraceae_sp.")

##LABELS Short-cut##
("","","","","","","")


install.packages("plotrix")

plotrix::thigmophobe.labels(1*mynmdsspe[c("ASV00051","ASV00061","ASV00062","ASV00136","ASV00202","ASV00210","ASV00216","ASV00225","ASV00237","ASV00286","ASV00347","ASV00355"),1], 1*mynmdsspe[c("ASV00051","ASV00061","ASV00062","ASV00136","ASV00202","ASV00210","ASV00216","ASV00225","ASV00237","ASV00286","ASV00347","ASV00355"),2], labels = a, cex = .6, font = 2, col = rgb(120,120,120,max = 255, alpha =200)) 

#simantiko# plotrix::thigmophobe.labels(1.2*mynmdsspe[1:15,1], 1.2*mynmdsspe[1:15,2], labels = mytax_tbl[row.names(mynmdsspe)[1:15],"forplt"], cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
# the following prints the legend
graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()

======================================================================
  
  MUST_Bacteria_Vidiano_2020_Merged_Vin_100  

myTaxa11_MUST_Bacteria_Vidiano_2020_Merged_Vin_100 <- names(sort(taxa_sums(MUST_Bacteria_Vidiano_2020_Merged_Vin_100), decreasing = TRUE)[1:12]) 

Top11_MUST_Bacteria_Vidiano_2020_Merged_Vin_100 <- prune_taxa(myTaxa11_MUST_Bacteria_Vidiano_2020_Merged_Vin_100, MUST_Bacteria_Vidiano_2020_Merged_Vin_100)


Top11_MUST_Fungi_Vidiano_2020_Merged_Vin_100

View(data.frame(otu_table(Top11_MUST_Bacteria_Vidiano_2020_Merged_Vin_100)))

View(data.frame(tax_table(Top11_MUST_Bacteria_Vidiano_2020_Merged_Vin_100)))



















##1st Approach ##stages replicates (2) are allready merged ## i will seperate stages (S1, S2, S3, S4, S5)## 

View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged)))

MUST_Bacteria_Vidiano_2020_Merged_Genus <- tax_glom(MUST_Bacteria_Vidiano_2020_Merged, taxrank = "Genus")

View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_Genus)))

###SEPERATE MERGED STAGES
MUST_Bacteria_Vidiano_2020_Merged_S1 <- subset_samples(MUST_Bacteria_Vidiano_2020_Merged_Genus, stage %in% (c("S1")))
MUST_Bacteria_Vidiano_2020_Merged_S1 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_Merged_S1)>0,MUST_Bacteria_Vidiano_2020_Merged_S1)
View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_S1)))


MUST_Bacteria_Vidiano_2020_Merged_S2 <- subset_samples(MUST_Bacteria_Vidiano_2020_Merged_Genus, stage %in% (c("S2")))
MUST_Bacteria_Vidiano_2020_Merged_S2 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_Merged_S2)>0,MUST_Bacteria_Vidiano_2020_Merged_S2)
View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_S2)))


MUST_Bacteria_Vidiano_2020_Merged_S3 <- subset_samples(MUST_Bacteria_Vidiano_2020_Merged_Genus, stage %in% (c("S3")))
MUST_Bacteria_Vidiano_2020_Merged_S3 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_Merged_S3)>0,MUST_Bacteria_Vidiano_2020_Merged_S3)
View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_S3)))


MUST_Bacteria_Vidiano_2020_Merged_S4 <- subset_samples(MUST_Bacteria_Vidiano_2020_Merged_Genus, stage %in% (c("S4")))
MUST_Bacteria_Vidiano_2020_Merged_S4 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_Merged_S4)>0,MUST_Bacteria_Vidiano_2020_Merged_S4)
View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_S4)))


MUST_Bacteria_Vidiano_2020_Merged_S5 <- subset_samples(MUST_Bacteria_Vidiano_2020_Merged_Genus, stage %in% (c("S5")))
MUST_Bacteria_Vidiano_2020_Merged_S5 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_Merged_S5)>0,MUST_Bacteria_Vidiano_2020_Merged_S5)
View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_S5)))

##SEPERATED AND CLEAN

##S1 Triplot
MUST_Bacteria_Vidiano_2020_Merged_S1
MUST_Bacteria_Vidiano_2020_Merged_S1_100 <- transform_sample_counts(MUST_Bacteria_Vidiano_2020_Merged_S1, function(OTU) 100*OTU/sum(OTU))

ord.nmds.brayS1 <- ordinate(MUST_Bacteria_Vidiano_2020_Merged_S1_100, method="CCA", distance="bray")

plot_ordination(MUST_Bacteria_Vidiano_2020_Merged_S1_100, ord.nmds.brayS1, color="vinification", shape ="", label = "", title=paste("S1")) + geom_point(size = 4)

View(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_S1_100)))
View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_S1_100)))
View(data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S1_100)))

Top_12_S1 <- names(sort(taxa_sums(MUST_Bacteria_Vidiano_2020_Merged_S1_100), decreasing = TRUE)[1:12])

##S1_Triplot_Start
VolatileVidiano2020_S1 <- read.table("VolatileVidiano2020_S1.txt", header=T,sep = "\t",row.names = 1)

fitVidiano_S1 <- envfit(ord.nmds.brayS1, VolatileVidiano2020_S1[,])
# prepare the names for the to be plotted taxa and save them in a new column
mytax_tbl <- data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S1_100))
# the following command pastes the Phylum column with the Genus one
# !!! make sure to get rid of the NAs with the for loop as we discussed before !!!

# For ITS - Remove letter from taxonomy
mytax_tbl$forplt <- for (i in c(1:nrow(mytax_tbl))) {
  for(j in c(1:ncol(mytax_tbl))) {
    mytax_tbl[i,j] <- gsub("[a-z]__","",mytax_tbl[i,j])
  }
}

# prepare the plotting colours
myplotcols <- RColorBrewer::brewer.pal(n = 6, name = 'RdBu')
# prep the plot points
mynmdssit <- scores(ord.nmds.brayS1, display = "site")
# prep the species points
mynmdsspe <- scores(ord.nmds.brayS1, display = "species")
# prep the variable of interest (e.g. vinification) has to be a factor
my_sel_var <- factor(MUST_Bacteria_Vidiano_2020_Merged_S1_100@sam_data$vinification)
# convert it also into numbers in orde to use it for colours etc
myterr_sel <- as.numeric(my_sel_var)
# get the arrow data that came out of the fitting of the environmental variables
arrowdata <- data.frame(fitVidiano_S1$vectors$arrows)

# start the graphics device
cairo_pdf("Bacteria_Vidiano_2020_Volatile_S1.pdf", height = 7, width = 7)

# draw and empty plot that you will start to populate with points, arrows, ellipses, names etc.
plot(mynmdssit, frame = F, cex = 0, pch = 21, xlim = c(min(1.5*mynmdssit[,1]),max(1.5*mynmdssit[,1])))

# first add the ellipses in order tto keep them in the background and prevent them to obstruct information
vegan::ordiellipse(mynmdssit, groups = my_sel_var, kind = "ehull", lty = 2, lwd=1)

# add the sample points and give them colors one for each treatment
points(mynmdssit, bg = myplotcols[myterr_sel], pch = 21, cex = 1.5
       ##gia keimeno sta simeia##
       , text(mynmdssit, labels = row.names(mynmdssit), font = 6, col = rgb(55,55,55, max = 255, alpha = 100)))

# write the stress value in the bottom left of the plot (you need to change the global "adj" parameter from 0.5 to 0 in order to add the text in the left instead of the middle (if you gave the value of 1 it would add the text in the rightmost side... I provide some text as example))
par(adj = 1)
title(sub = "Bacteria_Vidiano_2020_Volatile_S1")
par(adj = .5)

# add the env parameter arrows (in the rgb colouring mode we have the red, green and blue values from 0 to 255 and also the alpha value or transparency which also ranges between 0 and 255)
#small
arrows(0,0,0.3*arrowdata[,1] , 0.3*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)

#big
arrows(0,0,1.5*arrowdata[,1] , 1.5*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)


# also the labeling text
text(0.5*arrowdata[,1], 0.5*arrowdata[,2], labels = row.names(arrowdata), cex = 1.3, font = 2, col = rgb(20,20,20, max = 255, alpha = 255)) 

##simantiko species arrows and labels##
-----------------------------------------------------------------
  ##Epilego ta ASVs pou thelo ite aritmitika ite onomastika
  ##Aritmitika##
  arrows(0,0,mynmdsspe[1:12,1] , mynmdsspe[1:12,2], angle = 25, length = 0.15, col = rgb(40,40,40, max = 255, alpha = 100))

plotrix::thigmophobe.labels(1.2*mynmdsspe[1:12,1], 1.2*mynmdsspe[1:12,2], labels = mytax_tbl[row.names(mynmdsspe)[1:12],"forplt"], cex = .6, font = 2, col = rgb(120,120,120, max = 255, alpha =200))


##Onomastika anti gia diaforika afthona##
arrows(0,0,0.5*mynmdsspe[c("ASV00062","ASV00051","ASV00202","ASV00136","ASV00210","ASV00355","ASV00225","ASV00241","ASV00061","ASV00237","ASV00527","ASV00286"),1] , 0.5*mynmdsspe[c("ASV00062","ASV00051","ASV00202","ASV00136","ASV00210","ASV00355","ASV00225","ASV00241","ASV00061","ASV00237","ASV00527","ASV00286"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#big
arrows(0,0,1.5*mynmdsspe[c("ASV00062","ASV00051","ASV00202","ASV00136","ASV00210","ASV00355","ASV00225","ASV00241","ASV00061","ASV00237","ASV00527","ASV00286"),1] , 1.5*mynmdsspe[c("ASV00062","ASV00051","ASV00202","ASV00136","ASV00210","ASV00355","ASV00225","ASV00241","ASV00061","ASV00237","ASV00527","ASV00286"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

Top_12_S1
("ASV00062","ASV00051","ASV00202","ASV00136","ASV00210","ASV00355","ASV00225","ASV00241","ASV00061","ASV00237","ASV00527","ASV00286")

write.table(data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_Genus)), file="tax_table_MUST_Bacteria_Vidiano_2020_Merged_Genus.txt", quote = F,col.names = NA, sep="\t")

write.table(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_100)), file="otu_table_MUST_Bacteria_Vidiano_2020_Merged_100.txt", quote = F,col.names = NA, sep="\t")

##LABELS Short-cut##
a <- c("Burkholderia_sp.","Gluconobacter_sp.","Comamonadaceae_sp.","Bacillales_sp.","Asinibacterium_sp.","Obscuribacteraceae_sp.","Xanthobacteraceae_sp.","Chlamydiales_sp.","Komagataeibacter_sp.","Acetobacter_sp.","Reyranella_sp.","Sphingomonadaceae_sp.")

install.packages("plotrix")

plotrix::thigmophobe.labels(1*mynmdsspe[c("ASV00062","ASV00051","ASV00202","ASV00136","ASV00210","ASV00355","ASV00225","ASV00241","ASV00061","ASV00237","ASV00527","ASV00286"),1], 1*mynmdsspe[c("ASV00062","ASV00051","ASV00202","ASV00136","ASV00210","ASV00355","ASV00225","ASV00241","ASV00061","ASV00237","ASV00527","ASV00286"),2], labels = a, cex = .6, font = 2, col = rgb(120,120,120,max = 255, alpha =200)) 

#simantiko# plotrix::thigmophobe.labels(1.2*mynmdsspe[1:15,1], 1.2*mynmdsspe[1:15,2], labels = mytax_tbl[row.names(mynmdsspe)[1:15],"forplt"], cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
# the following prints the legend
graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()

##_.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.--..--.-..-.-.##

##S2 Triplot
MUST_Bacteria_Vidiano_2020_Merged_S2
MUST_Bacteria_Vidiano_2020_Merged_S2_100 <- transform_sample_counts(MUST_Bacteria_Vidiano_2020_Merged_S2, function(OTU) 100*OTU/sum(OTU))

ord.nmds.brayS2 <- ordinate(MUST_Bacteria_Vidiano_2020_Merged_S2_100, method="CCA", distance="bray")

plot_ordination(MUST_Bacteria_Vidiano_2020_Merged_S2_100, ord.nmds.brayS2, color="vinification", shape ="", label = "", title=paste("S2")) + geom_point(size = 4)

View(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_S2_100)))

View(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_S2_100)))
View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_S2_100)))
View(data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S2_100)))


Top_12_S2 <- names(sort(taxa_sums(MUST_Bacteria_Vidiano_2020_Merged_S2_100), decreasing = TRUE)[1:12])

##S2_Triplot_Start
VolatileVidiano2020_S2 <- read.table("VolatileVidiano2020_S2.txt", header=T,sep = "\t",row.names = 1)

fitVidiano_S2 <- envfit(ord.nmds.brayS2, VolatileVidiano2020_S2[,])

# prepare the names for the to be plotted taxa and save them in a new column
mytax_tbl <- data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S2_100))
# the following command pastes the Phylum column with the Genus one
# !!! make sure to get rid of the NAs with the for loop as we discussed before !!!

# For ITS - Remove letter from taxonomy
mytax_tbl$forplt <- for (i in c(1:nrow(mytax_tbl))) {
  for(j in c(1:ncol(mytax_tbl))) {
    mytax_tbl[i,j] <- gsub("[a-z]__","",mytax_tbl[i,j])
  }
}

# prepare the plotting colours
myplotcols <- RColorBrewer::brewer.pal(n = 6, name = 'RdBu')
# prep the plot points
mynmdssit <- scores(ord.nmds.brayS2, display = "site")
# prep the species points
mynmdsspe <- scores(ord.nmds.brayS2, display = "species")
# prep the variable of interest (e.g. vinification) has to be a factor
my_sel_var <- factor(MUST_Bacteria_Vidiano_2020_Merged_S2_100@sam_data$vinification)
# convert it also into numbers in orde to use it for colours etc
myterr_sel <- as.numeric(my_sel_var)
# get the arrow data that came out of the fitting of the environmental variables
arrowdata <- data.frame(fitVidiano_S2$vectors$arrows)

# start the graphics device
cairo_pdf("Bacteria_Vidiano_2020_Volatile_S2.pdf", height = 7, width = 7)

# draw and empty plot that you will start to populate with points, arrows, ellipses, names etc.
plot(mynmdssit, frame = F, cex = 0, pch = 21, xlim = c(min(1.5*mynmdssit[,1]),max(1.5*mynmdssit[,1])))

# first add the ellipses in order tto keep them in the background and prevent them to obstruct information
vegan::ordiellipse(mynmdssit, groups = my_sel_var, kind = "ehull", lty = 2, lwd=1)

# add the sample points and give them colors one for each treatment
points(mynmdssit, bg = myplotcols[myterr_sel], pch = 21, cex = 1.5
       ##gia keimeno sta simeia##
       , text(mynmdssit, labels = row.names(mynmdssit), font = 6, col = rgb(55,55,55, max = 255, alpha = 100)))

# write the stress value in the bottom left of the plot (you need to change the global "adj" parameter from 0.5 to 0 in order to add the text in the left instead of the middle (if you gave the value of 1 it would add the text in the rightmost side... I provide some text as example))
par(adj = 1)
title(sub = "Bacteria_Vidiano_2020_Volatile_S2")
par(adj = .5)

# add the env parameter arrows (in the rgb colouring mode we have the red, green and blue values from 0 to 255 and also the alpha value or transparency which also ranges between 0 and 255)
#small
arrows(0,0,0.3*arrowdata[,1] , 0.3*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)

#big
arrows(0,0,1.5*arrowdata[,1] , 1.5*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)


# also the labeling text
text(0.5*arrowdata[,1], 0.5*arrowdata[,2], labels = row.names(arrowdata), cex = 1.3, font = 2, col = rgb(20,20,20, max = 255, alpha = 255)) 

##simantiko species arrows and labels##
-----------------------------------------------------------------
  ##Epilego ta ASVs pou thelo ite aritmitika ite onomastika
  ##Aritmitika##
  arrows(0,0,mynmdsspe[1:12,1] , mynmdsspe[1:12,2], angle = 25, length = 0.15, col = rgb(40,40,40, max = 255, alpha = 100))

plotrix::thigmophobe.labels(1.2*mynmdsspe[1:12,1], 1.2*mynmdsspe[1:12,2], labels = mytax_tbl[row.names(mynmdsspe)[1:12],"forplt"], cex = .6, font = 2, col = rgb(120,120,120, max = 255, alpha =200))

##Onomastika anti gia diaforika afthona##
arrows(0,0,0.5*mynmdsspe[c("ASV00062","ASV00051","ASV00202","ASV00136","ASV00061","ASV00210","ASV00237","ASV00355","ASV00286","ASV00241","ASV00225","ASV00527"),1] , 0.5*mynmdsspe[c("ASV00062","ASV00051","ASV00202","ASV00136","ASV00061","ASV00210","ASV00237","ASV00355","ASV00286","ASV00241","ASV00225","ASV00527"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#big
arrows(0,0,2.5*mynmdsspe[c("ASV00062","ASV00051","ASV00202","ASV00136","ASV00061","ASV00210","ASV00237","ASV00355","ASV00286","ASV00241","ASV00225","ASV00527"),1] , 2.5*mynmdsspe[c("ASV00062","ASV00051","ASV00202","ASV00136","ASV00061","ASV00210","ASV00237","ASV00355","ASV00286","ASV00241","ASV00225","ASV00527"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))


Top_12_S2
("ASV00062","ASV00051","ASV00202","ASV00136","ASV00061","ASV00210","ASV00237","ASV00355","ASV00286","ASV00241","ASV00225","ASV00527")

write.table(data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged)), file="tax_table_MUST_Bacteria_Vidiano_2020_Merged.txt", quote = F,col.names = NA, sep="\t")

write.table(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_100)), file="otu_table_MUST_Bacteria_Vidiano_2020_Merged_100.txt", quote = F,col.names = NA, sep="\t")

##LABELS Short-cut##
a <- c("Burkholderia_sp.","Gluconobacter_sp.","Comamonadaceae_sp.","Bacillales_sp.","Komagataeibacter_sp.","Asinibacterium_sp.","Acetobacter_sp.","Obscuribacteraceae_sp.","Sphingomonadaceae_sp.","Chlamydiales_sp.","Xanthobacteraceae_sp.","Reyranella_sp.")

install.packages("plotrix")

plotrix::thigmophobe.labels(1*mynmdsspe[c("ASV00062","ASV00051","ASV00202","ASV00136","ASV00061","ASV00210","ASV00237","ASV00355","ASV00286","ASV00241","ASV00225","ASV00527"),1], 1*mynmdsspe[c("ASV00062","ASV00051","ASV00202","ASV00136","ASV00061","ASV00210","ASV00237","ASV00355","ASV00286","ASV00241","ASV00225","ASV00527"),2], labels = a, cex = .6, font = 2, col = rgb(120,120,120,max = 255, alpha =200)) 

#simantiko# plotrix::thigmophobe.labels(1.2*mynmdsspe[1:15,1], 1.2*mynmdsspe[1:15,2], labels = mytax_tbl[row.names(mynmdsspe)[1:15],"forplt"], cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
# the following prints the legend
graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()

##_.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.--..--.-..-.-.##

##S3 Triplot
MUST_Bacteria_Vidiano_2020_Merged_S3
MUST_Bacteria_Vidiano_2020_Merged_S3_100 <- transform_sample_counts(MUST_Bacteria_Vidiano_2020_Merged_S3, function(OTU) 100*OTU/sum(OTU))

ord.nmds.brayS3 <- ordinate(MUST_Bacteria_Vidiano_2020_Merged_S3_100, method="CCA", distance="bray")

plot_ordination(MUST_Bacteria_Vidiano_2020_Merged_S3_100, ord.nmds.brayS3, color="vinification", shape ="", label = "", title=paste("S3")) + geom_point(size = 4)

View(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_S3_100)))

View(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_S3_100)))
View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_S3_100)))
View(data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S3_100)))


Top_12_S3 <- names(sort(taxa_sums(MUST_Bacteria_Vidiano_2020_Merged_S3_100), decreasing = TRUE)[1:12])

##S3_Triplot_Start
VolatileVidiano2020_S3 <- read.table("VolatileVidiano2020_S3.txt", header=T,sep = "\t",row.names = 1)

fitVidiano_S3 <- envfit(ord.nmds.brayS3, VolatileVidiano2020_S3[,])

# prepare the names for the to be plotted taxa and save them in a new column
mytax_tbl <- data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S3_100))
# the following command pastes the Phylum column with the Genus one
# !!! make sure to get rid of the NAs with the for loop as we discussed before !!!

# For ITS - Remove letter from taxonomy
mytax_tbl$forplt <- for (i in c(1:nrow(mytax_tbl))) {
  for(j in c(1:ncol(mytax_tbl))) {
    mytax_tbl[i,j] <- gsub("[a-z]__","",mytax_tbl[i,j])
  }
}


# prepare the plotting colours
myplotcols <- RColorBrewer::brewer.pal(n = 6, name = 'RdBu')
# prep the plot points
mynmdssit <- scores(ord.nmds.brayS3, display = "site")
# prep the species points
mynmdsspe <- scores(ord.nmds.brayS3, display = "species")
# prep the variable of interest (e.g. vinification) has to be a factor
my_sel_var <- factor(MUST_Bacteria_Vidiano_2020_Merged_S3_100@sam_data$vinification)
# convert it also into numbers in orde to use it for colours etc
myterr_sel <- as.numeric(my_sel_var)
# get the arrow data that came out of the fitting of the environmental variables
arrowdata <- data.frame(fitVidiano_S3$vectors$arrows)

# start the graphics device
cairo_pdf("Bacteria_Vidiano_2020_Volatile_S3.pdf", height = 7, width = 7)

# draw and empty plot that you will start to populate with points, arrows, ellipses, names etc.
plot(mynmdssit, frame = F, cex = 0, pch = 21, xlim = c(min(1.5*mynmdssit[,1]),max(1.5*mynmdssit[,1])))

# first add the ellipses in order tto keep them in the background and prevent them to obstruct information
vegan::ordiellipse(mynmdssit, groups = my_sel_var, kind = "ehull", lty = 2, lwd=1)

# add the sample points and give them colors one for each treatment
points(mynmdssit, bg = myplotcols[myterr_sel], pch = 21, cex = 1.5
       ##gia keimeno sta simeia##
       , text(mynmdssit, labels = row.names(mynmdssit), font = 6, col = rgb(55,55,55, max = 255, alpha = 100)))

# write the stress value in the bottom left of the plot (you need to change the global "adj" parameter from 0.5 to 0 in order to add the text in the left instead of the middle (if you gave the value of 1 it would add the text in the rightmost side... I provide some text as example))
par(adj = 1)
title(sub = "Bacteria_Vidiano_2020_Volatile_S3")
par(adj = .5)

# add the env parameter arrows (in the rgb colouring mode we have the red, green and blue values from 0 to 255 and also the alpha value or transparency which also ranges between 0 and 255)
#small
arrows(0,0,0.3*arrowdata[,1] , 0.3*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)

#big
arrows(0,0,1.5*arrowdata[,1] , 1.5*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)


# also the labeling text
text(0.5*arrowdata[,1], 0.5*arrowdata[,2], labels = row.names(arrowdata), cex = 1.3, font = 2, col = rgb(20,20,20, max = 255, alpha = 255)) 

##simantiko species arrows and labels##
-----------------------------------------------------------------
  ##Epilego ta ASVs pou thelo ite aritmitika ite onomastika
  ##Aritmitika##
  arrows(0,0,mynmdsspe[1:12,1] , mynmdsspe[1:12,2], angle = 25, length = 0.15, col = rgb(40,40,40, max = 255, alpha = 100))

plotrix::thigmophobe.labels(1.2*mynmdsspe[1:12,1], 1.2*mynmdsspe[1:12,2], labels = mytax_tbl[row.names(mynmdsspe)[1:12],"forplt"], cex = .6, font = 2, col = rgb(120,120,120, max = 255, alpha =200))

##Onomastika anti gia diaforika afthona##
arrows(0,0,0.5*mynmdsspe[c("ASV00051","ASV00062","ASV00136","ASV00061","ASV00202","ASV00286","ASV00237","ASV00355","ASV00216","ASV00522","ASV00347","ASV00527"),1] , 0.5*mynmdsspe[c("ASV00051","ASV00062","ASV00136","ASV00061","ASV00202","ASV00286","ASV00237","ASV00355","ASV00216","ASV00522","ASV00347","ASV00527"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#big
arrows(0,0,2.5*mynmdsspe[c("ASV00051","ASV00062","ASV00136","ASV00061","ASV00202","ASV00286","ASV00237","ASV00355","ASV00216","ASV00522","ASV00347","ASV00527"),1] , 2.5*mynmdsspe[c("ASV00051","ASV00062","ASV00136","ASV00061","ASV00202","ASV00286","ASV00237","ASV00355","ASV00216","ASV00522","ASV00347","ASV00527"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

Top_12_S3
("ASV00051","ASV00062","ASV00136","ASV00061","ASV00202","ASV00286","ASV00237","ASV00355","ASV00216","ASV00522","ASV00347","ASV00527")

write.table(data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged)), file="tax_table_MUST_Bacteria_Vidiano_2020_Merged.txt", quote = F,col.names = NA, sep="\t")

write.table(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_100)), file="otu_table_MUST_Bacteria_Vidiano_2020_Merged_100.txt", quote = F,col.names = NA, sep="\t")

##LABELS Short-cut##
a <- c("Gluconobacter_sp.","Burkholderia_sp.","Bacillales_sp.","Komagataeibacter_sp.","Comamonadaceae_sp.","Sphingomonadaceae_sp.","Acetobacter_sp.","Obscuribacteraceae_sp.","Sphingomonas_sp.","Novosphingobium_sp.","Magnetospirillaceae_sp.","Reyranella_sp.")


install.packages("plotrix")

plotrix::thigmophobe.labels(1*mynmdsspe[c("ASV00051","ASV00062","ASV00136","ASV00061","ASV00202","ASV00286","ASV00237","ASV00355","ASV00216","ASV00522","ASV00347","ASV00527"),1], 1*mynmdsspe[c("ASV00051","ASV00062","ASV00136","ASV00061","ASV00202","ASV00286","ASV00237","ASV00355","ASV00216","ASV00522","ASV00347","ASV00527"),2], labels = a, cex = .6, font = 2, col = rgb(120,120,120,max = 255, alpha =200)) 

#simantiko# plotrix::thigmophobe.labels(1.2*mynmdsspe[1:15,1], 1.2*mynmdsspe[1:15,2], labels = mytax_tbl[row.names(mynmdsspe)[1:15],"forplt"], cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
# the following prints the legend
graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()

##_.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.--..--.-..-.-.##

##S4 Triplot
MUST_Bacteria_Vidiano_2020_Merged_S4
MUST_Bacteria_Vidiano_2020_Merged_S4_100 <- transform_sample_counts(MUST_Bacteria_Vidiano_2020_Merged_S4, function(OTU) 100*OTU/sum(OTU))

ord.nmds.brayS4 <- ordinate(MUST_Bacteria_Vidiano_2020_Merged_S4_100, method="CCA", distance="bray")

plot_ordination(MUST_Bacteria_Vidiano_2020_Merged_S4_100, ord.nmds.brayS4, color="vinification", shape ="", label = "", title=paste("S4")) + geom_point(size = 4)

View(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_S4_100)))


View(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_S4_100)))
View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_S4_100)))
View(data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S4_100)))


Top_12_S4 <- names(sort(taxa_sums(MUST_Bacteria_Vidiano_2020_Merged_S4_100), decreasing = TRUE)[1:12])

##S4_Triplot_Start
VolatileVidiano2020_S4 <- read.table("VolatileVidiano2020_S4.txt", header=T,sep = "\t",row.names = 1)

fitVidiano_S4 <- envfit(ord.nmds.brayS4, VolatileVidiano2020_S4[,])

# prepare the names for the to be plotted taxa and save them in a new column
mytax_tbl <- data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S4_100))
# the following command pastes the Phylum column with the Genus one
# !!! make sure to get rid of the NAs with the for loop as we discussed before !!!

# For ITS - Remove letter from taxonomy
mytax_tbl$forplt <- for (i in c(1:nrow(mytax_tbl))) {
  for(j in c(1:ncol(mytax_tbl))) {
    mytax_tbl[i,j] <- gsub("[a-z]__","",mytax_tbl[i,j])
  }
}


# prepare the plotting colours
myplotcols <- RColorBrewer::brewer.pal(n = 6, name = 'RdBu')
# prep the plot points
mynmdssit <- scores(ord.nmds.brayS4, display = "site")
# prep the species points
mynmdsspe <- scores(ord.nmds.brayS4, display = "species")
# prep the variable of interest (e.g. vinification) has to be a factor
my_sel_var <- factor(MUST_Bacteria_Vidiano_2020_Merged_S4_100@sam_data$vinification)
# convert it also into numbers in orde to use it for colours etc
myterr_sel <- as.numeric(my_sel_var)
# get the arrow data that came out of the fitting of the environmental variables
arrowdata <- data.frame(fitVidiano_S4$vectors$arrows)

# start the graphics device
cairo_pdf("Bacteria_Vidiano_2020_Volatile_S4.pdf", height = 7, width = 7)

# draw and empty plot that you will start to populate with points, arrows, ellipses, names etc.
plot(mynmdssit, frame = F, cex = 0, pch = 21, xlim = c(min(1.5*mynmdssit[,1]),max(1.5*mynmdssit[,1])))

# first add the ellipses in order tto keep them in the background and prevent them to obstruct information
vegan::ordiellipse(mynmdssit, groups = my_sel_var, kind = "ehull", lty = 2, lwd=1)

# add the sample points and give them colors one for each treatment
points(mynmdssit, bg = myplotcols[myterr_sel], pch = 21, cex = 1.5
       ##gia keimeno sta simeia##
       , text(mynmdssit, labels = row.names(mynmdssit), font = 6, col = rgb(55,55,55, max = 255, alpha = 100)))

# write the stress value in the bottom left of the plot (you need to change the global "adj" parameter from 0.5 to 0 in order to add the text in the left instead of the middle (if you gave the value of 1 it would add the text in the rightmost side... I provide some text as example))
par(adj = 1)
title(sub = "Bacteria_Vidiano_2020_Volatile_S4")
par(adj = .5)

# add the env parameter arrows (in the rgb colouring mode we have the red, green and blue values from 0 to 255 and also the alpha value or transparency which also ranges between 0 and 255)

#small
arrows(0,0,0.3*arrowdata[,1] , 0.3*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)

#big
arrows(0,0,1.5*arrowdata[,1] , 1.5*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)


# also the labeling text
text(0.5*arrowdata[,1], 0.5*arrowdata[,2], labels = row.names(arrowdata), cex = 1.3, font = 2, col = rgb(20,20,20, max = 255, alpha = 255)) 

##simantiko species arrows and labels##
-----------------------------------------------------------------
  ##Epilego ta ASVs pou thelo ite aritmitika ite onomastika
  ##Aritmitika##
  arrows(0,0,mynmdsspe[1:12,1] , mynmdsspe[1:12,2], angle = 25, length = 0.15, col = rgb(40,40,40, max = 255, alpha = 100))

plotrix::thigmophobe.labels(1.2*mynmdsspe[1:12,1], 1.2*mynmdsspe[1:12,2], labels = mytax_tbl[row.names(mynmdsspe)[1:12],"forplt"], cex = .6, font = 2, col = rgb(120,120,120, max = 255, alpha =200))

##Onomastika anti gia diaforika afthona##

#small
arrows(0,0,0.5*mynmdsspe[c("ASV00062","ASV00051","ASV00136","ASV00061","ASV00286","ASV00225","ASV00237","ASV00202","ASV00355","ASV00347","ASV00216","ASV00019"),1] , 0.5*mynmdsspe[c("ASV00062","ASV00051","ASV00136","ASV00061","ASV00286","ASV00225","ASV00237","ASV00202","ASV00355","ASV00347","ASV00216","ASV00019"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#big
arrows(0,0,2.5*mynmdsspe[c("ASV00062","ASV00051","ASV00136","ASV00061","ASV00286","ASV00225","ASV00237","ASV00202","ASV00355","ASV00347","ASV00216","ASV00019"),1] , 2.5*mynmdsspe[c("ASV00062","ASV00051","ASV00136","ASV00061","ASV00286","ASV00225","ASV00237","ASV00202","ASV00355","ASV00347","ASV00216","ASV00019"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

Top_12_S4
("ASV00062","ASV00051","ASV00136","ASV00061","ASV00286","ASV00225","ASV00237","ASV00202","ASV00355","ASV00347","ASV00216","ASV00019")

write.table(data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged)), file="tax_table_MUST_Bacteria_Vidiano_2020_Merged.txt", quote = F,col.names = NA, sep="\t")

write.table(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_100)), file="otu_table_MUST_Bacteria_Vidiano_2020_Merged_100.txt", quote = F,col.names = NA, sep="\t")

##LABELS Short-cut##
a <- c("Burkholderia_sp.","Gluconobacter_sp.","Bacillales_sp.","Komagataeibacter_sp.","Sphingomonadaceae_sp.","Xanthobacteraceae_sp.","Acetobacter_sp.","Comamonadaceae_sp.","Obscuribacteraceae_sp","Magnetospirillaceae","Sphingomonas_sp.","Oenococcus_sp.")


install.packages("plotrix")

plotrix::thigmophobe.labels(1*mynmdsspe[c("ASV00062","ASV00051","ASV00136","ASV00061","ASV00286","ASV00225","ASV00237","ASV00202","ASV00355","ASV00347","ASV00216","ASV00019"),1], 1*mynmdsspe[c("ASV00062","ASV00051","ASV00136","ASV00061","ASV00286","ASV00225","ASV00237","ASV00202","ASV00355","ASV00347","ASV00216","ASV00019"),2], labels = a, cex = .6, font = 2, col = rgb(120,120,120,max = 255, alpha =200)) 

#simantiko# plotrix::thigmophobe.labels(1.2*mynmdsspe[1:15,1], 1.2*mynmdsspe[1:15,2], labels = mytax_tbl[row.names(mynmdsspe)[1:15],"forplt"], cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
# the following prints the legend
graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()

##_.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.--..--.-..-.-.##

##S5 Triplot
MUST_Bacteria_Vidiano_2020_Merged_S5
MUST_Bacteria_Vidiano_2020_Merged_S5_100 <- transform_sample_counts(MUST_Bacteria_Vidiano_2020_Merged_S5, function(OTU) 100*OTU/sum(OTU))

ord.nmds.brayS5 <- ordinate(MUST_Bacteria_Vidiano_2020_Merged_S5_100, method="CCA", distance="bray")

plot_ordination(MUST_Bacteria_Vidiano_2020_Merged_S5_100, ord.nmds.brayS5, color="vinification", shape ="", label = "", title=paste("S5")) + geom_point(size = 4)

View(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_S5_100)))

View(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_S5_100)))
View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_Merged_S5_100)))
View(data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S5_100)))


Top_12_S5 <- names(sort(taxa_sums(MUST_Bacteria_Vidiano_2020_Merged_S5_100), decreasing = TRUE)[1:12])

##S5_Triplot_Start
VolatileVidiano2020_S5 <- read.table("VolatileVidiano2020_S5.txt", header=T,sep = "\t",row.names = 1)

fitVidiano_S5 <- envfit(ord.nmds.brayS5, VolatileVidiano2020_S5[,])

# prepare the names for the to be plotted taxa and save them in a new column
mytax_tbl <- data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S5_100))
# the following command pastes the Phylum column with the Genus one
# !!! make sure to get rid of the NAs with the for loop as we discussed before !!!

# For ITS - Remove letter from taxonomy
mytax_tbl$forplt <- for (i in c(1:nrow(mytax_tbl))) {
  for(j in c(1:ncol(mytax_tbl))) {
    mytax_tbl[i,j] <- gsub("[a-z]__","",mytax_tbl[i,j])
  }
}


# prepare the plotting colours
myplotcols <- RColorBrewer::brewer.pal(n = 6, name = 'RdBu')
# prep the plot points
mynmdssit <- scores(ord.nmds.brayS5, display = "site")
# prep the species points
mynmdsspe <- scores(ord.nmds.brayS5, display = "species")
# prep the variable of interest (e.g. vinification) has to be a factor
my_sel_var <- factor(MUST_Bacteria_Vidiano_2020_Merged_S5_100@sam_data$vinification)
# convert it also into numbers in orde to use it for colours etc
myterr_sel <- as.numeric(my_sel_var)
# get the arrow data that came out of the fitting of the environmental variables
arrowdata <- data.frame(fitVidiano_S5$vectors$arrows)


# start the graphics device
cairo_pdf("Bacteria_Vidiano_2020_Volatile_S5.pdf", height = 7, width = 7)

# draw and empty plot that you will start to populate with points, arrows, ellipses, names etc.
plot(mynmdssit, frame = F, cex = 0, pch = 21, xlim = c(min(1.5*mynmdssit[,1]),max(1.5*mynmdssit[,1])))

# first add the ellipses in order tto keep them in the background and prevent them to obstruct information
vegan::ordiellipse(mynmdssit, groups = my_sel_var, kind = "ehull", lty = 2, lwd=1)

# add the sample points and give them colors one for each treatment
points(mynmdssit, bg = myplotcols[myterr_sel], pch = 21, cex = 1.5
       ##gia keimeno sta simeia##
       , text(mynmdssit, labels = row.names(mynmdssit), font = 6, col = rgb(55,55,55, max = 255, alpha = 100)))

# write the stress value in the bottom left of the plot (you need to change the global "adj" parameter from 0.5 to 0 in order to add the text in the left instead of the middle (if you gave the value of 1 it would add the text in the rightmost side... I provide some text as example))

par(adj = 1)
title(sub = "Bacteria_Vidiano_2020_Volatile_S5")
par(adj = .5)

# add the env parameter arrows (in the rgb colouring mode we have the red, green and blue values from 0 to 255 and also the alpha value or transparency which also ranges between 0 and 255)

#small
arrows(0,0,0.3*arrowdata[,1] , 0.3*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)

#big
arrows(0,0,1.5*arrowdata[,1] , 1.5*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)


# also the labeling text
text(0.5*arrowdata[,1], 0.5*arrowdata[,2], labels = row.names(arrowdata), cex = 1.3, font = 2, col = rgb(20,20,20, max = 255, alpha = 255)) 

##simantiko species arrows and labels##
-----------------------------------------------------------------
  ##Epilego ta ASVs pou thelo ite aritmitika ite onomastika
  ##Aritmitika##
  arrows(0,0,mynmdsspe[1:12,1] , mynmdsspe[1:12,2], angle = 25, length = 0.15, col = rgb(40,40,40, max = 255, alpha = 100))

plotrix::thigmophobe.labels(1.2*mynmdsspe[1:12,1], 1.2*mynmdsspe[1:12,2], labels = mytax_tbl[row.names(mynmdsspe)[1:12],"forplt"], cex = .6, font = 2, col = rgb(120,120,120, max = 255, alpha =200))

##Onomastika anti gia diaforika afthona##

#small
arrows(0,0,0.5*mynmdsspe[c("ASV00062","ASV00051","ASV00216","ASV00286","ASV00347","ASV00340","ASV00202","ASV00225","ASV00491","ASV00061","ASV00136","ASV00237"),1] , 0.5*mynmdsspe[c("ASV00062","ASV00051","ASV00216","ASV00286","ASV00347","ASV00340","ASV00202","ASV00225","ASV00491","ASV00061","ASV00136","ASV00237"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#big
arrows(0,0,2.5*mynmdsspe[c("ASV00062","ASV00051","ASV00216","ASV00286","ASV00347","ASV00340","ASV00202","ASV00225","ASV00491","ASV00061","ASV00136","ASV00237"),1] , 2.5*mynmdsspe[c("ASV00062","ASV00051","ASV00216","ASV00286","ASV00347","ASV00340","ASV00202","ASV00225","ASV00491","ASV00061","ASV00136","ASV00237"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

Top_12_S5
("ASV00062","ASV00051","ASV00216","ASV00286","ASV00347","ASV00340","ASV00202","ASV00225","ASV00491","ASV00061","ASV00136","ASV00237")

write.table(data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged)), file="tax_table_MUST_Bacteria_Vidiano_2020_Merged.txt", quote = F,col.names = NA, sep="\t")

write.table(data.frame(otu_table(MUST_Bacteria_Vidiano_2020_Merged_100)), file="otu_table_MUST_Bacteria_Vidiano_2020_Merged_100.txt", quote = F,col.names = NA, sep="\t")

##LABELS Short-cut##
a <- c("Burkholderia_sp.","Gluconobacter_sp.","Sphingomonas_sp.","Sphingomonadaceae_sp.","Magnetospirillaceae_sp.","Methylobacterium_sp.","Comamonadaceae_sp.","Xanthobacteraceae_sp.","Alphaproteobacteria_sp.","Komagataeibacter_sp.","Bacillales_sp.","Acetobacter")


install.packages("plotrix")

plotrix::thigmophobe.labels(1*mynmdsspe[c("ASV00062","ASV00051","ASV00216","ASV00286","ASV00347","ASV00340","ASV00202","ASV00225","ASV00491","ASV00061","ASV00136","ASV00237"),1], 1*mynmdsspe[c("ASV00062","ASV00051","ASV00216","ASV00286","ASV00347","ASV00340","ASV00202","ASV00225","ASV00491","ASV00061","ASV00136","ASV00237"),2], labels = a, cex = .6, font = 2, col = rgb(120,120,120,max = 255, alpha =200)) 

#simantiko# plotrix::thigmophobe.labels(1.2*mynmdsspe[1:15,1], 1.2*mynmdsspe[1:15,2], labels = mytax_tbl[row.names(mynmdsspe)[1:15],"forplt"], cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
# the following prints the legend
graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()

##############################NEW ANALYSIS###########################
####Again TRIPLOT ANALYSIS for D.A Bacteria per Stage##
###Epeidi antimetopizo provlima sto Kruskal-Wallis gia tin euresi ton diaforikon afthonon per Stage exaitias tis apousias duplicate stin aythormiti oinopoiisi. Kano texnito duplicate diplasiazontas to mono deigma apo tin aythormiti kai kanontas merged ta phyloseq objects, gia na doylewei omos to script prepei na to allaxo elafros petontas exo merika taxa##
MUST_Bacteria_Vidiano_2020

write.table(data.frame(sample_data(MUST_Bacteria_Vidiano_2020)), file="sample_data_MUST_Bacteria_Vidiano_2020.txt", quote = F,col.names = NA, sep="\t")

SampleDataNew67 <- read.table("sample_data_MUST_Bacteria_Vidiano_2020.txt", header=T,sep = "\t",row.names = 1)

sample_data(MUST_Bacteria_Vidiano_2020) <- SampleDataNew67

View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020)))

MUST_Bacteria_Vidiano_2020_Genus <- tax_glom(MUST_Bacteria_Vidiano_2020, taxrank = "Genus")

MUST_Bacteria_Vidiano_2020_Genus <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_Genus)>0,MUST_Bacteria_Vidiano_2020_Genus)

##S1 Prepared for duplicate Spontaneous
MUST_Bacteria_Vidiano_2020_S1 <- subset_samples(MUST_Bacteria_Vidiano_2020_Genus, stage %in% (c("S1")))

MUST_Bacteria_Vidiano_2020_S1 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_S1)>0,MUST_Bacteria_Vidiano_2020_S1)

View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_S1)))

MUST_Bacteria_Vidiano_2020_S1 ##UNmerged S1 without spontaneous duplicate

##Duplicate Preperation
MUST_Bacteria_Vidiano_2020_S_S1 <- subset_samples(MUST_Bacteria_Vidiano_2020_Genus,replicate %in% (c("FS1")))

MUST_Bacteria_Vidiano_2020_S_S1 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_S_S1)>0,MUST_Bacteria_Vidiano_2020_S_S1)

MUST_Bacteria_Vidiano_2020_S_S1 ##Duplicate of spontaneous S1 (single sample) (97 taxa) (must be trimed so to be different from the original)

bS1 <- MUST_Bacteria_Vidiano_2020_S_S1 #97-1=96

myTaxa_bS1 <- names(sort(taxa_sums(bS1), decreasing = TRUE)[1:96]) 

Top_bS1 <- prune_taxa(myTaxa_bS1, bS1)

Top_bS1 <- prune_taxa(taxa_sums(Top_bS1)>0,Top_bS1)

sample_names(Top_bS1) <- paste("samp",sample_names(bS1), sep = "")

View(data.frame(sample_data(Top_bS1)))

##Merging Main Object with Spontaneous Duplicate#
c_S1 <- merge_phyloseq(MUST_Bacteria_Vidiano_2020_S1, Top_bS1)

View(data.frame(sample_data(c_S1)))

c_S1
##edo xekinaei to D.A part daneismeno apo to script ton heatmaps
c_S1 ##136 taxa

psdt <- c_S1

##transform phyloseq object raw counts to relative abundance##
psdt_ra <- transform_sample_counts(psdt, function(x) x / sum(x))

##Define the number of OTUs to be used in the analysis, 200 ASVs in our case## 
mynumOTUs <- 326 #326 taxa

## Prepare working object and select experimental variable of interest##
ps_htmp <- prune_taxa(names(taxa_sums(psdt_ra)[order(taxa_sums(psdt_ra), decreasing = T)][1:mynumOTUs]),psdt_ra)

## Prepare experimental variables table##
mydesign <- data.frame(sample_data(ps_htmp))
## Prepare the relative abundance table for the top 200 ASVs
dtcounts <- decostand(data.frame(otu_table(ps_htmp))[row.names(mydesign),], method = "range", MARGIN = 2)
## ASVs 
myotus <- colnames(dtcounts)

## Run Kruskal-Wallis, non-parametrix multiple comparison test##
mystatsout <- list()

for(myotu in myotus){
  if(sum(dtcounts[,myotu]) == 0){
    mystatsout[[myotu]]$krusk$statistics$p.chisq <- 1
  } else {
    mykrusk <- kruskal(dtcounts[,myotu], mydesign[rownames(dtcounts),]$vinification, group = T)
    mystatsout[[myotu]][["krusk"]] <- mykrusk
  }
}

## Extract p-values from the test results without p.adj##
mystatsoutkruskpvals <- unlist(lapply(myotus, function(x) mystatsout[[x]]$krusk$statistics$p.chisq))
names(mystatsoutkruskpvals) <- myotus

## Adjust p-values with "False Discovery Rate" method
mystatsoutkruskpvals.adj <- p.adjust(mystatsoutkruskpvals, method = "fdr")

# Check if there are any ASVs that significantly differ with p.adjust
for (i in c(1:length(mystatsoutkruskpvals.adj))) {
  if (mystatsoutkruskpvals.adj[i] < 0.05) {
    print(names(mystatsoutkruskpvals.adj[i]))
  }
}
###without p.adjust
for (i in c(1:length(mystatsoutkruskpvals))) {
  if (mystatsoutkruskpvals[i] < 0.06) {
    print(names(mystatsoutkruskpvals[i]))
  }
}
###without p.adjust finally gave results
[1] "ASV00470"
[1] "ASV00494"
[1] "ASV01141"
[1] "ASV01633"
[1] "ASV02088"
[1] "ASV02558"
[1] "ASV02791"
[1] "ASV03044"
[1] "ASV03428"
[1] "ASV04749"
[1] "ASV06136"
[1] "ASV06194"
[1] "ASV07177"
[1] "ASV08439"
[1] "ASV09167"
[1] "ASV09851"
[1] "ASV11684"
[1] "ASV13504"
[1] "ASV14049"
[1] "ASV15075"
[1] "ASV15162"
[1] "ASV15826"
[1] "ASV18925"
[1] "ASV27982"
[1] "ASV36703"

##S1_Triplot_Start_with D.A taxa
VolatileVidiano2020_S1 <- read.table("VolatileVidiano2020_S1.txt", header=T,sep = "\t",row.names = 1)

fitVidiano_S1 <- envfit(ord.nmds.brayS1, VolatileVidiano2020_S1[,])
# prepare the names for the to be plotted taxa and save them in a new column
mytax_tbl <- data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S1_100))
# prepare the plotting colours
myplotcols <- RColorBrewer::brewer.pal(n = 6, name = 'RdBu')
# prep the plot points
mynmdssit <- scores(ord.nmds.brayS1, display = "site")
# prep the species points
mynmdsspe <- scores(ord.nmds.brayS1, display = "species")
# prep the variable of interest (e.g. vinification) has to be a factor
my_sel_var <- factor(MUST_Bacteria_Vidiano_2020_Merged_S1_100@sam_data$vinification)
# convert it also into numbers in orde to use it for colours etc
myterr_sel <- as.numeric(my_sel_var)
# get the arrow data that came out of the fitting of the environmental variables
arrowdata <- data.frame(fitVidiano_S1$vectors$arrows)

# start the graphics device
cairo_pdf("Bacteria_Vidiano_2020_Volatile_S1.DA.pdf", height = 7, width = 7)

# draw and empty plot that you will start to populate with points, arrows, ellipses, names etc.
plot(mynmdssit, frame = F, cex = 0, pch = 21, xlim = c(min(1.5*mynmdssit[,1]),max(1.5*mynmdssit[,1])))

# first add the ellipses in order tto keep them in the background and prevent them to obstruct information
vegan::ordiellipse(mynmdssit, groups = my_sel_var, kind = "ehull", lty = 2, lwd=1)

# add the sample points and give them colors one for each treatment
points(mynmdssit, bg = myplotcols[myterr_sel], pch = 21, cex = 1.5
       ##gia keimeno sta simeia##
       , text(mynmdssit, labels = row.names(mynmdssit), font = 6, col = rgb(55,55,55, max = 255, alpha = 100)))

par(adj = 1)
title(sub = "Bacteria_Vidiano_2020_Volatile_S1.DA")
par(adj = .5)

# add the env parameter arrows (in the rgb colouring mode we have the red, green and blue values from 0 to 255 and also the alpha value or transparency which also ranges between 0 and 255)
#small
arrows(0,0,0.3*arrowdata[,1] , 0.3*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)
#big
arrows(0,0,1.5*arrowdata[,1] , 1.5*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)

# also the labeling text
text(0.5*arrowdata[,1], 0.5*arrowdata[,2], labels = row.names(arrowdata), cex = 1.3, font = 2, col = rgb(20,20,20, max = 255, alpha = 255)) 

##simantiko species arrows and labels##
-----------------------------------------------------------------
  ##diaforika afthona##
  #small
  arrows(0,0,0.5*mynmdsspe[c("ASV00470","ASV00494","ASV01141","ASV01633","ASV02088","ASV02558","ASV02791","ASV03044","ASV03428","ASV04749","ASV06136","ASV06194","ASV07177","ASV08439","ASV09167","ASV09851","ASV11684","ASV13504","ASV14049","ASV15075","ASV15162","ASV15826","ASV18925","ASV27982","ASV36703"),1] , 0.5*mynmdsspe[c("ASV00470","ASV00494","ASV01141","ASV01633","ASV02088","ASV02558","ASV02791","ASV03044","ASV03428","ASV04749","ASV06136","ASV06194","ASV07177","ASV08439","ASV09167","ASV09851","ASV11684","ASV13504","ASV14049","ASV15075","ASV15162","ASV15826","ASV18925","ASV27982","ASV36703"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#big
arrows(0,0,3.5*mynmdsspe[c("ASV00470","ASV00494","ASV01141","ASV01633","ASV02088","ASV02558","ASV02791","ASV03044","ASV03428","ASV04749","ASV06136","ASV06194","ASV07177","ASV08439","ASV09167","ASV09851","ASV11684","ASV13504","ASV14049","ASV15075","ASV15162","ASV15826","ASV18925","ASV27982","ASV36703"),1] , 3.5*mynmdsspe[c("ASV00470","ASV00494","ASV01141","ASV01633","ASV02088","ASV02558","ASV02791","ASV03044","ASV03428","ASV04749","ASV06136","ASV06194","ASV07177","ASV08439","ASV09167","ASV09851","ASV11684","ASV13504","ASV14049","ASV15075","ASV15162","ASV15826","ASV18925","ASV27982","ASV36703"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#D.A from Kruskal
("ASV00470","ASV00494","ASV01141","ASV01633","ASV02088","ASV02558","ASV02791","ASV03044","ASV03428","ASV04749","ASV06136","ASV06194","ASV07177","ASV08439","ASV09167","ASV09851","ASV11684","ASV13504","ASV14049","ASV15075","ASV15162","ASV15826","ASV18925","ASV27982","ASV36703")


##LABELS Short-cut##
a <- c("Fructobacillus_sp.","Cutibacterium_sp.","Microbacterium_sp.","Psychroglaciecola_sp.","Yersiniaceae_sp.","Modestobacter_sp.","Moraxellaceae_sp.","Shewanella_sp.","Methyloversatilis_sp.","Mucilaginibacter_sp.","Nubsella_sp.","Halobacteroidaceae_sp.","Lysobacter_sp.","Hymenobacter_sp.","Cellvibrio_sp.","Xanthomonadaceae_sp.","Thermoanaerobacterales_sp.","Amaricoccus_sp.","Gaiellales_sp.","Halomonas_sp.","Flectobacillus_sp.","Varibaculum_sp.","Actinobacteriota_sp.","Syntrophomonas_sp.","Defluviicoccus_sp.")

install.packages("plotrix")

plotrix::thigmophobe.labels(1*mynmdsspe[c("ASV00470","ASV00494","ASV01141","ASV01633","ASV02088","ASV02558","ASV02791","ASV03044","ASV03428","ASV04749","ASV06136","ASV06194","ASV07177","ASV08439","ASV09167","ASV09851","ASV11684","ASV13504","ASV14049","ASV15075","ASV15162","ASV15826","ASV18925","ASV27982","ASV36703"),1], 1*mynmdsspe[c("ASV00470","ASV00494","ASV01141","ASV01633","ASV02088","ASV02558","ASV02791","ASV03044","ASV03428","ASV04749","ASV06136","ASV06194","ASV07177","ASV08439","ASV09167","ASV09851","ASV11684","ASV13504","ASV14049","ASV15075","ASV15162","ASV15826","ASV18925","ASV27982","ASV36703"),2], labels = a, cex = .6, font = 2, col = rgb(120,120,120,max = 255, alpha =200)) 

#simantiko# plotrix::thigmophobe.labels(1.2*mynmdsspe[1:15,1], 1.2*mynmdsspe[1:15,2], labels = mytax_tbl[row.names(mynmdsspe)[1:15],"forplt"], cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
# the following prints the legend
graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()

##_.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.--..--.-..-.-.##

##S2 Prepared for duplicate Spontaneous

MUST_Bacteria_Vidiano_2020_S2 <- subset_samples(MUST_Bacteria_Vidiano_2020_Genus, stage %in% (c("S2")))

MUST_Bacteria_Vidiano_2020_S2 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_S2)>0,MUST_Bacteria_Vidiano_2020_S2)

View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_S2)))

MUST_Bacteria_Vidiano_2020_S2 ##UNmerged S2 without spontaneous duplicate

##Duplicate Preperation
MUST_Bacteria_Vidiano_2020_S_S2 <- subset_samples(MUST_Bacteria_Vidiano_2020_Genus,replicate %in% (c("FS2")))

MUST_Bacteria_Vidiano_2020_S_S2 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_S_S2)>0,MUST_Bacteria_Vidiano_2020_S_S2)

MUST_Bacteria_Vidiano_2020_S_S2 ##Duplicate of spontaneous S1 (single sample) (252 taxa) (must be trimed so to be different from the original)

bS2 <- MUST_Bacteria_Vidiano_2020_S_S2 #112-1=111

myTaxa_bS2 <- names(sort(taxa_sums(bS2), decreasing = TRUE)[1:111]) 

Top_bS2 <- prune_taxa(myTaxa_bS2, bS2)

Top_bS2 <- prune_taxa(taxa_sums(Top_bS2)>0,Top_bS2)

sample_names(Top_bS2) <- paste("samp",sample_names(bS2), sep = "")

View(data.frame(sample_data(Top_bS2)))

##Merging Main Object with Spontaneous Duplicate#
c_S2 <- merge_phyloseq(MUST_Bacteria_Vidiano_2020_S2, Top_bS2)

View(data.frame(sample_data(c_S2)))

c_S2
##edo xekinaei to D.A part daneismeno apo to script ton heatmaps
c_S2 ##456 taxa

psdt <- c_S2

##transform phyloseq object raw counts to relative abundance##
psdt_ra <- transform_sample_counts(psdt, function(x) x / sum(x))

##Define the number of OTUs to be used in the analysis, 200 ASVs in our case## 
mynumOTUs <-274 #274 taxa 

## Prepare working object and select experimental variable of interest##
ps_htmp <- prune_taxa(names(taxa_sums(psdt_ra)[order(taxa_sums(psdt_ra), decreasing = T)][1:mynumOTUs]),psdt_ra)

## Prepare experimental variables table##
mydesign <- data.frame(sample_data(ps_htmp))
## Prepare the relative abundance table for the top 200 ASVs
dtcounts <- decostand(data.frame(otu_table(ps_htmp))[row.names(mydesign),], method = "range", MARGIN = 2)
## ASVs 
myotus <- colnames(dtcounts)

## Run Kruskal-Wallis, non-parametrix multiple comparison test##
mystatsout <- list()

for(myotu in myotus){
  if(sum(dtcounts[,myotu]) == 0){
    mystatsout[[myotu]]$krusk$statistics$p.chisq <- 1
  } else {
    mykrusk <- kruskal(dtcounts[,myotu], mydesign[rownames(dtcounts),]$vinification, group = T)
    mystatsout[[myotu]][["krusk"]] <- mykrusk
  }
}

## Extract p-values from the test results without p.adj##
mystatsoutkruskpvals <- unlist(lapply(myotus, function(x) mystatsout[[x]]$krusk$statistics$p.chisq))
names(mystatsoutkruskpvals) <- myotus

## Adjust p-values with "False Discovery Rate" method
mystatsoutkruskpvals.adj <- p.adjust(mystatsoutkruskpvals, method = "fdr")

# Check if there are any ASVs that significantly differ with p.adjust
for (i in c(1:length(mystatsoutkruskpvals.adj))) {
  if (mystatsoutkruskpvals.adj[i] < 0.05) {
    print(names(mystatsoutkruskpvals.adj[i]))
  }
}
###without p.adjust
for (i in c(1:length(mystatsoutkruskpvals))) {
  if (mystatsoutkruskpvals[i] < 0.06) {
    print(names(mystatsoutkruskpvals[i]))
  }
}
###without p.adjust finally gave results
[1] "ASV02414"
[1] "ASV02626"
[1] "ASV03623"
[1] "ASV06332"
[1] "ASV08244"
[1] "ASV08341"
[1] "ASV08488"
[1] "ASV10318"
[1] "ASV13449"
[1] "ASV14724"
[1] "ASV17439"
[1] "ASV25895"
[1] "ASV31849"
[1] "ASV34685"
[1] "ASV40731"
[1] "ASV42794"

##S2_Triplot_Start_with D.A taxa
VolatileVidiano2020_S2 <- read.table("VolatileVidiano2020_S2.txt", header=T,sep = "\t",row.names = 1)

fitVidiano_S2 <- envfit(ord.nmds.brayS2, VolatileVidiano2020_S2[,])
# prepare the names for the to be plotted taxa and save them in a new column
mytax_tbl <- data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S2_100))
# prepare the plotting colours
myplotcols <- RColorBrewer::brewer.pal(n = 6, name = 'RdBu')
# prep the plot points
mynmdssit <- scores(ord.nmds.brayS2, display = "site")
# prep the species points
mynmdsspe <- scores(ord.nmds.brayS2, display = "species")
# prep the variable of interest (e.g. vinification) has to be a factor
my_sel_var <- factor(MUST_Bacteria_Vidiano_2020_Merged_S2_100@sam_data$vinification)
# convert it also into numbers in orde to use it for colours etc
myterr_sel <- as.numeric(my_sel_var)
# get the arrow data that came out of the fitting of the environmental variables
arrowdata <- data.frame(fitVidiano_S2$vectors$arrows)

# start the graphics device
cairo_pdf("Bacteria_Vidiano_2020_Volatile_S2.DA...pdf", height = 7, width = 7)

# draw and empty plot that you will start to populate with points, arrows, ellipses, names etc.
plot(mynmdssit, frame = F, cex = 0, pch = 21, xlim = c(min(1.5*mynmdssit[,1]),max(1.5*mynmdssit[,1])))

# first add the ellipses in order tto keep them in the background and prevent them to obstruct information
vegan::ordiellipse(mynmdssit, groups = my_sel_var, kind = "ehull", lty = 2, lwd=1)

# add the sample points and give them colors one for each treatment
points(mynmdssit, bg = myplotcols[myterr_sel], pch = 21, cex = 1.5
       ##gia keimeno sta simeia##
       , text(mynmdssit, labels = row.names(mynmdssit), font = 6, col = rgb(55,55,55, max = 255, alpha = 100)))

par(adj = 1)
title(sub = "Bacteria_Vidiano_2020_Volatile_S2.DA")
par(adj = .5)

# add the env parameter arrows (in the rgb colouring mode we have the red, green and blue values from 0 to 255 and also the alpha value or transparency which also ranges between 0 and 255)
#small
arrows(0,0,0.3*arrowdata[,1] , 0.3*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)
#big
arrows(0,0,1.5*arrowdata[,1] , 1.5*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)

# also the labeling text
text(0.5*arrowdata[,1], 0.5*arrowdata[,2], labels = row.names(arrowdata), cex = 1.3, font = 2, col = rgb(20,20,20, max = 255, alpha = 255)) 

##simantiko species arrows and labels##
-----------------------------------------------------------------
  ##diaforika afthona##
  #small
  arrows(0,0,0.5*mynmdsspe[c("ASV02414","ASV02626","ASV03623","ASV06332","ASV08244","ASV08341","ASV08488","ASV10318","ASV13449","ASV14724","ASV17439","ASV25895","ASV31849","ASV34685","ASV40731","ASV42794"),1] , 0.5*mynmdsspe[c("ASV02414","ASV02626","ASV03623","ASV06332","ASV08244","ASV08341","ASV08488","ASV10318","ASV13449","ASV14724","ASV17439","ASV25895","ASV31849","ASV34685","ASV40731","ASV42794"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#big
arrows(0,0,3.5*mynmdsspe[c("ASV02414","ASV02626","ASV03623","ASV06332","ASV08244","ASV08341","ASV08488","ASV10318","ASV13449","ASV14724","ASV17439","ASV25895","ASV31849","ASV34685","ASV40731","ASV42794"),1] , 3.5*mynmdsspe[c("ASV02414","ASV02626","ASV03623","ASV06332","ASV08244","ASV08341","ASV08488","ASV10318","ASV13449","ASV14724","ASV17439","ASV25895","ASV31849","ASV34685","ASV40731","ASV42794"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#D.A from Kruskal
("ASV02414","ASV02626","ASV03623","ASV06332","ASV08244","ASV08341","ASV08488","ASV10318","ASV13449","ASV14724","ASV17439","ASV25895","ASV31849","ASV34685","ASV40731","ASV42794")

[1] "ASV02414"
[1] "ASV02626"
[1] "ASV03623"
[1] "ASV06332"
[1] "ASV08244"
[1] "ASV08341"
[1] "ASV08488"
[1] "ASV10318"
[1] "ASV13449"
[1] "ASV14724"
[1] "ASV17439"
[1] "ASV25895"
[1] "ASV31849"
[1] "ASV34685"
[1] "ASV40731"
[1] "ASV42794"

##LABELS Short-cut##
a <- c("Deinococcus_sp.","Porphyromonas_sp.","Alkalicoccus_sp.","Anaerolineae_sp.","Veillonella_sp.","Pirellula_sp.","Rubinisphaeraceae_sp.","Ruminococcaceae_sp.","Pseudonocardia_sp.","Candidatus_Paracaedibacter_sp.","Sandaracinus_sp.","Longimicrobiaceae_sp.","Centipeda_sp.","Ammoniphilus_sp","Actinobacteriota_sp.","Olivibacter_sp.")

install.packages("plotrix")

plotrix::thigmophobe.labels(1*mynmdsspe[c("ASV02414","ASV02626","ASV03623","ASV06332","ASV08244","ASV08341","ASV08488","ASV10318","ASV13449","ASV14724","ASV17439","ASV25895","ASV31849","ASV34685","ASV40731","ASV42794"),1], 1*mynmdsspe[c("ASV02414","ASV02626","ASV03623","ASV06332","ASV08244","ASV08341","ASV08488","ASV10318","ASV13449","ASV14724","ASV17439","ASV25895","ASV31849","ASV34685","ASV40731","ASV42794"),2], labels = a, cex = .6, font = 2, col = rgb(120,120,120,max = 255, alpha =200)) 

#simantiko# plotrix::thigmophobe.labels(1.2*mynmdsspe[1:15,1], 1.2*mynmdsspe[1:15,2], labels = mytax_tbl[row.names(mynmdsspe)[1:15],"forplt"], cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
# the following prints the legend
graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()

##_.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.--..--.-..-.-.##

##S3 Prepared for duplicate Spontaneous

MUST_Bacteria_Vidiano_2020_S3 <- subset_samples(MUST_Bacteria_Vidiano_2020_Genus, stage %in% (c("S3")))

MUST_Bacteria_Vidiano_2020_S3 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_S3)>0,MUST_Bacteria_Vidiano_2020_S3)

View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_S3)))

MUST_Bacteria_Vidiano_2020_S3 ##UNmerged S2 without spontaneous duplicate

##Duplicate Preperation
MUST_Bacteria_Vidiano_2020_S_S3 <- subset_samples(MUST_Bacteria_Vidiano_2020_Genus,replicate %in% (c("FS3")))

MUST_Bacteria_Vidiano_2020_S_S3 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_S_S3)>0,MUST_Bacteria_Vidiano_2020_S_S3)

MUST_Bacteria_Vidiano_2020_S_S3 ##Duplicate of spontaneous S3 (single sample) (53 taxa) (must be trimed so to be different from the original)

bS3 <- MUST_Bacteria_Vidiano_2020_S_S3 #34-1=33

myTaxa_bS3 <- names(sort(taxa_sums(bS3), decreasing = TRUE)[1:33]) 

Top_bS3 <- prune_taxa(myTaxa_bS3, bS3)

Top_bS3 <- prune_taxa(taxa_sums(Top_bS3)>0,Top_bS3)

sample_names(Top_bS3) <- paste("samp",sample_names(bS3), sep = "")

View(data.frame(sample_data(Top_bS3)))

##Merging Main Object with Spontaneous Duplicate#
c_S3 <- merge_phyloseq(MUST_Bacteria_Vidiano_2020_S3, Top_bS3)

View(data.frame(sample_data(c_S3)))

c_S3
##edo xekinaei to D.A part daneismeno apo to script ton heatmaps
c_S3 ##56 taxa

psdt <- c_S3

##transform phyloseq object raw counts to relative abundance##
psdt_ra <- transform_sample_counts(psdt, function(x) x / sum(x))

##Define the number of OTUs to be used in the analysis, 200 ASVs in our case## 
mynumOTUs <- 249 #249 taxa

## Prepare working object and select experimental variable of interest##
ps_htmp <- prune_taxa(names(taxa_sums(psdt_ra)[order(taxa_sums(psdt_ra), decreasing = T)][1:mynumOTUs]),psdt_ra)

## Prepare experimental variables table##
mydesign <- data.frame(sample_data(ps_htmp))
## Prepare the relative abundance table for the top 200 ASVs
dtcounts <- decostand(data.frame(otu_table(ps_htmp))[row.names(mydesign),], method = "range", MARGIN = 2)
## ASVs 
myotus <- colnames(dtcounts)

## Run Kruskal-Wallis, non-parametrix multiple comparison test##
mystatsout <- list()

for(myotu in myotus){
  if(sum(dtcounts[,myotu]) == 0){
    mystatsout[[myotu]]$krusk$statistics$p.chisq <- 1
  } else {
    mykrusk <- kruskal(dtcounts[,myotu], mydesign[rownames(dtcounts),]$vinification, group = T)
    mystatsout[[myotu]][["krusk"]] <- mykrusk
  }
}

## Extract p-values from the test results without p.adj##
mystatsoutkruskpvals <- unlist(lapply(myotus, function(x) mystatsout[[x]]$krusk$statistics$p.chisq))
names(mystatsoutkruskpvals) <- myotus

## Adjust p-values with "False Discovery Rate" method
mystatsoutkruskpvals.adj <- p.adjust(mystatsoutkruskpvals, method = "fdr")

# Check if there are any ASVs that significantly differ with p.adjust
for (i in c(1:length(mystatsoutkruskpvals.adj))) {
  if (mystatsoutkruskpvals.adj[i] < 0.05) {
    print(names(mystatsoutkruskpvals.adj[i]))
  }
}
###without p.adjust
for (i in c(1:length(mystatsoutkruskpvals))) {
  if (mystatsoutkruskpvals[i] < 0.06) {
    print(names(mystatsoutkruskpvals[i]))
  }
}
###without p.adjust finally gave results
[1] "ASV00493"
[1] "ASV03432"
[1] "ASV04749"
[1] "ASV05606"
[1] "ASV06706"
[1] "ASV16590"


##S3_Triplot_Start_with D.A taxa
VolatileVidiano2020_S3 <- read.table("VolatileVidiano2020_S3.txt", header=T,sep = "\t",row.names = 1)

fitVidiano_S3 <- envfit(ord.nmds.brayS3, VolatileVidiano2020_S3[,])
# prepare the names for the to be plotted taxa and save them in a new column
mytax_tbl <- data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S3_100))
# prepare the plotting colours
myplotcols <- RColorBrewer::brewer.pal(n = 6, name = 'RdBu')
# prep the plot points
mynmdssit <- scores(ord.nmds.brayS3, display = "site")
# prep the species points
mynmdsspe <- scores(ord.nmds.brayS3, display = "species")
# prep the variable of interest (e.g. vinification) has to be a factor
my_sel_var <- factor(MUST_Bacteria_Vidiano_2020_Merged_S3_100@sam_data$vinification)
# convert it also into numbers in orde to use it for colours etc
myterr_sel <- as.numeric(my_sel_var)
# get the arrow data that came out of the fitting of the environmental variables
arrowdata <- data.frame(fitVidiano_S3$vectors$arrows)

# start the graphics device
cairo_pdf("Bacteria_Vidiano_2020_Volatile_S3.DA...pdf", height = 7, width = 7)

# draw and empty plot that you will start to populate with points, arrows, ellipses, names etc.
plot(mynmdssit, frame = F, cex = 0, pch = 21, xlim = c(min(1.5*mynmdssit[,1]),max(1.5*mynmdssit[,1])))

# first add the ellipses in order tto keep them in the background and prevent them to obstruct information
vegan::ordiellipse(mynmdssit, groups = my_sel_var, kind = "ehull", lty = 2, lwd=1)

# add the sample points and give them colors one for each treatment
points(mynmdssit, bg = myplotcols[myterr_sel], pch = 21, cex = 1.5
       ##gia keimeno sta simeia##
       , text(mynmdssit, labels = row.names(mynmdssit), font = 6, col = rgb(55,55,55, max = 255, alpha = 100)))

par(adj = 1)
title(sub = "Bacteria_Vidiano_2020_Volatile_S3.DA")
par(adj = .5)

# add the env parameter arrows (in the rgb colouring mode we have the red, green and blue values from 0 to 255 and also the alpha value or transparency which also ranges between 0 and 255)
#small
arrows(0,0,0.3*arrowdata[,1] , 0.3*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)
#big
arrows(0,0,1.5*arrowdata[,1] , 1.5*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)

# also the labeling text
text(0.5*arrowdata[,1], 0.5*arrowdata[,2], labels = row.names(arrowdata), cex = 1.3, font = 2, col = rgb(20,20,20, max = 255, alpha = 255)) 

##simantiko species arrows and labels##
-----------------------------------------------------------------
  ##diaforika afthona##
  #small
  arrows(0,0,0.5*mynmdsspe[c("ASV00493","ASV03432","ASV04749","ASV05606","ASV06706","ASV16590"),1] , 0.5*mynmdsspe[c("ASV00493","ASV03432","ASV04749","ASV05606","ASV06706","ASV16590"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#big
arrows(0,0,1.5*mynmdsspe[c("ASV00493","ASV03432","ASV04749","ASV05606","ASV06706","ASV16590"),1] , 1.5*mynmdsspe[c("ASV00493","ASV03432","ASV04749","ASV05606","ASV06706","ASV16590"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#D.A from Kruskal
("ASV00493","ASV03432","ASV04749","ASV05606","ASV06706","ASV16590")

[1] "ASV00493"
[1] "ASV03432"
[1] "ASV04749"
[1] "ASV05606"
[1] "ASV06706"
[1] "ASV16590"

##LABELS Short-cut##
a <- c("Clostridium_sp.","Stenotrophomonas_sp.","Mucilaginibacter_sp.","Syntrophomonadaceae_sp.","Roseomonas_sp.","Geodermatophilus_sp.")

install.packages("plotrix")

plotrix::thigmophobe.labels(1*mynmdsspe[c("ASV00493","ASV03432","ASV04749","ASV05606","ASV06706","ASV16590"),1], 1*mynmdsspe[c("ASV00493","ASV03432","ASV04749","ASV05606","ASV06706","ASV16590"),2], labels = a, cex = .6, font = 2, col = rgb(120,120,120,max = 255, alpha =200)) 

#simantiko# plotrix::thigmophobe.labels(1.2*mynmdsspe[1:15,1], 1.2*mynmdsspe[1:15,2], labels = mytax_tbl[row.names(mynmdsspe)[1:15],"forplt"], cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
# the following prints the legend
graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()


##_.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.--..--.-..-.-.##

##S4 Prepared for duplicate Spontaneous
MUST_Bacteria_Vidiano_2020_S4 <- subset_samples(MUST_Bacteria_Vidiano_2020_Genus, stage %in% (c("S4")))

MUST_Bacteria_Vidiano_2020_S4 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_S4)>0,MUST_Bacteria_Vidiano_2020_S4)

View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_S4)))

MUST_Bacteria_Vidiano_2020_S4 ##UNmerged S2 without spontaneous duplicate

##Duplicate Preperation
MUST_Bacteria_Vidiano_2020_S_S4 <- subset_samples(MUST_Bacteria_Vidiano_2020_Genus,replicate %in% (c("FS4")))

MUST_Bacteria_Vidiano_2020_S_S4 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_S_S4)>0,MUST_Bacteria_Vidiano_2020_S_S4)

MUST_Bacteria_Vidiano_2020_S_S4 ##Duplicate of spontaneous S4 (single sample) (119 taxa) (must be trimed so to be different from the original)

bS4 <- MUST_Bacteria_Vidiano_2020_S_S4 #119-1=118

myTaxa_bS4 <- names(sort(taxa_sums(bS4), decreasing = TRUE)[1:118]) 

Top_bS4 <- prune_taxa(myTaxa_bS4, bS4)

Top_bS4 <- prune_taxa(taxa_sums(Top_bS4)>0,Top_bS4)

sample_names(Top_bS4) <- paste("samp",sample_names(bS4), sep = "")

View(data.frame(sample_data(Top_bS4)))

##Merging Main Object with Spontaneous Duplicate#
c_S4 <- merge_phyloseq(MUST_Bacteria_Vidiano_2020_S4, Top_bS4)

View(data.frame(sample_data(c_S4)))

c_S4
##edo xekinaei to D.A part daneismeno apo to script ton heatmaps
c_S4 ##294 taxa

psdt <- c_S4

##transform phyloseq object raw counts to relative abundance##
psdt_ra <- transform_sample_counts(psdt, function(x) x / sum(x))

##Define the number of OTUs to be used in the analysis, 200 ASVs in our case## 
mynumOTUs <- 294 #294 taxa

## Prepare working object and select experimental variable of interest##
ps_htmp <- prune_taxa(names(taxa_sums(psdt_ra)[order(taxa_sums(psdt_ra), decreasing = T)][1:mynumOTUs]),psdt_ra)

## Prepare experimental variables table##
mydesign <- data.frame(sample_data(ps_htmp))
## Prepare the relative abundance table for the top 200 ASVs
dtcounts <- decostand(data.frame(otu_table(ps_htmp))[row.names(mydesign),], method = "range", MARGIN = 2)
## ASVs 
myotus <- colnames(dtcounts)

## Run Kruskal-Wallis, non-parametrix multiple comparison test##
mystatsout <- list()

for(myotu in myotus){
  if(sum(dtcounts[,myotu]) == 0){
    mystatsout[[myotu]]$krusk$statistics$p.chisq <- 1
  } else {
    mykrusk <- kruskal(dtcounts[,myotu], mydesign[rownames(dtcounts),]$vinification, group = T)
    mystatsout[[myotu]][["krusk"]] <- mykrusk
  }
}

## Extract p-values from the test results without p.adj##
mystatsoutkruskpvals <- unlist(lapply(myotus, function(x) mystatsout[[x]]$krusk$statistics$p.chisq))
names(mystatsoutkruskpvals) <- myotus

## Adjust p-values with "False Discovery Rate" method
mystatsoutkruskpvals.adj <- p.adjust(mystatsoutkruskpvals, method = "fdr")

# Check if there are any ASVs that significantly differ with p.adjust
for (i in c(1:length(mystatsoutkruskpvals.adj))) {
  if (mystatsoutkruskpvals.adj[i] < 0.05) {
    print(names(mystatsoutkruskpvals.adj[i]))
  }
}
###without p.adjust
for (i in c(1:length(mystatsoutkruskpvals))) {
  if (mystatsoutkruskpvals[i] < 0.06) {
    print(names(mystatsoutkruskpvals[i]))
  }
}
###without p.adjust finally gave results
[1] "ASV00175"
[1] "ASV01690"
[1] "ASV01996"
[1] "ASV02073"
[1] "ASV02606"
[1] "ASV02784"
[1] "ASV03295"
[1] "ASV04107"
[1] "ASV04485"
[1] "ASV06844"
[1] "ASV09857"
[1] "ASV10141"
[1] "ASV10145"
[1] "ASV12124"
[1] "ASV12630"
[1] "ASV16187"
[1] "ASV16188"
[1] "ASV21224"
[1] "ASV29178"
[1] "ASV40863"

##S4_Triplot_Start_with D.A taxa
VolatileVidiano2020_S4 <- read.table("VolatileVidiano2020_S4.txt", header=T,sep = "\t",row.names = 1)

fitVidiano_S4 <- envfit(ord.nmds.brayS4, VolatileVidiano2020_S4[,])
# prepare the names for the to be plotted taxa and save them in a new column
mytax_tbl <- data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S4_100))
# prepare the plotting colours
myplotcols <- RColorBrewer::brewer.pal(n = 6, name = 'RdBu')
# prep the plot points
mynmdssit <- scores(ord.nmds.brayS4, display = "site")
# prep the species points
mynmdsspe <- scores(ord.nmds.brayS4, display = "species")
# prep the variable of interest (e.g. vinification) has to be a factor
my_sel_var <- factor(MUST_Bacteria_Vidiano_2020_Merged_S4_100@sam_data$vinification)
# convert it also into numbers in orde to use it for colours etc
myterr_sel <- as.numeric(my_sel_var)
# get the arrow data that came out of the fitting of the environmental variables
arrowdata <- data.frame(fitVidiano_S4$vectors$arrows)

# start the graphics device
cairo_pdf("Bacteria_Vidiano_2020_Volatile_S4.DA.pdf", height = 7, width = 7)

# draw and empty plot that you will start to populate with points, arrows, ellipses, names etc.
plot(mynmdssit, frame = F, cex = 0, pch = 21, xlim = c(min(1.5*mynmdssit[,1]),max(1.5*mynmdssit[,1])))

# first add the ellipses in order tto keep them in the background and prevent them to obstruct information
vegan::ordiellipse(mynmdssit, groups = my_sel_var, kind = "ehull", lty = 2, lwd=1)

# add the sample points and give them colors one for each treatment
points(mynmdssit, bg = myplotcols[myterr_sel], pch = 21, cex = 1.5
       ##gia keimeno sta simeia##
       , text(mynmdssit, labels = row.names(mynmdssit), font = 6, col = rgb(55,55,55, max = 255, alpha = 100)))

par(adj = 1)
title(sub = "Bacteria_Vidiano_2020_Volatile_S4.DA")
par(adj = .5)

# add the env parameter arrows (in the rgb colouring mode we have the red, green and blue values from 0 to 255 and also the alpha value or transparency which also ranges between 0 and 255)
#small
arrows(0,0,0.3*arrowdata[,1] , 0.3*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)
#big
arrows(0,0,1.5*arrowdata[,1] , 1.5*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)

# also the labeling text
text(0.5*arrowdata[,1], 0.5*arrowdata[,2], labels = row.names(arrowdata), cex = 1.3, font = 2, col = rgb(20,20,20, max = 255, alpha = 255)) 

##simantiko species arrows and labels##
-----------------------------------------------------------------
  ##diaforika afthona##
  #small
  arrows(0,0,0.5*mynmdsspe[c("ASV00175","ASV01690","ASV01996","ASV02073","ASV02606","ASV02784","ASV03295","ASV04107","ASV04485","ASV06844","ASV09857","ASV10141","ASV10145","ASV12124","ASV12630","ASV16187","ASV16188","ASV21224","ASV29178","ASV40863"),1] , 0.5*mynmdsspe[c("ASV00175","ASV01690","ASV01996","ASV02073","ASV02606","ASV02784","ASV03295","ASV04107","ASV04485","ASV06844","ASV09857","ASV10141","ASV10145","ASV12124","ASV12630","ASV16187","ASV16188","ASV21224","ASV29178","ASV40863"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#big
arrows(0,0,1.5*mynmdsspe[c("ASV00175","ASV01690","ASV01996","ASV02073","ASV02606","ASV02784","ASV03295","ASV04107","ASV04485","ASV06844","ASV09857","ASV10141","ASV10145","ASV12124","ASV12630","ASV16187","ASV16188","ASV21224","ASV29178","ASV40863"),1] , 1.5*mynmdsspe[c("ASV00175","ASV01690","ASV01996","ASV02073","ASV02606","ASV02784","ASV03295","ASV04107","ASV04485","ASV06844","ASV09857","ASV10141","ASV10145","ASV12124","ASV12630","ASV16187","ASV16188","ASV21224","ASV29178","ASV40863"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#D.A from Kruskal
("ASV00175","ASV01690","ASV01996","ASV02073","ASV02606","ASV02784","ASV03295","ASV04107","ASV04485","ASV06844","ASV09857","ASV10141","ASV10145","ASV12124","ASV12630","ASV16187","ASV16188","ASV21224","ASV29178","ASV40863")

[1] "ASV00175"
[1] "ASV01690"
[1] "ASV01996"
[1] "ASV02073"
[1] "ASV02606"
[1] "ASV02784"
[1] "ASV03295"
[1] "ASV04107"
[1] "ASV04485"
[1] "ASV06844"
[1] "ASV09857"
[1] "ASV10141"
[1] "ASV10145"
[1] "ASV12124"
[1] "ASV12630"
[1] "ASV16187"
[1] "ASV16188"
[1] "ASV21224"
[1] "ASV29178"
[1] "ASV40863"

##LABELS Short-cut##
a <- c("Leuconostoc_sp.","Phenylobacterium_sp.","Ferrovibrionales_sp.","Rothia_sp.","Pluralibacter_sp.","Mycobacterium_sp.","Gallicola_sp.","Erwiniaceae_sp.","Actinomyces_sp.","Granulicatella_sp.","Sedimentibacter_sp.","Pseudactinotalea_sp.","Limnochordia_sp.","Hyphomonadaceae_sp.","Tepidisphaerales_sp.","Azospirillum_sp.","Rickettsiales_sp.","Dysgonomonadaceae_sp.","Nocardioidaceae_sp.","Flavobacteriaceae_sp.")

install.packages("plotrix")

plotrix::thigmophobe.labels(1*mynmdsspe[c("ASV00175","ASV01690","ASV01996","ASV02073","ASV02606","ASV02784","ASV03295","ASV04107","ASV04485","ASV06844","ASV09857","ASV10141","ASV10145","ASV12124","ASV12630","ASV16187","ASV16188","ASV21224","ASV29178","ASV40863"),1], 1*mynmdsspe[c("ASV00175","ASV01690","ASV01996","ASV02073","ASV02606","ASV02784","ASV03295","ASV04107","ASV04485","ASV06844","ASV09857","ASV10141","ASV10145","ASV12124","ASV12630","ASV16187","ASV16188","ASV21224","ASV29178","ASV40863"),2], labels = a, cex = .6, font = 2, col = rgb(120,120,120,max = 255, alpha =200)) 

#simantiko# plotrix::thigmophobe.labels(1.2*mynmdsspe[1:15,1], 1.2*mynmdsspe[1:15,2], labels = mytax_tbl[row.names(mynmdsspe)[1:15],"forplt"], cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
# the following prints the legend
graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()

##_.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.--..--.-..-.-.##

##S5 Prepared for duplicate Spontaneous
MUST_Bacteria_Vidiano_2020_S5 <- subset_samples(MUST_Bacteria_Vidiano_2020_Genus, stage %in% (c("S5")))

MUST_Bacteria_Vidiano_2020_S5 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_S5)>0,MUST_Bacteria_Vidiano_2020_S5)

View(data.frame(sample_data(MUST_Bacteria_Vidiano_2020_S5)))

MUST_Bacteria_Vidiano_2020_S5 ##UNmerged S5 without spontaneous duplicate

##Duplicate Preperation
MUST_Bacteria_Vidiano_2020_S_S5 <- subset_samples(MUST_Bacteria_Vidiano_2020_Genus,replicate %in% (c("FS5")))

MUST_Bacteria_Vidiano_2020_S_S5 <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2020_S_S5)>0,MUST_Bacteria_Vidiano_2020_S_S5)

MUST_Bacteria_Vidiano_2020_S_S5 ##Duplicate of spontaneous S4 (single sample) (142 taxa) (must be trimed so to be different from the original)

bS5 <- MUST_Bacteria_Vidiano_2020_S_S5 #142-1=141

myTaxa_bS5 <- names(sort(taxa_sums(bS5), decreasing = TRUE)[1:141]) 

Top_bS5 <- prune_taxa(myTaxa_bS5, bS5)

Top_bS5 <- prune_taxa(taxa_sums(Top_bS5)>0,Top_bS5)

sample_names(Top_bS5) <- paste("samp",sample_names(bS5), sep = "")

View(data.frame(sample_data(Top_bS5)))

##Merging Main Object with Spontaneous Duplicate#
c_S5 <- merge_phyloseq(MUST_Bacteria_Vidiano_2020_S5, Top_bS5)

View(data.frame(sample_data(c_S5)))

c_S5
##edo xekinaei to D.A part daneismeno apo to script ton heatmaps
c_S5 ##120 taxa

psdt <- c_S5

##transform phyloseq object raw counts to relative abundance##
psdt_ra <- transform_sample_counts(psdt, function(x) x / sum(x))

##Define the number of OTUs to be used in the analysis, 200 ASVs in our case## 
mynumOTUs <- 354 #354 taxa

## Prepare working object and select experimental variable of interest##
ps_htmp <- prune_taxa(names(taxa_sums(psdt_ra)[order(taxa_sums(psdt_ra), decreasing = T)][1:mynumOTUs]),psdt_ra)

## Prepare experimental variables table##
mydesign <- data.frame(sample_data(ps_htmp))
## Prepare the relative abundance table for the top 200 ASVs
dtcounts <- decostand(data.frame(otu_table(ps_htmp))[row.names(mydesign),], method = "range", MARGIN = 2)
## ASVs 
myotus <- colnames(dtcounts)

## Run Kruskal-Wallis, non-parametrix multiple comparison test##
mystatsout <- list()

for(myotu in myotus){
  if(sum(dtcounts[,myotu]) == 0){
    mystatsout[[myotu]]$krusk$statistics$p.chisq <- 1
  } else {
    mykrusk <- kruskal(dtcounts[,myotu], mydesign[rownames(dtcounts),]$vinification, group = T)
    mystatsout[[myotu]][["krusk"]] <- mykrusk
  }
}

## Extract p-values from the test results without p.adj##
mystatsoutkruskpvals <- unlist(lapply(myotus, function(x) mystatsout[[x]]$krusk$statistics$p.chisq))
names(mystatsoutkruskpvals) <- myotus

## Adjust p-values with "False Discovery Rate" method
mystatsoutkruskpvals.adj <- p.adjust(mystatsoutkruskpvals, method = "fdr")

# Check if there are any ASVs that significantly differ with p.adjust
for (i in c(1:length(mystatsoutkruskpvals.adj))) {
  if (mystatsoutkruskpvals.adj[i] < 0.05) {
    print(names(mystatsoutkruskpvals.adj[i]))
  }
}
###without p.adjust
for (i in c(1:length(mystatsoutkruskpvals))) {
  if (mystatsoutkruskpvals[i] < 0.06) {
    print(names(mystatsoutkruskpvals[i]))
  }
}
###without p.adjust finally gave results
[1] "ASV00568"
[1] "ASV00622"
[1] "ASV03753"
[1] "ASV04624"
[1] "ASV04644"
[1] "ASV05089"
[1] "ASV06572"
[1] "ASV06592"
[1] "ASV06829"
[1] "ASV07465"
[1] "ASV08341"
[1] "ASV15440"
[1] "ASV17466"
[1] "ASV17467"
[1] "ASV17944"
[1] "ASV18953"
[1] "ASV18954"
[1] "ASV19496"
[1] "ASV21907"
[1] "ASV21909"
[1] "ASV22630"
[1] "ASV23351"
[1] "ASV24172"
[1] "ASV30490"
[1] "ASV33428"


##S5_Triplot_Start_with D.A taxa
VolatileVidiano2020_S5 <- read.table("VolatileVidiano2020_S5.txt", header=T,sep = "\t",row.names = 1)

fitVidiano_S5 <- envfit(ord.nmds.brayS5, VolatileVidiano2020_S5[,])
# prepare the names for the to be plotted taxa and save them in a new column
mytax_tbl <- data.frame(tax_table(MUST_Bacteria_Vidiano_2020_Merged_S5_100))
# prepare the plotting colours
myplotcols <- RColorBrewer::brewer.pal(n = 6, name = 'RdBu')
# prep the plot points
mynmdssit <- scores(ord.nmds.brayS5, display = "site")
# prep the species points
mynmdsspe <- scores(ord.nmds.brayS5, display = "species")
# prep the variable of interest (e.g. vinification) has to be a factor
my_sel_var <- factor(MUST_Bacteria_Vidiano_2020_Merged_S5_100@sam_data$vinification)
# convert it also into numbers in orde to use it for colours etc
myterr_sel <- as.numeric(my_sel_var)
# get the arrow data that came out of the fitting of the environmental variables
arrowdata <- data.frame(fitVidiano_S5$vectors$arrows)

# start the graphics device
cairo_pdf("Bacteria_Vidiano_2020_Volatile_S5.DA.pdf", height = 7, width = 7)

# draw and empty plot that you will start to populate with points, arrows, ellipses, names etc.
plot(mynmdssit, frame = F, cex = 0, pch = 21, xlim = c(min(1.5*mynmdssit[,1]),max(1.5*mynmdssit[,1])))

# first add the ellipses in order tto keep them in the background and prevent them to obstruct information
vegan::ordiellipse(mynmdssit, groups = my_sel_var, kind = "ehull", lty = 2, lwd=1)

# add the sample points and give them colors one for each treatment
points(mynmdssit, bg = myplotcols[myterr_sel], pch = 21, cex = 1.5
       ##gia keimeno sta simeia##
       , text(mynmdssit, labels = row.names(mynmdssit), font = 6, col = rgb(55,55,55, max = 255, alpha = 100)))

par(adj = 1)
title(sub = "Bacteria_Vidiano_2020_Volatile_S5.DA")
par(adj = .5)

# add the env parameter arrows (in the rgb colouring mode we have the red, green and blue values from 0 to 255 and also the alpha value or transparency which also ranges between 0 and 255)
#small
arrows(0,0,0.3*arrowdata[,1] , 0.3*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)
#big
arrows(0,0,1.5*arrowdata[,1] , 1.5*arrowdata[,2], angle = 25, length = 0.15, col = rgb(20,20,20, max = 255, alpha = 100), lwd = 3)

# also the labeling text
text(0.5*arrowdata[,1], 0.5*arrowdata[,2], labels = row.names(arrowdata), cex = 1.3, font = 2, col = rgb(20,20,20, max = 255, alpha = 255)) 

##simantiko species arrows and labels##
-----------------------------------------------------------------
  ##diaforika afthona##
  #small
  arrows(0,0,0.5*mynmdsspe[c("ASV00568","ASV00622","ASV03753","ASV04624","ASV04644","ASV05089","ASV06572","ASV06592","ASV06829","ASV07465","ASV08341","ASV15440","ASV17466","ASV17467","ASV17944","ASV18953","ASV18954","ASV19496","ASV21907","ASV21909","ASV22630","ASV23351","ASV24172","ASV30490","ASV33428"),1] , 0.5*mynmdsspe[c("ASV00568","ASV00622","ASV03753","ASV04624","ASV04644","ASV05089","ASV06572","ASV06592","ASV06829","ASV07465","ASV08341","ASV15440","ASV17466","ASV17467","ASV17944","ASV18953","ASV18954","ASV19496","ASV21907","ASV21909","ASV22630","ASV23351","ASV24172","ASV30490","ASV33428"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#big
arrows(0,0,1.5*mynmdsspe[c("ASV00568","ASV00622","ASV03753","ASV04624","ASV04644","ASV05089","ASV06572","ASV06592","ASV06829","ASV07465","ASV08341","ASV15440","ASV17466","ASV17467","ASV17944","ASV18953","ASV18954","ASV19496","ASV21907","ASV21909","ASV22630","ASV23351","ASV24172","ASV30490","ASV33428"),1] , 1.5*mynmdsspe[c("ASV00568","ASV00622","ASV03753","ASV04624","ASV04644","ASV05089","ASV06572","ASV06592","ASV06829","ASV07465","ASV08341","ASV15440","ASV17466","ASV17467","ASV17944","ASV18953","ASV18954","ASV19496","ASV21907","ASV21909","ASV22630","ASV23351","ASV24172","ASV30490","ASV33428"),2], angle = 25, length = 0.20, col = rgb(40,40,40, max = 255, alpha = 100))

#D.A from Kruskal
("ASV00568","ASV00622","ASV03753","ASV04624","ASV04644","ASV05089","ASV06572","ASV06592","ASV06829","ASV07465","ASV08341","ASV15440","ASV17466","ASV17467","ASV17944","ASV18953","ASV18954","ASV19496","ASV21907","ASV21909","ASV22630","ASV23351","ASV24172","ASV30490","ASV33428")

[1] "ASV00568"
[1] "ASV00622"
[1] "ASV03753"
[1] "ASV04624"
[1] "ASV04644"
[1] "ASV05089"
[1] "ASV06572"
[1] "ASV06592"
[1] "ASV06829"
[1] "ASV07465"
[1] "ASV08341"
[1] "ASV15440"
[1] "ASV17466"
[1] "ASV17467"
[1] "ASV17944"
[1] "ASV18953"
[1] "ASV18954"
[1] "ASV19496"
[1] "ASV21907"
[1] "ASV21909"
[1] "ASV22630"
[1] "ASV23351"
[1] "ASV24172"
[1] "ASV30490"
[1] "ASV33428"


##LABELS Short-cut##
a <- c("Zymobacter_sp.","Romboutsia_sp.","Jatrophihabitans_sp.","Cyanobacteria_sp.","Blastococcus_sp.","Candidatus_Ovatusbacter_sp.","Rhodocyclaceae_sp.","Fusobacterium_sp.","Anaerocolumna_sp.","Polyangium_sp.","Pirellula_sp.","Christensenellaceae_sp.","Myxococcaceae_sp.","Roseibacillus_sp.","Rhizobiales_sp.","Hydrogenoanaerobacterium_sp.","Oscillospiraceae_sp.","Alkalibaculum_sp.","Fermentimonas_sp.","Candidatus_Nucleicultrix_sp.","Schlesneria_sp.","Actinotignum_sp.","Ruminococcaceae_sp.","Aureimonas_sp.","Undibacterium_sp.")

install.packages("plotrix")

plotrix::thigmophobe.labels(1*mynmdsspe[c("ASV00568","ASV00622","ASV03753","ASV04624","ASV04644","ASV05089","ASV06572","ASV06592","ASV06829","ASV07465","ASV08341","ASV15440","ASV17466","ASV17467","ASV17944","ASV18953","ASV18954","ASV19496","ASV21907","ASV21909","ASV22630","ASV23351","ASV24172","ASV30490","ASV33428"),1], 1*mynmdsspe[c("ASV00568","ASV00622","ASV03753","ASV04624","ASV04644","ASV05089","ASV06572","ASV06592","ASV06829","ASV07465","ASV08341","ASV15440","ASV17466","ASV17467","ASV17944","ASV18953","ASV18954","ASV19496","ASV21907","ASV21909","ASV22630","ASV23351","ASV24172","ASV30490","ASV33428"),2], labels = a, cex = .6, font = 2, col = rgb(120,120,120,max = 255, alpha =200)) 

#simantiko# plotrix::thigmophobe.labels(1.2*mynmdsspe[1:15,1], 1.2*mynmdsspe[1:15,2], labels = mytax_tbl[row.names(mynmdsspe)[1:15],"forplt"], cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
# the following prints the legend
graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()

















##S4 Prepared for duplicate Spontaneous

##S5 Prepared for duplicate Spontaneous


