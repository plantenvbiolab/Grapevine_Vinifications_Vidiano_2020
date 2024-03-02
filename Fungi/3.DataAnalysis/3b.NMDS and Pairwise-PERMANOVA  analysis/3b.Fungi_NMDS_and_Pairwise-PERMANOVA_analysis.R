##nmds##
MUST_Fungi_Vidiano_2020_Merged_100

ord.nmds.bray1 <- ordinate(MUST_Fungi_Vidiano_2020_Merged_100, method="NMDS", distance="bray")


plot_ordination(MUST_Fungi_Vidiano_2020_Merged_100, ord.nmds.bray1, color="vinification", shape ="", label = "", title=paste("VIDIANO Fungi 2020")) + geom_point(size = 4) + scale_color_manual(values = mycols) + theme_classic()

plot_ordination(MUST_Fungi_Vidiano_2020_Merged_100, ord.nmds.bray1, color="stage", shape ="", label = "", title=paste("VIDIANO Fungi 2020")) + geom_point(size = 4) + scale_color_manual(values = mycols) + theme_classic()




plot_ordination(MUST_Fungi_Vidiano_2020_Merged_100, ord.nmds.bray1, color="vinification", shape ="", label = "stage", title=paste("VIDIANO FUNGI 2020")) + geom_point(size = 4) + scale_color_manual(values = "#66A61E","#E6AB02","#A6761D","#666666","#E41A1C")

plot_ordination(MUST_Fungi_Vidiano_2020_Merged_100, ord.nmds.bray1, color="stage", shape ="vinification", label = "stage", title=paste("VIDIANO FUNGI 2020")) + geom_point(size = 4) + scale_color_manual(values = mycols)

##pairwise permanova (sto RA% object)##
library(pairwiseAdonis)

mycmpfactor_MUST_Fungi_Vidiano_2020_Merged_100 <- interaction(data.frame(MUST_Fungi_Vidiano_2020_Merged_100@sam_data)$vinification, data.frame(MUST_Fungi_Vidiano_2020_Merged_100@sam_data)$vinification) 

mympairwiseperm_MUST_Fungi_Vidiano_2020_Merged_100 <- pairwise.adonis(MUST_Fungi_Vidiano_2020_Merged_100@otu_table, mycmpfactor_MUST_Fungi_Vidiano_2020_Merged_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")##variety vs variety
