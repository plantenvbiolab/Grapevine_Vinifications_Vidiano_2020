#MERGED REPLICATE##
MUST_Fungi_Vidiano_2020_Merged <- merge_samples(fungi_vinification_Annotated, "replicate")

MUST_Fungi_Vidiano_2020_Merged <- prune_taxa(taxa_sums(MUST_Fungi_Vidiano_2020_Merged)>0,MUST_Fungi_Vidiano_2020_Merged)

View(data.frame(sample_data(MUST_Fungi_Vidiano_2020_Merged)))

write.table(data.frame(sample_data(MUST_Fungi_Vidiano_2020_Merged)), file="MUST_Fungi_Vidiano_2020_Merged.txt", quote = F,col.names = NA, sep="\t")

SampleDataNew67 <- read.table("MUST_Fungi_Vidiano_2020_Merged.txt", header=T,sep = "\t",row.names = 1)

sample_data(MUST_Fungi_Vidiano_2020_Merged) <- SampleDataNew67

MUST_Fungi_Vidiano_2020_Merged

View(data.frame(sample_data(MUST_Fungi_Vidiano_2020_Merged)))

MUST_Fungi_Vidiano_2020_Merged <- prune_taxa(taxa_sums(MUST_Fungi_Vidiano_2020_Merged)>0,MUST_Fungi_Vidiano_2020_Merged)

saveRDS(MUST_Fungi_Vidiano_2020_Merged, file = "MUST_Fungi_Vidiano_2020_Merged.RDS")

MUST_Fungi_Vidiano_2020_Merged_100 <- transform_sample_counts(MUST_Fungi_Vidiano_2020_Merged, function(OTU) 100*OTU/sum(OTU))

write.table(data.frame(otu_table(MUST_Fungi_Vidiano_2020_Merged_100)), file="MUST_Fungi_Vidiano_2020_Merged_100_otu.txt", quote = F,col.names = NA, sep="\t")

write.table(data.frame(tax_table(MUST_Fungi_Vidiano_2020_Merged_100)), file="MUST_Fungi_Vidiano_2020_Merged_100_tax.txt", quote = F,col.names = NA, sep="\t")

plot_bar(MUST_Fungi_Vidiano_2020_Merged_100, x="stage", fill="Genus", title = "VIDIANO FUNGI 2020") + facet_grid(cols = vars(vinification)) + geom_col()

##ASVS TOP 11
MUST_Fungi_Vidiano_2020_Merged_100

myTaxa11_MUST_Fungi_Vidiano_2020_Merged_100 <- names(sort(taxa_sums(MUST_Fungi_Vidiano_2020_Merged_100), decreasing = TRUE)[1:12]) 

Top11_MUST_Fungi_Vidiano_2020_Merged_100 <- prune_taxa(myTaxa11_MUST_Fungi_Vidiano_2020_Merged_100, MUST_Fungi_Vidiano_2020_Merged_100)

plot_bar(Top11_MUST_Fungi_Vidiano_2020_Merged_100, x="stage", fill="OTU", title = "VIDIANO FUNGI ASVs 2020") + facet_grid(cols = vars(vinification)) + geom_col() + scale_fill_manual(values = mycols)

***
  #Genus Glommed#
  MUST_Fungi_Vidiano_2020_Merged_100

MUST_Fungi_Vidiano_2020_Merged_100_Genus <- tax_glom(MUST_Fungi_Vidiano_2020_Merged_100, taxrank = "Genus")

myTaxa11_MUST_Fungi_Vidiano_2020_Merged_100_Genus <- names(sort(taxa_sums(MUST_Fungi_Vidiano_2020_Merged_100_Genus), decreasing = TRUE)[1:11]) 

Top11_MUST_Fungi_Vidiano_2020_Merged_100_Genus <- prune_taxa(myTaxa11_MUST_Fungi_Vidiano_2020_Merged_100_Genus, MUST_Fungi_Vidiano_2020_Merged_100_Genus)

plot_bar(Top11_MUST_Fungi_Vidiano_2020_Merged_100_Genus, x="stage", fill="Genus", title = "VIDIANO FUNGI 2020") + facet_grid(cols = vars(vinification), rows = vars(year)) + geom_col() 

