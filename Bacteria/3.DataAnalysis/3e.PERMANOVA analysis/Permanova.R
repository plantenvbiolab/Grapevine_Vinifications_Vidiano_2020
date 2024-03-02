PERMANOVA Variables (perm. 999)
_______________________________________________________________________
library(vegan)

##PERMANOVA Variables
mypermanova_MUST_Bacteria_Vidiano_2020_Merged_100 <- adonis(MUST_Bacteria_Vidiano_2020_Merged_100@otu_table ~ vinification + stage, method = "bray", data = data.frame(MUST_Bacteria_Vidiano_2020_Merged_100@sam_data))

write.table(data.frame(mypermanova_MUST_Bacteria_Vidiano_2020_Merged_100), file="mypermanova_MUST_Bacteria_Vidiano_2020_Merged_100F.txt", quote = F,col.names = NA, sep="\t")