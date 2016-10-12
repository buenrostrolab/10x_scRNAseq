# Goal of this script is to get hg19 coordinates for all the
# genes in 'genes.tsv'. When multiple options are present, 
# we retain the first. 

mart <- read.table("mart_export.txt", sep = "\t", header = TRUE)
names(mart) <- c("EGID", "Chrom", "Start", "Stop", "Strand")

genes <- read.table("genes.tsv", sep = "\t")
names(genes) <- c("EGID", "GENEID")

mm <- merge(genes, mart, by.x = "EGID", by.y = "EGID", sort = FALSE, all.x = TRUE)
unique_mm_chord <- mm[!duplicated(mm[,c(1,2)]), ]

sum(duplicated(unique_mm_chord[,1])) # unique IDs so keep these

dfout <- unique_mm_chord[,c(3,4,5,6,1)]
write.table(dfout, file = "hg19_coords.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
