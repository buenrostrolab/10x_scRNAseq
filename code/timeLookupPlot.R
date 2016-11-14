library(ggplot2)

pA <- c(seq(180,909,185)/909,1)
dA <- readRDS("../data/timeLookup/aml035_post_transplant_times.rds")

pB <- c(seq(400,1985,400)/1985,1)
dB <- readRDS("../data/timeLookup/frozen_bmmc_healthy_donor1_times.rds")

pC <- c(seq(800,3965,800)/3965, 1)
dC <- readRDS("../data/timeLookup/aml027_post_transplant_times.rds")

pD <- c(seq(1550,7783,1560)/7783, 1)
dD <- readRDS("../data/timeLookup/frozen_pbmc_donor_b_times.rds")

pE <- c(seq(2100,10085,2000)/10085, 1)
dE <- readRDS("../data/timeLookup/b_cells_times.rds")


df <- rbind(
      data.frame(fraction = pA, time = rowMeans(dA), ncells = as.character(909)),
      data.frame(fraction = pB, time = rowMeans(dB), ncells = as.character(1985)),
      data.frame(fraction = pC, time = rowMeans(dC), ncells = as.character(3965)),
      data.frame(fraction = pD, time = rowMeans(dD), ncells = as.character(7783)),
      data.frame(fraction = pE, time = rowMeans(dE), ncells = as.character(10085)))

df$ncells <- as.character(df$ncells)

ggplot(df, aes(x=fraction, y=time, color = ncells)) + theme_bw()+ labs(title = "SQLite Sparse Lookup Speeds", 
    x = "Proportion of Data", y = "System Time (seconds)") + geom_line()+labs(colour = "# Cells in Table")