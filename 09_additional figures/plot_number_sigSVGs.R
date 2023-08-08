###################################
# Plot barplots of significant SVGs detected
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(paletteer)
library(here)

load(here("data", "input_number_of_sigSVGs.RData"))
csv = read.csv(here("data", "dataset_metadata.csv"), head=T)

pal = paletteer::palettes_d$ggthemes$Jewel_Bright[-3]

# order datasets
o=order(csv$spatial_technology, csv$first_author)
csv = csv[o, ]

num <- rbind(giotto.km.num, giotto.r.num, meringue.num, moransi.num, nnSVG.num, somde.num, spark.num, spatialde.num)

fn = here("plots", "number_sigSVG_barplots.pdf")
pdf(fn, width=8, height=4)
barplot(num[,csv$dataset], beside = TRUE, las=2, ylab="Number of detected SVGs", col=pal)
dev.off()