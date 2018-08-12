options(stringsAsFactors = F)
setwd("//eristore03/Department/HBDS/xulong/CATS/GCTA") 

igap = read.table("./IGAP_stage_1.txt")
names(igap) = c("Chr", "Pos", "Marker", "EA", "NEA", "b", "se", "p")
igap$N = 70000

g1k = read.table("../../1000G/plink.frq", header = T)
g1k = g1k[match(igap$Marker, g1k$SNP), ]

igap$freq = g1k$MAF
igap = igap[! is.na(igap$freq), ]

igap$b = exp(igap$b)
igap = igap[c("Marker", "EA", "NEA", "freq", "b", "se", "p", "N")]

write.table(igap, file = "input.ma", quote = F, sep = "\t", row.names = F)

rs10421385 = c(0.05, 0.4117)
rs12459419 = c(-0.09, -1.133)

x = as.data.frame(rbind(rs10421385, rs12459419))
names(x) = c("AD", "CD33")

plot(x$CD33, x$AD, 
     xlab = "CD33 protein changes", ylab = "AD susceptibility", col = "red", 
     xlim = c(-1.5, 1.5),
     ylim = c(-0.1, 0.1))
abline(v = 0)
abline(h = 0)

load("//eristore03/Department/HBDS/IGAP/IGAP_stage_1.rdt") 

options(stringsAsFactors = F)
pqtl = read.table("//eristore03/Department/HBDS/xulong/CATS/Cytokine/3166-92_1_one.out")

# pqtl = pqtl[pqtl$V17 < 5e-8, ]

stage1 = stage1[stage1$MarkerName %in% pqtl$V10, ]
pqtl = pqtl[match(stage1$MarkerName, pqtl$V10), ]

head(pqtl)
head(stage1)

.libPaths("C:/Users/esi16740/R/library")
library(ggplot2)
library(ggrepel)

y = data.frame(marker = stage1$MarkerName, 
               AD = stage1$Beta, 
               p = stage1$Pvalue, 
               CD33 = pqtl$V15)

y = y[y$marker %in% c("rs12459419",
                      "rs12985029",
                      "rs10421385",
                      "rs2455069"), ]

y = rbind(y, data.frame(marker = "rs1710354", AD = -0.0127, p = 0.4538, CD33 = -0.3874))

ggplot(y, aes(x = CD33, y = AD)) + 
  geom_point(aes(size = -log10(p))) + 
  geom_text_repel(aes(label = marker)) +
  theme_bw() + xlab("Plasma CD33 level (Suhre et al., 2017)") + ylab("AD susceptibility (Lambert et al., 2013)")


