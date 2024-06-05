
# 1. Reads mapping 

# 2. Variant calling

# 3. Population genetics

## 3.1 Admixture and phylogenetic tree
``` R
library("tidyverse")
library("reshape2")
library("ggtree")
library("treeio")
library("ggstance")

admixture_order <- read.table("./Admixture/09005/add1913.nosex") %>% 
  mutate(Sample = V1) %>% 
  select(Sample) 

admixture_dir <- dir("./Admixture/09005/", pattern = ".Q", full.names = T)

admixture_qmatrix <- purrr::map(
    admixture_dir, ~read.table(., header = F) %>% 
    cbind(admixture_order, .) %>% 
    melt()
    ) %>%
  setNames(paste0("K",c(10:12,2:9)))

mltree <- treeio::read.newick("./MLTree/09005_hapAB.treefile", node.label = "label")

mltree_plot <-  ggtree(mltree, size = 0.2) + 
  geom_tiplab(size = 1, offset = 1, as_ylab = T) 

mltree_plot_combine <- mltree_plot

for (i in 2:8){
  Kn <- paste0("K",i)
  print(Kn)
  mltree_plot_combine <- facet_plot(
    mltree_plot_combine, 
    panel = Kn, 
    data = admixture_qmatrix[[Kn]], 
    geom = geom_barh, 
    width = 1, 
    mapping = aes(x = value, fill = as.factor(variable)), 
    color = "white", stat = "identity")  + 
    theme(legend.position = 'none') 
}

mltree_plot_combine

ggsave(mltree_plot_combine, file = "./Admixture/admixture_tree.pdf", width = 8, height =8)

```
## 3.2PCA
```R
col_L6 <- c("L1" = "#cc79a7", "L2" = "#4caf49", "L3" = "#fdbf70", "L4" = "#e66434", "L5" = "#4dbbd6", "L6" = "#974ea2")
popInfoDNA <- read.table("./Data/popInfo.txt", sep = "\t", header = T) 

pca09005_pc12 <- ggplot(data = pca_plink_09005, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill = Lineage), shape = 21, color = "black", size = 3) +
  scale_fill_manual(values = col_L6)+
  ylab("PC2 (20.91%)") +
  xlab("PC1 (34.08%)") + 
  theme_bw()+
  theme(axis.line = element_line(colour = "black"), 
                   axis.text.x = element_text(colour = "black"), 
                   axis.text.y = element_text(colour = "black"), 
                   axis.ticks = element_line(colour = "black"), 
                   axis.title.x = element_text(colour = "black"), 
                   axis.title.y = element_text(colour = "black")) 
ggsave(pca09005_pc12, file= "./PCA/hapA_09005pc12.pdf", width = 6, height = 4)
```
## 3.3 fst
```R

pairwise_fst <- read.table("./selection/fst/pairwise_fst_longdata", header = F)
pairwise_fst_tileplot <- ggplot(pairwise_fst, aes(x = V1, y = V2, fill = V3)) + 
  geom_tile(color = rgb(45, 67,121, 0, maxColorValue = 255))+
  geom_text(aes(label = V3), color = "black")+
  scale_y_discrete(limits = paste0("L", c(2:6)))+
  coord_fixed(expand = FALSE) +
  scale_fill_distiller(
    palette = "Reds",
    direction = 1, limits = c(0.1, 0.5),
    name = "Fst")+
  theme_classic()+  
  labs(x = NULL, y = NULL) +
  theme(panel.border = element_rect(color = NA, fill = NA),
        legend.position = c(1.2, 1.2),
        axis.ticks = element_blank(),
        axis.line = element_blank())

ggsave(pairwise_fst_tileplot, file = "./selection/fst/pairwise_fst_tileplot.pdf", width = 5, height = 5)

highfst_num <- read.table("./selection/fst/pairwise_highfst_num_longdata", header = F)
highfst_num_tileplot <- ggplot(highfst_num, aes(x = V1, y = V2, fill = V3)) + 
  geom_tile(color = rgb(45, 67,121, 0, maxColorValue = 255))+
  geom_text(aes(label = V3), color = "black")+
  scale_y_discrete(limits = paste0("L", c(2:6)))+
  coord_fixed(expand = FALSE) +
  scale_fill_distiller(
    palette = "Blues",
    direction = 1, limits = c(0, 25352),
    name = "Fst")+
  theme_classic()+  
  labs(x = NULL, y = NULL) +
  theme(panel.border = element_rect(color = NA, fill = NA),
        legend.position = c(1.2,1.2),
        axis.ticks = element_blank(),
        axis.line = element_blank())+
  coord_flip()
highfst_num_tileplot
```
# Selection - raisd
```R
## raisd ----------------------------------------------------------------------
library(magrittr)
### data prcessing ---------------------------------------------------------------
raisd_dir <- "./selection/raisd/data/"

raisd_data <- dir(raisd_dir, pattern = "Ppz_chr", full.names = T) %>%
  purrr::map(
    ., ~read.table(., header = F) %>%
      set_colnames(c("Pos", "Start", "End", "Var", "Sfs", "Ld", "mu")) %>%
      mutate(Haploid = paste0("hap", substr(.x, nchar(.x), nchar(.x)))) %>%
      mutate(Chr =substr(.x, str_locate(.x, "Ppz_chr.*")[1], str_locate(.x, "Ppz_chr.*")[2])) %>%
      mutate(Lineage =substr(.x, str_locate(.x, "L\\d")[1], str_locate(.x, "L\\d")[2]))
  )%>% 
  bind_rows() 

chr_length <- read.table("./selection/raisd/GD1913_hapAB.fasta.fai") %>%
  reshape2::melt(id.vars = "V1") %>% 
  select(-2) %>%
  set_colnames(c("Chr","Pos")) %>%
  mutate(Start = Pos-1, End = Pos, Var = 0, Sfs = 0, Ld = 0, mu = -0.02) %>%
  mutate(Haploid = ifelse(substr(Chr,nchar(Chr),nchar(Chr)) == "A", "hapA",ifelse(substr(Chr,nchar(Chr),nchar(Chr)) == "B", "hapB", NA)))

raisd_data <- bind_rows(raisd_data, chr_length)

raisd_threshold <- raisd_data %>% 
  group_by(Lineage, Haploid) %>% 
  summarise(Threshold = quantile(mu, probs = 0.99)) %>% 
  drop_na()

raisd_data <- raisd_data %>% 
  left_join(raisd_threshold, by = c("Lineage","Haploid")) %>% 
  mutate(Top = mu > Threshold) %>%
  mutate(length = End - Start)

raisd_data_selected <- raisd_data %>% 
  filter(Top == TRUE) 

### circle plot -------------------------------------------------------------------

raisd_data_selected_circ <- raisd_data_selected %>% 
  select(Chr, Start, End, mu, Lineage) %>% 
  mutate(Chr = substr(Chr, 8,10)) 

library(GenomicRanges)

gr <- GRanges(raisd_data_selected_circ)
gr_list <- split(gr, mcols(gr)$Lineage)
gr_reduced_list <- lapply(gr_list, IRanges::reduce)

gr_reduced <- unlist(gr_reduced_list)
gr_reduced_matrix <- list()

for(i in 1:length(gr_reduced)){
  df <- as.data.frame(gr_reduced[[i]])
  df$Lineage <- deparse(names(gr_reduced[i]))
  gr_reduced_matrix[[i]] <- df
}

gr_reduced_matrix <- do.call(rbind, gr_reduced_matrix) %>% 
  dplyr::rename("Chromosome" = "seqnames")
write.table(gr_reduced_matrix, file = "./selection/raisd/raisd_regions_combined.txt", quote = F, row.names = F, sep = "\t")

total_length <- gr_reduced_matrix %>% 
  group_by(Lineage) %>% 
  summarise(total_length = sum(width)) %>% 
  mutate(percent = total_length / 1700000000)

te_density <- read.table("selection/raisd/te_density2", header = T) %>% 
  mutate(chr = str_sub(chr, 8,11), start = start + 1, value = value/max(value)) %>% 
  filter(value != 0)

library("circlize")

pdf("./selection/raisd/plot/circle990-te.pdf", height = 5, width = 5)
circos.clear()
circos.par("clock.wise" = T, start.degree = 80, gap.after = c(rep(1,17), 1, rep(1,17), 20), cell.padding = c(0,1,0,1))
cytoband.df <- read.table(file = "./selection/raisd/GD1913_hapAB.fasta.fai", header = F, colClasses = c("character", "numeric", "numeric"),sep = "\t") %>% 
  mutate(V1 = str_replace_all(V1, "Ppz_chr", ""))
circos.genomicInitialize(
  cytoband.df,major.by = 50000000,
  axis.labels.cex = 0.25*par("cex"),
  labels.cex = 0.6*par("cex"))
chrcolor = c(rep("#fda277",18), rep("#89b6ff", 18))
circos.genomicTrackPlotRegion(cytoband.df,ylim=c(0,0.1),track.height=0.03,bg.col=chrcolor,bg.lwd=0.4,cell.padding=c(0.01, 0.5, 0.01, 0.5))
circos.genomicDensity(te_density, col = c("#FF8800"), track.height = 0.05, area = F, border = NA, type = "l")

for (lineage in names(gr_reduced_list)){
  gr_sub <- as.data.frame(gr_reduced_list[[lineage]])
  circos.genomicTrack(
    gr_sub, 
    track.height=0.05, 
    ylim = c(0,1), 
    bg.border = NA, bg.col = "#fdecdd",
    panel.fun = function(region, value, ...){
      circos.rect(
        xleft = region$start,
        xright = region$end, 
        ytop = rep(1,nrow(region)), 
        ybottom = rep(0,nrow(region)), 
        col = col_L6[lineage], border = NA)
    }
  )
}
dev.off()

```
