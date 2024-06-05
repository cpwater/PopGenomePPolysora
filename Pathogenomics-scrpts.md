
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