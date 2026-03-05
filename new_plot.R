### tree plot
library(tidyverse)
library(treeio)
library(tidytree)
library(ggtree)
library(ggtreeExtra)
library(ggsci)
library(ggstar)
library(ape)
library(phytools)
library(geiger)

### configuration
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

### data locations
iqTreeFile<-"/lisc/data/scratch/molevo/jmontenegro/gaba_evo/results/iqtree/pLGIC.einsi.contree"
annoFile<-"/lisc/data/scratch/molevo/jmontenegro/gaba_evo/data/annot_table.tsv"

### load data
ml_tree<-read.tree(iqTreeFile)
annot<-read_tsv(annoFile)

# correct names
ml_tree$tip.label <- gsub("'", "", ml_tree$tip.label)

# now fix the family of nematostella GLUCLR family
annot$family <- ifelse(annot$species == "N. vectensis", "Unknown", annot$family)

# now fix the family of nematostella GLUCLR family
annot$family <- ifelse(grepl(":GluCl", annot$primaryName, ignore.case=TRUE), "GLUCLR", annot$family)

### now subset the tree to drop the cationic branch
tips_to_keep <- tips(ml_tree, node=481)
tips_to_drop <- ml_tree$tip.label[!ml_tree$tip.label %in% tips_to_keep]
subtree <- drop.tip(ml_tree, tips_to_drop)

### convert ml_tree into a treedata object
ml_treeTab <- as_tibble(subtree)
# add extra meta information
ml_treeTab <- ml_treeTab %>% left_join(annot, by="label")
ml_treeTab$family <- ifelse(is.na(ml_treeTab$family), "Unknown", ml_treeTab$family)
# convert the annotated tree back into plotable object treedata
ml_treeData <- as.treedata(ml_treeTab)

### automatic propagation of annotations into internal nodes
propagate_annotations <- function(treedata, annotation_column) {
  require(treeio)
  require(dplyr)
  require(tidytree)
  
  node_data <- as_tibble(treedata)
  n_tips <- Ntip(treedata@phylo)
  
  # Only modify internal nodes
  internal_nodes <- (n_tips + 1):(n_tips + Nnode(treedata@phylo))
  
  for (node_id in internal_nodes) {
    tips <- tidytree::offspring(treedata, node_id, type = "tips") %>% unlist()
    
    if (length(tips) > 0) {
      known_annotations <- node_data %>% 
        filter(node %in% tips) %>% 
        pull(!!sym(annotation_column)) %>% 
        .[!is.na(.)]  # Still exclude NAs from propagation
      
      if (length(known_annotations) > 0) {
        unique_annots <- unique(known_annotations)
        if (length(unique_annots) == 1) {
          node_data[node_data$node == node_id, annotation_column] <- unique_annots[1]
        } else {
          node_data[node_data$node == node_id, annotation_column] <- "Unknown"
        }
      } else {
        node_data[node_data$node == node_id, annotation_column] <- "Unknown"
      }
    }
  }
  
  treedata <- as.treedata(node_data)
  return(treedata)
}

### now time to root the tree
rootNodes <- ml_treeData@data %>% filter(grepl("5HT", entryName)) %>% select(node) %>% unlist() %>% unname()
ml_treeData <- root(ml_treeData, node=MRCA(ml_treeData, rootNodes))

### Now we can propagate the annotations where needed
ml_treeData <- propagate_annotations(ml_treeData, "family")
ml_treeData <- propagate_annotations(ml_treeData, "species")
ml_treeData <- propagate_annotations(ml_treeData, "Phylum")
ml_treeData <- propagate_annotations(ml_treeData, "ligand")
ml_treeData@data$ligand <- ifelse(is.na(ml_treeData@data$ligand), "Unknown", ml_treeData@data$ligand)
ml_treeData@data$species <- factor (ml_treeData@data$species, 
  levels=c("N. vectensis", "S. pistillata", "A. coerulea", "H. vulgaris", "C. elegans",
  "D. melanogaster", "P. dumerilii", "B. floridae", "H. sapiens", "Unknown"))
ml_treeData@data$ligand <- factor(ml_treeData@data$ligand,
  levels=c("Unknown", "Glutamate", "GABA", "Glycine", "ACh", 
           "Proton", "Serotonine", "Dopamine", "Tyramine", "Histamine"))

### reattached the correct labels to the phylo element:
ml_treeData@phylo$node.label <- subtree$node.label

# collor palette
colPal=c("#003c67", "#e73e26", "#009a67", "#0069b4", "#efc000", 
         "#9f569d", "#ea781d", "#f051a6", "#3a751d", "#bc9787", 
         "#2c797d", "#adaa4d", "#94c5c7")

# now plot
anionTreePlot <- ggtree(ml_treeData, size=1, aes(col=ligand)) %>% 
  scaleClade(node=332, scale = 0.3) %>% 
  scaleClade(node=342, scale = 0.3) %>% 
  scaleClade(node=356, scale = 0.3) %>% 
  scaleClade(node=201, scale = 0.3) %>% 
  scaleClade(node=259, scale = 0.3) %>% 
  scaleClade(node=275, scale = 0.3) %>% 
  scaleClade(node=286, scale = 0.3) %>% 
  scaleClade(node=347, scale = 0.3) %>% 
  scaleClade(node=345, scale = 0.3) %>% 
  scaleClade(node=324, scale = 0.3) %>% 
  collapse(node=332, mode="min", alpha=0.5, fill="grey70") %>%
  collapse(node=342, mode="min", alpha=0.5, fill="grey70") %>%
  collapse(node=356, mode="min", alpha=0.5, fill="grey70") %>%
  collapse(node=201, mode="min", alpha=0.5, fill="grey70") %>%
  collapse(node=259, mode="min", alpha=0.5, fill="grey70") %>%
  collapse(node=275, mode="min", alpha=0.5, fill="grey70") %>%
  collapse(node=286, mode="min", alpha=0.5, fill="grey70") %>%
  collapse(node=347, mode="min", alpha=0.5, fill="grey70") %>%
  collapse(node=345, mode="min", alpha=0.5, fill="grey70") %>%
  collapse(node=324, mode="min", alpha=0.5, fill="grey70")
  

anionTreePlot +
  geom_text2(aes(subset=!isTip, label=label), size=3, nudge_x=-0.5, nudge_y=0.5) + 
  geom_tiplab(align=TRUE, size=4, offset=1.5, show.legend=FALSE, 
#    aes(label=ifelse((family != "Unknown" | species %in% c("N. vectensis", "H. sapiens", "D. melanogaster", "C. elegans")), primaryName, NA))) + 
    aes(label=primaryName)) + 
  geom_fruit(geom="geom_star", offset=0.02, starstroke=0.2, size=4, aes(starshape=species, fill=ligand)) + 
  theme(legend.position=c(.05, .85), legend.text=element_text(size=12), 
    legend.title=element_text(size=14, face="bold")) +
  coord_fixed(clip="off", ratio=0.5) +
  labs(starshape="Species", col="Ligand") + 
  scale_color_manual(values=colPal) + 
  scale_fill_manual(values=colPal) + 
  guides(fill="none",
    starshape = guide_legend( 
      theme = theme(
        legend.text=element_text(face="italic")),
      override.aes = list(
        color=c(rep("#7f5d13", 4), rep("#00b1eb", 5)),
        starstroke=0.5)
    )) +
  geom_highlight(node=268, alpha=0.3, fill="#e73e26", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=294, alpha=0.3, fill="#003c67", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=242, alpha=0.3, fill="#003c67", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=280, alpha=0.3, fill="#003c67", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=258, alpha=0.3, fill="#003c67", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=275, alpha=0.3, fill="#003c67", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=231, alpha=0.3, fill="#003c67", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=236, alpha=0.3, fill="#003c67", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=303, alpha=0.3, fill="#003c67", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=316, alpha=0.3, fill="#003c67", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=345, alpha=0.3, fill="#003c67", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=347, alpha=0.3, fill="#003c67", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=354, alpha=0.3, fill="#003c67", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=196, alpha=0.3, fill="#009a67", type="gradient", gradient.direction="tr") + 
  geom_highlight(node=341, alpha=0.3, fill="#bc9787", type="gradient", gradient.direction="tr") +
  geom_highlight(node=328, alpha=0.3, fill="#0069b4", type="gradient", gradient.direction="tr") +
  geom_highlight(node=311, alpha=0.3, fill="#efc000", type="gradient", gradient.direction="tr") +
  geom_highlight(node=322, alpha=0.3, fill="#e73e26", type="gradient", gradient.direction="tr")


