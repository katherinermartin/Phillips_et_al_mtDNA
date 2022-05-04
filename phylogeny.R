## Phillips et al. 2022: Mapping life stage and ocean basin information onto the mtDNA long fragment phylogeny

library(treeio) # for read.mrbayes function
library(ggtree) # for phylogeny rendering
library(viridis) # for viridis palette
library(tidyverse) # for general tidying of data
library(ape) # for various phylo functions
library(ggnewscale) # for heatmap
library(picante) # for various phylo functions
library(emojifont) # for triangles on beast tree

setwd("/Users/KatieMartin/Documents/UCF/Research/mtDNA_longfragment/analysis/phylo_render/")

df <- read.csv("all_haps_places_stages_31Mar22.csv")

# make the haplotype column the name of the rows
df <- df %>% remove_rownames %>% column_to_rownames(var="haplotype")

# turn the matrix of zeroes and ones into "Atlantic" or N/A for the Basins

df <- df %>% mutate_at("Atlantic", str_replace, "1", "Atlantic")
df <- df %>% mutate_at("IndoPacific", str_replace, "1", "IndoPacific")
df <- df %>% mutate_at("Mediterranean", str_replace, "1", "Mediterranean")

# remove underscore from columns

df <- df %>%
  rename_with(~ gsub('_', ' ', .x))

# rename oceanic juveniles to dispersal stage juveniles
df <- rename(df, "dispersal-stage juveniles" = "oceanic juveniles")

# replace any time there is a record with a '1'

df <- df %>% mutate_at("dispersal-stage juveniles", str_replace, "1", "dispersal-stage juveniles")
df <- df %>% mutate_at("inwater juveniles", str_replace, "1", "inwater juveniles")
df <- df %>% mutate_at("rookery", str_replace, "1", "rookery")
df <- df %>% mutate_at("inwater adults", str_replace, "1", "inwater adults")
df <- df %>% mutate_at("inwater juveniles and adults", str_replace, "1", "inwater juveniles and adults")


# make all zeroes NAs in the dataframe
df[df == 0] <- NA


lifestages <- df %>% select(c("dispersal-stage juveniles", "inwater juveniles", "inwater juveniles and adults", "inwater adults", "rookery"))


basins <- df %>% select(c("Atlantic", "Mediterranean", "IndoPacific"))

tree <- read.mrbayes("/Users/KatieMartin/Documents/UCF/Research/mtDNA_longfragment/analysis/MrBayes/less_freq_sampling/25 Jan 22 corrected mtdna; mtDNA long frag bayes 27Jan22/CIPRES/infile.nex.con.tre")
tree # this has 716 tips.

# drop the outgroup tips:
tree <- treeio::drop.tip(tree, c("EF071948.1_trimmed", "EF122793.1_trimmed"))

# Initialize tree: tree aesthetics
p1 <- ggtree(tree, # tree read in
             ladderize = TRUE,
             right = TRUE, 
             size = 0.25) + 
  geom_treescale(fontsize=2,
                 linesize=1,
                 offset=4)


p2_rotated <- ggtree::rotate(p1, 715) %>% # rotate Dc to bottom
  ggtree::rotate(716) %>% # rotate Dc about its axis
  ggtree::rotate(749) %>% # rotate Cm+Nd down
  ggtree::rotate(1091) %>% # rotate Ei down
  ggtree::rotate(1092) %>% # rotate Ei about its axis
  ggtree::rotate(1249) %>% # rotate Cc down
  ggtree::rotate(1250) %>% # rotate Lk down
  ggtree::rotate(1251) %>% # rotate Lk about its axis
  ggtree::rotate(1260) # rotate Lo about its axis

p2_rotated

# tree with internal nodes colored by posterior probability
p2 <- p2_rotated + # previous tree
  geom_nodepoint(color = "black", # node border
                 fill = "#000000", # node fill
                 aes(subset = prob >= 0.995), # aesthetics: subset those with a probability greater than or equal to 0.995
                 size = 0.5, # size of nodes
                 shape = 21) # shape of nodes (circle)
p2

p3 <- gheatmap(p2, # previous tree to use
               basins, # variable to create a heatmap from
               width=0.17, # width of bars
               color = "NA", # color of bar borders
               colnames=FALSE) + # no column names
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3)) +
  scale_fill_manual(values = c("Atlantic" = "#dc180a",
                               "Mediterranean" = "#00BA38",
                               "IndoPacific" = "#619CFF"),
                    na.value = "white", name = "ocean basin") # have to put name of legend here
p3

# add another matrix; add ggnewscale since ggplot2 doesn't allow multiple fill scales

p4 <- p3 + new_scale_fill()
large_phylo_heatmap <- gheatmap(p4, # previous tree to use
                                lifestages, # variable to create a heatmap from
                                color = NA, # color of border around bins
                                offset=0.06, # how close 2 matrices are; width = width of bins
                                legend_title = "life stages", # legend title
                                width=0.3, # width of bars
                                font.size = 1, # font size for labels
                                colnames_position = "top", # If there were column labels, where they'd be
                                colnames = FALSE, colnames_offset_y = 1) + # column names
  scale_fill_manual(values = c("dispersal-stage juveniles" = "#fa8116",
                               "inwater juveniles" = "#e16462",
                               "inwater juveniles and adults" = "#b12a90",
                               "inwater adults" = "#5c01a6",
                               "rookery" = "#0b076f"),
                    na.value = "white", name = "life stage")
# dispersal-stage juveniles: fca636, fb9004
# rookery:0d0887
# inwater: e16462


large_phylo_heatmap
# plasma colors from viridis package; could try going a little darker monochromatically.
# switch colnames = FALSE to = 45 if I want them there.

#ggsave("figures/figure_phylo_mrbayes_19April22.svg")


### Zoom in on specific parts of tree for species-level phylogeny

tree <- read.mrbayes("/Users/KatieMartin/Documents/UCF/Research/mtDNA_longfragment/analysis/MrBayes/less_freq_sampling/25 Jan 22 corrected mtdna; mtDNA long frag bayes 27Jan22/CIPRES/infile.nex.con.tre")
tree # this has 716 tips.

# drop the outgroup tips:
tree <- treeio::drop.tip(tree, c("EF071948.1_trimmed", "EF122793.1_trimmed"))

# Initialize tree: tree aesthetics
p1 <- ggtree(tree, # tree read in
             ladderize = TRUE,
             right = TRUE, 
             size = 0.25) + 
  geom_treescale(fontsize=2,
                 linesize=1,
                 offset=4)
p1

# Call tree p5 to distinguish from previous phylogeny
p5 <- p1 + geom_tiplab(size = 1) + geom_treescale(fontsize=2, linesize=1, offset=1)
p5

p6 <- p5 + 
  geom_nodepoint(color = "black",
                 fill = "#000000",
                 aes(subset = prob >= 0.995),
                 size = 0.5,
                 shape = 21) + #posterior prob between 0.995- 1 (0.995 rounds up to 1)
  geom_nodepoint(color = "black",
                 fill = "#525252",
                 aes(subset = prob < 0.995 & prob >= 0.945),
                 size = 0.5,
                 shape = 21) + #posterior prob between 0.945-0.994 (0.945 rounds to 0.95) 
  geom_nodepoint(color = "black",
                 fill = "#cccccc",
                 aes(subset = prob < 0.945 & prob >= 0.895),
                 size = 0.5,
                 shape = 21) + # posterior prob between 0.895-0.944 (0.895 rounds up to 0.90)
  geom_nodepoint(color = "black",
                 fill = "#f7f7f7", aes(subset = prob < 0.895 & prob >= 0.795),
                 size = 0.5,
                 shape = 21) # posterior prob between 0.795-0.894 (0.795 rounds up to 0.80)

p7 <- gheatmap(p6,
               basins,
               width=0.3,
               colnames=FALSE,
               color = NA) +
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3)) +
  scale_fill_manual(values = c("Atlantic" = "#dc180a",
                               "Mediterranean" = "#00BA38",
                               "IndoPacific" = "#619CFF"),
                    na.value = "white", name = "ocean basin")

p7
# add another matrix; add ggnewscale since ggplot2 doesn't allow multiple fill scales

p8 <- p7 + new_scale_fill()
p9 <- gheatmap(p8,
               lifestages,
               color = NA,
               offset=0.1,
               width=0.5, # offset = how close 2 matrices are; width = width of bins
               font.size = 1, colnames_position = "top", colnames = FALSE, colnames_offset_y = 1) +
  scale_fill_manual(values = c("dispersal-stage juveniles" = "#fa8116",
                               "inwater juveniles" = "#e16462",
                               "inwater juveniles and adults" = "#b12a90",
                               "inwater adults" = "#5c01a6",
                               "rookery" = "#0b076f"),
                    na.value = "white", name = "life stage")
p9
# annotate a phylogenetic tree with insets: https://guangchuangyu.github.io/2016/01/annotate-a-phylogenetic-tree-with-insets/

phylo_Dc <- viewClade(p9, MRCA(p9, "Dc1.4", "Dc10.1"))
phylo_Dc
#ggsave("figures/phylo_Dc_19Apr22.svg")

phylo_Ei <- viewClade(p9, MRCA(p9, "EiIP_14", "Moz114"))
phylo_Ei
#ggsave("figures/phylo_Ei_19Apr22.svg")

phylo_Lk <- viewClade(p9, MRCA(p9, "Lk1.1", "Lk7.1"))
phylo_Lk
#ggsave("figures/phylo_Lk_19Apr22.svg")

phylo_Lo <- viewClade(p9, MRCA(p9, "Lo50", "Lo100"))
phylo_Lo
#ggsave("figures/phylo_Lo_19Apr22.svg")

phylo_Cc <- viewClade(p9, MRCA(p9, "CcA4.3", "CcA65.1"))
phylo_Cc
#ggsave("figures/phylo_Cc_19Apr22.svg")

phylo_Nd <- viewClade(p9, MRCA(p9, "Nd05", "Nd14"))
phylo_Nd
#ggsave("figures/phylo_Nd_19Apr22.svg")

phylo_Cm <- viewClade(p9, MRCA(p9, "CmP152.1", "CmA38.1"))
phylo_Cm
#ggsave("figures/phylo_Cm_19Apr22.svg")



#### BEAST2 tree

tree <- read.beast("/Users/KatieMartin/Documents/UCF/Research/mtDNA_longfragment/analysis/BEAST_beauti/treeannotator/13Mar22 treeannotator low mem/15Mar22_lowmem_treeannotator.tre")

tree #716 tips as expected; don't drop outgroups this time

# Initialize tree: tree aesthetics
p1 <- ggtree(tree, # tree read in
             ladderize = TRUE,
             right = TRUE, 
             size = 0.25)
p1

# rotate the clades to match previous phylo topology in literature
p1_rotated <- ggtree::rotate(p1, 717) %>% # rotate outgroups down
  ggtree::rotate(719) %>% # rotate Cm+Nd down
  ggtree::rotate(1067) %>% # rotate Ei
  ggtree::rotate(1068) %>% # rotate Cc down
  ggtree::rotate(1106) %>% # rotate Cc IA up
  ggtree::rotate(1172) %>% # rotate Lk down
  ggtree::rotate(1173) %>% # rotate Lo Indo up
  ggtree::rotate(1175) %>% #rotate Lo Atl/Pac down
  ggtree::rotate(721) # rotate Cm clades
p1_rotated

## collapse the nodes

p1_collapsed <- p1_rotated %>%
  ggtree::collapse(node=1431) %>% # collapse outgroups
  ggtree::collapse(node=1398) %>% # collapse Dc
  ggtree::collapse(node=1036) %>% # collapse Nd
  ggtree::collapse(node=1241) %>% # collapse Ei
  ggtree::collapse(node=1232) %>% # collapse Lk
  ggtree::collapse(node=722) %>% # collapse Cm Group 1/clade I and II
  ggtree::collapse(node=816) %>% # collapse Cm group 2/III-XI
  ggtree::collapse(node=1070) %>% # collapse Cc Atlantic Group 1/clade IB
  ggtree::collapse(node=1106) %>% # collapse Cc Atlantic Group 2/clade II
  ggtree::collapse(node=1160) %>% # collapse Cc Pacific Group/clade IA
  ggtree::collapse(node=1228) %>% # collapse Lo indo group/Indian ocean clade
  ggtree::collapse(node=1213) %>% # collapse Lo Pacific/Australia and E. pacific clade
  ggtree::collapse(node=1175) # collapse all other Lo/Atlantic and Pacific
p1_collapsed


#### Add triangles manually

p1_collapsed_triangle <- p1_collapsed +
  geom_text2(aes(subset=(node == 1431)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust=.42) +
  geom_text2(aes(subset=(node == 1398)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1036)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1241)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1232)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 722)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 816)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1070)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1106)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1160)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1228)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1213)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1175)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42)

p1_collapsed_triangle

# label the tree nodes
p1_collapsed_labeled <- p1_collapsed_triangle +
  geom_cladelabel(node=1431, "outgroups", fontsize = 3, hjust = -0.5) +
  geom_cladelabel(node=1398, "D. coriacea", fontsize = 3, hjust = -0.5) +
  geom_cladelabel(node=1036, "N. depressus", fontsize = 3, hjust = -0.5) +
  geom_cladelabel(node=1241, "E. imbricata", fontsize = 3, hjust = -0.5) +
  geom_cladelabel(node=1232, "L. kempii", fontsize = 3,hjust = -0.5) +
  geom_cladelabel(node=722, "C. mydas clades I & II", fontsize = 3, hjust = -0.5) +
  geom_cladelabel(node=816, "C. mydas clades III-XI", fontsize = 3, hjust = -0.5) +
  geom_cladelabel(node=1070, "C. caretta clade IB", fontsize = 3, hjust = -0.5) +
  geom_cladelabel(node=1106, "C. caretta clade II", fontsize = 3, hjust = -0.5) +
  geom_cladelabel(node=1160, "C. caretta clade IA", fontsize = 3, hjust = -0.5) +
  geom_cladelabel(node=1228, "L. olivacea Indian Ocean clade", fontsize = 3, hjust = -0.5) +
  geom_cladelabel(node=1213, "L. olivacea Australia and E. Pacific clade", fontsize = 3, hjust = -0.5) +
  geom_cladelabel(node=1175, "L. olivacea Atlantic and Pacific clade", fontsize = 3, hjust = -0.5)

# tree with internal nodes colored by posterior probability
p2 <- p1_collapsed_triangle + # previous tree
  geom_nodepoint(color = "black", # node border
                 fill = "#000000", # node fill
                 aes(subset = posterior >= 0.999), # aesthetics: subset those with a probability greater than or equal to 0.999
                 size = 0.5, # size of nodes
                 shape = 21) # shape of nodes (circle)
p2

# put some color on the ranges
p3 <- p2 + geom_range("CAheight_0.95_HPD", color = "#0000FF", size = 3, alpha = 0.4)
p3

# add the confidence intervals around height; I have tried so hard to round these numbers here but will have to do post-hoc in InkScape/Illustrator
p4 <- p3 + geom_text(aes(x=branch, label=CAheight_0.95_HPD), vjust = -1.5, hjust = 0.5, color='black', size = 2)
p4

p5 <- revts(p4) + theme_tree2() # theme_tree2 is needed for the revts
p5

p6 <- p5 + scale_x_continuous(breaks=c(-150, -125, -100, -75, -50, -25, -10, -5, 0), labels = abs)
p6

# another option: p5 <- revts(p4) + theme_tree2() + scale_x_continuous(labels = abs)

#ggsave("/Users/KatieMartin/Documents/UCF/Research/mtDNA_longfragment/analysis/phylo_render/figures/beast/beast_9Apr22.svg")

# Maximum likelihood tree

ml_mtdna <- read.tree("../../analysis/IQtree_29Apr22/seaturtle_mtDNA_longfragment_alm_final_25Jan22.fasta.treefile")

# drop the outgroup tips:
ml_mtdna <- treeio::drop.tip(ml_mtdna, c("EF071948.1_trimmed", "EF122793.1_trimmed"))

# Initialize tree: tree aesthetics
ml_mtdna <- ggtree(ml_mtdna,
                   ladderize = TRUE,
                   right = TRUE, 
                   size = 0.25) + 
  geom_treescale(fontsize=2,
                 linesize=1,
                 offset=4)

# rotate nodes to match previous presentations
ml_mtdna_rotated <- ggtree::rotate(ml_mtdna, 715) %>% # rotate Dc to bottom
  ggtree::rotate(716) %>% # rotate Dc about its axis
  ggtree::rotate(749) %>% # rotate Cm+Nd down
  ggtree::rotate(1097) %>% # rotate Ei down
  ggtree::rotate(1098) %>% # rotate Ei about its axis
  ggtree::rotate(1255) %>% # rotate Cc down
  ggtree::rotate(1256) %>% # rotate Lk down
  ggtree::rotate(1266) %>% # rotate Lo about its axis
  ggtree::rotate(1257) # rotate Lk about its axis
ml_mtdna_rotated <- ml_mtdna_rotated + geom_tiplab(size = 0.5)
ml_mtdna_rotated

ML_mtdna_plot <- ml_mtdna_rotated +
  geom_nodepoint(color = "black", fill = "#000000", aes(subset = label >= 100), size = 1, shape = 21) + # bootstrap greater than 100--it won't plot it for some reason if I specify only >=90 below
  geom_nodepoint(color = "black", fill = "#000000", aes(subset = label >= 90), size = 1, shape = 21) + # bootstrap greater than 90
  geom_nodepoint(color = "black", fill = "#525252", aes(subset = label < 89.5 & label >= 80), size = 1, shape = 21) + # bootstrap between 80 and 90
  geom_nodepoint(color = "black", fill = "#cccccc", aes(subset = label < 79.5 & label >= 70), size = 1, shape = 21) # bootstrap between 70 and 80
ML_mtdna_plot # bootstrap labels above 90 labeled with solid black dot


#ggsave(filename = "../analysis/phylo_render/figures/ML_mtdna_plot.svg", plot = ML_mtdna_plot)
