## Phillips et al. 2022: Mapping life stage and ocean basin information onto the mtDNA long fragment phylogeny

library(treeio) # for read.mrbayes function
library(ggtree) # for phylogeny rendering
library(viridis) # for viridis palette
library(tidyverse) # for general tidying of data
library(ape) # for various phylo functions
library(ggnewscale) # for heatmap
library(picante) # for various phylo functions
library(emojifont) # for triangles on beast tree

df <- read.csv("all_haps_places_stages_25May22_final.csv")

# make the haplotype column the name of the rows
df <- df %>% remove_rownames %>% column_to_rownames(var="haplotype")

# turn the matrix of zeroes and ones into "Atlantic" or N/A for the Basins

df <- df %>% mutate_at("Atlantic", str_replace, "1", "Atlantic")
df <- df %>% mutate_at("IndoPacific", str_replace, "1", "IndoPacific")
df <- df %>% mutate_at("Mediterranean", str_replace, "1", "Mediterranean")

# rename column
df <- rename(df, `Post-dispersal_juveniles` = `Post.dispersal_juveniles`)

# remove underscore from columns

df <- df %>%
  rename_with(~ gsub('_', ' ', .x))


# replace any time there is a record with a '1'

df <- df %>% mutate_at("Dispersing juveniles", str_replace, "1", "Dispersing juveniles")
df <- df %>% mutate_at("Post-dispersal juveniles", str_replace, "1", "Post-dispersal juveniles")
df <- df %>% mutate_at("Rookery", str_replace, "1", "Rookery")
df <- df %>% mutate_at("Inwater adults", str_replace, "1", "Inwater adults")
df <- df %>% mutate_at("Mixed juveniles and inwater adults", str_replace, "1", "Mixed juveniles and inwater adults")


# make all zeroes NAs in the dataframe
df[df == 0] <- NA


lifestages <- df %>% select(c("Dispersing juveniles", "Post-dispersal juveniles", "Mixed juveniles and inwater adults", "Inwater adults", "Rookery"))


basins <- df %>% select(c("Atlantic", "Mediterranean", "IndoPacific"))

tree <- read.mrbayes("infile.nex.con.tre")
tree # this has 709 tips.

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

#p1 <- p1 + geom_text(aes(label=node), hjust=-.3) # add only if node labels desired


p2_rotated <- ggtree::rotate(p1, 708) %>% # rotate Dc to bottom
  ggtree::rotate(710) %>% # rotate Dc about its axis
  ggtree::rotate(735) %>% # rotate Cm+Nd down
  ggtree::rotate(1075) %>% # rotate Ei + Cc + Lo + lk
  ggtree::rotate(1076) %>% # rotate Ei about its axis
  ggtree::rotate(1233) %>% # rotate Cc + Lk + Lo
  ggtree::rotate(1234) %>% # rotate Lo + Lk
  ggtree::rotate(1235) %>% # rotate Lk about its axis
  ggtree::rotate(1244) # rotate Lo about its axis

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
  scale_fill_manual(values = c("Dispersing juveniles" = "#fa8116",
                               "Post-dispersal juveniles" = "#e16462",
                               "Mixed juveniles and inwater adults" = "#b12a90",
                               "Inwater adults" = "#5c01a6",
                               "Rookery" = "#0b076f"),
                    na.value = "white", name = "life stage")

large_phylo_heatmap


ggsave("figures/MrBayes/figure_phylo_mrbayes_01Sept22.svg")


### Zoom in on specific parts of tree for species-level phylogeny

tree <- read.mrbayes("infile.nex.con.tre")
tree # this has 709 tips.

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
  scale_fill_manual(values = c("Dispersing juveniles" = "#fa8116",
                               "Post-dispersal juveniles" = "#e16462",
                               "Mixed juveniles and inwater adults" = "#b12a90",
                               "Inwater adults" = "#5c01a6",
                               "Rookery" = "#0b076f"),
                    na.value = "white", name = "life stage")
p9
# annotate a phylogenetic tree with insets: https://guangchuangyu.github.io/2016/01/annotate-a-phylogenetic-tree-with-insets/

phylo_Dc <- viewClade(p9, MRCA(p9, "Dc12.1", "Dc1.2"))
phylo_Dc
ggsave("figures/phylo_Dc_25May22.svg")

phylo_Ei <- viewClade(p9, MRCA(p9, "EiIP_14", "EiIP_49"))
phylo_Ei
ggsave("figures/phylo_Ei_25May22.svg")

phylo_Lk <- viewClade(p9, MRCA(p9, "Lk1.1", "Lk3.1"))
phylo_Lk
ggsave("figures/phylo_Lk_25May22.svg")

phylo_Lo <- viewClade(p9, MRCA(p9, "Lo50", "Lo20"))
phylo_Lo
ggsave("figures/phylo_Lo_25May22.svg")

phylo_Cc <- viewClade(p9, MRCA(p9, "CcA4.3", "CcA66.1"))
phylo_Cc
ggsave("figures/phylo_Cc_25May22.svg")

phylo_Nd <- viewClade(p9, MRCA(p9, "Nd01", "Nd32"))
phylo_Nd
ggsave("figures/phylo_Nd_25May22.svg")

phylo_Cm <- viewClade(p9, MRCA(p9, "CmA1.1", "JF926556.1"))
phylo_Cm
ggsave("figures/phylo_Cm_25May22.svg")

#### BEAST2 tree

tree <- read.beast("run2_treeannotator_Yule_uniform.txt")

p <- ggtree(tree) 

p1_rotated <- ggtree::rotate(p, 712) %>% # rotate Cm down next to Dc
  ggtree::rotate(1060) %>% # rotate Ei above Cm
  ggtree::rotate(1061) %>% # rotate Cc above Ei
  ggtree::rotate(1165) %>% # rotate Lo up
  ggtree::rotate(1166) # rotate Lo Indian Ocean clade to top

p1_collapsed <- p1_rotated %>%
  ggtree::collapse(node=1234) %>% # collapse Ei
  ggtree::collapse(node=1225) %>% # collapse Lk
  ggtree::collapse(node=715) %>% # collapse Cm Group 1/clade I and II
  ggtree::collapse(node=809) %>% # collapse Cm group 2/III-XI
  ggtree::collapse(node=1063) %>% # collapse Cc Atlantic Group 1/clade IB
  ggtree::collapse(node=1099) %>% # collapse Cc Atlantic Group 2/clade II
  ggtree::collapse(node=1153) %>% # collapse Cc Pacific Group/clade IA
  ggtree::collapse(node=1221) %>% # collapse Lo indo group/Indian ocean clade
  ggtree::collapse(node=1206) %>% # collapse Lo Pacific/Australia and E. pacific clade
  ggtree::collapse(node=1168) %>% # collapse Lo Atlantic/Pacific clade
  ggtree::collapse(node=1391) %>% # collapse Dc
  ggtree::collapse(node=1029) %>% # collapse Nd
  ggtree::collapse(node=1417) # collapse outgroup


p1_collapsed_triangle <- p1_collapsed +
  geom_text2(aes(subset=(node == 1234)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust=.42) +
  geom_text2(aes(subset=(node == 1225)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 715)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 809)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1063)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1099)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1153)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1221)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1206)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1168)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1391)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1029)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42) +
  geom_text2(aes(subset=(node == 1417)), cex=15, alpha=0.2,label=intToUtf8(9668), hjust =.2,vjust = 0.42)
p1_collapsed_triangle

# tree with internal nodes colored by posterior probability
p2 <- p1_collapsed_triangle + # previous tree with triangles
  geom_nodepoint(color = "black", # node border
                 fill = "#000000", # node fill
                 aes(subset = posterior >= 0.999), # aesthetics: subset those with a probability greater than or equal to 0.999
                 size = 1.5, # size of nodes
                 shape = 21) # shape of nodes (circle)
p2

# put some color on the ranges
p3 <- p2 + geom_range("CAheight_0.95_HPD", color = "#0000FF", size = 3, alpha = 0.2)
p3

#Time-Scale bar
p5 <- revts(p3) + theme_tree2() # theme_tree2 is needed for the revts
p5

ggsave(filename = "new_Beast_plot_Run2_wide.svg", plot = p5, width = 12, height = 8)

# Maximum likelihood tree

ml_mtdna <- read.tree("seaturtle_mtDNA_longfragment_alm.fasta.treefile")

# drop the outgroup tips:
ml_mtdna <- treeio::drop.tip(ml_mtdna, c("EF071948.1_trimmed", "EF122793.1_trimmed"))

# Initialize tree: tree aesthetics
ml_p1 <- ggtree(ml_mtdna,
                   ladderize = TRUE,
                   right = TRUE, 
                   size = 0.25) + 
  geom_treescale(fontsize=2,
                 linesize=1,
                 offset=4)
#ml_p1 <- ml_p1 + geom_text(aes(label=node), hjust=-.3, size = 2) # add only if node labels desired
ml_p1

# rotate nodes to match previous presentations
ml_p1_rotated <- ggtree::rotate(ml_p1, 708) %>% # rotate Dc to bottom
  ggtree::rotate(735) %>% # rotate Chelonidae
  ggtree::rotate(1083) %>% # rotate Ei + Cc + Lo + Lk
  ggtree::rotate(1241) %>% # rotate Cc + Lo + Lk
  ggtree::rotate(1242) %>% # rotate Lo + Lk
  ggtree::rotate(1084) %>% # rotate Ei about its axis
  ggtree::rotate(1252) %>% # rotate Lo about its axis
  ggtree::rotate(1243) %>% # rotate Lk about its axis
  ggtree::rotate(709) # rotate Dc about axis

ml_p1_rotated

ml_mtdna_rotated <- ml_p1_rotated + geom_tiplab(size = 0.5) # add tip labels
ml_mtdna_rotated

ML_mtdna_plot <- ml_mtdna_rotated +
  geom_nodepoint(color = "black", fill = "#000000", aes(subset = label >= 100), size = 1, shape = 21) + # bootstrap greater than 100--it won't plot it for some reason if I specify only >=90 below
  geom_nodepoint(color = "black", fill = "#000000", aes(subset = label >= 90), size = 1, shape = 21) + # bootstrap greater than 90
  geom_nodepoint(color = "black", fill = "#525252", aes(subset = label < 89.5 & label >= 80), size = 1, shape = 21) + # bootstrap between 80 and 90
  geom_nodepoint(color = "black", fill = "#cccccc", aes(subset = label < 79.5 & label >= 70), size = 1, shape = 21) # bootstrap between 70 and 80
ML_mtdna_plot # bootstrap labels above 90 labeled with solid black dot

ggsave(filename = "figure_phylo_IQtree_21May22.svg", plot = ML_mtdna_plot)
