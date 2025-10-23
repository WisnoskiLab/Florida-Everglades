# ---- Core data handling & manipulation ----
library(dplyr)
library(tidyr)
library(tidyverse)
library(plyr)
library(magrittr)
library(parallel)
library(reshape2)

# ---- Microbiome data processing ----
library(phyloseq)
library(biomformat)
library(file2meco)
library(microeco)
library(MicrobiomeStat)
library(meconetcomp)
library(WGCNA)
library(ggClusterNet)
library(ape)
library(picante)
library(Biostrings)

# ---- Differential abundance & compositional analysis ----
library(metagenomeSeq)
library(ALDEx2)
library(ANCOMBC)

# ---- Visualization ----
library(ggplot2)
library(ggpubr)
library(ggtree)
library(tidygraph)
library(paletteer)
library(colorspace)
library(ComplexHeatmap)
library(circlize)
library(vegan)

# ---- Statistical modeling ----
library(lme4)
library(lmerTest)
library(multcomp)
library(emmeans)
library(multcompView)
library(dplyr)
library(usethis)

# ---- Package manager ----
library(BiocManager)


### Importing data ####

q2_biom = import_biom("Processed_data/fcelter-q2.biom")
metadata_fce = import_qiime_sample_data("Processed_data/metadata-q2fce.txt")
tree = read_tree("Processed_data/unrooed-tree.nwk")
rep_fasta = readDNAStringSet("Processed_data/dna-sequences.fasta", format = "fasta")
str(metadata_fce)
#merging the imported components

q2_fcelterbiom  = merge_phyloseq(q2_biom, metadata_fce, tree, rep_fasta)

colnames(tax_table(q2_fcelterbiom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

q2_fcelterbiom

# converting phyloseq object to microeco onject
meco_dataset <- phyloseq2meco(q2_fcelterbiom)
meco_dataset
class(meco_dataset)

meco_dataset$cal_abund()

# phylum based

abun = trans_abund$new(dataset = meco_dataset, taxrank = "Phylum", ntaxa = 15)
abundance_plot = abun$plot_bar(others_color = "grey70", facet = c("env_broad_scale", "depth"), xtext_keep = FALSE, legend_text_italic = FALSE)
abundance_plot
ggsave("abundance.pdf", plot = abundance_plot, width =20, height = 12, dpi = 1000)


abundance_broad <- trans_abund$new(dataset = meco_dataset, taxrank = "Phylum", ntaxa = 20, groupmean = "env_broad_scale", group_morestats = T)

plot_broad <- abundance_broad$plot_bar(others_color = "grey70", legend_text_italic = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

plot_broad

abundance_depth <- trans_abund$new(dataset = meco_dataset, taxrank = "Phylum", ntaxa = 20, groupmean = "depth")

plot_depth <- abundance_depth$plot_bar(others_color = "grey70", legend_text_italic = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

plot_depth

ggsave("plot_depth.pdf", plot = plot_depth, width =10, height = 7, dpi = 1000)

# show 40 taxa at Genus level
genus <- trans_abund$new(dataset = meco_dataset, taxrank = "Genus", ntaxa = 40)

genus_heat <- genus$plot_heatmap(facet = "depth", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))

genus_heat
l= genus_heat + theme(axis.text.y = element_text(face = 'italic'))
ggsave("genus_heatmap_depth.pdf", plot = l, width =12, height = 7, dpi = 1000)

genus_heat_broad <- genus$plot_heatmap(facet = "env_broad_scale", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
m = genus_heat_broad + theme(axis.text.y = element_text(face = 'italic'))
ggsave("genus_heat_broad.pdf", plot = m, width =12, height = 7, dpi = 1000)

# taxa change along the gradient, space, time

genus_line <- trans_abund$new(dataset = meco_dataset, taxrank = "Phylum", ntaxa = 5, groupmean = "depth")

line_depth = genus_line$plot_line()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

line_depth
ggsave("Output/PDFs/line_depth.pdf", plot = line_depth, width =8, height = 7, dpi = 1000)

genus_line_broad <- trans_abund$new(dataset = meco_dataset, taxrank = "Phylum", ntaxa = 5, groupmean = "env_broad_scale")
line_broad = genus_line_broad$plot_line()
ggsave("line_broad.pdf", plot = line_broad, width =8, height = 7, dpi = 1000)


##################################
###### Alpha diversity ###### 

meco_dataset$cal_abund()
meco_dataset$cal_alphadiv()
meco_dataset$cal_betadiv()

water = c("dodgerblue","#f54260")

alpha <- trans_alpha$new(dataset = meco_dataset, group = "depth") 

alpha$cal_diff(method = "t.test")

alpha_table_q2fcelter = alpha[["data_alpha"]]

write.csv(alpha_table_q2fcelter, "alpha_table_q2fcelter.csv")

dat_alpha = readxl::read_xlsx("alpha_results_q2fcelter.xlsx")

# ---- Ensure factors ----
# 1) Build clean, ordered depth labels
depth_levels <- c("surface","09mm","27mm","45mm","63mm","81mm")

dat_alpha <- dat_alpha %>%
  mutate(
    depth_chr = as.character(depth),
    depth_lab = case_when(
      grepl("surface", depth_chr, ignore.case = TRUE) ~ "surface",
      TRUE ~ sprintf("%02dmm", as.integer(parse_number(depth_chr)))
    ),
    depth = factor(depth_lab, levels = depth_levels)  # ordered!
  )

dat_alpha <- dat_alpha %>% mutate(depth = fct_rev(depth))
# Example for Shannon index
mod_shannon <- lmer(Shannon ~ env_broad_scale * depth + (1 | Sitename), data = dat_alpha)
mod_chao1 <- lmer(Chao1 ~ env_broad_scale * depth + (1 | Sitename), data = dat_alpha)
mod_shannon

anova(mod_shannon)
anova(mod_chao1)

# Post-hoc pairwise
# Shannon
emm_shannon <- emmeans(mod_shannon, pairwise ~ env_broad_scale * depth)
cld_shannon <- multcomp::cld(emm_shannon$emmeans, Letters = letters, adjust = "tukey")

# Chao1
emm_chao1 <- emmeans(mod_chao1, pairwise ~ env_broad_scale * depth)
cld_chao1 <- multcomp::cld(emm_chao1$emmeans, Letters = letters, adjust = "tukey")# Shannon

# --- Shannon ---
shannon_depth = ggplot(dat_alpha, aes(x = env_broad_scale, y = Shannon, fill = depth)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_text(
    data = cld_shannon,
    aes(x = env_broad_scale, y = emmean + 0.05, label = .group, group = depth),
    position = position_dodge(width = 0.8),
    color = "black",
    size = 4
  ) +
  theme_bw(base_size = 14) +
  labs(y = "Shannon Diversity", x = "Environment Type") +
  scale_fill_brewer(palette = "Set2")

shannon_depth

# A small offset so letters float just above the group mean in the plot scale
offset_shan <- 0.03 * diff(range(dat_alpha$Shannon, na.rm = TRUE))

# --- 4) Plot Shannon depth profiles ---
shannon_depth <- ggplot(dat_alpha, aes(x = depth, y = Shannon, fill = env_broad_scale)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8), outlier.shape = 21) +
  # model means/SE + letters as you had:
  geom_errorbar(data = cld_shannon, aes(x = depth, ymin = emmean - SE, ymax = emmean + SE, group = env_broad_scale), 
                position = position_dodge(width = 0.8), width = 0.2, inherit.aes = FALSE) +
  geom_point(data = cld_shannon, aes(x = depth, y = emmean, group = env_broad_scale, color = env_broad_scale),
             position = position_dodge(width = 0.8), size = 2.7, inherit.aes = FALSE) +
  geom_text(data = cld_shannon, aes(x = depth, y = emmean + 0.03*diff(range(dat_alpha$Shannon, na.rm=TRUE)), label = .group, group = env_broad_scale),
            position = position_dodge(width = 0.8), size = 4, color = "black", inherit.aes = FALSE) +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(x = "Depth", y = "Shannon diversity", fill = "Biome", color = "Biome")


shannon_depth

# --- Chao1 ---
chao1_depth = ggplot(dat_alpha, aes(x = env_broad_scale, y = Chao1, fill = depth)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_text(
    data = cld_chao1,
    aes(x = env_broad_scale, y = emmean + 10, label = .group, group = depth),
    position = position_dodge(width = 0.8),
    color = "black",
    size = 4
  ) +
  theme_bw(base_size = 14) +
  labs(y = "Chao1 Richness", x = "Environment Type") +
  scale_fill_brewer(palette = "Set2")

chao1_depth_v <- ggplot(dat_alpha, aes(x = depth, y = Chao1, fill = env_broad_scale)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8), outlier.shape = 21) +
  # model means/SE + letters as you had:
  geom_errorbar(data = cld_chao1, aes(x = depth, ymin = emmean - SE, ymax = emmean + SE, group = env_broad_scale), 
                position = position_dodge(width = 0.8), width = 0.2, inherit.aes = FALSE) +
  geom_point(data = cld_chao1, aes(x = depth, y = emmean, group = env_broad_scale, color = env_broad_scale),
             position = position_dodge(width = 0.8), size = 2.7, inherit.aes = FALSE) +
  geom_text(data = cld_chao1, aes(x = depth, y = emmean + 0.03*diff(range(dat_alpha$Chao1, na.rm=TRUE)), label = .group, group = env_broad_scale),
            position = position_dodge(width = 0.8), size = 4, color = "black", inherit.aes = FALSE) +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(x = "Depth", y = "Chao1 index", fill = "Biome", color = "Biome")

chao1_depth_v


ggsave("shannon depth.pdf", plot = shannon_depth, width =8, height =5, dpi = 1000)
ggsave("chao1 depth.pdf", plot = chao1_depth, width =8, height =5, dpi = 1000)


ggsave("shannon depth ver.pdf", plot = shannon_depth, width =8, height =5, dpi = 1000)
ggsave("chao1 depth ver.pdf", plot = chao1_depth_v, width =8, height =5, dpi = 1000)
##################################
###### beta diversity ###### 


beta_fce <- trans_beta$new(dataset = meco_dataset, group = "env_broad_scale", measure = "bray")

beta_fce$cal_ordination(method = "PCoA")

beta_fce_loc = beta_fce$plot_ordination( plot_shape = "depth", plot_color = "env_broad_scale", point_size = 2, plot_type = c("point", "ellipse")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

beta_fce_loc

ggsave("beta_fce_brod.pdf", plot = beta_fce_loc, width =9, height =6, dpi = 1000)

## comparing the group distances

broad_beta_grp <- trans_beta$new(dataset = meco_dataset, group = "env_broad_scale", measure = "bray")
depth_beta_grp <- trans_beta$new(dataset = meco_dataset, group = "depth", measure = "bray")

# calculate and plot sample distances within groups
broad_beta_grp$cal_group_distance(within_group = TRUE)
depth_beta_grp$cal_group_distance(within_group = TRUE)

# return t1$res_group_distance
# perform Wilcoxon Rank Sum and Signed Rank Tests
broad_beta_grp$cal_group_distance_diff(method = "wilcox")
depth_beta_grp$cal_group_distance_diff(method = "wilcox")
# plot_group_order parameter can be used to adjust orders in x axis
broad_bray_plot = broad_beta_grp$plot_group_distance(add = "mean") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

depth_bray_plot = depth_beta_grp$plot_group_distance(add = "mean")  + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

ggsave("broad_bray_plot.pdf", plot = broad_bray_plot, width = 7, height =6, dpi = 1000)
ggsave("depth_bray_plot.pdf", plot = depth_bray_plot, width = 7, height =6, dpi = 1000)

## PERMANOVA can be applied to the differential test of distances among groups via the cal_manova
#function developed based on the adonis2 function of vegan package

panova_broad <- trans_beta$new(dataset = meco_dataset, group = "env_broad_scale", measure = "bray")

# manova for all groups when manova_all = TRUE
panova_broad$cal_manova(manova_set = "env_broad_scale + depth")
panova_broad$res_manova

## Community assembly mechanisms

nullmodel = trans_nullmodel$new(meco_dataset, filter_thres = 0.0005)

# see null.model parameter for other null models
# null model run 500 times for the example

#beta NRI
nullmodel$cal_ses_betampd(runs = 500, abundance.weighted = TRUE)
nullmodel$res_ses_betampd

#beta NTI
nullmodel$cal_ses_betamntd(runs = 500, abundance.weighted = TRUE, null.model = "taxa.labels")
nullmodel$res_ses_betamntd

# add betaNRI matrix to beta_diversity list
meco_dataset$beta_diversity[["betaNRI"]] <- nullmodel$res_ses_betampd
meco_dataset$beta_diversity[["betaNTI"]] <- nullmodel$res_ses_betamntd

# create trans_beta class, use measure "betaNRI"
plot_betaNTRI <- trans_beta$new(dataset = meco_dataset, group = "env_broad_scale", measure = "betaNRI")
plot_betaNTRI_depth <- trans_beta$new(dataset = meco_dataset, group = "depth", measure = "betaNRI")

plot_betaNTI <- trans_beta$new(dataset = meco_dataset, group = "env_broad_scale", measure = "betaNTI")
plot_betaNTI_depth <- trans_beta$new(dataset = meco_dataset, group = "depth", measure = "betaNTI")

# transform the distance for each group
plot_betaNTRI$cal_group_distance()
plot_betaNTRI_depth$cal_group_distance()

plot_betaNTI $cal_group_distance()
plot_betaNTI_depth$cal_group_distance()


# see the help document for more methods, e.g. "anova" and "KW_dunn"
plot_betaNTRI$cal_group_distance_diff(method = "wilcox")
plot_betaNTRI_depth$cal_group_distance_diff(method = "wilcox")

plot_betaNTI$cal_group_distance_diff(method = "wilcox")
plot_betaNTI_depth$cal_group_distance_diff(method = "wilcox")

# plot the results
betaNTRI <- plot_betaNTRI$plot_group_distance(add = "mean") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)  # Add border

betaNTRI_depth <- plot_betaNTRI_depth$plot_group_distance(add = "mean") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)  # Add border

betaNTI_broad <- plot_betaNTI$plot_group_distance(add = "mean") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)  # Add border

betaNTI_depth <- plot_betaNTI_depth$plot_group_distance(add = "mean") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)  # Add border

ggsave("betaNTRI_broad.pdf", plot = betaNTRI, width = 8, height =7, dpi = 1000)
ggsave("betaNTRI_depth.pdf", plot = betaNTRI_depth, width = 8, height =7, dpi = 1000)

ggsave("betaNTI_broad.pdf", plot = betaNTI_broad, width = 8, height =7, dpi = 1000)
ggsave("betaNTI_depth.pdf", plot = betaNTI_depth, width = 8, height =7, dpi = 1000)

# RC bray (Bray-Curtis-based Raup-Crick)
nullmodel$cal_rcbray(runs = 1000)
nullmodel$res_rcbray

# use betaNTI and rcbray to evaluate processes
broad_assem = nullmodel$cal_process(use_betamntd = TRUE, group = "env_broad_scale") 
depth_assem = nullmodel$cal_process(use_betamntd = TRUE, group = "depth")

assem = nullmodel$cal_process(use_betamntd = T)

nullmodel$res_process

# require NST package to be installed
nst_broad = nullmodel$cal_NST(method = "tNST", group = "env_broad_scale", dist.method = "bray", abundance.weighted = TRUE, output.rand = TRUE, SES = TRUE)
nst_broad$res_NST$index.grp
# test the NST difference between each pair of groups
nullmodel$cal_NST_test(method = "nst.boot")

# convert long format table to square matrix
# the 10th column: MST.ij.bray in t1$res_NST$index.pair
test <- nullmodel$cal_NST_convert(10)

# for pNST method, phylogenetic tree is needed
nullmodel$cal_NST(method = "pNST", group = "env_broad_scale", output.rand = TRUE, SES = TRUE)
nullmodel$cal_NST_test(method = "nst.boot")

# For nearest Taxon Index (NTI) and nearest Relative Index (NRI), please use cal_NTI and cal_NRI, respectively.

nullmodel$cal_NRI(null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)
nullmodel$cal_NTI(null.model = "taxa.labels", abundance.weighted = TRUE, runs = 999)

nst_depth = nullmodel$cal_NST(method = "tNST", group = "depth", dist.method = "bray", abundance.weighted = TRUE, output.rand = TRUE, SES = TRUE)
nst_depth$res_NST$index.grp
