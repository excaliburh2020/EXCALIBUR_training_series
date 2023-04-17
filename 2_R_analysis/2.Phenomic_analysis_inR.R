###############################################################################
# Author: Francesco Vitali                                                    #
# Affiliation: Research Centre for Agriculture and Environment, Council       #
# for Agricultural Research and Economics (CREA-AA)                           #
# e-mail: francesco.vitali@crea.gov.it                                        #
# twitter: @svito_fi                                                          #
# github: FrancescoVit                                                        #
# Date: 15/04/2023                                                            #
#                                                                             #
# Licence: MIT licence                                                        #
# Copyright (c) 2023 Excalibur H2020 Project                                  #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a copy#
# of this software and associated documentation files (the "Software"), to    #
# deal in the Software without restriction, including without limitation the  #
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or #
# sell copies of the Software, and to permit persons to whom the Software is  #
# furnished to do so, subject to the following conditions:                    #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     # 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
###############################################################################
###############################################################################
# Title of the script: Phenomic data analysis with R                          #
# prepared as supporting material for the workshop "EXCALIBUR Training Series:#
# Addressing microbial metabolic profile by means of Phenotype Microarray     #
# technology (BIOLOG)                                                         #
#                                                                             #
# Purpose of the script: use the activity values (AV) obtained as result from #
# the DuctApe workflow to perform some data analysis and visualization        #
###############################################################################

# Loading needed libraries

library(tidyverse)
library(ggsci)
library(ggpubr)
library(reshape2)
library(ggside)
library(ggdist)
library(FactoMineR)
library(factoextra)
library(pathview)


# Import in AV data, from DuctApe. We will use the data obtained by the first
# script used in the workshop. The script assumes that the entire repository
# was downloaded from GitHub, for this reason the AV data will be in the 
# "1_DuctApe_workflow" folder which is one "step" ahead of the "2_R_analysis"
# folder, from which this script was launched. This is the meaning of 
# "../1_DuctApe_workflow/phenome_combined.tsv" the "../" indicate the folder
# before (at least under linux environment). If needed, edit this part with the
# correct file path for the file in your computer.

ductape_data <- read.table(file = "../1_DuctApe_workflow/phenome_combined.tsv", 
                           header = F,
                           sep = "\t")

# Set column names

colnames(ductape_data) <- c("plate_id",	
                            "well_id",	
                            "chemical",	
                            "category",	
                            "moa",	
                            "co_id",	
                            "replica",	
                            "CoInoculum",	
                            "StrainA",	
                            "StrainB")

# Compare the obtained AV in the different compound category for the CoInoculum 
# strain. The same plot could be obtained for the StrainA and StrainB by 
# editing the "CoInoculum" part. We will also add a general Kruskal wallis test 
# to know if AVs in the different compound categories are different for the 
# selected strain

ductape_data %>%
  ggboxplot(x = "category", y = "CoInoculum", add = "jitter") +
  ylab("Activity Value (AV)") + 
  xlab("Compound category") + 
  stat_compare_means(method = "kruskal.test", label.x = 2)

# The above plot is a asic visualization. To demonstrate potentiality of use of 
# R software, we could obtain a more complex visualization with a 
# Raincloud like plot for cpd category 
# credit https://www.cedricscherer.com/2021/06/06/visualizing-distributions-with-raincloud-plots-and-how-to-create-them-with-ggplot2/


ductape_data %>% 
  ggplot(aes(x = category, y = CoInoculum)) + 
  stat_halfeye(adjust = 0.5, 
               width = 0.5, 
               .width = 0, 
               justification = -0.4, 
               point_colour = NA, 
               fill = "grey70") + 
  geom_boxplot(width = .25,
               outlier.shape = NA, 
               fill = "white") +
  geom_point(aes(fill = category), 
             shape = 21, 
             size = 1.3,
             alpha = .5,
             position = position_jitter(seed = 1, width = .1)) +
  theme_pubclean() + 
  ylab("Activity Value (AV)") + 
  xlab("Compound category") + 
  scale_fill_jama() +
  stat_compare_means(method = "kruskal.test", label.x = 2) +
  theme(legend.position="top",
        legend.title = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(0.8, 'cm'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) 

# The same analysis, could be performed on specific compound categories, 
# comparing in this case the AVs obtained with the three different strains

# Here we test the Carbon sources
ductape_data %>%
  select(-replica) %>% 
  melt() %>% 
  filter(category == "carbon") %>% 
  ggboxplot(x = "variable", y = "value", add = "jitter") +
  ylab("Activity Value (AV) on carbon sources") + 
  xlab("Compound category") + 
  stat_compare_means(method = "kruskal.test", label.x = 2, label.y = 11)

# Here we test the Nitrogen sources
ductape_data %>%
  select(-replica) %>% 
  melt() %>% 
  filter(category == "nitrogen") %>% 
  ggboxplot(x = "variable", y = "value", add = "jitter") +
  ylab("Activity Value (AV) on nitrogen sources") + 
  xlab("Compound category") + 
  stat_compare_means(method = "kruskal.test", label.x = 2, label.y = 11)

# Here we obtain a more elaborated visualization for the Nitrogen sources
ductape_data %>% 
  select(-replica) %>% 
  melt() %>% 
  filter(category == "nitrogen") %>% 
  ggplot(aes(x = variable, y = value)) + 
  stat_halfeye(adjust = 0.5, 
               width = 0.5, 
               .width = 0, 
               justification = -0.4, 
               point_colour = NA, 
               fill = "grey70") + 
  geom_boxplot(width = .25,
               outlier.shape = NA, 
               fill = "white") +
  geom_point(aes(fill = variable), 
             shape = 21, 
             size = 1.3,
             alpha = .5, 
             position = position_jitter(seed = 1, width = .1)) +
  theme_pubclean() + 
  ylab("Activity Value (AV) on nitrogen sources") + 
  xlab("Strain") + 
  scale_fill_jama() +
  stat_compare_means(method = "kruskal.test", label.x = 2, label.y = 11) +
  theme(legend.position="top",
        legend.title = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(0.8, 'cm'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) 

# Above analysis suggests no difference in how the two strains and the co-inoc
# experiment utilize Carbon sources, but a significant difference regarding the
# Nitrogen sources.
# We may futher explore the difference on Nitrogen sources with multivariate 
# analysis. In this case, we use a Principal Component Analysis (PCA)

ductape_data %>% 
  filter(category == "nitrogen") %>% 
  select(c(3,8,9,10)) %>% 
  column_to_rownames(var = "chemical") %>% 
  t() %>% 
  PCA() -> PCA_nitrogen

# The PCA function in the FactomineR package already produces plot. This can 
# be avoided by adding graph = F in the PCA() command.
# Plotting is usually better done by using the "companion" package, which is
# factoextra. We will produce a biplot for the PCA ordination, which visualize
# in the same graphics both the samples (black points) and the variables (blue
# arrows). This is useful for interpretation of the results, as we can know
# which variable is more strongly contributing to the difference that we see in
# the samples. The arrows (which are usually named vectors) represent the 
# direction in which each variable determine the placement of samples. Each
# variable "pull" the ordination in some direction

# Basic biplot
fviz_pca_biplot(X = PCA_nitrogen)

# The argument "repel = T" allows to avoid overlapping labels of variables, but 
# some label are lost.
fviz_pca_biplot(X = PCA_nitrogen, 
                repel = T)


# The fviz_contrib() fuction can be used to inspect variable contribution to 
# each axis. 
fviz_contrib(X = PCA_nitrogen, 
             xtickslab.rt = 90,
             choice = "var", 
             axes = 1) # change with 2 for second axis

fviz_contrib(X = PCA_nitrogen, 
             xtickslab.rt = 90,
             choice = "var",
             top = 30, # set the number of variables to plot
             axes = 1) # change with 2 for second axis

# If genomic data are available for the strain, we may explore the genomic
# basis for the phenomic data observed by using Kegg map. This approach
# takes a specific pathway of interest; in this example the map00910 for the
# Nitrogen metabolism, as we have evidences of different use of nitrogen 
# sources between the two strains. Compound from the PM plate which are included
# in the kegg pathway will be colored based on observed AV. Orthologs genes 
# from the strain genome which are included in the kegg pathway will be colored 

# Import gene annotation data. This contains mapping of protein sequences to 
# the Kegg ortholog. Can be obtained using the KAAS service 
# (https://www.genome.jp/kegg/kaas/) while protein sequences can be obtained
# by any pipeline for assembled genome annotation. As the function for plotting
# also allow to record the copy numbers for each orthologs in the genome, we 
# include a column with all values = 1
StrainA_gene <- read.table(file = "./StrainA_KAAS.csv", 
                           header = T,
                           sep = ",")

# we need to build a "named vector" for input in the function below
gene_data_strainA <- StrainA_gene$Presence
names(gene_data_strainA) <- StrainA_gene$Ortholog


# Same as above for StrainB
StrainB_gene <- read.table(file = "./StrainB_KAAS.csv", 
                           header = T,
                           sep = ",")
gene_data_strainB <- StrainB_gene$Presence
names(gene_data_strainB) <- StrainB_gene$Ortholog

# From the DuctApe output we can obtain the AV of each compound, and name it 
# using the column "co_id" which is the compound code in Kegg
cpd_data_strainA <- ductape_data$StrainA
names(cpd_data_strainA) <- ductape_data$co_id

# Same as above for StrainB
cpd_data_strainB <- ductape_data$StrainB
names(cpd_data_strainB) <- ductape_data$co_id

# Here we obtain the annotated Kegg map.
pv.out.N <- pathview(gene.data = gene_data_strainA, 
                     cpd.data = cpd_data_strainA,
                     both.dirs = list(gene = FALSE, cpd = FALSE),
                     bins = list(gene = 1, cpd = 15),
                     discrete = list(gene = TRUE, cpd = FALSE),
                     limit = list(gene = 1, cpd = 9),
                     species = "ko",
                     cpd.idtype = "kegg",
                     gene.idtype = "KEGG", 
                     pathway.id = "00910", 
                     out.suffix = "strainA.N",
                     keys.align = "y", 
                     kegg.native = T, 
                     key.pos = "topright")

pv.out.N <- pathview(gene.data = gene_data_strainB, 
                     cpd.data = cpd_data_strainB,
                     both.dirs = list(gene = FALSE, cpd = FALSE),
                     bins = list(gene = 1, cpd = 15),
                     discrete = list(gene = TRUE, cpd = FALSE),
                     limit = list(gene = 1, cpd = 9),
                     species = "ko",
                     cpd.idtype = "kegg",
                     gene.idtype = "KEGG", 
                     pathway.id = "00910", 
                     out.suffix = "strainB.N",
                     keys.align = "y", 
                     kegg.native = T, 
                     key.pos = "topright")
