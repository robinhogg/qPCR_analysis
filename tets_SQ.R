install.packages("gridExtra")
library("gridExtra")

install.packages("cowplot")
library("cowplot")

install.packages('dplyr', type ="source")
library('dplyr')

install.packages('tidyr', type ="source")
library('tidyr')

install.packages('ggplot2', type ="source")
library('ggplot2')

setwd('/pasteur/homes/rhogg/Desktop/tpdroso/qPCR_analysis/')

a <-qPCR_analysis_SQ(filename = "bin4_Sin3A_EH.csv", row_hk =   "E", row_test   = "H", rep ="1")
b <-qPCR_analysis_SQ(filename = "bin3_sin3a_DF.csv", row_hk =   "D", row_test   = "F", rep= "2")
sin3a <- gene_graph(a,b, prim ="sin3a")
sin3a

c <-qPCR_analysis_SQ(filename = "bin1_bam_AB.csv", row_hk =   "A", row_test   = "B", rep ="1")
d <-qPCR_analysis_SQ(filename = "bin2_bam_CD.csv", row_hk =   "C", row_test   = "D", rep= "2")
bam <- gene_graph(c,d, prim ="bam")
bam

e <-qPCR_analysis_SQ(filename = "bin1_sf_AB.csv", row_hk =   "A", row_test   = "B", rep ="1")
f <-qPCR_analysis_SQ(filename = "bin2_sf_CD.csv", row_hk =   "C", row_test   = "D", rep= "2")
sf <- gene_graph(e,f, prim ="SF")
sf

g <-qPCR_analysis_SQ(filename = "bin3_dhd_DF.csv", row_hk =   "D", row_test   = "F", rep ="1")
h <-qPCR_analysis_SQ(filename = "bin4_Dhd_EH.csv", row_hk =   "E", row_test   = "H", rep= "2")
dhd <- gene_graph(g,h, prim ="dhd")
dhd

j <-qPCR_analysis_SQ_bin5(filename = "bin5_rhi_FB_SQ.csv", row_hk =   "F", row_test   = "B", rep ="1")
k <-qPCR_analysis_SQ_bin6_rhi(filename = "bin6_rhi_GD_SQ.csv", row_hk =   "G", row_test   = "D", rep= "2")
rhi <- gene_graph(j,k, prim ="rhi")
rhi

l <-qPCR_analysis_SQ_bin7_jhmd(filename = "bin7_jhdm2_BF_SQ.csv", row_hk =   "B", row_test   = "F", rep ="1")
m <-qPCR_analysis_SQ(filename = "bin8_jhdm2_CH.csv", row_hk =   "C", row_test   = "H", rep= "2")
jh <- gene_graph(l,m, prim ="jhdm")
jh

grid.arrange(sin3a, sf, bam, dhd, rhi,jh, ncol = 3, nrow =2)
