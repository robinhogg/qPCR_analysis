#qPCR Analaysis by Robin HOGG (robinhogg60600@gmail.com)

#Script for the analysis of qPCR data using delta delta Cq method
#Your data should be like this :
# row;Sample.Name;Cq
# A0;Name;(a number)

#Defining the function
qPCR_analysis <- function(filename, row_hk, row_test, rep){

###########################################################################
#######IMPORTING DATA FROM EXCEL TO R AND SETUP ###########################
###########################################################################

qpcr_data <- read.csv(file = filename, stringsAsFactors = FALSE, sep=';')
head(qpcr_data)

primer_key <- data.frame(row = c(row_hk,row_test), #Give the the row number of your housekeeping gene and test gene
                         primers = c(rep("hk_gene", 1), rep("test_gene", 1)))

tidy_data <- separate(qpcr_data, Well, into = c("row", "column"), #Create a column row
                      sep = 1, convert = TRUE)
head(tidy_data)

tidy_data <- left_join(tidy_data, primer_key, by = "row") #Use Left_join by dplyr to past both tabs using row as matrix

#We control the sheet disposition
ggplot(tidy_data, aes(x = column, y = row, fill = primers, label = Sample.Name)) +
  geom_tile(colour = "black") +
  geom_text() +
  scale_y_discrete(limits = c("row_hk", "row_test")) +
  scale_x_continuous(breaks = 1:12)

###########################################################################
##############################ANALYZE DATA#################################
###########################################################################

#Calculate average Cq values

head(tidy_data$Cq)
class(tidy_data$Cq)

summarised_data <- tidy_data %>%
  mutate(Cq = as.numeric(as.character(Cq))) %>%
  group_by(Sample.Name, primers) %>%
  summarise(mean_Cq = mean(Cq))

head(summarised_data)

#Look the homogeneity of our result and allow us to eliminate aberrant data
ggplot(summarised_data, aes(x = Sample.Name, y = mean_Cq, colour = primers)) +
  geom_point()


###########################################################################
#############################DeltaDelta Cq ################################
###########################################################################

test_data <- summarised_data %>%
  filter(primers == "test_gene")

ref_data <- summarised_data %>%
  filter(primers == "hk_gene") %>%
  rename("ref_Cq" = "mean_Cq")

combined_data <- left_join(test_data, ref_data, by = c("Sample.Name"))

combined_data <- mutate(combined_data, delta_Cq = ref_Cq - mean_Cq)
ggplot(combined_data, aes(x = Sample.Name, y = delta_Cq)) +
  geom_point()

treatment_summary <- combined_data %>%
  group_by(Sample.Name) %>%
  summarise(mean_delta_Cq = mean(delta_Cq))

mean_control <- filter(treatment_summary, Sample.Name == "WT1" |  Sample.Name == "WT2"|  Sample.Name == "WT3") %>% pull(mean_delta_Cq)
mean_control <- mean(mean_control)
mean_control

combined_data <- combined_data %>% 
  mutate(delta_delta_Cq = mean_control - delta_Cq)

ggplot(combined_data, aes(x = Sample.Name, y = delta_delta_Cq)) +
  geom_point()

###########################################################################
############################Relative DNA C ################################
###########################################################################

combined_data_c <- combined_data %>%
  mutate(rel_conc = 2^-delta_delta_Cq)

ggplot(combined_data_c, aes(x = Sample.Name, y = rel_conc)) +
  geom_point() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.3))

combined_data_mc <- combined_data_c %>%
  group_by(Sample.Name) %>%
  mutate(mean_rel_conc = mean(rel_conc))

mean_rc_wt <- filter(combined_data_c, Sample.Name == "WT1" |  Sample.Name == "WT2"|  Sample.Name == "WT3") %>% pull(rel_conc)
mean_rc_wt2 <- mean(mean_rc_wt)
mean_rc_wt2

mean_rc_sf <- filter(combined_data_c, Sample.Name == "SF1" |  Sample.Name == "SF2"|  Sample.Name == "SF3") %>% pull(rel_conc)
mean_rc_sf2 <- mean(mean_rc_sf)
mean_rc_sf2

primer <- c("rpl32")
mean_rc <- c(mean_rc_wt)
wt_summary <- data.frame(x = primer, y = mean_rc)
colnames(wt_summary) <- c("primers", "mean_rc")

primer <- c("test")
mean_rc <- c(mean_rc_sf)
sf_summary <- data.frame(x = primer, y = mean_rc)
colnames(sf_summary) <- c("primers", "mean_rc")

sd_rc_wt <- sd(wt_summary$mean_rc)
sd_rc_wt

sd_rc_sf <- sd(sf_summary$mean_rc)
sd_rc_sf

sd_sum <- c(sd_rc_wt, sd_rc_sf)

primer <- c("ctrl", "SF")
mean_rc <- c(mean_rc_wt2, mean_rc_sf2)

replicatt <- c(rep,rep)

final_summary <- data.frame(x = primer, y = mean_rc, w = sd_sum, z= replicatt)
colnames(final_summary) <- c("primers", "mean_rc", "sd", "replicat")
final_summary
}

qPCR_analysis_SQ_bin5 <- function(filename, row_hk, row_test, rep){
  
  ###########################################################################
  #######IMPORTING DATA FROM EXCEL TO R AND SETUP ###########################
  ###########################################################################
  
  qpcr_data <- read.csv(file = filename, stringsAsFactors = FALSE, sep=';')
  head(qpcr_data)
  
  primer_key <- data.frame(row = c(row_hk,row_test), #Give the the row number of your housekeeping gene and test gene
                           primers = c(rep("hk_gene", 1), rep("test_gene", 1)))
  
  tidy_data <- separate(qpcr_data, Well, into = c("row", "column"), #Create a column row
                        sep = 1, convert = TRUE)
  head(tidy_data)
  
  tidy_data <- left_join(tidy_data, primer_key, by = "row") #Use Left_join by dplyr to past both tabs using row as matrix
  
  #We control the sheet disposition
  ggplot(tidy_data, aes(x = column, y = row, fill = primers, label = Sample.Name)) +
    geom_tile(colour = "black") +
    geom_text() +
    scale_y_discrete(limits = c("row_hk", "row_test")) +
    scale_x_continuous(breaks = 1:12)
  
  ###########################################################################
  ##############################ANALYZE DATA#################################
  ###########################################################################
  
  #Calculate average SQ values
  
  head(tidy_data$SQ)
  class(tidy_data$SQ)
  
  summarised_data <- tidy_data %>%
    mutate(Cq = as.numeric(as.character(SQ))) %>%
    group_by(Sample.Name, primers) %>%
    summarise(mean_SQ = mean(SQ))
  
  head(summarised_data)
  
  #Look the homogeneity of our result and allow us to eliminate aberrant data
  ggplot(summarised_data, aes(x = Sample.Name, y = mean_SQ, colour = primers)) +
    geom_point()
  
  #for i in [1,lenght(summarized_data[3])], j in [2,lenght(summarized_data)+1];
  #ratio_wt1 <- summarised_data[i,3] /  summarised_data[j,3]
  #j=j+2
  #i=i+2
  
  ratio_sf1 <- summarised_data[2,3] /  summarised_data[1,3]
  ratio_sf2 <- summarised_data[4,3] /  summarised_data[3,3]
  ratio_sf3 <- summarised_data[6,3] /  summarised_data[5,3]
  ratio_wt1 <- summarised_data[8,3] /  summarised_data[7,3]
  ratio_wt2 <- summarised_data[10,3] /  summarised_data[9,3]
  ratio_wt3 <- summarised_data[12,3] /  summarised_data[11,3]
  
  WT1 <- "WT1"
  WT2 <- "WT2"
  WT3 <- "WT3"
  SF1 <- "SF1"
  SF2 <- "SF2"
  SF3 <- "SF3"
  
  df_wt1 <- data.frame(x=WT1, y=ratio_wt1)
  colnames(df_wt1) <- c("primers", "ratio")
  
  df_wt2 <- data.frame(x=WT2, y=ratio_wt2)
  colnames(df_wt2) <- c("primers", "ratio")
  
  df_wt3 <- data.frame(x=WT3, y=ratio_wt3)
  colnames(df_wt3) <- c("primers", "ratio")
  
  df_sf1 <- data.frame(x=SF1, y=ratio_sf1)
  colnames(df_sf1) <- c("primers", "ratio")
  
  df_sf2 <- data.frame(x=SF2, y=ratio_sf2)
  colnames(df_sf2) <- c("primers", "ratio")
  
  df_sf3 <- data.frame(x=SF3, y=ratio_sf3)
  colnames(df_sf3) <- c("primers", "ratio")
  
  inter_1 <- union(df_wt1, df_wt2)
  inter_2 <- union(inter_1, df_wt3)
  inter_3 <- union(inter_2, df_sf1)
  inter_4 <- union(inter_3, df_sf2)
  treatment_summary <- union(inter_4, df_sf3)
  
  mean_sq_wt <- filter(treatment_summary, primers == "WT1") %>% pull(ratio)
  mean_SQ_wt_ <- mean(mean_sq_wt)
  mean_SQ_wt_
  
  mean_sq_sf <- filter(treatment_summary,  primers == "SF2"|  primers == "SF3") %>% pull(ratio)
  mean_SQ_sf_ <- mean(mean_sq_sf)
  mean_SQ_sf_
  
  primer <- c("rpl32")
  mean_sq <- c(mean_sq_wt)
  wt_summary <- data.frame(x = primer, y = mean_sq)
  colnames(wt_summary) <- c("primers", "mean_sq")
  
  primer <- c("test")
  mean_sq <- c(mean_sq_sf)
  sf_summary <- data.frame(x = primer, y = mean_sq)
  colnames(sf_summary) <- c("primers", "mean_sq")
  
  sd_rc_wt <- sd(wt_summary$mean_sq)
  sd_rc_wt
  
  sd_rc_sf <- sd(sf_summary$mean_sq)
  sd_rc_sf
  
  sd_sum <- c(sd_rc_wt, sd_rc_sf)
  
  primer <- c("ctrl", "SF")
  mean_sq <- c(mean_SQ_wt_, mean_SQ_sf_)
  
  replicatt <- c(rep,rep)
  
  final_summary <- data.frame(x = primer, y = mean_sq, w = sd_sum, z= replicatt)
  colnames(final_summary) <- c("primers", "mean_rc", "sd", "replicat")
  final_summary
  
}


qPCR_analysis_SQ_bin6_rhi <- function(filename, row_hk, row_test, rep){
  
  ###########################################################################
  #######IMPORTING DATA FROM EXCEL TO R AND SETUP ###########################
  ###########################################################################
  
  qpcr_data <- read.csv(file = filename, stringsAsFactors = FALSE, sep=';')
  head(qpcr_data)
  
  primer_key <- data.frame(row = c(row_hk,row_test), #Give the the row number of your housekeeping gene and test gene
                           primers = c(rep("hk_gene", 1), rep("test_gene", 1)))
  
  tidy_data <- separate(qpcr_data, Well, into = c("row", "column"), #Create a column row
                        sep = 1, convert = TRUE)
  head(tidy_data)
  
  tidy_data <- left_join(tidy_data, primer_key, by = "row") #Use Left_join by dplyr to past both tabs using row as matrix
  
  #We control the sheet disposition
  ggplot(tidy_data, aes(x = column, y = row, fill = primers, label = Sample.Name)) +
    geom_tile(colour = "black") +
    geom_text() +
    scale_y_discrete(limits = c("row_hk", "row_test")) +
    scale_x_continuous(breaks = 1:12)
  
  ###########################################################################
  ##############################ANALYZE DATA#################################
  ###########################################################################
  
  #Calculate average SQ values
  
  head(tidy_data$SQ)
  class(tidy_data$SQ)
  
  summarised_data <- tidy_data %>%
    mutate(Cq = as.numeric(as.character(SQ))) %>%
    group_by(Sample.Name, primers) %>%
    summarise(mean_SQ = mean(SQ))
  
  head(summarised_data)
  
  #Look the homogeneity of our result and allow us to eliminate aberrant data
  ggplot(summarised_data, aes(x = Sample.Name, y = mean_SQ, colour = primers)) +
    geom_point()
  
  #for i in [1,lenght(summarized_data[3])], j in [2,lenght(summarized_data)+1];
  #ratio_wt1 <- summarised_data[i,3] /  summarised_data[j,3]
  #j=j+2
  #i=i+2
  
  ratio_sf1 <- summarised_data[2,3] /  summarised_data[1,3]
  ratio_sf2 <- summarised_data[4,3] /  summarised_data[3,3]
  ratio_sf3 <- summarised_data[6,3] /  summarised_data[5,3]
  ratio_wt1 <- summarised_data[8,3] /  summarised_data[7,3]
  ratio_wt2 <- summarised_data[10,3] /  summarised_data[9,3]
  ratio_wt3 <- summarised_data[12,3] /  summarised_data[11,3]
  
  WT1 <- "WT1"
  WT2 <- "WT2"
  WT3 <- "WT3"
  SF1 <- "SF1"
  SF2 <- "SF2"
  SF3 <- "SF3"
  
  df_wt1 <- data.frame(x=WT1, y=ratio_wt1)
  colnames(df_wt1) <- c("primers", "ratio")
  
  df_wt2 <- data.frame(x=WT2, y=ratio_wt2)
  colnames(df_wt2) <- c("primers", "ratio")
  
  df_wt3 <- data.frame(x=WT3, y=ratio_wt3)
  colnames(df_wt3) <- c("primers", "ratio")
  
  df_sf1 <- data.frame(x=SF1, y=ratio_sf1)
  colnames(df_sf1) <- c("primers", "ratio")
  
  df_sf2 <- data.frame(x=SF2, y=ratio_sf2)
  colnames(df_sf2) <- c("primers", "ratio")
  
  df_sf3 <- data.frame(x=SF3, y=ratio_sf3)
  colnames(df_sf3) <- c("primers", "ratio")
  
  inter_1 <- union(df_wt1, df_wt2)
  inter_2 <- union(inter_1, df_wt3)
  inter_3 <- union(inter_2, df_sf1)
  inter_4 <- union(inter_3, df_sf2)
  treatment_summary <- union(inter_4, df_sf3)
  
  mean_sq_wt <- filter(treatment_summary, primers == "WT2" |primers == "WT3") %>% pull(ratio)
  mean_SQ_wt_ <- mean(mean_sq_wt)
  mean_SQ_wt_
  
  mean_sq_sf <- filter(treatment_summary, primers == "SF1"| primers == "SF2"|  primers == "SF3") %>% pull(ratio)
  mean_SQ_sf_ <- mean(mean_sq_sf)
  mean_SQ_sf_
  
  primer <- c("rpl32")
  mean_sq <- c(mean_sq_wt)
  wt_summary <- data.frame(x = primer, y = mean_sq)
  colnames(wt_summary) <- c("primers", "mean_sq")
  
  primer <- c("test")
  mean_sq <- c(mean_sq_sf)
  sf_summary <- data.frame(x = primer, y = mean_sq)
  colnames(sf_summary) <- c("primers", "mean_sq")
  
  sd_rc_wt <- sd(wt_summary$mean_sq)
  sd_rc_wt
  
  sd_rc_sf <- sd(sf_summary$mean_sq)
  sd_rc_sf
  
  sd_sum <- c(sd_rc_wt, sd_rc_sf)
  
  primer <- c("ctrl", "SF")
  mean_sq <- c(mean_SQ_wt_, mean_SQ_sf_)
  
  replicatt <- c(rep,rep)
  
  final_summary <- data.frame(x = primer, y = mean_sq, w = sd_sum, z= replicatt)
  colnames(final_summary) <- c("primers", "mean_rc", "sd", "replicat")
  final_summary
  
}

qPCR_analysis_SQ_bin7_jhmd <- function(filename, row_hk, row_test, rep){
  
  ###########################################################################
  #######IMPORTING DATA FROM EXCEL TO R AND SETUP ###########################
  ###########################################################################
  
  qpcr_data <- read.csv(file = filename, stringsAsFactors = FALSE, sep=';')
  head(qpcr_data)
  
  primer_key <- data.frame(row = c(row_hk,row_test), #Give the the row number of your housekeeping gene and test gene
                           primers = c(rep("hk_gene", 1), rep("test_gene", 1)))
  
  tidy_data <- separate(qpcr_data, Well, into = c("row", "column"), #Create a column row
                        sep = 1, convert = TRUE)
  head(tidy_data)
  
  tidy_data <- left_join(tidy_data, primer_key, by = "row") #Use Left_join by dplyr to past both tabs using row as matrix
  
  #We control the sheet disposition
  ggplot(tidy_data, aes(x = column, y = row, fill = primers, label = Sample.Name)) +
    geom_tile(colour = "black") +
    geom_text() +
    scale_y_discrete(limits = c("row_hk", "row_test")) +
    scale_x_continuous(breaks = 1:12)
  
  ###########################################################################
  ##############################ANALYZE DATA#################################
  ###########################################################################
  
  #Calculate average SQ values
  
  head(tidy_data$SQ)
  class(tidy_data$SQ)
  
  summarised_data <- tidy_data %>%
    mutate(Cq = as.numeric(as.character(SQ))) %>%
    group_by(Sample.Name, primers) %>%
    summarise(mean_SQ = mean(SQ))
  
  head(summarised_data)
  
  #Look the homogeneity of our result and allow us to eliminate aberrant data
  ggplot(summarised_data, aes(x = Sample.Name, y = mean_SQ, colour = primers)) +
    geom_point()
  
  #for i in [1,lenght(summarized_data[3])], j in [2,lenght(summarized_data)+1];
  #ratio_wt1 <- summarised_data[i,3] /  summarised_data[j,3]
  #j=j+2
  #i=i+2
  
  ratio_sf1 <- summarised_data[2,3] /  summarised_data[1,3]
  ratio_sf2 <- summarised_data[4,3] /  summarised_data[3,3]
  ratio_sf3 <- summarised_data[6,3] /  summarised_data[5,3]
  ratio_wt1 <- summarised_data[8,3] /  summarised_data[7,3]
  ratio_wt2 <- summarised_data[10,3] /  summarised_data[9,3]
  ratio_wt3 <- summarised_data[12,3] /  summarised_data[11,3]
  
  WT1 <- "WT1"
  WT2 <- "WT2"
  WT3 <- "WT3"
  SF1 <- "SF1"
  SF2 <- "SF2"
  SF3 <- "SF3"
  
  df_wt1 <- data.frame(x=WT1, y=ratio_wt1)
  colnames(df_wt1) <- c("primers", "ratio")
  
  df_wt2 <- data.frame(x=WT2, y=ratio_wt2)
  colnames(df_wt2) <- c("primers", "ratio")
  
  df_wt3 <- data.frame(x=WT3, y=ratio_wt3)
  colnames(df_wt3) <- c("primers", "ratio")
  
  df_sf1 <- data.frame(x=SF1, y=ratio_sf1)
  colnames(df_sf1) <- c("primers", "ratio")
  
  df_sf2 <- data.frame(x=SF2, y=ratio_sf2)
  colnames(df_sf2) <- c("primers", "ratio")
  
  df_sf3 <- data.frame(x=SF3, y=ratio_sf3)
  colnames(df_sf3) <- c("primers", "ratio")
  
  inter_1 <- union(df_wt1, df_wt2)
  inter_2 <- union(inter_1, df_wt3)
  inter_3 <- union(inter_2, df_sf1)
  inter_4 <- union(inter_3, df_sf2)
  treatment_summary <- union(inter_4, df_sf3)
  
  mean_sq_wt <- filter(treatment_summary, primers == "WT1" |primers == "WT3") %>% pull(ratio)
  mean_SQ_wt_ <- mean(mean_sq_wt)
  mean_SQ_wt_
  
  mean_sq_sf <- filter(treatment_summary, primers == "SF1"| primers == "SF2"|  primers == "SF3") %>% pull(ratio)
  mean_SQ_sf_ <- mean(mean_sq_sf)
  mean_SQ_sf_
  
  primer <- c("rpl32")
  mean_sq <- c(mean_sq_wt)
  wt_summary <- data.frame(x = primer, y = mean_sq)
  colnames(wt_summary) <- c("primers", "mean_sq")
  
  primer <- c("test")
  mean_sq <- c(mean_sq_sf)
  sf_summary <- data.frame(x = primer, y = mean_sq)
  colnames(sf_summary) <- c("primers", "mean_sq")
  
  sd_rc_wt <- sd(wt_summary$mean_sq)
  sd_rc_wt
  
  sd_rc_sf <- sd(sf_summary$mean_sq)
  sd_rc_sf
  
  sd_sum <- c(sd_rc_wt, sd_rc_sf)
  
  primer <- c("ctrl", "SF")
  mean_sq <- c(mean_SQ_wt_, mean_SQ_sf_)
  
  replicatt <- c(rep,rep)
  
  final_summary <- data.frame(x = primer, y = mean_sq, w = sd_sum, z= replicatt)
  colnames(final_summary) <- c("primers", "mean_rc", "sd", "replicat")
  final_summary
  
}


gene_graph <- function(tab1, tab2, prim) {
  
  final_summary <-bind_rows(tab1, tab2)
  
  
  final_summary$primers <- as.character(final_summary$primers)
  final_summary$primers <- factor(final_summary$primers)
  
  title <- paste("           ", prim)
  
  ggplot(final_summary, aes(x = replicat, y = mean_rc, fill = primers)) +
    geom_bar(stat="identity", position=position_dodge())+
    ggtitle(prim) +
    geom_errorbar(aes(ymin=mean_rc-sd, ymax=mean_rc+sd), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_brewer(palette="Reds")+
    xlab("Replicates")+
    ylab("ratio mean SQ")+
    coord_cartesian(ylim=c(0,3))+
    labs(fill = "Conditions")+
    theme_classic()+
    theme(
      plot.title = element_text(color="black", size=16, face="bold", hjust =0.5),
      axis.title.x = element_text(color="blue", size=14, face="bold"),
      axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) 
}


val_abe <- function(filename, hk_gene, test_gene, row_hk, row_test, plot_name){
  
  ###########################################################################
  #######IMPORTING DATA FROM EXCEL TO R AND SETUP ###########################
  ###########################################################################
  
  qpcr_data <- read.csv(file = filename, stringsAsFactors = FALSE, sep=';')
  head(qpcr_data)
  
  primer_key <- data.frame(row = c(row_hk,row_test), #Donne le nom de la ligne (A01,A02,.. =A) contenant le HK gene et l'autre
                           primers = c(rep("hk_gene", 1), rep("test_gene", 1)))
  
  tidy_data <- separate(qpcr_data, Well, into = c("row", "column"), #Cr?er une colonne "Row" qui contien B ou F
                        sep = 1, convert = TRUE)
  head(tidy_data)
  
  tidy_data <- left_join(tidy_data, primer_key, by = "row") #Tulise Left_join de dplyr pour coller les deux tableaux en utilisant row comme matrice
  
  #On visualise la plaque pour contr?ler que notre tableau est correcte
  ggplot(tidy_data, aes(x = column, y = row, fill = primers, label = Sample.Name)) +
    geom_tile(colour = "black") +
    geom_text() +
    scale_y_discrete(limits = c("row_hk", "row_test")) +
    scale_x_continuous(breaks = 1:12)
  
  ###########################################################################
  ##############################ANALYZE DATA#################################
  ###########################################################################
  
  #Calculate average Cq values#
  
  head(tidy_data$Cq)
  class(tidy_data$Cq)
  
  summarised_data <- tidy_data %>%
    mutate(Cq = as.numeric(as.character(Cq))) %>%
    group_by(Sample.Name, primers) %>%
    summarise(mean_Cq = mean(Cq))
  
  head(summarised_data)
  
  #Permet de controler nos data et ?liminer des r?sultats ?
  ggplot(summarised_data, aes(x = Sample.Name, y = mean_Cq, colour = primers)) +
    geom_point()
}
qPCR_analysis_SQ <- function(filename, row_hk, row_test, rep){
  
  ###########################################################################
  #######IMPORTING DATA FROM EXCEL TO R AND SETUP ###########################
  ###########################################################################
  
  qpcr_data <- read.csv(file = filename, stringsAsFactors = FALSE, sep=';')
  head(qpcr_data)
  
  primer_key <- data.frame(row = c(row_hk,row_test), #Give the the row number of your housekeeping gene and test gene
                           primers = c(rep("hk_gene", 1), rep("test_gene", 1)))
  
  tidy_data <- separate(qpcr_data, Well, into = c("row", "column"), #Create a column row
                        sep = 1, convert = TRUE)
  head(tidy_data)
  
  tidy_data <- left_join(tidy_data, primer_key, by = "row") #Use Left_join by dplyr to past both tabs using row as matrix
  
  #We control the sheet disposition
  ggplot(tidy_data, aes(x = column, y = row, fill = primers, label = Sample.Name)) +
    geom_tile(colour = "black") +
    geom_text() +
    scale_y_discrete(limits = c("row_hk", "row_test")) +
    scale_x_continuous(breaks = 1:12)
  
  ###########################################################################
  ##############################ANALYZE DATA#################################
  ###########################################################################
  
  #Calculate average SQ values
  
  head(tidy_data$SQ)
  class(tidy_data$SQ)
  
  summarised_data <- tidy_data %>%
    mutate(Cq = as.numeric(as.character(SQ))) %>%
    group_by(Sample.Name, primers) %>%
    summarise(mean_SQ = mean(SQ))
  
  head(summarised_data)
  
  #Look the homogeneity of our result and allow us to eliminate aberrant data
  ggplot(summarised_data, aes(x = Sample.Name, y = mean_SQ, colour = primers)) +
    geom_point()
  
  #for i in [1,lenght(summarized_data[3])], j in [2,lenght(summarized_data)+1];
  #ratio_wt1 <- summarised_data[i,3] /  summarised_data[j,3]
  #j=j+2
  #i=i+2
  
  ratio_wt1 <- summarised_data[2,3] /  summarised_data[1,3]
  ratio_wt2 <- summarised_data[4,3] /  summarised_data[3,3]
  ratio_wt3 <- summarised_data[6,3] /  summarised_data[5,3]
  ratio_sf1 <- summarised_data[8,3] /  summarised_data[7,3]
  ratio_sf2 <- summarised_data[10,3] /  summarised_data[9,3]
  ratio_sf3 <- summarised_data[12,3] /  summarised_data[11,3]
  
  WT1 <- "WT1"
  WT2 <- "WT2"
  WT3 <- "WT3"
  SF1 <- "SF1"
  SF2 <- "SF2"
  SF3 <- "SF3"
  
  df_wt1 <- data.frame(x=WT1, y=ratio_wt1)
  colnames(df_wt1) <- c("primers", "ratio")
  
  df_wt2 <- data.frame(x=WT2, y=ratio_wt2)
  colnames(df_wt2) <- c("primers", "ratio")
  
  df_wt3 <- data.frame(x=WT3, y=ratio_wt3)
  colnames(df_wt3) <- c("primers", "ratio")
  
  df_sf1 <- data.frame(x=SF1, y=ratio_sf1)
  colnames(df_sf1) <- c("primers", "ratio")
  
  df_sf2 <- data.frame(x=SF2, y=ratio_sf2)
  colnames(df_sf2) <- c("primers", "ratio")
  
  df_sf3 <- data.frame(x=SF3, y=ratio_sf3)
  colnames(df_sf3) <- c("primers", "ratio")
  
  inter_1 <- union(df_wt1, df_wt2)
  inter_2 <- union(inter_1, df_wt3)
  inter_3 <- union(inter_2, df_sf1)
  inter_4 <- union(inter_3, df_sf2)
  treatment_summary <- union(inter_4, df_sf3)
  
  mean_sq_wt <- filter(treatment_summary, primers == "WT1" |  primers == "WT2"|  primers == "WT3") %>% pull(ratio)
  mean_SQ_wt_ <- mean(mean_sq_wt)
  mean_SQ_wt_
  
  mean_sq_sf <- filter(treatment_summary, primers == "SF1" |  primers == "SF2"|  primers == "SF3") %>% pull(ratio)
  mean_SQ_sf_ <- mean(mean_sq_sf)
  mean_SQ_sf_
  
  primer <- c("rpl32")
  mean_sq <- c(mean_sq_wt)
  wt_summary <- data.frame(x = primer, y = mean_sq)
  colnames(wt_summary) <- c("primers", "mean_sq")
  
  primer <- c("test")
  mean_sq <- c(mean_sq_sf)
  sf_summary <- data.frame(x = primer, y = mean_sq)
  colnames(sf_summary) <- c("primers", "mean_sq")
  
  sd_rc_wt <- sd(wt_summary$mean_sq)
  sd_rc_wt
  
  sd_rc_sf <- sd(sf_summary$mean_sq)
  sd_rc_sf
  
  sd_sum <- c(sd_rc_wt, sd_rc_sf)
  
  primer <- c("SF", "ctrl")
  mean_sq <- c(mean_SQ_wt_, mean_SQ_sf_)
  
  replicatt <- c(rep,rep)
  
  final_summary <- data.frame(x = primer, y = mean_sq, w = sd_sum, z= replicatt)
  colnames(final_summary) <- c("primers", "mean_rc", "sd", "replicat")
  final_summary
  
}