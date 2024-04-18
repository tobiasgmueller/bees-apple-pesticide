# code copied from pesticide_analysis.R
# but without ld50 adjustments for bees
# running this code will create and save all graphs and tables used for analysis

rm(list = ls()) # clean environment to start


# load packages ####
library(reshape2)
library(gridExtra)
library(tidyverse) # for tidyness
library(vegan) # for ordinations
library(ggpubr) # for multipanel graphs
library(cowplot) # new way to multi graph
library(rstatix) # pipe friendly stats
library(plotly) # 3d NMDS
library(ggrepel) # for graphy pretty
library(rcompanion) # auto letters after pairwise tests
library(ggnewscale) #allows for multiple color palletes

# for color patterns
remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

# for pairwise permanova 
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis",force = TRUE)
library(pairwiseAdonis)

#a few packages for graph colors
library(scales)
library(ggpattern)
library(RColorBrewer)

# non CRAN package for nice NMDS graphs in ggplot
install.packages("remotes")
remotes::install_github("jfq3/ggordiplots")
library(ggordiplots)


# a few custom functions
#for pretty log tick marks on graphs
base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

# for calculating SE
se <- function(x) sqrt(var(x)/length(x))





# DATA PREP -------------------
### read in data --------------------------------------------------



# read in data
alldf<- read.csv("input/all_pesticide.csv", stringsAsFactors = TRUE) # field collected data
ld50<- read.csv("input/ld50_per_bee.csv", stringsAsFactors = TRUE) # list of ld50 values for pesticides
pesticidenames<-read.csv("input/pesticide_type.csv", stringsAsFactors = TRUE) # some information on what each pesticide is 
pesticidenames$class<-as.factor(pesticidenames$class)
nesting_type <- read.csv("input/nesting_type.csv", stringsAsFactors = TRUE) # information on bee nesting type                                                            


#all column of short names
shortname <- read.csv("input/short_name.csv", stringsAsFactors = TRUE)
alldf1 <- merge(alldf,shortname, by="Sample.type")
alldf1$Sample.type <- alldf1$bee.name


# drop columns that I dont use
alldf1 <- alldf1 %>%
  select(!c("Unique.ID", "Filename", "bee.name")) 


# add a column of summed pesticide quantity
alldf2<-alldf1 %>% 
  mutate(total = select(.,!c("Sample.type",
                                   "Sample.ID", 
                                   "Site", "mass_g"))%>% 
  rowSums(na.rm = TRUE))

# delete colums that sum to zero (never detected pesticides)
alldf<- alldf2 %>% select(c("Sample.type","Sample.ID", "Site", "mass_g", where( ~ is.numeric(.x) && sum(.x) != 0)))



#adjusting all df to HQ ####


# first we need to adjust exposure values from PPB to ug/bee

# PPB = ug/kg so to get to ug/bee we do PPB*bee weight (IN KG)

#drop identifier colummns
pesticide_per_bee <- alldf%>%
  select(!c("Sample.ID", "Sample.type", "Site", "total"))%>%
  mutate(across(.fns = ~.*(mass_g/1000)))%>%
  select(!"mass_g")
  


#split LD50 into contact and oral
ld50contact <- ld50%>%
  filter(Common.name == "Honey bee contact LD50 (ug/bee)")

ld50oral <- ld50%>%
  filter(Common.name == "Honey bee oral LD50 (ug/bee)")


# first find common columns 
common_cols <- intersect(colnames(pesticide_per_bee), colnames(ld50))

# drop extra columns from ld50
common_oral<-  ld50oral %>% select({common_cols})
common_contact<-  ld50contact %>% select({common_cols})

# and then divide pesticides by oral and contact ld50 values
oralhq<-pesticide_per_bee / common_oral[col(pesticide_per_bee)]
contacthq<-pesticide_per_bee / common_contact[col(pesticide_per_bee)]

# now add back in identifying information
info <- alldf%>%
  select("Sample.type","Site", "Sample.ID", "mass_g")

contactdf<- cbind(info,contacthq)
oraldf<- cbind(info,oralhq)


# now lets add a column to see total additive ld50 
# (this assumes pesticides are purely additive)
df_oral<-oraldf %>% 
  mutate(chemsum = select(.,!c("Sample.type",
                                   "Sample.ID", 
                                   "Site", "mass_g"))%>% 
  rowSums(na.rm = TRUE))


df_contact<-contactdf %>% 
  mutate(chemsum = select(.,!c("Sample.type",
                                   "Sample.ID", 
                                   "Site", "mass_g"))%>% 
  rowSums(na.rm = TRUE))





### combining oral and contact ####
oralsum <- df_oral %>%
  select("Sample.type", "Site","Sample.ID", "chemsum")%>%
  rename(oralhq = chemsum)

contactsum <- df_contact %>%
  select("Sample.ID", "chemsum")%>%
  rename(contacthq = chemsum)

oral_contacthq <- left_join(oralsum, contactsum, 
                           by = "Sample.ID")

melt_hq <- reshape2::melt(oral_contacthq, 
                         id.vars=c("Sample.type", "Site","Sample.ID"),
                         value.name = "hq",
                         variable.name="exposure.type")







### melting for graphing ####
# now lets melt that data for graphing and playing with


# melt all data
df_melt <- alldf %>%
  pivot_longer(!c(Site, Sample.type, Sample.ID, mass_g), names_to = "ai", values_to = "quantity")
df_melt$ai <- as.factor(df_melt$ai)


oralmelt1 <- df_oral %>%
  pivot_longer(!c(Site, Sample.type, Sample.ID, mass_g), names_to = "pesticide", values_to = "oralhq")
oralmelt1$pesticide <- as.factor(oralmelt1$pesticide)

contactmelt1 <- df_contact %>%
  pivot_longer(!c(Site, Sample.type, Sample.ID, mass_g), names_to = "pesticide", values_to = "contacthq")
contactmelt1$pesticide <- as.factor(contactmelt1$pesticide)


#now add information about the pesticide type
melt_oral <- left_join(oralmelt1, pesticidenames, by = "pesticide")
melt_contact <- left_join(contactmelt1, pesticidenames, by = "pesticide")
melt_all <- left_join(df_melt, pesticidenames, by = c("ai" = "pesticide"))




# cleaning ####

# then rename 2 things for clarity and delete the rest
df_all <- alldf

rm(list = c("alldf1",
            "alldf2",
            "contacthq",
            "oralhq",
            "ld50",
            "ld50contact",
            "ld50oral",
            "common_cols",
            "pesticidenames",
            "pesticide_per_bee",
            "shortname",
            "df_melt",
            "oraldf",
            "contactdf",
            "oralmelt1",
            "contactmelt1",
            "alldf",
            "oral_contacthq",
            "oralsum",
            "contactsum",
            "info",
            "common_contact",
            "common_oral"
            ))


# END OF DATA PREP ------------------ 




# proportion of bees over  limit ####
sample_count<-df_all %>%  
  group_by(Sample.type)%>%
   summarize(n  = n())

prop_over_03<- melt_hq %>% 
  filter(exposure.type == "oralhq", hq >= .03) %>%
  distinct(Sample.ID, .keep_all = TRUE)%>%
  group_by(Sample.type)%>%
  mutate(count  = n()) %>%
  left_join(sample_count)%>%
  distinct(Sample.type, .keep_all = TRUE)%>%
  mutate(percent_over_03 = (count/n)*100)%>%
  select("Sample.type","percent_over_03") 
  # then do the same with other thresholdes
  
prop_over_2<- melt_hq %>% 
  filter(exposure.type == "contacthq", hq >= .2) %>%
  distinct(Sample.ID, .keep_all = TRUE)%>%
  group_by(Sample.type)%>%
  mutate(count  = n()) %>%
  left_join(sample_count)%>%
  distinct(Sample.type, .keep_all = TRUE)%>%
  mutate(percent_over_2 = (count/n)*100)%>%
  select("Sample.type","percent_over_2")
  
prop_over_4<- melt_hq %>% 
  filter(exposure.type == "contacthq", hq >= .4) %>%
  distinct(Sample.ID, .keep_all = TRUE)%>%
  group_by(Sample.type)%>%
  mutate(count  = n()) %>%
  left_join(sample_count)%>%
  distinct(Sample.type, .keep_all = TRUE)%>%
  mutate(percent_over_4 = (count/n)*100)%>%
  select("Sample.type","percent_over_4")
  
  
  
  # then merge those 3 df into large on and make a table

prop_over_limit <- left_join(prop_over_03, prop_over_2) %>%
                left_join(., prop_over_4) 



prop_over_limit


# FIGURE - oral and contact HQ ####


# create thresholds for lines

thresholds <- data.frame (agency  = c("epa_acute",
                                        "efsa_acute",
                                        "efsa_chronic"),
                            limit = c(.4, .2,  .03))%>% 
    mutate(agency=fct_relevel(agency, c("epa_acute",
                                        "efsa_acute",
                                        "efsa_chronic")))
    
  
rq_unadjusted_boxplot<- 
   melt_hq %>%
      filter(
           !Sample.type %in% c("Apple flower","Dandelion flower") # comment out if you want all samples. otherwise this removed flowers
           )%>% 
    mutate(exposure.type=fct_relevel(exposure.type, c("contacthq","oralhq"))) %>%
    
    
    
    
    ggplot( 
         aes(
           x=Sample.type, y=hq))+
    geom_hline(data=thresholds, aes(yintercept = limit, linetype = agency),
               size=1.3)+
    geom_hline(data=thresholds, aes(yintercept = limit),
               size=1.3, 
               color = "black",
               alpha = .25)+
    scale_linetype_manual(name = "Levels of Concern",
                          labels = c("EPA acute",
                                     "EFSA acute",
                                     "EFSA chronic"
                                     ),
                          values = c(2, 4, 3),
     ) +
    
    scale_color_brewer(palette="Set1", labels = c("Contact", "Oral"), direction = -1)+
    scale_fill_brewer(palette="Set1", labels = c("Contact", "Oral"), direction = -1)+
  
    geom_boxplot(aes(color=exposure.type), linewidth = 1, fill = NA, outlier.shape=NA)+
    geom_point(aes(fill=exposure.type), position = position_jitterdodge(jitter.width = 0.1, 
                                               jitter.height = 0, 
                                               dodge.width = 0.75),
               size=3,
               stroke=1,
               pch=21,
               color="black")+
    scale_y_continuous(
      trans = log_trans(),
                     breaks = base_breaks(),
                     labels = function(x) paste0(x*100, "%")  # change scientific notation to percent
  ) +
    ylab("Risk Quotient (percent of LD50)")+
    xlab("Sample Type")+
    labs(fill="Exposure Type")+
    expand_limits(y=15)+
    guides(color="none")+
    theme_bw()+
    theme(axis.text.x = element_text(face = "italic", angle = 90, vjust = 0.5, hjust=1),
          text = element_text(size = 20),
          axis.title.x = element_blank(),
          legend.key.width=unit(2.1,"cm"))
  
  
  
  rq_unadjusted_boxplot

ggsave(plot = rq_unadjusted_boxplot, file = "output/graphs/unadjusted/RQ_unadjusted.jpeg", width  =10, height =8, dpi = 300)

ggsave(plot = rq_unadjusted_boxplot, file = "output/graphs/unadjusted/RQ_unadjusted.svg", width  =10, height =8)



# HQ across bees ####

# kw of hq across just bees
melt_hq %>%
   filter(!Sample.type %in% c("Apple flower","Dandelion flower"))%>% 
  group_by(exposure.type)%>%
  kruskal_test(hq ~ Sample.type)%>%
  ungroup()%>%
  write_csv(file="output/tables/unadjusted/kw_hq_across_bees_UNADJUSTED.csv")

# dunn test of hq across just bees
dunnhq<- melt_hq %>%
   filter(!Sample.type %in% c("Apple flower","Dandelion flower"))%>% 
  group_by(exposure.type)%>%
  dunn_test(hq ~ Sample.type, p.adjust.method = "holm")%>%
  ungroup()%>%
  write_csv(file="output/tables/unadjusted/dunn_hq_across_bees_UNADJUSTED.csv")

# auto generate sig letters
letter_oral<- dunnhq%>%
  filter(exposure.type == "oralhq")%>%
  unite(comparison, c(group1, group2), sep = "-")%>%
  cldList(p.adj ~ comparison, data = .)

# auto generate sig letters
letter_contact<-dunnhq%>%
  filter(exposure.type == "contacthq")%>%
  unite(comparison, c(group1, group2), sep = "-")%>%
  cldList(p.adj ~ comparison, data = .)



#permanova RQ ####
# overall permanova of unadjusted contact RQ across bees

perm_df_hq<- df_contact %>%
  arrange(Sample.ID)%>%
  filter(!chemsum == 0)

perm_matrix_hq <- perm_df_hq%>% 
  select(-c(Site, 
            Sample.type,
            Sample.ID,
            chemsum,
            mass_g
            ))%>%
  as.matrix()

perm_bee_hq <- (perm_df_hq[,"Sample.type"]) # make vector of bees for later grouping

bee.dis.hq <- vegdist(perm_matrix_hq, method = 'bray')
beepermanhq<- adonis2(bee.dis.hq~Sample.type, data=perm_df_hq, permutations=999, method="bray")
beepermanhq

write_csv(beepermanhq, "output/tables/unadjusted/permanova_contact_bees_UNADJUSTED.csv")


## pairwise permanova RQ ####

# contact HQ
pair.contact<-pairwise.adonis(bee.dis.hq, factors = perm_bee_hq)
write_csv(pair.contact, "output/tables/unadjusted/pairwise_permanova_contact_UNADJUSTED.csv")



# simper ####
hq_sparse<- perm_df_hq%>% 
  select(-c(Site, 
            Sample.type,
            chemsum,
            mass_g
            ))%>%
 column_to_rownames('Sample.ID') # get rid of descriptors and make Sample ID the row names


#  overall simper
sim_hq_all<- simper(hq_sparse, permutations = 999)
#and store output
simper_hq_all_output <- summary(sim_hq_all, ordered = TRUE,
    digits = max(3,getOption("digits") - 3))
capture.output(summary(sim_hq_all), file = "output/tables/unadjusted/sim_hq_all_UNADJUSTED.txt")



## pairwise Simper ####
sim_hq<- simper(hq_sparse, group = perm_bee_hq, permutations = 999)

# now put output into txt file
simperoutput_hq <- summary(sim_hq, ordered = TRUE,
    digits = max(3,getOption("digits") - 3))
capture.output(summary(sim_hq), file = "output/tables/unadjusted/sim_hq_UNADJUSTED.txt")