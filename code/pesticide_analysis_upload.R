# Tobias Mueller
# code for apple pesticide manuscript


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

#a few packages for graph colors
library(scales)
library(ggpattern)
library(RColorBrewer)
# and for color patterns
remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

# nice NMDS graphs in ggplot
install.packages("remotes")
remotes::install_github("jfq3/ggordiplots")
library(ggordiplots)

# pairwise permanova 
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis",force = TRUE)
library(pairwiseAdonis)

# then a few custom functions
# for pretty log tick marks on graphs
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
ld50<- read.csv("input/ld50.csv", stringsAsFactors = TRUE) # list of ld50 values for pesticides
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

                          
#drop identifier colummns
justpesticides <- alldf%>%
  select(!c("Sample.ID", "Sample.type", "Site", "total", "mass_g"))

#split LD50 into contact and oral
ld50contact <- ld50%>%
  filter(Common.name == "Honey bee contact LD50 (ppb)")

ld50oral <- ld50%>%
  filter(Common.name == "Honey bee oral LD50 (ppb)")


# first find common columns 
common_cols <- intersect(colnames(justpesticides), colnames(ld50))

# drop extra columns from ld50
common_oral<-  ld50oral %>% select({common_cols})
common_contact<-  ld50contact %>% select({common_cols})

# and then divide pesticides by oral and contact ld50 values
oralhq<-justpesticides / common_oral[col(justpesticides)]
contacthq<-justpesticides / common_contact[col(justpesticides)]



## adjust ld50 to bees ####

# now I want to take all ld50 values and then adjust them based on the weight of the bee (mass_g)
# currently I have everything as percent of ld50 of honeybee

# this is done by calculating a ratio of bee ld50 based on the samples mass
# vs average mass honeybee ld 50 
# then applying this ratio to all LD50 findings for each bee

# first to make bee specific ratio number

# save lm variables
# averaged by species and compound
#these numbers are taken from R code in LD50 adjust
m <- 1.218592
b <- -3.176504



ratio <- alldf%>%
  select("Sample.type","Site", "Sample.ID", "mass_g")%>%
  mutate(ratio = 
     (10^((m * (log10(100))) + b))/   #honey bee ld50(per 100mg bee) calculated based on log weight
     (10^((m * (log10(mass_g*1000))) + b))  # sample bee ld50 based on log weight adjusted to mg
      )



# thencombine with honeybee ld50 
ratio_contact<- cbind(ratio,contacthq)
ratio_oral<- cbind(ratio,oralhq)


# and then divide  across by the ratio
# and  drop flower samples since these make no sense adjusted to weight
# and drop the ratio column

contactdf<- ratio_contact %>%
  mutate(
    across(!c("Sample.type","Site", "Sample.ID", "mass_g"),
           .fns = ~.*ratio)) %>%
  select(!c("ratio")) %>%
  filter(!Sample.type %in% c("Apple flower","Dandelion flower"))


oraldf<- ratio_oral %>%
  mutate(
    across(!c("Sample.type","Site", "Sample.ID", "mass_g"),
           .fns = ~.*ratio))%>%
  select(!c("ratio")) %>%
  filter(!Sample.type %in% c("Apple flower","Dandelion flower"))



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

# then rename for clarity and delete the rest
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
            "justpesticides",
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
            "ratio",
            "ratio_contact",
            "ratio_oral",
            "common_contact",
            "common_oral"
            ))




# END OF DATA PREP ------------------ 


## proportion of clean samples #
df_all %>%
  group_by(Sample.type)%>%
  summarize(prop_clean = mean(!total))%>%
  write_csv("output/tables/proportion_clean_samples.csv")


## list of samples that exceed thresholds #
df_oral %>%
  filter(chemsum >= .03)%>%
  write_csv("output/tables/bees_over_efsa_oral_ADJUSTED.csv")

df_contact %>%
  filter(chemsum >= .2)%>%
  write_csv("output/tables/bees_over_efsa_contact_ADJUSTED.csv")

df_contact %>%
  filter(chemsum >= .4)%>%
  write_csv("output/tables/bees_over_epa_contact_ADJUSTED.csv")


# bees over LOC ####
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
  
prop_over_2contact<- melt_hq %>% 
  filter(exposure.type == "contacthq", hq >= .2) %>%
  distinct(Sample.ID, .keep_all = TRUE)%>%
  group_by(Sample.type)%>%
  mutate(count  = n()) %>%
  left_join(sample_count)%>%
  distinct(Sample.type, .keep_all = TRUE)%>%
  mutate(percent_over_2contact = (count/n)*100)%>%
  select("Sample.type","percent_over_2contact")
  
prop_over_2oral<- melt_hq %>% 
  filter(exposure.type == "oralhq", hq >= .2) %>%
  distinct(Sample.ID, .keep_all = TRUE)%>%
  group_by(Sample.type)%>%
  mutate(count  = n()) %>%
  left_join(sample_count)%>%
  distinct(Sample.type, .keep_all = TRUE)%>%
  mutate(percent_over_2oral = (count/n)*100)%>%
  select("Sample.type","percent_over_2oral")

prop_over_4contact<- melt_hq %>% 
  filter(exposure.type == "contacthq", hq >= .4) %>%
  distinct(Sample.ID, .keep_all = TRUE)%>%
  group_by(Sample.type)%>%
  mutate(count  = n()) %>%
  left_join(sample_count)%>%
  distinct(Sample.type, .keep_all = TRUE)%>%
  mutate(percent_over_4contact = (count/n)*100)%>%
  select("Sample.type","percent_over_4contact")
 
prop_over_4oral<- melt_hq %>% 
  filter(exposure.type == "oralhq", hq >= .4) %>%
  distinct(Sample.ID, .keep_all = TRUE)%>%
  group_by(Sample.type)%>%
  mutate(count  = n()) %>%
  left_join(sample_count)%>%
  distinct(Sample.type, .keep_all = TRUE)%>%
  mutate(percent_over_4oral = (count/n)*100)%>%
  select("Sample.type","percent_over_4oral") 
  
  
  # then merge those dfs and make a table

loc_list = list(prop_over_03,
                             prop_over_2contact,
                             prop_over_2oral,
                             prop_over_4contact,
                             prop_over_4oral)

prop_over_limit <- reduce(loc_list,left_join, by = 'Sample.type')%>%
  write_csv("output/tables/bees_over_loc_ADJUSTED.csv")





## percent of hq from pesticides ####

percent_hq_contact<- melt_contact %>%
  filter(pesticide != "chemsum")%>%
  group_by(Sample.ID)%>%
  mutate(percent = contacthq/sum(contacthq))

percent_hq_oral<- melt_oral %>%
  filter(pesticide != "chemsum")%>%
  group_by(Sample.ID)%>%
  mutate(percent = oralhq/sum(oralhq))


## mean PPB of bee pesticide ####
mean_exposure <- melt_all %>%
  filter(! Sample.type %in% c("Apple flower", "Dandelion flower")) %>%
  group_by(Sample.type, ai)%>%
  summarise(mean = mean(quantity))


mean_exposure_by_bee <- df_all %>%
  filter(! Sample.type %in% c("Apple flower", "Dandelion flower")) %>%
  group_by(Sample.type)%>%
  summarise(mean = mean(total), se = se(total))

mean_exposure_all_bee <- df_all %>%
  filter(! Sample.type %in% c("Apple flower", "Dandelion flower")) %>%
  summarise(mean = mean(total), se = se(total))


mean_contact_hq <- melt_contact %>%
  filter(! Sample.type %in% c("Apple flower", "Dandelion flower")) %>%
  group_by(Sample.type, pesticide)%>%
  summarise(mean = mean(contacthq))



## pesticide count across bee ####

df_count <- df_all %>%
  mutate(count  = rowSums((select(.,!c("Sample.type",
                                   "Sample.ID", 
                                   "Site", 
                              "total",
                              "mass_g")))!=0)) 
           
                  
# average number of pesticide findings
df_count %>%
  group_by(Sample.type) %>% 
  summarise(pest_num = mean(count), se = se(count))%>%
  ungroup()%>%
  write_csv(file="output/tables/pesticide_count.csv")



# FIGURE - heat map ####


# this creates 2 heat maps, one for bees and one for flowers
# and then stacks them 

heatmap_exp_bees <- melt_all %>%
   # relevel so total is at the bottom
  mutate(ai=fct_relevel(ai,"total", after = Inf),
         # then also create new column to facet total seperately
         total = ifelse(ai == "total", "yes", "no"),
         # and relevel bees so flowers are all on one side
         Sample.type = fct_relevel(Sample.type, c("Apple flower", "Dandelion flower"), after = Inf)
         )%>%
  filter(!Sample.type %in% c("Apple flower", "Dandelion flower"))%>%
  ggplot(aes(Sample.ID,ai, fill = log10(quantity)))+
  #add border white colour of line thickness 0.25
  geom_tile(colour="white", size=.3)+
  #remove x and y axis labels
  labs(x="", y="")+
  #remove extra space
  scale_y_discrete(expand=c(0, 0))+
  #set a base size for all fonts
  theme_grey(base_size=20)+
#coord_equal()+
  #theme options
  theme(
    #bold font for legend text
    legend.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.4),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank(),
    #remove x axis ticks
    axis.text.x=element_blank(),
    # removed y axis facet labels
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 18, face="bold.italic"))+
  scale_fill_gradientn(values=c(1, .6, .5, .4, 0),
                       breaks=c(-1,0,1,2,3,4), 
                       labels=c(10^-1, 
                                10^0, 
                                10^1, 
                                10^2, 
                                10^3,
                                10^4),
                       colours=c("black", 
                                                               "darkred", 
                                                               "red", 
                                                               "orange", 
                                                               "yellow2"))+
  guides(fill = guide_colourbar(barwidth = 1.2,
                                barheight = 28,
                                title = "ppb"))+
  facet_grid(type~Sample.type, scale = 'free', space = "free", switch="both")

heatmap_exp_bees


# and now for flowers
heatmap_exp_flowers <- melt_all %>%
   # relevel so total is at the bottom
  mutate(ai=fct_relevel(ai,"total", after = Inf),
         # then also create new column to facet total seperately
         total = ifelse(ai == "total", "yes", "no"),
         # and relevel bees so flowers are all on one side
         Sample.type = fct_relevel(Sample.type, c("Apple flower", "Dandelion flower"), after = Inf)
         )%>%
  filter(Sample.type %in% c("Apple flower", "Dandelion flower"))%>%
  ggplot(aes(Sample.ID,ai, fill = log10(quantity)))+
  #add border white colour of line thickness 0.25
  geom_tile(colour="white", size=.3)+
  #remove x and y axis labels
  labs(x="", y="")+
  #remove extra space
  scale_y_discrete(expand=c(0, 0))+
  #set a base size for all fonts
  theme_grey(base_size=20)+
#coord_equal()+
  #theme options
  theme(
    #bold font for legend text
    legend.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.4),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank(),
    #remove x axis ticks
    axis.text.x=element_blank(),
    # removed y axis facet labels
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 18, face="bold.italic"))+
  scale_fill_gradientn(values=c(1, .6, .5, .4, 0),
                       breaks=c(-1,0,1,2,3,4), 
                       labels=c(10^-1, 
                                10^0, 
                                10^1, 
                                10^2, 
                                10^3,
                                10^4),
                       colours=c("black", 
                                                               "darkred", 
                                                               "red", 
                                                               "orange", 
                                                               "yellow2"))+
  guides(fill = guide_colourbar(barwidth = 1.2,
                                barheight = 28,
                                title = "ppb"))+
  facet_grid(type~Sample.type, scale = 'free', space = "free", switch="both")

heatmap_exp_flowers


# and then combine heatmaps

heatmap_combo <- ggarrange(heatmap_exp_bees, heatmap_exp_flowers, common.legend = TRUE,
          ncol=1,
          legend = "right")

ggsave(plot=heatmap_combo, file = "output/graphs/heatmap_combo.png", width = 20, height = 17, bg = "#FFFFFF")

ggsave(plot=heatmap_combo, file = "output/graphs/heatmap_combo.svg", width = 20, height = 17, bg = "#FFFFFF")




# FIGURE - oral and contact HQ ####

# create LOC thresholds for lines
thresholds <- data.frame (agency  = c("epa_acute",
                                        "efsa_acute",
                                        "efsa_chronic"),
                            limit = c(.4, .2,  .03))%>% 
    mutate(agency=fct_relevel(agency, c("epa_acute",
                                        "efsa_acute",
                                        "efsa_chronic")))
    
  
  rq_boxplot<- 
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
  
  
  
  rq_boxplot

ggsave(plot = rq_boxplot, file = "output/graphs/RQ_adjusted.png", width  =10, height =8, dpi = 300)

ggsave(plot = rq_boxplot, file = "output/graphs/RQ_adjusted.svg", width  =10, height =8)



# FIGURE - counts by bee ####

df_count <- df_all %>%
  mutate(count  = rowSums((select(.,!c("Sample.type",
                                   "Sample.ID", 
                                   "Site", 
                              "total",
                              "mass_g")))!=0)) 
           

pest_count <- df_count %>%
  filter(!Sample.type %in% c("Apple flower","Dandelion flower"))%>%
  ggplot(aes(x=Sample.type, y=count, color = Sample.type))+
  geom_boxplot(outlier.shape=NA, linewidth = 1)+
  geom_point(aes(fill = Sample.type), position = position_jitter(width=.08, 
                                                                 height = 0),
             size=3,
             stroke = 1,
             pch=21,
             color="black")+
  ylab("Pesticide Count")+
  theme_bw()+
  theme(text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "italic", angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")+
  scale_y_continuous(limits = c(0, 15))+ # make y axis longer 
    # adjust colors for bees
  scale_color_manual(values=c("#FFC543",
                              "#5FAEDA",
                              "#487761",
                              "#3A27CE",
                              "#C14510",
                              "#5F1C3B"))+
  scale_fill_manual(values=c("#FFC543",
                             "#5FAEDA",
                             "#487761",
                             "#3A27CE",
                             "#C14510",
                             "#5F1C3B"))

  pest_count
  
ggsave(plot = pest_count, file = "output/graphs/pesticide_count.png", width  =10, height =8, dpi = 300)
ggsave(plot = pest_count, file = "output/graphs/pesticide_count.svg", width  =10, height =8, dpi = 300)
   


# FIGURE - ppb bees ####
pest_concentration <- df_all %>%
  filter(!Sample.type %in% c("Apple flower","Dandelion flower"))%>%
  ggplot(aes(x=Sample.type, y=total, color = Sample.type))+
  geom_boxplot(outlier.shape=NA, linewidth = 1)+
  geom_point(aes(fill = Sample.type), position = position_jitter(width=.08, 
                                                                 height = 0), 
             stroke = 1,
             size=3,
             pch=21,
             color="black")+
  ylab("Total Pesticide Concentration (Log PPB)")+
  theme_bw()+
  theme(axis.text.x = element_text(face = "italic", angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 20),
        axis.title.x = element_blank(), 
        legend.position = "none")+
  scale_y_log10(
    breaks = base_breaks(),
    limits=c(1,10000))+ 
    # adjust colors for bees
  scale_color_manual(values=c("#FFC543",
                              "#5FAEDA",
                              "#487761",
                              "#3A27CE",
                              "#C14510",
                              "#5F1C3B"))+
  scale_fill_manual(values=c("#FFC543",
                             "#5FAEDA",
                             "#487761",
                             "#3A27CE",
                             "#C14510",
                             "#5F1C3B"))





  pest_concentration


  
ggsave(plot = pest_concentration, file = "output/graphs/pesticide_cocentration_log.png", width  =10, height =8, dpi = 300)
ggsave(plot = pest_concentration, file = "output/graphs/pesticide_cocentration_log.svg", width  =10, height =8, dpi = 300)
   

#then combine both
combo<- ggarrange(pest_count, pest_concentration)
ggsave(plot = combo, file = "output/graphs/pesticide_count_conc_combo.png", width  =20, height =8, dpi = 300)
ggsave(plot = combo, file = "output/graphs/pesticide_count_conc_combo.svg", width  =20, height =8, dpi = 300)




# SUPP FIGURE - correlation graph ####

df_thia_corr_wide<- df_all %>%
  filter(Sample.type %in% c("Apple flower",
                            "Dandelion flower",
                            "A. mellifera"))%>%
  select(Site, Sample.type, Thiamethoxam)%>%
  pivot_wider(names_from = Sample.type, values_from = Thiamethoxam,
              values_fn = mean) %>% # this takes the mean flower for that site
  filter(!is.na(`A. mellifera`))

  
df_thia_corr_long<-df_thia_corr_wide%>% 
  pivot_longer(cols=c(`Apple flower`,`Dandelion flower`))
  
 corr_graph<- df_thia_corr_long%>%
  ggplot(aes(y=`A. mellifera`, x=value, color=name, shape=name))+
    geom_point(size=4)+
    geom_smooth(method=lm)+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    theme_bw()+
    labs(y="Thiamethoxam in A. mellifera ",
         x="Thiamethoxam in Flowers",
         color="Orchard flower",
         shape="Orchard flower")+
    scale_shape_manual(values = c(17, 20))+
    scale_color_manual(values = c("#FFC543",
                               "#5FAEDA"))


ggsave(corr_graph, file = "output/graphs/correlation_graph.png", width  =10, height =8, dpi = 300)  
ggsave(corr_graph, file = "output/graphs/correlation_graph.svg", width  =10, height =8, dpi = 300)  



cor.test(df_thia_corr_wide$`A. mellifera`, df_thia_corr_wide$`Apple flower`)
cor.test(df_thia_corr_wide$`A. mellifera`, df_thia_corr_wide$`Dandelion flower`)




#KW AND DUNN ####
df_count <- df_all %>%
  mutate(count  = rowSums((select(.,!c("Sample.type",
                                   "Sample.ID", 
                                   "Site", 
                              "total",
                              "mass_g")))!=0)) 
           

# kw of count across just bees
df_count%>%
  filter(!Sample.type %in% c("Apple flower","Dandelion flower"))%>% 
  kruskal_test(count~ Sample.type)%>%
    write_csv(file="output/tables/kw_count_across_bees.csv")
  
# dunn test of count across just bees
df_count%>%
  filter(!Sample.type %in% c("Apple flower","Dandelion flower"))%>% 
  dunn_test(count~ Sample.type, p.adjust.method = "holm")%>%
  write_csv(file="output/tables/dunn_count_across_bees.csv")


# kw of concentration across just bees
df_all%>%
  filter(!Sample.type %in% c("Apple flower","Dandelion flower"))%>% 
  kruskal_test(total~ Sample.type)%>%
    write_csv(file="output/tables/kw_conc_across_bees.csv")
  
# dunn test of concentration across bees
df_all%>%
  filter(!Sample.type %in% c("Apple flower","Dandelion flower"))%>% 
  dunn_test(total ~ Sample.type, p.adjust.method = "holm")




## HQ across bees ####

# kw of hq across just bees
melt_hq %>%
   filter(!Sample.type %in% c("Apple flower","Dandelion flower"))%>% 
  group_by(exposure.type)%>%
  kruskal_test(hq ~ Sample.type)%>%
  ungroup()%>%
  write_csv(file="output/tables/kw_hq_across_bees_ADJUSTED.csv")


# dunn test of hq across just bees
melt_hq %>%
   filter(!Sample.type %in% c("Apple flower","Dandelion flower"))%>% 
  group_by(exposure.type)%>%
  dunn_test(hq ~ Sample.type, p.adjust.method = "holm")%>%
  ungroup()%>%
  write_csv(file="output/tables/dunn_hq_across_bees_ADJUSTED.csv")

  



#NMDS    -------------------------------------
# now lets make an NMDS plot of pesticides by species

## HQ NMDS ####
### data prep hq ####
nmds_df_hq<- df_contact %>%
  arrange(Sample.ID)%>% 
  filter(!chemsum == 0, # get rid of bees with no pesticides as it makes nmds sad
         !Sample.type %in% c("Apple flower","Dandelion flower") # comment out if you want all samples. otherwise this removed flowers
         )%>%
  select(c("Sample.type","Sample.ID", "Site", "mass_g", "chemsum", where( ~ is.numeric(.x) && sum(.x) != 0))) %>% #drop pesticides that arnt in any samples
  left_join(nesting_type, by = c("Sample.type" = "bee"))%>%# add information about nesting type
  droplevels() 


nmds_hq_bee <- (nmds_df_hq[,"Sample.type"]) # make vector of bees for later grouping
nmds_hq_nesting <- (nmds_df_hq[,"nesting"]) # make vector of nesting type for later grouping


nmds_hq_sparse<- nmds_df_hq%>% 
  select(-c(Site, 
            Sample.type,
            chemsum,
            nesting,
            mass_g
            ))%>%
 column_to_rownames('Sample.ID') # get rid of descriptors and make Sample ID the row names



### create nmds ####
nmds_hq <- metaMDS(nmds_hq_sparse,  
                        distance = "bray",       
                        k = 3,
                        maxit = 999,
                        trymax = 500,
                        wascores = TRUE,
                        autotransform = FALSE)

 
stressplot(nmds_hq)

#color scheme for BEES ONLY in HQ
beecolor_hq= c("#00BA38", "#B79F00", "#00BFC4",   "#619CFF", "#F8766D", "#F564E3") # create colors for bees
nmds_df_hq$Sample.type <- as.factor(nmds_df_hq$Sample.type)
beecolor_hq_loop = c("#00BA38", "#B79F00", "#00BFC4",   "#619CFF", "#F8766D", "#F564E3")[nmds_df_hq$Sample.type]






## hq nmds ggplot ####

# first we create the plot with arrows

choices_hq = c(1,2) # first select which axis to plot

hq_nmds_plot<- gg_envfit(ord = nmds_hq,
          env =nmds_hq_sparse,
          groups = nmds_hq_bee,
          scaling = 1,
          choices = choices_hq,
          perm = 999,
          alpha = 1, 
          angle = 20, 
          len = 0.5, 
          arrow.col = "black", 
          pt.size = 3,
          plot = TRUE)

# add in arrows that are top drivers from simper
hq_arrows<- hq_nmds_plot$df_arrows %>%
  filter(var %in% c("Thiamethoxam",
                 "Difenoconazole",
                 "Indoxacarb",
                 "Carbaryl",
                 "Phosmet"))

# then create the same plot with ellipses
hq_ellipse_plot<- gg_ordiplot(ord = nmds_hq,
            groups = nmds_hq_bee,
            scaling = 1,
            choices = choices_hq,
            kind = "sd",
            conf = NULL,
            show.groups = "all",
            ellipse = TRUE,
            label = FALSE,
            hull = FALSE,
            spiders = FALSE,
            plot = FALSE)

#and now put those together


# then make a base plot
 hq_ord_gg <- ggplot() + theme_bw(16) + 
      labs(x="NMDS Axis 1", 
           y="NMDS Axis 2",
           color = "Bee Group")+
      geom_vline(xintercept = 0,  color="gray", linetype = "longdash",
                 size=.7) +
      geom_hline(yintercept = 0, color="gray", linetype = "longdash",
                 size=.7) +
      theme(panel.grid=element_blank(),
            legend.text = element_text(face = "italic"),
            text = element_text(size = 30))+
      guides(colour = guide_legend(override.aes = list(size=8)))
 
 
 # below should increase legend circle size
# + guides(colour = guide_legend(override.aes = list(size=10)))
 
hq12<-hq_ord_gg + 
    #add ellipses
  geom_polygon(data = hq_ellipse_plot$df_ellipse,
            aes(x=x, y=y, 
                fill=Group,
                color=Group),
            show.legend = FALSE,
            size=1,
            alpha=.3)+
  #add in points
  geom_point(data=hq_nmds_plot$df_ord, 
             aes(x=x, y=y, 
                 colour=Group),
             size=3)+
  geom_point(data=hq_nmds_plot$df_ord,
             aes(x=x, y=y),
             colour="black",
             size=3, 
             stroke=1,
             shape=1)+

  #add arrows
  geom_segment(data=hq_arrows, 
                   aes(x=0, y=0, 
                       xend=x, 
                       yend=y),
                   arrow=arrow(length = unit(0.03, "npc")),
                   lwd=.9,
               color = "black") +
  #add arrow labels
  geom_text_repel(data=hq_arrows,
                aes(x=x*1.2, 
                    y=y*1.2, 
                    label=var), 
      #          nudge_x = 0.06, 
     #           nudge_y =-0.05,
                color = "white", 
                size = 7,
                fontface = "bold",
                bg.colour = "black",
                bg.r = .2,
                force = .7)+
  # adjust colors for bees
  scale_color_manual(values=c("#FFB91E",
                              "#87D4FF",
                              "#28885B",
                              "#4531DE",
                              "#D55E00",
                              "#690C37"))+
  scale_fill_manual(values=c("#FFB91E",
                             "#87D4FF",
                             "#28885B",
                             "#4531DE",
                             "#D55E00",
                             "#690C37"))
                                         
       
# now do it again for 1 and 3 axis


choices_hq = c(1,3) # first select which axis to plot

hq_nmds_plot<- gg_envfit(ord = nmds_hq,
          env =nmds_hq_sparse,
          groups = nmds_hq_bee,
          scaling = 1,
          choices = choices_hq,
          perm = 999,
          alpha = 1, 
          angle = 20, 
          len = 0.5, 
          arrow.col = "black", 
          pt.size = 3,
          plot = TRUE)

hq_arrows<- hq_nmds_plot$df_arrows %>%
  filter(var %in% c("Thiamethoxam",
                 "Difenoconazole",
                 "Indoxacarb",
                 "Carbaryl",
                 "Phosmet"))

# then create the same plot with ellipses
hq_ellipse_plot<- gg_ordiplot(ord = nmds_hq,
            groups = nmds_hq_bee,
            scaling = 1,
            choices = choices_hq,
            kind = "sd",
            conf = NULL,
            show.groups = "all",
            ellipse = TRUE,
            label = FALSE,
            hull = FALSE,
            spiders = FALSE,
            plot = FALSE)

#and now put those together


# then make a base plot
 hq_ord_gg <- ggplot() + theme_bw(16) + 
      labs(x="NMDS Axis 1", 
           y="NMDS Axis 3",
           color = "Bee Group")+
      geom_vline(xintercept = 0, lty=3, color="darkgrey") +
      geom_hline(yintercept = 0, lty=3, color="darkgrey") +
      theme(panel.grid=element_blank(),
            legend.text = element_text(face = "italic"),
            text = element_text(size = 30)
            )
 
 
hq13<-hq_ord_gg + 
  #add ellipses
  geom_polygon(data = hq_ellipse_plot$df_ellipse,
            aes(x=x, y=y, 
                fill=Group,
                color=Group),
            show.legend = FALSE,
            size=1,
            alpha=.3)+
  #add in points
  geom_point(data=hq_nmds_plot$df_ord, 
             aes(x=x, y=y, 
                 colour=Group),
             size=3)+
  geom_point(data=hq_nmds_plot$df_ord,
             aes(x=x, y=y),
             colour="black",
             size=3, 
             stroke=1,
             shape=1)+

  #add arrows
  geom_segment(data=hq_arrows, 
                   aes(x=0, y=0, 
                       xend=x, 
                       yend=y),
                   arrow=arrow(length = unit(0.03, "npc")),
                   lwd=.9,
               color = "black") +
  #add arrow labels
  geom_text_repel(data=hq_arrows,
                aes(x=x*1.2, 
                    y=y*1.2, 
                    label=var), 
      #          nudge_x = 0.06, 
     #           nudge_y =-0.05,
                color = "white", 
                size = 7,
                fontface = "bold",
                bg.colour = "black",
                bg.r = .2,
                force = .7)+
  # adjust colors for bees
  scale_color_manual(values=c("#FFB91E",
                              "#87D4FF",
                              "#28885B",
                              "#4531DE",
                              "#D55E00",
                              "#690C37"))+
  scale_fill_manual(values=c("#FFB91E",
                             "#87D4FF",
                             "#28885B",
                             "#4531DE",
                             "#D55E00",
                             "#690C37"))
                                         
       


## combine hq nmds and save ####

                                  


 nmds_hq_combined<- ggarrange(hq12, hq13, 
                              common.legend = TRUE, 
                              legend="bottom")
nmds_hq_combined
ggsave(plot =nmds_hq_combined, file = "output/graphs/nmds_hq_adj_combo.png",
       height  = 9, width = 15,
       bg = "white")                 





## exposure NMDS ####
### data prep exposure ####
nmds_df<- df_all %>%
  arrange(Sample.ID)%>% 
  filter(!total == 0, # get rid of bees with no pesticides as it makes nmds sad
         !Sample.type %in% c("Apple flower","Dandelion flower") 
         )%>%
  select(c("Sample.type","Sample.ID", "Site", "mass_g", "total", where( ~ is.numeric(.x) && sum(.x) != 0))) %>% #drop pesticides that arnt in any samples
  left_join(nesting_type, by = c("Sample.type" = "bee"))%>%# add information about nesting type
  droplevels() 


nmds_df_bee <- (nmds_df[,"Sample.type"]) # make vector of bees for later grouping
nmds_df_nesting <- (nmds_df[,"nesting"]) # make vector of nesting type for later grouping


nmds_df_sparse<- nmds_df%>% 
  select(-c(Site, 
            Sample.type,
            total,
            mass_g,
            nesting
            ))%>%
 column_to_rownames('Sample.ID') # get rid of descriptors and make Sample ID the row names


### create nmds ####
nmds_exp <- metaMDS(nmds_df_sparse,  
                        distance = "bray",       
                        k = 3,
                        maxit = 999,
                        trymax = 500,
                        wascores = TRUE,
                        autotransform = FALSE)


stressplot(nmds_exp)

#color scheme for BEES ONLY in presence absence
beecolor= c("#00BA38", "#B79F00", "#00BFC4",   "#619CFF", "#F8766D", "#F564E3") # create colors for bees
nmds_df$Sample.type <- as.factor(nmds_df$Sample.type)
beecolor_loop = c("#00BA38", "#B79F00", "#00BFC4",   "#619CFF", "#F8766D", "#F564E3")[nmds_df$Sample.type]



## ggplot nmds ####

# create plot for axis 1 and 2
choices_exp = c(1,2) # first select which axis to plot

exp_nmds_plot<- gg_envfit(ord = nmds_exp,
          env = nmds_df_sparse,
          groups = nmds_df_bee,
          scaling = 1,
          choices = choices_exp,
          perm = 999,
          alpha = 1, 
          angle = 20, 
          len = 0.5, 
          arrow.col = "black", 
          pt.size = 3,
          plot = TRUE)
# after making arrows, subset to 5 main drivers from simper

exp_arrows<- exp_nmds_plot$df_arrows %>%
  filter(var %in% c("Difenoconazole",
                 "Fluxapyroxad",
                 "Thiamethoxam",
                 "Trifloxystrobin",
                 "Pyraclostrobin"))


# then create the same plot with ellipses
exp_ellipse_plot<- gg_ordiplot(ord = nmds_exp,
            groups = nmds_df_bee,
            scaling = 1,
            choices = choices_exp,
            kind = "sd",
            conf = NULL,
            show.groups = "all",
            ellipse = TRUE,
            label = FALSE,
            hull = FALSE,
            spiders = FALSE,
            plot = FALSE)

#and now put those together

# first make a base plot
 ord_gg <- ggplot() + theme_bw(16) + 
      labs(x="NMDS Axis 1", 
           y="NMDS Axis 2",
           color="Bee Group") +
      geom_vline(xintercept = 0,  color="gray", linetype = "longdash",
                 size=.7) +
      geom_hline(yintercept = 0, color="gray", linetype = "longdash",
                 size=.7) +
      theme(panel.grid=element_blank(),
            legend.text = element_text(face = "italic"),
            text = element_text(size = 30))+
      guides(colour = guide_legend(override.aes = list(size=8)))
      
     
         
exp12<-ord_gg + 
    #add ellipses
  geom_polygon(data = exp_ellipse_plot$df_ellipse,
            aes(x=x, y=y, 
                fill=Group,
                color=Group),
            show.legend = FALSE,
            size=.75,
            alpha=.3,
            color="black")+
  #add in points
  geom_point(data=exp_nmds_plot$df_ord, 
             aes(x=x, y=y, 
                 colour=Group),
             size=3)+
  geom_point(data=exp_nmds_plot$df_ord,
             aes(x=x, y=y),
             colour="black",
             size=3, 
             stroke=1,
             shape=1)+
  #add arrows
  geom_segment(data=exp_arrows, 
                   aes(x=0, y=0, 
                       xend=x, 
                       yend=y),
                   arrow=arrow(length = unit(0.03, "npc")),
                   lwd=.9,
               color = "black") +
  #add arrow labels
  geom_text_repel(data=exp_arrows,
                aes(x=x*1.2, 
                    y=y*1.2, 
                    label=var), 
   #             nudge_x = 0.06, 
    #            nudge_y =-0.05,
                color = "white", 
                size = 7,
                fontface = "bold",
                bg.colour = "black",
                bg.r = .2,
                force = .7)+
  # adjust colors for bees
  scale_color_manual(values=c("#FFB91E",
                              "#87D4FF",
                              "#28885B",
                              "#4531DE",
                              "#D55E00",
                              "#690C37"))+
  scale_fill_manual(values=c("#FFB91E",
                             "#87D4FF",
                             "#28885B",
                             "#4531DE",
                             "#D55E00",
                             "#690C37"))

#now repeat for nmds 1 and 3
choices_exp = c(1,3) # first select which axis to plot

exp_nmds_plot<- gg_envfit(ord = nmds_exp,
          env =nmds_df_sparse,
          groups = nmds_df_bee,
          scaling = 1,
          choices = choices_exp,
          perm = 999,
          alpha = 1, 
          angle = 20, 
          len = 0.5, 
          arrow.col = "black", 
          pt.size = 3,
          plot = TRUE)

# remake arrows for new axis
exp_arrows<- exp_nmds_plot$df_arrows %>%
  filter(var %in% c("Difenoconazole",
                 "Fluxapyroxad",
                 "Thiamethoxam",
                 "Trifloxystrobin",
                 "Pyraclostrobin"))


# then create the same plot with ellipses
exp_ellipse_plot<- gg_ordiplot(ord = nmds_exp,
            groups = nmds_df_bee,
            scaling = 1,
            choices = choices_exp,
            kind = "sd",
            conf = NULL,
            show.groups = "all",
            ellipse = TRUE,
            label = FALSE,
            hull = FALSE,
            spiders = FALSE,
            plot = FALSE)

#and now put those together

# first make a base plot
 ord_gg <- ggplot() + theme_bw(16) + 
      labs(x="NMDS Axis 1", 
           y="NMDS Axis 3",
           color="Bee Group") +
      geom_vline(xintercept = 0,  color="gray", linetype = "longdash",
                 size=.7) +
      geom_hline(yintercept = 0, color="gray", linetype = "longdash",
                 size=.7) +
      theme(panel.grid=element_blank(),
            legend.text = element_text(face = "italic"),
            text = element_text(size = 30))+
      guides(colour = guide_legend(override.aes = list(size=8)))

exp13<-ord_gg + 
   #add ellipses
  geom_polygon(data = exp_ellipse_plot$df_ellipse,
            aes(x=x, y=y, 
                fill=Group,
                color=Group),
            show.legend = FALSE,
            size=.75,
            alpha=.3,
            color="black")+
  #add in points
  geom_point(data=exp_nmds_plot$df_ord, 
             aes(x=x, y=y, 
                 colour=Group),
             size=3)+
  geom_point(data=exp_nmds_plot$df_ord,
             aes(x=x, y=y),
             colour="black",
             size=3, 
             stroke=1,
             shape=1)+
  #add arrows
  geom_segment(data=exp_arrows, 
                   aes(x=0, y=0, 
                       xend=x, 
                       yend=y),
                   arrow=arrow(length = unit(0.03, "npc")),
                   lwd=.9,
               color = "black") +
  #add arrow labels
  geom_text_repel(data=exp_arrows,
                aes(x=x*1.2, 
                    y=y*1.2, 
                    label=var), 
   #             nudge_x = 0.06, 
    #            nudge_y =-0.05,
                color = "white", 
                size = 7,
                fontface = "bold",
                bg.colour = "black",
                bg.r = .2,
                force = .7)+
  # adjust colors for bees
  scale_color_manual(values=c("#FFB91E",
                              "#87D4FF",
                              "#28885B",
                              "#4531DE",
                              "#D55E00",
                              "#690C37"))+
  scale_fill_manual(values=c("#FFB91E",
                             "#87D4FF",
                             "#28885B",
                             "#4531DE",
                             "#D55E00",
                             "#690C37"))
                                 

## combine exp nmds and save ####
 
nmds_exp_combined<- ggarrange(exp12, exp13, 
                              common.legend = TRUE, 
                              legend="bottom")
nmds_exp_combined

 ggsave(plot =nmds_exp_combined, file = "output/graphs/nmds_exp_combo.png", height = 9, width = 15,
       bg = "white")                 


 

## flower NMDS ####
 
nmds_df_flowers<- df_all %>%
  mutate(Sample.type=fct_relevel(Sample.type,"Dandelion flower", after = Inf))%>%
  mutate(Sample.type=fct_relevel(Sample.type,"Apple flower", after = Inf))%>%
  arrange(Sample.ID)%>% 
  filter(!total == 0, # get rid of bees with no pesticides as it makes nmds sad
         )%>%
  select(c("Sample.type","Sample.ID", "Site", "mass_g", "total", where( ~ is.numeric(.x) && sum(.x) != 0))) %>% #drop pesticides that arnt in any samples
  left_join(nesting_type, by = c("Sample.type" = "bee"))%>%# add information about nesting type
  droplevels() 


nmds_df_flower_bee <- (nmds_df_flowers[,"Sample.type"]) # make vector of bees for later grouping
nmds_df__flower_nesting <- (nmds_df_flowers[,"nesting"]) # make vector of nesting type for later grouping


nmds_df_sparse_flower<- nmds_df_flowers%>% 
  select(-c(Site, 
            Sample.type,
            total,
            mass_g,
            nesting
            ))%>%
 column_to_rownames('Sample.ID') # get rid of descriptors and make Sample ID the row names


### create nmds ####
nmds_exp_flower <- metaMDS(nmds_df_sparse_flower,  
                        distance = "bray",       
                        k = 3,
                        maxit = 999,
                        trymax = 500,
                        wascores = TRUE,
                        autotransform = FALSE)


stressplot(nmds_exp_flower)



## ggplot flower nmds ####



# create plot for axis 1 and 2
choices_flowers = c(1,2) # first select which axis to plot

flower_nmds_plot<- gg_envfit(ord = nmds_exp_flower,
          env =nmds_df_sparse_flower,
          groups = nmds_df_flower_bee,
          scaling = 1,
          choices = choices_flowers,
          perm = 999,
          alpha = 1, 
          angle = 20, 
          len = 0.5, 
          arrow.col = "black", 
          pt.size = 3,
          plot = TRUE)

flower_exp_arrows<- flower_nmds_plot$df_arrows %>%
  filter(var %in% c("Difenoconazole",
                 "Fluxapyroxad",
                 "Thiamethoxam",
                 "Trifloxystrobin",
                 "Cyprodinil"))

# then create the same plot with ellipses
flower_ellipse_plot<- gg_ordiplot(ord = nmds_exp_flower,
            groups = nmds_df_flower_bee,
            scaling = 1,
            choices = choices_flowers,
            kind = "sd",
            conf = NULL,
            show.groups = "all",
            ellipse = TRUE,
            label = FALSE,
            hull = FALSE,
            spiders = FALSE,
            plot = FALSE)

#and now put those together

# first make a base plot
 ord_gg <- ggplot() + theme_bw(16) + 
      labs(x="NMDS Axis 1", 
           y="NMDS Axis 2",
           color="Bee Group") +
      ggtitle("Exposure")+
      geom_vline(xintercept = 0,  color="gray", linetype = "longdash",
                 size=.7) +
      geom_hline(yintercept = 0, color="gray", linetype = "longdash",
                 size=.7) +
      theme(panel.grid=element_blank(),
            legend.text = element_text(face = "italic"),
            text = element_text(size = 30),
            legend.title=element_blank())+
      guides(colour = guide_legend(override.aes = list(size=8)))

flower12<-ord_gg + 
    #add ellipses 
  #split by apples and bees so different aesthetics
  geom_polygon(data = flower_ellipse_plot$df_ellipse %>%
                         filter(!Group %in% c("Apple flower", "Dandelion flower")),
            aes(x=x, y=y, 
                fill=Group,
                color=Group),
            show.legend = FALSE,
            size=1,
            alpha=.3)+
  geom_polygon_pattern(data = flower_ellipse_plot$df_ellipse %>% # add pattern just for flowers
                         filter(Group %in% c("Apple flower", "Dandelion flower"))%>%
                         rename(`Flower Type` = Group), #rename it to flowers for legend
            aes(x=x, y=y, 
                pattern = `Flower Type`,
                pattern_angle = `Flower Type`),
            size=.75,
            alpha=.3,
            color="black")+
  #add in points for bees
  geom_point(data=flower_nmds_plot$df_ord%>%
                         filter(!Group %in% c("Apple flower", "Dandelion flower")) , 
             aes(x=x, y=y, 
                 colour=Group),
             size=3)+
  geom_point(data=flower_nmds_plot$df_ord%>%
                         filter(!Group %in% c("Apple flower", "Dandelion flower")),
             aes(x=x, y=y),
             colour="black",
             size=3, 
             stroke=1,
             shape=1)+
  #add arrows
  geom_segment(data=flower_exp_arrows, 
                   aes(x=0, y=0, 
                       xend=x, 
                       yend=y),
                   arrow=arrow(length = unit(0.03, "npc")),
                   lwd=.9,
               color = "black") +
  #add arrow labels
  geom_text_repel(data=flower_exp_arrows,
                aes(x=x*1.2, 
                    y=y*1.2, 
                    label=var), 
   #             nudge_x = 0.06, 
    #            nudge_y =-0.05,
                color = "white", 
                size = 7,
                fontface = "bold",
                bg.colour = "black",
                bg.r = .2,
                force = .1)+
  # adjust colors for bees
  scale_color_manual(values=c("#FFB91E",
                              "#87D4FF",
                              "#28885B",
                              "#4531DE",
                              "#D55E00",
                              "#690C37"))+
  scale_fill_manual(values=c("#FFB91E",
                             "#87D4FF",
                             "#28885B",
                             "#4531DE",
                             "#D55E00",
                             "#690C37"))



#now repeat for nmds 1 and 3
choices_flowers = c(1,3) # first select which axis to plot

flower_nmds_plot<- gg_envfit(ord = nmds_exp_flower,
          env =nmds_df_sparse_flower,
          groups = nmds_df_flower_bee,
          scaling = 1,
          choices = choices_flowers,
          perm = 999,
          alpha = 1, 
          angle = 20, 
          len = 0.5, 
          arrow.col = "black", 
          pt.size = 3,
          plot = TRUE)

flower_exp_arrows<- flower_nmds_plot$df_arrows %>%
  filter(var %in% c("Difenoconazole",
                 "Fluxapyroxad",
                 "Thiamethoxam",
                 "Trifloxystrobin",
                 "Cyprodinil"))

# then create the same plot with ellipses
flower_ellipse_plot<- gg_ordiplot(ord = nmds_exp_flower,
            groups = nmds_df_flower_bee,
            scaling = 1,
            choices = choices_flowers,
            kind = "sd",
            conf = NULL,
            show.groups = "all",
            ellipse = TRUE,
            label = FALSE,
            hull = FALSE,
            spiders = FALSE,
            plot = FALSE)

#and now put those together

# first make a base plot
 ord_gg <- ggplot() + theme_bw(16) + 
      labs(x="NMDS Axis 1", 
           y="NMDS Axis 3",
           color="Bee Group") +
      ggtitle("Exposure")+
      geom_vline(xintercept = 0,  color="gray", linetype = "longdash",
                 size=.7) +
      geom_hline(yintercept = 0, color="gray", linetype = "longdash",
                 size=.7) +
      theme(panel.grid=element_blank(),
            legend.text = element_text(face = "italic"),
            text = element_text(size = 30),
            legend.title=element_blank())+
      guides(colour = guide_legend(override.aes = list(size=8)))

flower13<-ord_gg + 
    #add ellipses 
  #split by apples and bees so different aesthetics
  geom_polygon(data = flower_ellipse_plot$df_ellipse %>%
                         filter(!Group %in% c("Apple flower", "Dandelion flower")),
            aes(x=x, y=y, 
                fill=Group,
                color=Group),
            show.legend = FALSE,
            size=1,
            alpha=.3)+
  geom_polygon_pattern(data = flower_ellipse_plot$df_ellipse %>% # add pattern just for flowers
                         filter(Group %in% c("Apple flower", "Dandelion flower"))%>%
                         rename(`Flower Type` = Group), #rename it to flowers for legend
            aes(x=x, y=y, 
                pattern = `Flower Type`,
                pattern_angle = `Flower Type`),
            size=.75,
            alpha=.3,
            color="black")+
  #add in points for bees
  geom_point(data=flower_nmds_plot$df_ord%>%
                         filter(!Group %in% c("Apple flower", "Dandelion flower")) , 
             aes(x=x, y=y, 
                 colour=Group),
             size=3)+
  geom_point(data=flower_nmds_plot$df_ord%>%
                         filter(!Group %in% c("Apple flower", "Dandelion flower")),
             aes(x=x, y=y),
             colour="black",
             size=3, 
             stroke=1,
             shape=1)+
  #add arrows
  geom_segment(data=flower_exp_arrows, 
                   aes(x=0, y=0, 
                       xend=x, 
                       yend=y),
                   arrow=arrow(length = unit(0.03, "npc")),
                   lwd=.9,
               color = "black") +
  #add arrow labels
  geom_text_repel(data=flower_exp_arrows,
                aes(x=x*1.2, 
                    y=y*1.2, 
                    label=var), 
   #             nudge_x = 0.06, 
    #            nudge_y =-0.05,
                color = "white", 
                size = 7,
                fontface = "bold",
                bg.colour = "black",
                bg.r = .2,
                force = .7)+
  # adjust colors for bees
  scale_color_manual(values=c("#FFB91E",
                              "#87D4FF",
                              "#28885B",
                              "#4531DE",
                              "#D55E00",
                              "#690C37"))+
  scale_fill_manual(values=c("#FFB91E",
                             "#87D4FF",
                             "#28885B",
                             "#4531DE",
                             "#D55E00",
                             "#690C37"))

## combine flower  nmds and save ####
 
nmds_flower_combined<- ggarrange(flower12, flower13, 
                              common.legend = TRUE, 
                              legend="bottom")
nmds_flower_combined

 ggsave(plot =nmds_flower_combined, file = "output/graphs/nmds_flower_combo_pattern.jpeg", 
        height = 9, width = 16,
        bg = "white")                 



 ggsave(plot =nmds_flower_combined, file = "output/graphs/nmds_flower_combo_pattern.svg", 
        height = 9, width = 16,
        bg = "white")                 



# PERMANOVA  ####

## hq permanova ####
 
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

write_csv(beepermanhq, "output/tables/permanova_contact_bees_ADJUSTED.csv")



## exposure permanova with flowers ####
perm_df<- df_all %>%
  arrange(Sample.ID)%>%
  filter(!total == 0 )%>% 
  left_join(nesting_type, by = c("Sample.type" = "bee")) # add information about nesting type

perm_matrix <- perm_df%>% 
  select(-c(Site, 
            Sample.type,
            Sample.ID,
            total,
            nesting,
            mass_g
            ))%>%
  as.matrix()

 
perm_flower <- (perm_df[,"Sample.type"]) # make vector of bees for later grouping

bee_flower_dis <- vegdist(perm_matrix, method = 'bray')
bee_flower_perman <- adonis2(bee_flower_dis~Sample.type, data=perm_df, permutations=999, method="bray")
bee_flower_perman
write_csv(bee_flower_perman, "output/tables/permanova_exposure_all.csv")


## exposure permanova bees ####

perm_df_bee<- df_all %>%
  arrange(Sample.ID)%>%
  filter(!total == 0,  # get rid of bees with no pesticides as it makes nmds and permanova sad
         !Sample.type %in% c("Apple flower","Dandelion flower") # comment out if you want all samples. otherwise this removes flowers
    )%>% 
  left_join(nesting_type, by = c("Sample.type" = "bee")) # add information about nesting type

perm_matrix_bee <- perm_df_bee%>% 
  select(-c(Site, 
            Sample.type,
            Sample.ID,
            total,
            nesting,
            mass_g
            ))%>%
  as.matrix()

 
perm_bee <- (perm_df_bee[,"Sample.type"]) # make vector of bees for later grouping

bee_dis <- vegdist(perm_matrix_bee, method = 'bray')
beeperman <- adonis2(bee_dis~Sample.type, data=perm_df_bee, permutations=999, method="bray")
beeperman
write_csv(beeperman, "output/tables/permanova_exposure_all.csv")




## pairwise permanova ####

#exposure across all samples
pair.exposure<-pairwise.adonis(bee_flower_dis, factors = perm_flower)
write_csv(pair.exposure, "output/tables/pairwise_permanova_exposure.csv")

#exposure across bees
pair_bee_exposure<-pairwise.adonis(bee_dis, factors = perm_bee)
write_csv(pair_bee_exposure, "output/tables/pairwise_permanova_exposure_bees.csv")

# contact HQ
pair.contact<-pairwise.adonis(bee.dis.hq, factors = perm_bee_hq)
write_csv(pair.contact, "output/tables/pairwise_permanova_contact_ADJUSTED.csv")



# SIMPER exposure ####

# now run a simper model to assess how different pesticides
# are contributing to dissimilarity between species 


## pairwiseexposure with flowers ####
sim_flower<- simper(nmds_df_sparse_flower, group = nmds_df_flower_bee, permutations = 999)
## S3 method for class 'simper'
simper_flower_output <- summary(sim_flower, ordered = TRUE,
    digits = max(3,getOption("digits") - 3))
capture.output(summary(sim_flower), file = "output/tables/sim_exp_withflowers.txt")



## all exposure with flowers ####
sim_flower_all<- simper(nmds_df_sparse_flower, permutations = 999)
## S3 method for class 'simper'
simper_flower_all_output <- summary(sim_flower_all, ordered = TRUE,
    digits = max(3,getOption("digits") - 3))
capture.output(summary(sim_flower_all), file = "output/tables/sim_exp_withflowers_all.txt")


## pairwise exp bee####
sim<- simper(nmds_df_sparse, group = nmds_df_bee, permutations = 999)
## S3 method for class 'simper'
simperoutput <- summary(sim, ordered = TRUE,
    digits = max(3,getOption("digits") - 3))
capture.output(summary(sim), file = "output/tables/sim_exp_bees.txt")



# ungrouped exposure (bee) simper
sim_all<- simper(nmds_df_sparse, permutations = 999)

simper_all_output <- summary(sim_all, ordered = TRUE,
    digits = max(3,getOption("digits") - 3))
capture.output(summary(sim_all), file = "output/tables/sim_exp_all.txt")



# SIMPER hq ####

#pairwise simper
sim_hq<- simper(nmds_hq_sparse, group = nmds_hq_bee, permutations = 999)

## S3 method for class 'simper'
simperoutput_hq <- summary(sim_hq, ordered = TRUE,
    digits = max(3,getOption("digits") - 3))
capture.output(summary(sim_hq), file = "output/tables/sim_hq_ADJUSTED.txt")


# ungrouped simper
sim_hq_all<- simper(nmds_hq_sparse, permutations = 999)

simper_hq_all_output <- summary(sim_hq_all, ordered = TRUE,
    digits = max(3,getOption("digits") - 3))
capture.output(summary(sim_hq_all), file = "output/tables/sim_hq_all_ADJUSTED.txt")