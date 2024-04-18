# script for adjusting ld50 values by weight to other bees
# based on tobias pammingers paper
# Extrapolating Acute Contact Bee Sensitivity to Insecticides Based on Body Weight Using a Phylogenetically Informed Interspecies Scaling Framework
# https://doi.org/10.1002/etc.5045

# load packages ####
library(tidyverse)
library(broom)
library(psych) # for calculating geometric mean 
library(ggrepel)

# read in data from publication ####
rm(list = ls()) # clean environment

df <-read.csv("ld50_adjust/2020_BeeTox_database_acute_contact_publication_final_R1.csv")%>%
  mutate_if(is.character, factor) # make all characters factors


##  averaging within paper, then within chemical group, then within species ####
species_avg_df<-df %>%  
  group_by(tox_ref, chemical_group, Scientific.Name) %>%  
  mutate(paper_avg = geometric.mean(LD50_bee))%>% 
  distinct(tox_ref, .keep_all = TRUE) %>% 
  ungroup()%>%
  group_by(chemical_group, Scientific.Name,) %>% 
  mutate(overall_avg = geometric.mean(paper_avg))%>% 
  distinct(Scientific.Name, .keep_all = TRUE)%>%
  group_by(Scientific.Name,) %>% 
  mutate(species_avg = geometric.mean(overall_avg))%>% 
  distinct(Scientific.Name, .keep_all = TRUE)%>%
  mutate(log_weight = log10(weight), log_ld50 = log10(species_avg)) # then long 10 transform weight and ld50


## linear model ####

#create lm that is forced through A. mellifera
regression <- lm((log_ld50 - -0.73931903) ~ I(log_weight - 2) - 1, 
              data = species_avg_df)

print(summary(regression), digits=20)

# then calculate out the intercept
intercept <- with(species_avg_df,
    -0.73931903 - 2 * coef(regression)[1])
intercept


## graph ####
adjust_plot <- species_avg_df %>%
  ggplot(aes(x= log_weight, y=log_ld50))+
#  geom_abline(slope= 0  , intercept = -0.73931903, size = 1.3, color = "blue3")+
#  geom_abline(slope= 1.34036  , intercept = -3.420038, size= 1.3, color = "")+
    geom_abline(slope= 0  , intercept = -0.73931903, size= 1.3, color = "blue")+
    geom_abline(slope = 1.219, intercept = -3.176504, size = 1.3, color = "orange")+
  geom_point(aes(fill = Scientific.Name), color= "black", shape = 21, size = 2.8) + 
  xlab("Log(bee weight) (mg)")+
  ylab("Log(LD50)")+
  geom_text_repel(aes(label= Scientific.Name),  # add point labels
                  fontface = "italic",
                  box.padding = 0.4,
                  force = 1.5,
                  force_pull = .5)+
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(size = 20))

adjust_plot

ggsave(plot = adjust_plot, file = "output/graphs/adjustplot.jpeg", width  =10, height =8)
