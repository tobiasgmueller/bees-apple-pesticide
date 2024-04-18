# script for trying to adjust ld50 values by weight to other bees


# based on tobias pammingers paper
# Extrapolating Acute Contact Bee Sensitivity to Insecticides Based on Body Weight Using a Phylogenetically Informed Interspecies Scaling Framework


# load packages ####
library(tidyverse)
library(broom)
library(psych) # for calculating geometric mean 
library(ggrepel)

# read in data from publication ####
rm(list = ls()) # clean environment

df <-read.csv("ld50_adjust/2020_BeeTox_database_acute_contact_publication_final_R1.csv")

df <- df %>% mutate_if(is.character, factor) # make all characters factors


#  create a number of publications for each species
# this is  used as a weight in some models
df<- df %>%
  group_by(Scientific.Name, chemical_group) %>%
  mutate(pub_num_group = n_distinct(tox_ref)) %>% 
  ungroup() %>% 
  group_by(Scientific.Name, substance) %>%
  mutate(pub_num_compound = n_distinct(tox_ref)) %>% 
  ungroup()



# the code will be structured into 4 sections each using a different way of averaging data

# each section will be structured identically
# containing 1) data prep (averaging), 2) lm creation, 3) graph of that data with lm plotted

# there is a bonus section facetted by moa




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


### linear model ####
# this would be simple lm
species_avg_lm<- species_avg_df %>%
  lm(log_ld50 ~ log_weight, data=.)

species_avg_df$residuals <- species_avg_lm$residuals

print(summary(species_avg_lm), digits=20)



#create lm that is forced through A. mellifera
regression <- lm((log_ld50 - -0.73931903) ~ I(log_weight - 2) - 1, 
              data = species_avg_df)

print(summary(regression), digits=20)

# then calculate out the intercept
intercept <- with(species_avg_df,
    -0.73931903 - 2 * coef(regression)[1])
intercept





#------------#

### graph ####
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






# below is not used analyses ####

























## avg within paper, then within chem group #### ----------------
avg_df<-df %>%  
  group_by(tox_ref, chemical_group, Scientific.Name) %>%  
  mutate(paper_avg = geometric.mean(LD50_bee))%>% 
  distinct(tox_ref, .keep_all = TRUE) %>% 
  ungroup()%>%
  group_by(chemical_group, Scientific.Name,) %>% 
  mutate(overall_avg = geometric.mean(paper_avg))%>% 
  distinct(Scientific.Name, .keep_all = TRUE)%>%
  mutate(log_weight = log10(weight), log_ld50 = log10(overall_avg)) # then long 10 transform weight and ld50


# abbreviate species names for easy to read labeled points
avg_df$name <- abbreviate(avg_df$Scientific.Name, minlength = 4,
                          method = "both")


### linear model ####
lm_moa_avg<-lm(log_ld50 ~ log_weight, weights = pub_num_group, data = avg_df)
print(summary(lm_moa_avg), digits=20)




### graph ####
avg_df %>%
  ggplot(aes(x= log_weight, y=log_ld50, label=Scientific.Name
             ))+
  geom_point(aes(color=Scientific.Name)) + 
#geom_text()+  # this labels points. remove for cleaner graph
  geom_abline(slope=1, intercept = -2.964, size= 1) +
  geom_abline(slope= 1.106, intercept = -2.964, size= 1, color = "blue")+ 
  xlab("Log(bee weight) (g)")+
  ylab("Log(LD50)")+
  labs(color="Bee species")

# annotate("text", x = .8, y = 1.5, label = "Adjusted R2 = 0.16")











##  averaging within paper, then within chemical compound #### -------------
compound_avg_df<-df %>%  
  group_by(tox_ref, chemical_group, Scientific.Name)%>% 
  mutate(paper_avg = geometric.mean(LD50_bee))%>%
  distinct(tox_ref, .keep_all = TRUE) %>%
  ungroup() %>%
  group_by(substance, Scientific.Name) %>% 
  mutate(overall_avg = geometric.mean(paper_avg))%>% 
  distinct(Scientific.Name, .keep_all = TRUE)%>%
  mutate(log_weight = log10(weight), log_ld50 = log10(overall_avg))


### linear model ####
lm_comp<- compound_avg_df %>%
  lm(log_ld50 ~ log_weight, weights = pub_num_compound, data=.)
print(summary(lm_comp), digits=20)



### graph ####
compound_avg_df %>%
  ggplot(aes(x= log_weight, y=log_ld50, label=Scientific.Name
             ))+
  geom_point(aes(color=Scientific.Name)) + 
  geom_abline(slope= 1  , intercept = -2.936, size =1) +
  geom_abline(slope= 1.0888 , intercept = -2.936, size =1, color = "blue")+ 
  xlab("Log(bee weight) (g)")+
  ylab("Log(LD50)")+
  labs(color="Bee species")




## no averaging at all (plotting all studies) ####  --------------
df <- df %>%
  mutate(log_weight = log10(weight), log_ld50 = log10(LD50_bee))

### linear model ####

lm_no_avg<- df %>%
  lm(log_ld50 ~ log_weight, data=.)
print(summary(lm_no_avg), digits=20)



### graph ####
df %>%
  ggplot(aes(x= log_weight, y=log_ld50, label=Scientific.Name
             ))+
  geom_point(aes(color=Scientific.Name)) + 
  geom_abline(slope= 1  , intercept = -3.15427, size=1)+
  geom_abline(slope= 1.18856  , intercept = -3.15427, size=1, color="blue")+ 
  xlab("Log(bee weight) (g)")+
  ylab("Log(LD50)")+
  labs(color="Bee species")




















# lm split by group ####

# this corresponds to the facetted plot below
all_regress <-  avg_df %>% group_by(chemical_group) %>%
  do(mod1 = lm(log_ld50 ~ log_weight, weight = pub_num, data = .))%>% ungroup()


# use broom the extract the slope and rsq per group
glance <-all_regress %>% mutate(tidy = map(mod1, broom::tidy),
                                   glance = map(mod1, broom::glance),
                                   augment = map(mod1, broom::augment),
                                   rsq = glance %>% map_dbl('r.squared'),
                                   slope = tidy %>% map_dbl(function(x)
                                     x$estimate[2]))

glance




## facetted by ai moa ####
avg_df %>%
  ggplot(aes(x= log_weight, y=log_ld50, label=name
             ))+
  geom_point() + 
geom_text()+   # this labels points. remove for cleaner graph
facet_wrap(~chemical_group, scales = "free")+
  geom_smooth(method = "lm", mapping = aes(weight = pub_num_group))
  





comp_avg %>%
  ggplot(aes(x= log_weight, y=log_ld50
             ))+
  geom_point(aes(color=Scientific.Name)) + 
  geom_smooth(method = "lm", mapping = aes(weight = pub_num))










# graveyard ####


# old r2 calculation ####


#----------- #

# calculating R2 for fitting lm by hand
# residuals from your prediction
RSS = sum((species_avg_lm$residuals)^2)

RSS = sum((species_avg_lm$fitted.values - species_avg_lm$model$log_ld50)^2)

# total residuals from data
TSS = sum((species_avg_lm$model$log_ld50 - mean(species_avg_lm$model$log_ld50))^2)    

1-RSS/TSS

### linear model ####
species_avg_lm<- species_avg_df %>%
  lm(log_ld50 ~ log_weight, data=.)

species_avg_df$residuals <- species_avg_lm$residuals

print(summary(species_avg_lm), digits=20)



# below is the formula for the unadjusted 1:1 line if plotting on log log scale
unadjusted<-(lm(log_ld50-0.037119 ~ 1 + offset(log_weight), data = species_avg_df))



# calculating R2 for fitting lm by hand
# residuals from your prediction
RSS = sum((species_avg_lm$residuals)^2)

RSS = sum((species_avg_lm$fitted.values - species_avg_lm$model$log_ld50)^2)
# total residuals from data
TSS = sum((species_avg_lm$model$log_ld50 - mean(species_avg_lm$model$log_ld50))^2)    

1-RSS/TSS

####
# calculating R2 for 1:1 lm by hand

# residuals from your prediction
RSS_test = sum((unadjusted$residuals)^2)

RSS_test = sum((unadjusted$fitted.values - (unadjusted$model$`log_ld50 - 0.037119`))^2)

# total residuals from data
TSS_test = sum(((unadjusted$model$`log_ld50 - 0.037119`) - mean((unadjusted$model$`log_ld50 - 0.037119`)))^2)    

1-RSS_test/TSS_test







# calculate residuals
species_avg_df %>%
  filter(Scientific.Name == "Apis mellifera") %>%
  select(residuals)


# double checking a few discrepencies

#their ld50 value for A. cerana organochlorine is ~0.01 and im getting 0.09
# they have also removed osmia cornifrons from neonicotinoids (nitro)
# other things are wrong.. hmm...


#doble checking the LD50 value for A. cerana organochlorine
# im getting a different value than is plotted in manuscript... 
# manuscript is getting ~0.01
# im getting 0.09

x <-c(1.58000,1.40000,1.60000,1.23000,1.15000,1.34000)
y<- geometric.mean(x)

xprime<- c(0.00250,0.01600)
yprime<- geometric.mean(xtwo)

z<- c(y,ytwo)

geometric.mean(z)




## old graphs ####

## all data graph ####

# this graph would correspond to the first lm of geometric means
avg_df %>%
  ggplot(aes(x= log_weight, y=log_ld50, label=Scientific.Name
             ))+
  geom_point(aes(color=Scientific.Name)) + 
#geom_text()+  # this labels points. remove for cleaner graph
  geom_smooth(method = "lm", mapping = aes(weight = pub_num))+
  geom_abline(slope=1, intercept = -2.964)
  

# lets also look at non averaged numbers
df %>%
  ggplot(aes(x= log_weight, y=log_ld50, label=Scientific.Name
             ))+
  geom_point(aes(color=Scientific.Name)) + 
#geom_text()+  # this labels points. remove for cleaner graph
  geom_smooth(method = "lm", mapping = aes(weight = pub_num))+
  geom_abline(slope= 1  , intercept = -3.154)


# then averaged by bee

species_avg_df %>%
  ggplot(aes(x= log_weight, y=log_ld50, label=name
             ))+
  geom_point(aes(color=Scientific.Name)) + 
geom_text()+  # this labels points. remove for cleaner graph
  geom_smooth(method = "lm", mapping = aes(weight = pub_num))+
  geom_abline(slope= 1  , intercept = -3.324)



