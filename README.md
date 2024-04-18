# bees-apple-pesticide


Pesticide risk during commercial apple pollination is greater for honeybees than other managed and wild bees


This hosts all code for the analysis. The files are organized as below. All data (the "input" folder) is also hosted on zenodo at


input - contains the main data files.

all_pesticide2019.csv – contains the pesticide residue data for all samples

LD50.csv – contains pesticide LD50s in PPB

ld50_per_bee.csv – contains pesticide LD50s in ug / honeybee

nesting_type.csv – lists nesting types of the different bee samples

pesticide_type.csv – lists pesticides by their type

short_name.csv – lists short bee names

folder “gis”

hive_distance_matrix.csv – distances from each orchard to other sites
sitelocation.csv – contains coordinates of sites

code - contains the scripts to analyze the input files

pesticide_analysis.R – is the main anaylsis for the paper
unadjusted_pesticide_analysis.R – is a copy of most of the code above but without LD50 weight adjustements

ld50 adjust - contains the inputs and code for the ld50 adjustments I ran

ld50_adjust.R – is the code to calculate our the adjustements
2020_BeeTox_database_acute_contact_publication_final_R1.csv – is the data file taken from pamminger publication

