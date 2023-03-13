library(crestr)

Bambili1 <- read.csv("dataset_Bambili1.csv", header=T, sep = ",")

# creating empty pse taxa list
list_of_taxa <- Bambili1[-(1:6),1]
createPSE(list_of_taxa, loc = "pse_Bambili.xlsx")

# filling in pse list manually
getTaxonomy(family = "", genus = "leonardoxa", species = "",
            taxaType = 1, depth.out = 6, dbname = "gbif4crest_02" )

# loading in pse list
Bambili1_pse <- read.csv("pse_Bambili.csv", header=T, sep = ";")

Bambili_rec <- crest.get_modern_data(df=Bambili1,
                                  pse=Bambili1_pse,
                                  taxaType = 1,
                                  climate= c("bio1"),
                                  # selectedTaxa = crest_ex_selection, # taxa used for reconstruction
                                  # dbname= "crest_example", # db to extract data from
                                  # site_info = c(7.5, 7.5), 
                                  # continents = "Africa",
                                  # countries = c("Cameroon"),
                                  verbose=F)

Bambili_rec <- crest.calibrate(Bambili_rec, # prev created crestObj
                            climateSpaceWeighting = T, # apply correction for modern climate
                            bin_width = c(2), # bin sized used for correction
                            shape = c("normal"), # shape of the species pdfs
                            verbose = F)
