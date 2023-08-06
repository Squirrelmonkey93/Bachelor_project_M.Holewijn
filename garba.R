library(crestr)
library(readxl)

data_path <- "Datasets/Garba Guracha/"

# Loading in complete original data set
Garba_org <- read.csv(paste0(data_path, "GarbaGuracha_site.csv"), header=T, sep = ",")

# Loading in cleaned up data sets
# All terrestrial taxa
Garba <- read.csv(paste0(data_path, "GarbaGuracha_100%.csv"), header=T, sep = ";")

# Terrestrial taxa that are at least 1% per sample
Garba_1p <- read.csv(paste0(data_path, "GarbaGuracha_1%.csv"), header=T, sep = ";")

# Terrestrial taxa that are at least 0,5% per sample
Garba_05p <- read.csv(paste0(data_path, "GarbaGuracha_05%.csv"), header=T, sep = ";")


# Creating empty pse taxa list
taxa_list_Garba <- colnames(Garba[,-1])
createPSE(taxa_list_Garba, loc = "pse_Garba.xlsx")

# Loading in pse list
Garba_pse <- read_excel(paste0(data_path, "pse_Garba.xlsx"), sheet = 1, col_names = TRUE)
Garba_05p_pse <- read_excel(paste0(data_path, "pse_Garba.xlsx"), sheet = 2, col_names = TRUE)

# Getting modern day data for the area looked at
Garba_rec <- crest.get_modern_data(df=Garba_05p, pse=Garba_05p_pse,
                                  taxaType = 1, climate= c("bio1", "bio12", "ai"),
                                  dbname= "subset_bandAfrica.sqlite3",
                                  continents = "Africa",
                                  countries = countries,
                                  site_info = c(39.861517, 6.877983),
                                  verbose=T)

# Fitting the species and proxy PDFs
Garba_rec <- crest.calibrate(Garba_rec,
                            climateSpaceWeighting = T,
                            bin_width = c(0.5, 25, 0.25),
                            shape = c("normal", "lognormal", "normal"),
                            npoints = 2000,
                            verbose = T)

# Plotting general overview of the space
plot_climateSpace(Garba_rec,
                  climate = c('bio1', 'bio12', 'ai'), bin_width = c(0.5, 25, 0.25),
                  save = T, filename = "Garba_05p_ClimateSpace",
                  as.png = T, png.res = 500, width = 7.48, height = 6, y0 = 0.5,
                  add_modern = T)

# removing taxa in families "Amaranthaceae", "Asteraceae", "Poaceae"
Garba_rec <- excludeTaxa(Garba_rec, c("Chenopodiaceae", 
                                      "Artemisia", "Asteraceae", "Aster", "Carduus", "Senecio",
                                      "Poaceae"), c("bio1", "bio12", "ai"))
# print(Garba_rec$inputs$selectedTaxa)

# Actual reconstruction function
Garba_rec <- crest.reconstruct(Garba_rec, verbose = T)

# Plotting of annual temperature (bio1) and precipitation (bio12) of ancient data
plot(Garba_rec, climate = 'bio1', add_modern=T, 
     xlim= c(0,20200), ylim=c(12,30),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = F, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Garba_Temperature_05p.pdf",)

plot(Garba_rec, climate = 'bio12', add_modern=T, 
     xlim= c(0,20200), ylim=c(200, 2500), 
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = T, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Garba_Precip_05p.pdf",)

plot(Garba_rec, climate = 'ai', add_modern=T,
     xlim= c(0,20200), ylim=c(0,2), 
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5,
     save = T, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Garba_Aridity_05p.pdf",)

# Leave One Out analysis and plots
Garba_rec <- loo(Garba_rec)
plot_loo(Garba_rec, save = T, filename = "Garba_LOO_05p", 
         as.png=T, png.res = 500, width=7.5, height = 10, xlim= c(0,20200),
         sort='decr', col_pos = 'blue', col_neg = 'red')