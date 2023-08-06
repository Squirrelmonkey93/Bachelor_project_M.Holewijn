library(crestr)
library(readxl)

data_path <- "Datasets/Bambili/"

# Loading in complete original data set
Bambili_org <- read.csv(paste0(data_path, "Bambili_site23662.csv"), header=T, sep = ",")

# Loading in cleaned up data sets
# All terrestrial taxa
Bambili <- read.csv(paste0(data_path, "Bambili_100%.csv"), header=T, sep = ";")[,-1]
rownames(Bambili) <- read.csv(paste0(data_path, "Bambili_100%.csv"), header=T, sep = ";")[,1]
Bambili[is.na(Bambili)] <- 0

# Terrestrial taxa that are at least 1% per sample
Bambili_1p <- read.csv(paste0(data_path, "Bambili_1%.csv"), header=T, sep = ";")[,-1]
rownames(Bambili_1p) <- read.csv(paste0(data_path, "Bambili_1%.csv"), header=T, sep = ";")[,1]
Bambili_1p[is.na(Bambili_1p)] <- 0

# Terrestrial taxa that are at least 0,5% per sample
Bambili_05p <- read.csv(paste0(data_path, "Bambili_05%.csv"), header=T, sep = ";")[,-1]
rownames(Bambili_05p) <- read.csv(paste0(data_path, "Bambili_05%.csv"), header=T, sep = ";")[,1]
Bambili_05p[is.na(Bambili_05p)] <- 0


# Creating empty pse taxa list
taxa_list_Bambili <- colnames(Bambili[,-1])
createPSE(taxa_list_Bambili, loc = "pse_Bambili1.xlsx")

# Loading in pse list
Bambili_pse <- read_excel(paste0(data_path, "pse_Bambili.xlsx"), sheet = 1, col_names = TRUE)
Bambili_05p_pse <- read_excel(paste0(data_path, "pse_Bambili.xlsx"), sheet = 2, col_names = TRUE)

# Getting modern day data for the area looked at
Bambili_rec <- crest.get_modern_data(df=Bambili_05p, pse=Bambili_05p_pse,
                                     taxaType = 1, climate= c("bio1", "bio12", "ai"),
                                     dbname= "subset_bandAfrica.sqlite3",
                                     continents = "Africa",
                                     countries = countries,
                                     site_info = c(10.24383, 5.93487),
                                     verbose=T)

# Fitting the species and proxy PDFs
Bambili_rec <- crest.calibrate(Bambili_rec,
                               climateSpaceWeighting = T,
                               bin_width = c(0.5, 25, 0.25),
                               shape = c("normal", "lognormal", "normal"),
                               npoints = 2000, verbose = T)

# Plotting general overview of the space
plot_climateSpace(Bambili_rec,
                  climate = c('bio1', 'bio12', 'ai'), bin_width = c(0.5, 25, 0.25),
                  save = F, filename = "Bambili_05p_ClimateSpace",
                  as.png = T, png.res = 500, width = 7.48, height = 6, y0 = 0.5,
                  add_modern = T)

# removing taxa in families "Amaranthaceae", "Asteraceae", "Poaceae"
Bambili_rec <- excludeTaxa(Bambili_rec, c("Amaranthaceae.undiff.", "Artemisia",
                                          "Asteraceae.undiff.", "Centaurea.type",
                                          "Poaceae.undiff.")
                                      , c("bio1", "bio12", "ai"))
print(Bambili_rec$inputs$selectedTaxa)

# Actual reconstruction function
Bambili_rec <- crest.reconstruct(Bambili_rec, verbose = T)

# Plotting of annual temperature (bio1) and precipitation (bio12) of ancient data
plot(Bambili_rec, climate = 'bio1', add_modern=T, 
     xlim= c(0,20200), ylim=c(12,30),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5,
     save = F, width = 10, height = 7.5,  as.png = T,
     png.res = 500, filename = "Bambili_Temperature_05p.pdf",)

plot(Bambili_rec, climate = 'bio12', add_modern=T, 
     xlim= c(0,20200),ylim=c(200, 2500),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     ave = F, width = 10, height = 7.5,  as.png = T,
     png.res = 500, filename = "Bambili_Precip_05p.pdf",)

plot(Bambili_rec, climate = 'ai', add_modern=T,
     xlim= c(0,20200), ylim=c(0,2), 
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = F, width = 10, height = 7.5,  as.png = T,
     png.res = 500, filename = "Bambili_Aridity_05p.pdf",)

# Leave One Out analysis and plots
Bambili_rec <- loo(Bambili_rec)
plot_loo(Bambili_rec, save = T, filename = "Bambili_LOO_05p", 
         as.png=T, png.res = 500, width=7.5, height = 10, xlim= c(0,20200),
         sort='decr', col_pos = 'blue', col_neg = 'red')