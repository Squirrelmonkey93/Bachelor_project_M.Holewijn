library(crestr)
library(readxl)

data_path <- "Datasets/Rumuiku/"

# Loading in complete original data set
Rumuiku_org <- read.csv(paste0(data_path, "Rumuiku_site26667.csv"), header=T, sep = ",")

# Loading in cleaned up data sets
# All terrestrial taxa
Rumuiku <- read.csv(paste0(data_path, "Rumuiku_100%.csv"), header=T, sep = ";")[,-1]
rownames(Rumuiku) <- read.csv(paste0(data_path, "Rumuiku_100%.csv"), header=T, sep = ";")[,1]
Rumuiku[is.na(Rumuiku)] <- 0

# Terrestrial taxa that are at least 1% per sample
Rumuiku_1p <- read.csv(paste0(data_path, "Rumuiku_1%.csv"), header=T, sep = ";")[,-1]
rownames(Rumuiku_1p) <- read.csv(paste0(data_path, "Rumuiku_1%.csv"), header=T, sep = ";")[,1]
Rumuiku_1p[is.na(Rumuiku_1p)] <- 0

# Terrestrial taxa that are at least 0,5% per sample
Rumuiku_05p <- read.csv(paste0(data_path, "Rumuiku_05%.csv"), header=T, sep = ";")[,-1]
rownames(Rumuiku_05p) <- read.csv(paste0(data_path, "Rumuiku_05%.csv"), header=T, sep = ";")[,1]
Rumuiku_05p[is.na(Rumuiku_05p)] <- 0


# Creating empty pse taxa list
taxa_list_Rumuiku <- colnames(Rumuiku[,-1])
createPSE(taxa_list_Rumuiku, loc = "pse_Rumuiku1.xlsx")

# Loading in pse list
Rumuiku_pse <- read_excel(paste0(data_path, "pse_Rumuiku.xlsx"), sheet = 1, col_names = TRUE)
Rumuiku_05p_pse <- read_excel(paste0(data_path, "pse_Rumuiku.xlsx"), sheet = 2, col_names = TRUE)

# Getting modern day data for the area looked at
Rumuiku_rec <- crest.get_modern_data(df=Rumuiku_05p, pse=Rumuiku_05p_pse,
                                     taxaType = 1, climate= c("bio1", "bio12", "ai"),
                                     dbname= "subset_bandAfrica.sqlite3",
                                     continents = "Africa",
                                     countries = countries,
                                     site_info = c(37.535985, -0.46234),
                                     verbose=T)

# Fitting the species and proxy PDFs
Rumuiku_rec <- crest.calibrate(Rumuiku_rec,
                               climateSpaceWeighting = T,
                               bin_width = c(0.5, 25, 0.25),
                               shape = c("normal", "lognormal", "normal"),
                               npoints = 2000,
                               verbose = T)

# Plotting general overview of the space
plot_climateSpace(Rumuiku_rec,
                  climate = c('bio1', 'bio12', 'ai'), bin_width = c(0.5, 25, 0.25),
                  save = T, filename = "Rumuika_05p_ClimateSpace",
                  as.png = T, png.res = 500, width = 7.48, height = 6, y0 = 0.5,
                  add_modern = T)


# removing taxa in families "Amaranthaceae", "Asteraceae", "Poaceae"
Rumuiku_rec <- excludeTaxa(Rumuiku_rec, c("Amaranthaceae", "Artemisia", 
                                          "Asteraceae", "Stoebe", "Vernonia",
                                          "Poaceae")
                           , c("bio1", "bio12", "ai"))
# print(Rumuiku_rec$inputs$selectedTaxa)

# Actual reconstruction function
Rumuiku_rec <- crest.reconstruct(Rumuiku_rec, verbose = T)

# Plotting of annual temperature (bio1) and precipitation (bio12) of ancient data
plot(Rumuiku_rec, climate = 'bio1', add_modern=T, 
     xlim= c(0,20200), ylim=c(12,30),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = F, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Rumuiku_Temperature_05p.pdf",)

plot(Rumuiku_rec, climate = 'bio12', add_modern=T, 
     xlim= c(0,20200), ylim=c(200, 2500),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = T, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Rumuiku_Precip_05p.pdf",)

plot(Rumuiku_rec, climate = 'ai', add_modern=T,
     xlim= c(0,20200), ylim=c(0,2), 
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5,
     save = T, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Rumuiku_Aridity_05p.pdf",)

# Leave One Out analysis and plots
Rumuiku_rec <- loo(Rumuiku_rec)
plot_loo(Rumuiku_rec, save = T, filename = "Rumuiku_LOO_05p", 
         as.png=T, png.res = 500, width=7.5, height = 10, xlim= c(0,20200),
         sort='decr', col_pos = 'blue', col_neg = 'red')