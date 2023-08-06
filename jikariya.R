library(crestr)
library(readxl)

data_path <- "Datasets/Jikariya/"

# Loading in complete original data set
Jikariya_org <- read.csv(paste0(data_path, "Jikariya_site28404.csv"), header=T, sep = ",")

# Loading in cleaned up data sets
# All terrestrial taxa
Jikariya <- read.csv(paste0(data_path, "Jikariya_100%.csv"), header=T, sep = ";")[,-1]
rownames(Jikariya) <- read.csv(paste0(data_path, "Jikariya_100%.csv"), header=T, sep = ";")[,1]
Jikariya[is.na(Jikariya)] <- 0

# Terrestrial taxa that are at least 1% per sample
Jikariya_1p <- read.csv(paste0(data_path, "Jikariya_1%.csv"), header=T, sep = ";")[,-1]
rownames(Jikariya_1p) <- read.csv(paste0(data_path, "Jikariya_1%.csv"), header=T, sep = ";")[,1]
Jikariya_1p[is.na(Jikariya_1p)] <- 0

# Terrestrial taxa that are at least 0,5% per sample
Jikariya_05p <- read.csv(paste0(data_path, "Jikariya_05%.csv"), header=T, sep = ";")[,-1]
rownames(Jikariya_05p) <- read.csv(paste0(data_path, "Jikariya_05%.csv"), header=T, sep = ";")[,1]
Jikariya_05p[is.na(Jikariya_05p)] <- 0


# Creating empty pse taxa list
taxa_list_Jikariya <- colnames(Jikariya[,-1])
createPSE(taxa_list_Jikariya, loc = "pse_Jikariya1.xlsx")

# Loading in pse list
Jikariya_pse <- read_excel(paste0(data_path, "pse_Jikariya.xlsx"), sheet = 1, col_names = TRUE)
Jikariya_05p_pse <- read_excel(paste0(data_path, "pse_Jikariya.xlsx"), sheet = 2, col_names = TRUE)

# Getting modern day data for the area looked at
Jikariya_rec <- crest.get_modern_data(df=Jikariya_05p, pse=Jikariya_05p_pse,
                                     taxaType = 1, climate= c("bio1", "bio12", "ai"),
                                     dbname= "subset_bandAfrica.sqlite3",
                                     continents = "Africa",
                                     countries = countries,
                                     site_info = c(11.077, 13.313667),
                                     verbose=T)

# Fitting the species and proxy PDFs
Jikariya_rec <- crest.calibrate(Jikariya_rec,
                               climateSpaceWeighting = T,
                               bin_width = c(0.5, 25, 0.25),
                               shape = c("normal", "lognormal", "normal"),
                               npoints = 2000,
                               verbose = T)

# Plotting general overview of the space
plot_climateSpace(Jikariya_rec,
                  climate = c('bio1', 'bio12', 'ai'), bin_width = c(0.5, 25, 0.25),
                  save = T, filename = "Jikariya_05p_ClimateSpace",
                  as.png = T, png.res = 500, width = 7.48, height = 6, y0 = 0.5,
                  add_modern = T)


# removing taxa in families "Amaranthaceae", "Asteraceae", "Poaceae"
Jikariya_rec <- excludeTaxa(Jikariya_rec, c("Amaranthaceae", "Asteroideae", "Poaceae")
                           , c("bio1", "bio12", "ai"))
# print(Jikariya_rec$inputs$selectedTaxa)

# Actual reconstruction function
Jikariya_rec <- crest.reconstruct(Jikariya_rec, verbose = T)

# Plotting of annual temperature (bio1) and precipitation (bio12) of ancient data
plot(Jikariya_rec, climate = 'bio1', add_modern=T, 
     xlim= c(0,20200), ylim=c(12,30),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5,
     save = F, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Jikariya_Temperature_05p",)

plot(Jikariya_rec, climate = 'bio12', add_modern=T, 
     xlim= c(0,20200), ylim=c(200, 2500),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = T, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Jikariya_Precip_05p",)

plot(Jikariya_rec, climate = 'ai', add_modern=T,
     xlim= c(0,20200), ylim=c(0,2), 
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5,
     save = T, width = 10, height = 7.5,  as.png = T,
     png.res = 500, filename = "Jikariya_Aridity_05p",)

# Leave One Out analysis and plots
Jikariya_rec <- loo(Jikariya_rec)
plot_loo(Jikariya_rec, save = T, filename = "Jikariya_LOO_05p", 
         as.png=T, png.res = 500, width=7.5, height = 10, xlim= c(0,20200),
         sort='decr', col_pos = 'blue', col_neg = 'red')