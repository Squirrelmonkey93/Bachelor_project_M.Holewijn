library(crestr)
library(readxl)

data_path <- "Datasets/Deva-Deva/"

# Loading in complete original data set
Deva_org <- read.csv(paste0(data_path, "Deva-Deva_site26662.csv"), header=T, sep = ",")

# Loading in cleaned up data sets
# All terrestrial taxa
Deva <- read.csv(paste0(data_path, "Deva-Deva_100%.csv"), header=T, sep = ";")[,-1]
rownames(Deva) <- read.csv(paste0(data_path, "Deva-Deva_100%.csv"), header=T, sep = ";")[,1]
Deva[is.na(Deva)] <- 0

# Terrestrial taxa that are at least 1% per sample
Deva_1p <- read.csv(paste0(data_path, "Deva-Deva_1%.csv"), header=T, sep = ";")[,-1]
rownames(Deva_1p) <- read.csv(paste0(data_path, "Deva-Deva_1%.csv"), header=T, sep = ";")[,1]
Deva_1p[is.na(Deva_1p)] <- 0

# Terrestrial taxa that are at least 0,5% per sample
Deva_05p <- read.csv(paste0(data_path, "Deva-Deva_05%.csv"), header=T, sep = ";")[,-1]
rownames(Deva_05p) <- read.csv(paste0(data_path, "Deva-Deva_05%.csv"), header=T, sep = ";")[,1]
Deva_05p[is.na(Deva_05p)] <- 0


# Creating empty pse taxa list
taxa_list_Deva <- colnames(Deva[,-1])
createPSE(taxa_list_Deva, loc = "pse_Deva-Deva1.xlsx")

# Loading in pse list
Deva_pse <- read_excel(paste0(data_path, "pse_Deva-Deva.xlsx"), sheet = 1, col_names = TRUE)
Deva_05p_pse <- read_excel(paste0(data_path, "pse_Deva-Deva.xlsx"), sheet = 2, col_names = TRUE)

# Getting modern day data for the area looked at
Deva_rec <- crest.get_modern_data(df=Deva_05p, pse=Deva_05p_pse,
                                   taxaType = 1, climate= c("bio1", "bio12", "ai"),
                                   dbname= "subset_bandAfrica.sqlite3",
                                   continents = "Africa",
                                   countries = countries,
                                   site_info = c(37.62533, -7.1222),
                                   verbose=T)

# Fitting the species and proxy PDFs
Deva_rec <- crest.calibrate(Deva_rec,
                             climateSpaceWeighting = T,
                             bin_width = c(0.5, 25, 0.25),
                             shape = c("normal", "lognormal", "normal"),
                             npoints = 2000,
                             verbose = T)

# Plotting general overview of the space
plot_climateSpace(Deva_rec,
                  climate = c('bio1', 'bio12', 'ai'), bin_width = c(0.5, 25, 0.25),
                  save = T, filename = "Deva_05p_ClimateSpace",
                  as.png = T, png.res = 500, width = 7.48, height = 6, y0 = 0.5,
                  add_modern = T)

Deva_rec <- excludeTaxa(Deva_rec, c("Asteroideae", "Carduus.type",
                                    "Crassocephalum.type..C..montuosum.", "Vernonia.type",
                                    "Poaceae.undiff"), c("bio1", "bio12", "ai"))
# print(Deva_rec$inputs$selectedTaxa)

# Actual reconstruction function
Deva_rec <- crest.reconstruct(Deva_rec, verbose = T)

# Plotting of annual temperature (bio1) and precipitation (bio12) of ancient data
plot(Deva_rec, climate = 'bio1', add_modern=T, 
     xlim= c(0,20200), ylim=c(12,30),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     ave = F, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Deva_Temperature_05p.pdf",)

plot(Deva_rec, climate = 'bio12', add_modern=T, 
     xlim= c(0,20200), ylim=c(200, 2500), 
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = T, width = 10, height = 7.5, as.png = T, 
     png.res = 500, filename = "Deva_Precip_05p.pdf",)

plot(Deva_rec, climate = 'ai', add_modern=T,
     xlim= c(0,20200), ylim=c(0,2), 
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = T, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Deva_Aridity_05p.pdf",)

# Leave One Out analysis and plots
Deva_rec <- loo(Deva_rec)
plot_loo(Deva_rec, save = T, filename = "Deva_LOO_05p", 
         as.png=T, png.res = 500, width=7.5, height = 10, xlim= c(0,20200),
         sort='decr', col_pos = 'blue', col_neg = 'red')