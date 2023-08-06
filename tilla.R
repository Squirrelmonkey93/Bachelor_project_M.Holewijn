library(crestr)
library(readxl)

data_path <- "Datasets/Tilla/"

# Loading in complete original data set
Tilla_org <- read.csv(paste0(data_path, "Tilla_site27424.csv"), header=T, sep = ",")

# Loading in cleaned up data sets
# All terrestrial taxa
Tilla <- read.csv(paste0(data_path, "Tilla_100%.csv"), header=T, sep = ";")[,-1]
rownames(Tilla) <- read.csv(paste0(data_path, "Tilla_100%.csv"), header=T, sep = ";")[,1]
Tilla[is.na(Tilla)] <- 0

# Terrestrial taxa that are at least 1% per sample
Tilla_1p <- read.csv(paste0(data_path, "Tilla_1%.csv"), header=T, sep = ";")[,-1]
rownames(Tilla_1p) <- read.csv(paste0(data_path, "Tilla_1%.csv"), header=T, sep = ";")[,1]
Tilla_1p[is.na(Tilla_1p)] <- 0

# Terrestrial taxa that are at least 0,5% per sample
Tilla_05p <- read.csv(paste0(data_path, "Tilla_05%.csv"), header=T, sep = ";")[,-1]
rownames(Tilla_05p) <- read.csv(paste0(data_path, "Tilla_05%.csv"), header=T, sep = ";")[,1]
Tilla_05p[is.na(Tilla_05p)] <- 0


# Creating empty pse taxa list
taxa_list_Tilla <- colnames(Tilla[,-1])
createPSE(taxa_list_Tilla, loc = "pse_Tilla1.xlsx")

# Loading in pse list
Tilla_pse <- read_excel(paste0(data_path, "pse_Tilla.xlsx"), sheet = 1, col_names = TRUE)
Tilla_05p_pse <- read_excel(paste0(data_path, "pse_Tilla.xlsx"), sheet = 2, col_names = TRUE)

# Getting modern day data for the area looked at
Tilla_rec <- crest.get_modern_data(df=Tilla_05p, pse=Tilla_05p_pse,
                                   taxaType = 1, climate= c("bio1", "bio12", "ai"),
                                   dbname= "subset_bandAfrica.sqlite3",
                                   continents = "Africa",
                                   countries = countries,
                                   site_info = c(12.129722, 10.395556),
                                   verbose=T)

# Fitting the species and proxy PDFs
Tilla_rec <- crest.calibrate(Tilla_rec,
                             climateSpaceWeighting = T,
                             bin_width = c(0.5, 25, 0.25),
                             shape = c("normal", "lognormal", "normal"),
                             npoints = 2000,
                             verbose = T)

# Plotting general overview of the space
plot_climateSpace(Tilla_rec,
                  climate = c('bio1', 'bio12', 'ai'), bin_width = c(0.5, 25, 0.25),
                  save = T, filename = "Tilla_05p_ClimateSpace",
                  as.png = T, png.res = 500, width = 7.48, height = 6, y0 = 0.5,
                  add_modern = T)

Tilla_rec <- excludeTaxa(Tilla_rec, c("Asteraceae", "Amaranthaceae", "Poaceae"), c("bio1", "bio12", "ai"))
# print(Tilla_rec$inputs$selectedTaxa)

# Actual reconstruction function
Tilla_rec <- crest.reconstruct(Tilla_rec, verbose = T)

# Plotting of annual temperature (bio1) and precipitation (bio12) of ancient data
plot(Tilla_rec, climate = 'bio1', add_modern=T, 
     xlim= c(0,20200), ylim=c(12,30),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = F, width = 10, height = 7.5,  as.png = T,
     png.res = 500, filename = "Tilla_Temperature_05p.pdf",)

plot(Tilla_rec, climate = 'bio12', add_modern=T, 
     xlim= c(0,20200),ylim=c(200, 2500),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = T, width = 10, height = 7.5,  as.png = T,
     png.res = 500, filename = "Tilla_Precip_05p.pdf",)

plot(Tilla_rec, climate = 'ai', add_modern=T,
     xlim= c(0,20200), ylim=c(0,2), 
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = T, width = 10, height = 7.5,  as.png = T,
     png.res = 500, filename = "Tilla_Aridity_05p.pdf",)

# Leave One Out analysis and plots
Tilla_rec <- loo(Tilla_rec)
plot_loo(Tilla_rec, save = T, filename = "Tilla_LOO_05p", 
         as.png=T, png.res = 500, width=7.5, height = 10, xlim= c(0,20200),
         sort='decr', col_pos = 'blue', col_neg = 'red')