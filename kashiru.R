library(crestr)
library(readxl)

data_path <- "Datasets/Kashiru/"

# Loading in complete original data set
Kashiru595_org <- read.csv(paste0(data_path, "Kashiru(47595)_site27021.csv"), header=T, sep = ",")

# Loading in cleaned up data sets
# All terrestrial taxa
Kashiru595 <- read.csv(paste0(data_path, "kashiru595_100%.csv"), header=T, sep = ";")[,-1]
rownames(Kashiru595) <- read.csv(paste0(data_path, "kashiru595_100%.csv"), header=T, sep = ";")[,1]
Kashiru595[is.na(Kashiru595)] <- 0

# Terrestrial taxa that are at least 1% per sample
Kashiru595_1p <- read.csv(paste0(data_path, "kashiru595_1%.csv"), header=T, sep = ";")[,-1]
rownames(Kashiru595_1p) <- read.csv(paste0(data_path, "kashiru595_1%.csv"), header=T, sep = ";")[,1]
Kashiru595_1p[is.na(Kashiru595_1p)] <- 0

# Terrestrial taxa that are at least 0,5% per sample
Kashiru595_05p <- read.csv(paste0(data_path, "kashiru595_05%.csv"), header=T, sep = ";")[,-1]
rownames(Kashiru595_05p) <- read.csv(paste0(data_path, "kashiru595_05%.csv"), header=T, sep = ";")[,1]
Kashiru595_05p[is.na(Kashiru595_05p)] <- 0

# Creating empty pse taxa list
taxa_list_Kashiru595 <- colnames(Kashiru595[,-1])
createPSE(taxa_list_Kashiru595, loc = "pse_Kashiru595.xlsx")

# Loading in pse list
Kashiru595_pse <- read_excel(paste0(data_path, "pse_Kashiru595.xlsx"), sheet = 1, col_names = TRUE)
Kashiru595_05p_pse <- read_excel(paste0(data_path, "pse_Kashiru595.xlsx"), sheet = 2, col_names = TRUE)

# Getting modern day data for the area looked at
Kashiru595_rec <- crest.get_modern_data(df=Kashiru595_05p, pse=Kashiru595_05p_pse,
                                  taxaType = 1, climate= c("bio1", "bio12", "ai"),
                                  dbname= "subset_bandAfrica.sqlite3",
                                  continents = "Africa",
                                  countries = countries,
                                  site_info = c(29.53333, -3.45),
                                  verbose=T)

# Fitting the species and proxy PDFs
Kashiru595_rec <- crest.calibrate(Kashiru595_rec,
                            climateSpaceWeighting = T,
                            bin_width = c(0.5, 25, 0.25),
                            shape = c("normal", "lognormal", "normal"), # 
                            npoints = 2000,
                            verbose = T)

# Plotting general overview of the space
plot_climateSpace(Kashiru595_rec,
                  climate = c('bio1', 'bio12', 'ai'), bin_width = c(0.5, 25, 0.25),
                  save = T, filename = "Kashiru_05p_ClimateSpace",
                  as.png = T, png.res = 500, width = 7.48, height = 6, y0 = 0.5,
                  add_modern = T)

# removing taxa in families "Asteraceae", "Amaranthaceae", "Poaceae"
Kashiru595_rec <- excludeTaxa(Kashiru595_rec, 
                              c("Amaranthaceae",
                                "Artemisia",  "Asteroideae", "Carduus.type", "Crassocephalum.type..C..montuosum.",
                                "Poaceae"
                              ), 
                              c("bio1", "bio12", "ai"))
# print(kashiru595_rec$inputs$selectedTaxa)

# Actual reconstruction function
Kashiru595_rec <- crest.reconstruct(Kashiru595_rec, verbose = T)

# Plotting of annual temperature (bio1) and precipitation (bio12) of ancient data
plot(Kashiru595_rec, climate = 'bio1', add_modern=T, 
     xlim= c(0,20200), ylim=c(12,30),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = F, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Kashiru595_Temperature_05p.pdf",)

plot(Kashiru595_rec, climate = 'bio12', add_modern=T, 
     xlim= c(0,20200), ylim=c(200, 2500),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = T, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Kashiru595_Precip_05p.pdf",)

plot(Kashiru595_rec, climate = 'ai', add_modern=T,
     xlim= c(0,20200), ylim=c(0,2), 
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5,
     save = T, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Kashiru595_Aridity_05p.pdf",)

# Leave One Out analysis and plots
Kashiru595_rec <- loo(Kashiru595_rec)
plot_loo(Kashiru595_rec, save = T, filename = "Kashiru595_LOO_05p", 
         as.png=T, png.res = 500, width=7.5, height = 10, xlim= c(0,20200),
         sort='decr', col_pos = 'blue', col_neg = 'red')