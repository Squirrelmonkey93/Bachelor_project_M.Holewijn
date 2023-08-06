library(crestr)
library(readxl)

data_path <- "Datasets/Kitumbako/"

# Loading in complete original data set
Kitumbako_org <- read.csv(paste0(data_path, "Kitumbako_site26663.csv"), header=T, sep = ",")

# Loading in cleaned up data sets
# All terrestrial taxa
Kitumbako <- read.csv(paste0(data_path, "Kitumbako_100%.csv"), header=T, sep = ";")[,-1]
rownames(Kitumbako) <- read.csv(paste0(data_path, "Kitumbako_100%.csv"), header=T, sep = ";")[,1]
Kitumbako[is.na(Kitumbako)] <- 0

# Terrestrial taxa that are at least 1% per sample
Kitumbako_1p <- read.csv(paste0(data_path, "Kitumbako_1%.csv"), header=T, sep = ";")[,-1]
rownames(Kitumbako_1p) <- read.csv(paste0(data_path, "Kitumbako_1%.csv"), header=T, sep = ";")[,1]
Kitumbako_1p[is.na(Kitumbako_1p)] <- 0

# Terrestrial taxa that are at least 0,5% per sample
Kitumbako_05p <- read.csv(paste0(data_path, "Kitumbako_05%.csv"), header=T, sep = ";")[,-1]
rownames(Kitumbako_05p) <- read.csv(paste0(data_path, "Kitumbako_05%.csv"), header=T, sep = ";")[,1]
Kitumbako_05p[is.na(Kitumbako_05p)] <- 0


# Creating empty pse taxa list
taxa_list_Kitumbako <- colnames(Kitumbako[,-1])
createPSE(taxa_list_Kitumbako, loc = "pse_Kitumbako1.xlsx")

# Loading in pse list
Kitumbako_pse <- read_excel(paste0(data_path, "pse_Kitumbako.xlsx"), sheet = 1, col_names = TRUE)
Kitumbako_05p_pse <- read_excel(paste0(data_path, "pse_Kitumbako.xlsx"), sheet = 2, col_names = TRUE)

# Getting modern day data for the area looked at
Kitumbako_rec <- crest.get_modern_data(df=Kitumbako_05p, pse=Kitumbako_05p_pse,
                                  taxaType = 1, climate= c("bio1", "bio12", "ai"),
                                  dbname= "subset_bandAfrica.sqlite3",
                                  continents = "Africa",
                                  countries = countries,
                                  site_info = c(37.631333, -7.145333),
                                  verbose=T)

# Fitting the species and proxy PDFs
Kitumbako_rec <- crest.calibrate(Kitumbako_rec,
                            climateSpaceWeighting = T,
                            bin_width = c(0.5, 25, 0.25),
                            shape = c("normal", "lognormal", "normal"),
                            npoints = 2000,
                            verbose = T)

# Plotting general overview of the space
plot_climateSpace(Kitumbako_rec,
                  climate = c('bio1', 'bio12', 'ai'), bin_width = c(0.5, 25, 0.25),
                  save = T, filename = "Kitumbako_05p_ClimateSpace",
                  as.png = T, png.res = 500, width = 7.48, height = 6, y0 = 0.5,
                  add_modern = T)

# removing taxa in families "Asteraceae", "Poaceae"
Kitumbako_rec <- excludeTaxa(Kitumbako_rec, 
                             c("Asteroideae.undiff.", "Carduus.type", "Vernonia",
                               "Poaceae.undiff."), c("bio1", "bio12", "ai"))
print(Kitumbako_rec$inputs$selectedTaxa)

# Actual reconstruction function
Kitumbako_rec <- crest.reconstruct(Kitumbako_rec, verbose = T)

# Plotting of annual temperature (bio1) and precipitation (bio12) of ancient data
plot(Kitumbako_rec, climate = 'bio1', add_modern=T, 
     xlim= c(0,20200), ylim=c(12,30),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = F, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Kitumbako_Temperature_05p.pdf")

plot(Kitumbako_rec, climate = 'bio12', add_modern=T, 
     xlim= c(0,20200), ylim=c(200, 2500),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = T, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Kitumbako_Precip_05p.pdf")

plot(Kitumbako_rec, climate = 'ai', add_modern=T,
     xlim= c(0,20200), ylim=c(0,2), 
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = T, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Kitumbako_Aridity_05p.pdf")

# Leave One Out analysis and plots
Kitumbako_rec <- loo(Kitumbako_rec)
plot_loo(Kitumbako_rec, save = T, filename = "Kitumbako_LOO_05p",
         as.png=T, png.res = 500, width=7.5, height = 10, xlim= c(0,20200),
         sort='decr', col_pos = 'blue', col_neg = 'red')