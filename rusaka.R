library(crestr)
library(readxl)

data_path <- "Datasets/Rusaka/"

# Loading in complete original data set
Rusaka_org <- read.csv(paste0(data_path, "Rusaka_site2233.csv"), header=T, sep = ",")

# Loading in cleaned up data sets
# All terrestrial taxa
Rusaka <- read.csv(paste0(data_path, "Rusaka_100%.csv"), header=T, sep = ";")[,-1]
rownames(Rusaka) <- read.csv(paste0(data_path, "Rusaka_100%.csv"), header=T, sep = ";")[,1]
Rusaka[is.na(Rusaka)] <- 0

# Terrestrial taxa that are at least 1% per sample
Rusaka_1p <- read.csv(paste0(data_path, "Rusaka_1%.csv"), header=T, sep = ";")[,-1]
rownames(Rusaka_1p) <- read.csv(paste0(data_path, "Rusaka_1%.csv"), header=T, sep = ";")[,1]
Rusaka_1p[is.na(Rusaka_1p)] <- 0

# Terrestrial taxa that are at least 0,5% per sample
Rusaka_05p <- read.csv(paste0(data_path, "Rusaka_05%.csv"), header=T, sep = ";")[,-1]
rownames(Rusaka_05p) <- read.csv(paste0(data_path, "Rusaka_05%.csv"), header=T, sep = ";")[,1]
Rusaka_05p[is.na(Rusaka_05p)] <- 0


# Creating empty pse taxa list
taxa_list_Rusaka <- colnames(Rusaka[,-1])
createPSE(taxa_list_Rusaka, loc = "pse_Rusaka.xlsx")

# Loading in pse list
Rusaka_pse <- read_excel(paste0(data_path, "pse_Rusaka.xlsx"), sheet = 1, col_names = TRUE)
Rusaka_05p_pse <- read_excel(paste0(data_path, "pse_Rusaka.xlsx"), sheet = 2, col_names = TRUE)

# Getting modern day data for the area looked at
Rusaka_rec <- crest.get_modern_data(df=Rusaka_05p, pse=Rusaka_05p_pse,
                                     taxaType = 1, climate= c("bio1", "bio12", "ai"),
                                     dbname= "subset_bandAfrica.sqlite3",
                                     continents = "Africa",
                                     countries = countries,
                                     site_info = c(29.61667, -3.43333),
                                     verbose=T)

# Fitting the species and proxy PDFs
Rusaka_rec <- crest.calibrate(Rusaka_rec,
                               climateSpaceWeighting = T,
                               bin_width = c(0.5, 25, 0.25),
                               shape = c("normal", "lognormal", "normal"),
                               npoints = 2000,
                               verbose = T)

# Plotting general overview of the space
plot_climateSpace(Rusaka_rec,
                  climate = c('bio1', 'bio12', 'ai'), bin_width = c(0.5, 25, 0.25),
                  save = T, filename = "Rusaka_05p_ClimateSpace",
                  as.png = T, png.res = 500, width = 7.48, height = 6, y0 = 0.5,
                  add_modern = T)


# removing taxa in families "Amaranthaceae", "Asteraceae", "Poaceae"
Rusaka_rec <- excludeTaxa(Rusaka_rec, c("Achyranthes.aspera.type", 
                                        "Amaranthaceae", "Sericostachys.scandens.type",
                                        "Artemisia", "Asteroideae", "Carduus", 
                                        "Crassocephalum.montuosum.type", "Solanecio.mannii.type",
                                        "Poaceae")
                           , c("bio1", "bio12", "ai"))
# print(Rusaka_rec$inputs$selectedTaxa)

# Actual reconstruction function
Rusaka_rec <- crest.reconstruct(Rusaka_rec, verbose = T)

# Plotting of annual temperature (bio1) and precipitation (bio12) of ancient data
plot(Rusaka_rec, climate = 'bio1', add_modern=T, 
     xlim= c(0,20200), ylim=c(12,30),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5,
     save = F, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Rusaka_Temperature_05p.pdf",)

plot(Rusaka_rec, climate = 'bio12', add_modern=T, 
     xlim= c(0,20200), ylim=c(200, 2500),
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5, 
     save = T, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Rusaka_Precip_05p.pdf",)

plot(Rusaka_rec, climate = 'ai', add_modern=T,
     xlim= c(0,20200), ylim=c(0,2), 
     simplify = T, uncertainties = 0.5, pt.cex = 1, pt.lwd = 1.5,
     save = T, width = 10, height = 7.5, as.png = T,
     png.res = 500, filename = "Rusaka_Aridity_05p.pdf",)

# Leave One Out analysis and plots
Rusaka_rec <- loo(Rusaka_rec)
plot_loo(Rusaka_rec, save = T, filename = "Rusaka_LOO_05p", 
         as.png=T, png.res = 500, width=7.5, height = 10, xlim= c(0,20200),
         sort='decr', col_pos = 'blue', col_neg = 'red')