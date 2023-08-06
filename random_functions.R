library(crestr)

folder_name <- "Datasets/"
dataset_name <- "Tilla/"
data_path <- paste0(folder_name, dataset_name)

# Creating subset of the data as to reduce processing times
dbSubset(1, xmn=-20, xmx=55, ymn=-15, ymx=15, out="subset_bandAfrica")

# Creating selection of countries
countries <- accCountryNames()$Africa[
                -which(accCountryNames()$Africa %in% 
                    c('Seychelles', "Sao Tome and Principe", "Madagascar", "Mauritius", "Mayotte", "Réunion"))]

# Filling in pse list manually by putting in taxa from pse
getTaxonomy(family = "", genus = "rubia", species = "",
            taxaType = 1, depth.out = 6, dbname = "gbif4crest_02" )

plot_loo(Rusaka_rec, 
         save = T, filename = "Rusaka_LOO_05p", as.png=T, png.res = 500,
         width=7.5, height = 10, xlim= c(0,20200),
         sort='decr', col_pos = 'blue', col_neg = 'red')

plot(Rusaka_rec, add_modern=T,
     xlim= c(0,20200), simplify = T, uncertainties = 0.5,
     pt.cex = 1, pt.lwd = 1.5, save = T, width = 10, height = 7.5, 
     as.png=T, png.res = 500, filename = "Rusaka_05p.pdf",)

# Plotting general overview of the space
plot_climateSpace(Bambili_rec,
climate = c('bio1', 'bio12', 'ai'), bin_width = c(0.5, 25, 0.25),
save = T, filename = "Bambili_05p_ClimateSpace",
as.png = T, png.res = 500, width = 7.48, height = 6, y0 = 0.5,
add_modern = T)


debugging_get_modern_data <- function (pse, taxaType, climate, df = NA, ai.sqrt = FALSE,
          xmn = NA, xmx = NA, ymn = NA, ymx = NA, continents = NA,
          countries = NA, basins = NA, sectors = NA, realms = NA,
          biomes = NA, ecoregions = NA, minGridCells = 20, elev_min = NA,
          elev_max = NA, elev_range = NA, year_min = 1900, year_max = 2021,
          nodate = TRUE, type_of_obs = c(1, 2, 3, 8, 9), selectedTaxa = NA,
          site_info = c(NA, NA), site_name = NA, dbname = "gbif4crest_02",
          verbose = TRUE)
{
  if (base::missing(pse))
    pse
  if (base::missing(taxaType))
    taxaType
  if (base::missing(climate))
    climate
  if (methods::is(pse, "tbl"))
    pse <- as.data.frame(pse)
  if (methods::is(df, "tbl"))
    df <- as.data.frame(df)
  if (methods::is(selectedTaxa, "tbl"))
    selectedTaxa <- as.data.frame(selectedTaxa)
  if (verbose)
    cat("\n## Prepping data for database extraction\n")
  if (verbose)
    cat("  <> Checking database connection .......... ")
  if (dbname == "crest_example") {
    dbname <- .exampleDB()
  }
  if (!testConnection(dbname))
    return(NA)
  if (verbose)
    cat("[OK]\n  <> Checking pse .......................... ")
  if (!is.data.frame(pse)) {
    cat("[FAILED]\n")
    stop("The 'pse' variable (proxy_species_equivalency) must be a data frame.\n\n")
  }
  pse <- pse[!is.na((pse[, "ProxyName"])), ]
  pse <- pse[(pse[, "ProxyName"] != ""), ]
  taxa.name <- unique(as.character(pse[, "ProxyName"]))
  if (is.data.frame(df))
    taxa.name <- unique(c(taxa.name, colnames(df)[-1]))
  taxa_to_ignore = c()
  for (tax in taxa.name) {
    if (!tax %in% pse[, "ProxyName"])
      taxa_to_ignore = c(taxa_to_ignore, tax)
  }
  if (verbose)
    cat("[OK]\n  <> Checking climate variables ............ ")
  for (clim in 1:length(climate)) {
    climVar <- accClimateVariables()
    new_clim <- climate
    if (!(climate[clim] %in% climVar[, 1] | climate[clim] %in%
          climVar[, 2])) {
      cat("[FAILED]\n\n")
      stop(paste0("The variable '", climate[clim], "' is not an accepted value. Check the list of accepted values using 'accClimateVariables()'.\n"))
    }
    else {
      if (suppressWarnings(!is.na(as.numeric(climate[clim])))) {
        new_clim[clim] <- as.character(climVar[which(climVar[,
                                                             1] == as.numeric(climate[clim])), 2])
      }
    }
  }
  climate <- new_clim
  if (verbose)
    cat("[OK]\n  <> Checking taxaType ..................... ")
  if (taxaType > 6 | taxaType < 0) {
    cat("[FAILED]\n\n")
    stop("'taxaType' should be an integer between 0 and 6. See ?crest.get_modern_data for more information.\n")
  }
  if (verbose)
    cat("[OK]\n  <> Checking coordinates .................. ")
  coords <- check_coordinates(xmn, xmx, ymn, ymx)
  xmn <- coords[1]
  xmx <- coords[2]
  ymn <- coords[3]
  ymx <- coords[4]
  estimate_xlim <- coords[5]
  estimate_ylim <- coords[6]
  if (!is.na(elev_min) & !is.na(elev_max) & elev_min >= elev_max) {
    warning("elev_min was larger than elev_max. The two values were inverted.\n")
    tmp <- elev_min
    elev_min <- elev_max
    elev_max <- tmp
  }
  if (!is.na(elev_range) & elev_range <= 0) {
    stop("elev_range should be a positive integer.\n")
    tmp <- elev_min
    elev_min <- elev_max
    elev_max <- tmp
  }
  if (!is.na(year_min) & !is.na(year_max) & year_min >= year_max) {
    warning("year_min was larger than year_max. The two values were inverted.\n")
    tmp <- year_min
    year_min <- year_max
    year_max <- tmp
  }
  if (taxaType %in% c(1, 2, 3, 6)) {
    if (verbose)
      cat("[OK]\n  <> Checking continent and country names .. ")
    cont.list <- accCountryNames()
    if (!is.na(continents[1])) {
      for (cont in continents) {
        if (!cont %in% names(cont.list)) {
          stop(paste0("The continent '", cont, "' is not an accepted value. Please select a name from this list: '",
                      paste(names(cont.list), collapse = "', '"),
                      "'.\n"))
        }
      }
    }
    if (!is.na(countries)[1]) {
      for (country in countries) {
        if (!country %in% unlist(cont.list)) {
          acc_vals <- ifelse(is.na(continents[1]), "",
                             paste0("c('", paste(continents, collapse = "', '"),
                                    "')"))
          stop(paste0("The country '", country, "' is not an accepted value. Get the list of accepted values using 'accCountryNames(",
                      acc_vals, ")'.\n"))
        }
      }
    }
    if (!is.na(countries[1]) | !is.na(continents[1])) {
      res <- dbRequest(paste0("SELECT DISTINCT continent, name, count(*) FROM geopolitical_units WHERE ",
                              ifelse(is.na(continents)[1], "", paste0("continent IN ('",
                                                                      paste(continents, collapse = "', '"), "') ")),
                              ifelse(is.na(continents)[1] | is.na(countries)[1],
                                     "", "AND "), ifelse(is.na(countries)[1], "",
                                                         paste0("name IN ('", paste(countries, collapse = "', '"),
                                                                "') ")), " GROUP BY continent, name"), dbname)
      if (length(res) == 0) {
        cat(paste("Problem here. No result for any of the combination continent x country.\n",
                  sep = ""))
      }
      else {
        res
      }
    }
  }
  else {
    if (verbose)
      cat("[OK]\n  <> Checking basin and sector names ....... ")
    basin.list <- accBasinNames()
    if (!is.na(basins[1])) {
      for (bas in basins) {
        if (!bas %in% names(basin.list)) {
          stop(paste0("The basin '", bas, "' is not an accepted value. Please select a name from this list: '",
                      paste(names(basin.list), collapse = "', '"),
                      "'.\n"))
        }
      }
    }
    if (!is.na(sectors)[1]) {
      for (sector in sectors) {
        if (!sector %in% unlist(basin.list)) {
          acc_vals <- ifelse(is.na(sectors[1]), "",
                             paste0("c('", paste(basins, collapse = "', '"),
                                    "')"))
          stop(paste0("The sector '", sector, "' is not an accepted value. Get the list of accepted values using 'accBasinNames(",
                      acc_vals, ")'.\n"))
        }
      }
    }
    if (!is.na(basins[1]) | !is.na(sectors[1])) {
      res <- dbRequest(paste0("SELECT DISTINCT basin, name, count(*) FROM geopolitical_units WHERE ",
                              ifelse(is.na(basins)[1], "", paste0("basin IN ('",
                                                                  paste(basins, collapse = "', '"), "') ")),
                              ifelse(is.na(basins)[1] | is.na(sectors)[1],
                                     "", "AND "), ifelse(is.na(sectors)[1], "",
                                                         paste0("name IN ('", paste(sectors, collapse = "', '"),
                                                                "') ")), " GROUP BY basin, name"), dbname)
      if (length(res) == 0) {
        cat(paste("Problem here. No result for any of the combination basin x sector.\n",
                  sep = ""))
      }
      else {
        res
      }
    }
  }
  if (verbose)
    cat("[OK]\n  <> Checking realm/biome/ecoregion names .. ")
  realm.list <- accRealmNames()
  if (!is.na(realms[1])) {
    for (realm in realms) {
      if (!realm %in% names(realm.list)) {
        stop(paste0("The realm '", realm, "' is not an accepted value. Please select a name from this list: '",
                    paste(names(realm.list), collapse = "', '"),
                    "'.\n"))
      }
    }
  }
  if (taxaType %in% c(1, 2, 3, 6)) {
    if (!is.na(biomes)[1]) {
      for (biome in biomes) {
        if (!biome %in% unique(unlist(lapply(realm.list,
                                             function(x) return(unique(x[, 1])))))) {
          acc_vals <- ifelse(is.na(realms[1]), "", paste0("c('",
                                                          paste(realms, collapse = "', '"), "')"))
          stop(paste0("The realm '", realm, "' is not an accepted value. Get the list of accepted values using 'accRealmNames(",
                      acc_vals, ")'.\n"))
        }
      }
    }
    if (!is.na(ecoregions)[1]) {
      for (ecoregion in ecoregions) {
        if (!ecoregion %in% unique(unlist(lapply(realm.list,
                                                 function(x) return(unique(x[, 2])))))) {
          acc_vals <- ifelse(is.na(realms[1]), "", paste0("c('",
                                                          paste(realms, collapse = "', '"), "')"))
          stop(paste0("The ecoregion '", ecoregion,
                      "' is not an accepted value. Get the list of accepted values using 'accRealmNames(",
                      acc_vals, ")'.\n"))
        }
      }
    }
    if (!is.na(realms[1]) | !is.na(biomes[1]) | !is.na(ecoregions[1])) {
      s_realms <- ifelse(is.na(realms)[1], "", paste0("realm IN ('",
                                                      paste(realms, collapse = "', '"), "') "))
      s_biomes <- ifelse(is.na(biomes)[1], "", paste0("biome IN ('",
                                                      paste(biomes, collapse = "', '"), "') "))
      s_ecoregions <- ifelse(is.na(ecoregions)[1], "",
                             paste0("ecoregion IN ('", paste(ecoregions,
                                                             collapse = "', '"), "') "))
      res <- dbRequest(paste0("SELECT DISTINCT realm, biome, ecoregion, count(*) FROM biogeography WHERE ",
                              s_realms, ifelse(s_realms != "" & (s_biomes !=
                                                                   "" | s_ecoregions != ""), " AND ", ""), s_biomes,
                              ifelse(s_biomes != "" & s_ecoregions != "",
                                     " AND ", ""), s_ecoregions, " GROUP BY realm, biome,ecoregion"),
                       dbname)
      if (length(res) == 0) {
        cat(paste("Problem here. No result for any of the combination realm x biome x ecoregion .\n",
                  sep = ""))
      }
      else {
        res
      }
    }
  }
  if (verbose)
    cat("[OK]\n  <> Checking/Defining selectedTaxa ........ ")
  if (is.na(as.vector(t(selectedTaxa))[1])) {
    selectedTaxa <- data.frame(matrix(rep(1, length(climate) *
                                            length(taxa.name)), ncol = length(climate)))
    rownames(selectedTaxa) <- taxa.name
    colnames(selectedTaxa) <- climate
  }
  else {
    if (length(which(!rownames(selectedTaxa) %in% taxa.name)) >
        0) {
      selectedTaxa[which(!rownames(selectedTaxa) %in%
                           taxa.name), ] <- -2
      warning("One or more taxa recorded in the selectedTaxa were not recorded in either PSE or df. They are excluded for the rest of the study (their value is set to -2 is `x$inputs$selectedTaxa`.)\n")
    }
  }
  sendWarning <- FALSE
  for (clim in climate) {
    if (!clim %in% colnames(selectedTaxa)) {
      selectedTaxa <- cbind(selectedTaxa, rep(1, nrow(selectedTaxa)))
      colnames(selectedTaxa)[ncol(selectedTaxa)] <- clim
      sendWarning <- TRUE
    }
  }
  if (sendWarning)
    warning("One or more of the selected variables were not in the selectedTaxa table (check for typos?). Missing columns have been added with a default value of 1.\n")
  taxa_notes <- list()
  for (tax in taxa_to_ignore) {
    message <- "Taxon not in the proxy_species_equivalency table."
    if (!message %in% names(taxa_notes)) {
      taxa_notes[[message]] <- c()
      warning(paste0("One or more taxa were are not in the proxy-species equivalence table and have been ignored. Check `x$misc$taxa_notes` for details."))
    }
    taxa_notes[[message]] <- append(taxa_notes[[message]],
                                    tax)
    if (tax %in% rownames(selectedTaxa)) {
      selectedTaxa[tax, ] <- rep(-1, length(climate))
    }
    else {
      selectedTaxa <- rbind(selectedTaxa, rep(-1, length(climate)))
      rownames(selectedTaxa)[nrow(selectedTaxa)] <- tax
    }
  }
  taxa.name <- taxa.name[taxa.name %in% rownames(selectedTaxa)[apply(selectedTaxa,
                                                                     1, sum) >= 0]]
  w <- !(taxa.name %in% rownames(selectedTaxa))
  print("first w")
  print(w)
  if (sum(w) > 0) {
    for (w in which(!(taxa.name %in% rownames(selectedTaxa)))) {
      selectedTaxa <- rbind(selectedTaxa, rep(1, length(climate)))
      rownames(selectedTaxa)[nrow(selectedTaxa)] <- taxa.name[w]
      for (tax in taxa_to_ignore) {
        message <- "Not present in the original selectedTaxa table. Added by default as 1s."
        if (!message %in% names(taxa_notes)) {
          taxa_notes[[message]] <- c()
          warning(paste0("One or more taxa were are not in the selectedTaxa table. They have been added but are not selected for any variable. Check `x$misc$taxa_notes` for details."))
        }
        taxa_notes[[message]] <- append(taxa_notes[[message]],
                                        tax)
        selectedTaxa[tax, climate] <- rep(0, length(climate))
      }
    }
  }
  if (verbose)
    cat("[OK]\n  <> Checking the pse table ................ ")
  if (!("Level" %in% colnames(pse) & "Family" %in% colnames(pse) &
        "Genus" %in% colnames(pse) & "Species" %in% colnames(pse) &
        "ProxyName" %in% colnames(pse))) {
    stop(paste0("\nThe PSE data frame should contain columns with the following 5 names: 'Level', 'Family', 'Genus', 'Species' and 'ProxyName' .\n\n"))
  }
  else {
    pse <- pse[, c("Level", "Family", "Genus", "Species",
                   "ProxyName")]
  }
  w <- (pse$Level == 4)
  print("2nd w")
  print(w)
  if (sum(w) > 0) {
    for (tax in unique(pse$ProxyName[w])) {
      if (tax %in% rownames(selectedTaxa)) {
        selectedTaxa[tax, ] <- rep(-1, length(climate))
      }
      else {
        selectedTaxa <- rbind(selectedTaxa, rep(-1,
                                                length(climate)))
        rownames(selectedTaxa)[nrow(selectedTaxa)] <- tax
      }
      message <- "No association between the proxy names and species"
      if (!message %in% names(taxa_notes)) {
        taxa_notes[[message]] <- c()
        warning(paste0("One or more taxa were not associated with species. Check `x$misc$taxa_notes` for details."))
      }
      taxa_notes[[message]] <- append(taxa_notes[[message]],
                                      tax)
    }
    pse <- pse[!w, ]
  }
  taxa.name <- taxa.name[taxa.name %in% rownames(selectedTaxa)[apply(selectedTaxa,
                                                                     1, sum) >= 0]]
  if (verbose) {
    cat("[OK]\n  <> Extracting taxon species .............. \r")
  }
  crest <- crestObj(taxa.name, pse = pse, taxaType = taxaType,
                    climate = climate, xmn = xmn, xmx = xmx, ymn = ymn,
                    ymx = ymx, continents = continents, countries = countries,
                    basins = basins, sectors = sectors, realms = realms,
                    biomes = biomes, ecoregions = ecoregions, elev_min = elev_min,
                    elev_max = elev_max, elev_range = elev_range, year_min = year_min,
                    year_max = year_max, nodate = nodate, type_of_obs = type_of_obs,
                    selectedTaxa = selectedTaxa, dbname = dbname)
  crest$misc[["taxa_notes"]] <- taxa_notes
  crest$misc$site_info <- list()
  crest$misc$site_info[["long"]] <- site_info[1]
  crest$misc$site_info[["lat"]] <- site_info[2]
  crest$misc$site_info[["site_name"]] <- site_name
  if ((!is.na(crest$misc$site_info[["long"]])) & (!is.na(crest$misc$site_info[["lat"]]))) {
    resol <- ifelse(.ifExampleDB(dbname), 0.5, 0.25)
    crest$misc$site_info[["climate"]] <- climate_from_xy(crest$misc$site_info[["long"]],
                                                         crest$misc$site_info[["lat"]], crest$parameters$climate,
                                                         resol = resol, dbname = crest$misc$dbname)
  }
  if (is.data.frame(df)) {
    df[is.na(df)] <- 0
    crest$inputs$x <- df[, 1]
    crest$inputs$x.name <- colnames(df)[1]
    crest$inputs$taxa.name <- taxa.name
    crest$inputs$df <- df[, -1]
    if (unique(is.numeric(crest$inputs$x))) {
      crest$inputs$df <- crest$inputs$df[order(crest$inputs$x),
      ]
      crest$inputs$x <- crest$inputs$x[order(crest$inputs$x)]
      crest$inputs$x[is.na(crest$inputs$x)] <- 0
    }
    w <- (apply(crest$inputs$df, 2, sum) == 0)
    print("3rd w")
    print(w)
    if (sum(w) > 0) {
      for (tax in colnames(crest$inputs$df)[w]) {
        crest$inputs$selectedTaxa[tax, ] <- rep(-1,
                                                length(climate))
        message <- "All percentages equal to 0."
        if (!message %in% names(crest$misc[["taxa_notes"]])) {
          crest$misc[["taxa_notes"]][[message]] <- c()
          warning(paste0("The percentages of one or more taxa were always 0 and have been removed accordingly. Check `x$misc$taxa_notes` for details."))
        }
        crest$misc[["taxa_notes"]][[message]] <- append(crest$misc[["taxa_notes"]][[message]],
                                                        tax)
      }
    }
    w <- (!taxa.name %in% colnames(df)[-1])
    print("4th w")
    print(w)
    if (sum(w) > 0) {
      for (tax in taxa.name[w]) {
        crest$inputs$selectedTaxa[tax, ] <- rep(-1,
                                                length(climate))
        message <- "Taxon not recorded in the data file."
        if (!message %in% names(crest$misc[["taxa_notes"]])) {
          crest$misc[["taxa_notes"]][[message]] <- c()
          warning(paste0("One or more taxa were are not recorded in the data file. Check `x$misc$taxa_notes` for details."))
        }
        crest$misc[["taxa_notes"]][[message]] <- append(crest$misc[["taxa_notes"]][[message]],
                                                        tax)
      }
    }
  }
  taxonID2proxy <- data.frame(taxonID = NA, proxyName = NA,
                              stringsAsFactors = FALSE)
  pse$Level <- as.numeric(as.character(pse$Level))
  pse$Family <- as.character(pse$Family)
  pse$Genus <- as.character(pse$Genus)
  pse$Species <- as.character(pse$Species)
  pse$ProxyName <- as.character(pse$ProxyName)
  pbi <- 100
  for (taxLevel in 1:3) {
    for (tax in unique(pse$ProxyName[pse$Level == taxLevel])) {
      if (verbose) {
        cat(paste0("  <> Extracting taxon species .............. ",
                   stringr::str_pad(paste0(round(pbi/length(pse$ProxyName)),
                                           "%\r"), width = 4, side = "left")))
        utils::flush.console()
      }
      for (w in which(pse$ProxyName == tax & pse$Level ==
                      taxLevel)) {
        taxonIDs <- getTaxonID(pse$Family[w], pse$Genus[w],
                               pse$Species[w], taxaType, dbname)
        if (length(taxonIDs) > 0) {
          existingTaxa <- taxonIDs %in% taxonID2proxy[,
                                                      "taxonID"]
          if (sum(existingTaxa) > 0) {
            taxonID2proxy[taxonID2proxy[, "taxonID"] %in%
                            taxonIDs, "proxyName"] <- tax
          }
          if (sum(existingTaxa) != length(taxonIDs)) {
            taxonID2proxy <- rbind(taxonID2proxy, data.frame(taxonID = taxonIDs[!existingTaxa],
                                                             proxyName = rep(tax, sum(!existingTaxa)),
                                                             stringsAsFactors = FALSE))
          }
        }
        else {
          if (tax %in% crest$inputs$taxa.name) {
            crest$inputs$selectedTaxa[tax, ] <- rep(-1,
                                                    length(climate))
            message <- "No correspondance with specific species"
            if (!message %in% names(crest$misc[["taxa_notes"]])) {
              crest$misc[["taxa_notes"]][[message]] <- as.data.frame(matrix(0,
                                                                            ncol = 5, nrow = 0))
              colnames(crest$misc[["taxa_notes"]][[message]]) <- colnames(pse)
              warning(paste0("The classification of one or more taxa into species was not successful. Check `x$misc$taxa_notes` for details."))
            }
            crest$misc[["taxa_notes"]][[message]] <- rbind(crest$misc[["taxa_notes"]][[message]],
                                                           pse[w, ])
          }
        }
      }
      pbi <- pbi + 100
    }
  }
  crest$inputs$taxa.name <- crest$inputs$taxa.name[crest$inputs$taxa.name %in%
                                                     rownames(crest$inputs$selectedTaxa)[apply(crest$inputs$selectedTaxa,
                                                                                               1, sum) >= 0]]
  if (verbose) {
    cat("  <> Extracting taxon species .............. [OK]\n  <> Extracting species distributions ...... \r")
  }
  pbi <- 100
  taxonID2proxy <- taxonID2proxy[-1, ]
  taxonID2proxy <- taxonID2proxy[order(taxonID2proxy[, "proxyName"]),
  ]
  crest$modelling$taxonID2proxy <- taxonID2proxy
  distributions <- list()
  for (tax in crest$inputs$taxa.name) {
    taxIDs <- taxonID2proxy[taxonID2proxy[, "proxyName"] ==
                              tax, 1]
    if (length(taxIDs) == 0) {
      crest$inputs$selectedTaxa[tax, ] <- rep(-1, length(climate))
      message <- "No species remained associated with the proxy name at the end of the classification."
      if (!message %in% names(crest$misc[["taxa_notes"]])) {
        crest$misc[["taxa_notes"]][[message]] <- c()
        warning(paste0("For one or more taxa, no species remained associated with the proxy name at the end of the classification. Check `x$misc$taxa_notes` for details."))
      }
      crest$misc[["taxa_notes"]][[message]] <- append(crest$misc[["taxa_notes"]][[message]],
                                                      tax)
    }
    if (verbose) {
      cat(paste0("  <> Extracting species distributions ...... ",
                 stringr::str_pad(paste0(round(pbi/length(crest$inputs$taxa.name)),
                                         "%\r"), width = 4, side = "left")))
      utils::flush.console()
    }
    if (sum(crest$inputs$selectedTaxa[tax, climate] >= 0) >
        0) {
      distributions[[tax]] <- getDistribTaxa(taxIDs, climate = climate,
                                             xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx,
                                             continents = continents, countries = countries,
                                             basins = basins, sectors = sectors, realms = realms,
                                             biomes = biomes, ecoregions = ecoregions, elev_min = elev_min,
                                             elev_max = elev_max, elev_range = elev_range,
                                             year_min = year_min, year_max = year_max, nodate = nodate,
                                             type_of_obs = type_of_obs, dbname = dbname)
      if (nrow(distributions[[tax]]) == 0) {
        distributions[[tax]] <- NA
        crest$inputs$selectedTaxa[tax, ] <- rep(-1,
                                                length(climate))
        message <- "No data point available in the study area."
        if (!message %in% names(crest$misc[["taxa_notes"]])) {
          crest$misc[["taxa_notes"]][[message]] <- c()
          warning(paste0("No data were available within the study area for one or more taxa. Check `x$misc$taxa_notes` for details."))
        }
        crest$misc[["taxa_notes"]][[message]] <- append(crest$misc[["taxa_notes"]][[message]],
                                                        tax)
      }
      else {
        extent_taxa <- table(distributions[[tax]][,
                                                  "taxonid"])
        extent_taxa_id <- as.numeric(names(extent_taxa)[extent_taxa >=
                                                          minGridCells])
        distributions[[tax]] <- distributions[[tax]][distributions[[tax]][,
                                                                          "taxonid"] %in% extent_taxa_id, ]
        if (nrow(distributions[[tax]]) == 0) {
          distributions[[tax]] <- NA
          crest$inputs$selectedTaxa[tax, ] <- rep(-1,
                                                  length(climate))
          message <- "Present but insufficient data in the study area to fit a pdf"
          if (!message %in% names(crest$misc[["taxa_notes"]])) {
            crest$misc[["taxa_notes"]][[message]] <- c()
            warning(paste0("An insufficient amount of calibration data points was available within the study area for one or more taxa. Consider reducing 'minGridCells'. Check `x$misc$taxa_notes` for details."))
          }
          crest$misc[["taxa_notes"]][[message]] <- append(crest$misc[["taxa_notes"]][[message]],
                                                          tax)
        }
      }
    }
    pbi <- pbi + 100
  }
  crest$inputs$taxa.name <- crest$inputs$taxa.name[crest$inputs$taxa.name %in%
                                                     rownames(crest$inputs$selectedTaxa)[apply(crest$inputs$selectedTaxa,
                                                                                               1, sum) >= 0]]
  if (verbose) {
    cat("  <> Extracting species distributions ...... [OK]\n")
  }
  class_names <- rep(NA, nrow(crest$inputs$pse))
  if (crest$parameters$taxaType == 1) {
    pbi <- 100
    for (tax in crest$inputs$taxa.name) {
      if (verbose)
        cat(paste0("  <> Postprocessing plant data ............. ",
                   stringr::str_pad(paste0(round(pbi/length(crest$inputs$taxa.name)),
                                           "%\r"), width = 4, side = "left")))
      utils::flush.console()
      for (w in which(crest$inputs$pse[, "ProxyName"] ==
                      tax)) {
        taxonomy <- getTaxonomy(family = crest$inputs$pse[w,
                                                          "Family"], genus = crest$inputs$pse[w, "Genus"],
                                species = crest$inputs$pse[w, "Species"],
                                taxaType = crest$parameters$taxaType, depth.out = 3,
                                dbname = dbname)
        class_names[w] <- taxonomy[1, "class_name"]
      }
      pbi <- pbi + 100
    }
    if (verbose)
      cat("  <> Postprocessing plant data ............. [OK]\n")
  }
  crest$inputs$pse <- cbind(crest$inputs$pse, Class_name = class_names)
  crest$modelling$distributions <- distributions
  if (verbose) {
    cat("  <> Extracting climate space .............. ")
  }
  climate_space <- getClimateSpace(climate = crest$parameters$climate,
                                   xmn = crest$parameters$xmn, xmx = crest$parameters$xmx,
                                   ymn = crest$parameters$ymn, ymx = crest$parameters$ymx,
                                   continents = crest$parameters$continents, countries = crest$parameters$countries,
                                   basins = crest$parameters$basins, sectors = crest$parameters$sectors,
                                   realms = crest$parameters$realms, biomes = crest$parameters$biomes,
                                   ecoregions = crest$parameters$ecoregions, elev_min = elev_min,
                                   elev_max = elev_max, elev_range = elev_range, dbname)
  if (nrow(climate_space) == 0) {
    stop(paste0("No climate values available in the defined study area N: ",
                crest$parameters$ymx, " S: ", crest$parameters$ymn,
                " W: ", crest$parameters$xmn, " E: ", crest$parameters$xmx,
                ".\n\n"))
  }
  colnames(climate_space)[-c(1, 2)] <- crest$parameters$climate
  crest$modelling$climate_space <- climate_space
  if (ai.sqrt & "ai" %in% crest$parameters$climate) {
    crest$modelling$climate_space[, "ai"] <- sqrt(crest$modelling$climate_space[,
                                                                                "ai"])
    for (tax in crest$inputs$taxa.name) {
      crest$modelling$distributions[[tax]][, "ai"] <- sqrt(crest$modelling$distributions[[tax]][,
                                                                                                "ai"])
    }
    if ((!is.na(crest$misc$site_info[["long"]])) & (!is.na(crest$misc$site_info[["lat"]]))) {
      crest$misc$site_info$climate$ai <- sqrt(crest$misc$site_info$climate$ai)
    }
  }
  resol <- sort(unique(diff(sort(unique(crest$modelling$climate_space[,
                                                                      1])))))[1]/2
  xx <- range(climate_space[, 1])
  if (estimate_xlim) {
    crest$parameters$xmn <- xx[1] - resol
    crest$parameters$xmx <- xx[2] + resol
  }
  else {
    if (crest$parameters$xmn > xx[1] - resol)
      crest$parameters$xmn <- xx[1] - resol
    if (crest$parameters$xmx < xx[2] + resol)
      crest$parameters$xmx <- xx[2] + resol
  }
  resol <- sort(unique(diff(sort(unique(crest$modelling$climate_space[,
                                                                      2])))))[1]/2
  yy <- range(climate_space[, 2])
  if (estimate_ylim) {
    crest$parameters$ymn <- yy[1] - resol
    crest$parameters$ymx <- yy[2] + resol
  }
  else {
    if (crest$parameters$ymn > yy[1] - resol)
      crest$parameters$ymn <- yy[1] - resol
    if (crest$parameters$ymx < yy[2] + resol)
      crest$parameters$ymx <- yy[2] + resol
  }
  if (verbose) {
    cat("[OK]\n")
    cat(paste0("## Data extraction completed.\n\n"))
  }
  if (length(distributions) == 0) {
    warning(paste0("No distributions available in the defined study area N: ",
                   crest$parameters$ymx, " S: ", crest$parameters$ymn,
                   " W: ", crest$parameters$xmn, " E: ", crest$parameters$xmx,
                   ".\n\n"))
  }
  crest$misc$stage <- "data_extracted"
  crest
}


wfunc <- function (x, climate = x$parameters$climate, bin_width = x$parameters$bin_width, 
                   save = FALSE, filename = "Climate_space.pdf", as.png = FALSE, 
                   png.res = 300, width = 7.48, height = min(9, 3.5 * length(climate)), 
                   y0 = 0.5, add_modern = FALSE, resol = 0.25) 
{
  print("start")
  if (base::missing(x)) 
    x
  if (is.crestObj(x)) {
    if (add_modern) {
      if (length(x$misc$site_info) <= 3) {
        add_modern <- FALSE
      }
    }
    if (x$misc$stage == "data_extracted") {
      while (length(bin_width) < length(climate)) bin_width <- c(bin_width, 
                                                                 1)
      bin_width <- bin_width[1:length(climate)]
      ccs <- list()
      x$modelling$xrange <- list()
      bin_width <- as.data.frame(matrix(bin_width, ncol = 1))
      rownames(bin_width) <- climate
      x$parameters$bin_width <- bin_width
      for (clim in x$parameters$climate) {
        print(which(is.na(x$modelling$climate_space[, 
                                                    clim])))
        ccs[[clim]] <- calib_clim_space(x$modelling$climate_space[, 
                                                                  clim], x$parameters$bin_width[clim, ])
        x$modelling$xrange[[clim]] <- fit_xrange(ccs[[clim]], 
                                                 x$parameters$shape[clim, ], x$parameters$bin_width[clim, 
                                                 ], x$parameters$npoints)
      }
      x$modelling$ccs <- ccs
    }
    print("check1")
    ext <- c(x$parameters$xmn, x$parameters$xmx, x$parameters$ymn, 
             x$parameters$ymx)
    ext_eqearth <- eqearth_get_ext(ext)
    xy_ratio <- diff(ext_eqearth[1:2])/diff(ext_eqearth[3:4])
    y1 <- (height - length(climate) * y0)/length(climate)
    x1 <- min(c(width/3, xy_ratio * y1))
    y1 <- x1/xy_ratio
    x2 <- width - 2 * x1
    if (save) {
      opt_height <- round(length(climate) * (y0 + y1), 
                          3)
      if ((opt_height/height) >= 1.05 | (opt_height/height) <= 
          0.95) {
        cat("SUGGEST: Using height =", opt_height, "would get rid of all the white spaces.\n")
      }
      if (x1/width < 0.25) {
        opt_height <- round(length(climate) * (y0 + 
                                                 width/3/xy_ratio), 3)
        cat("SUGGEST: Using height =", opt_height, "would increase the width of the maps to a more optimal size.\n")
      }
    }
    if (save) {
      if (as.png) {
        grDevices::png(paste0(strsplit(filename, ".png")[[1]], 
                              ".png"), width = width, height = height, units = "in", 
                       res = png.res)
      }
      else {
        grDevices::pdf(filename, width = width, height = height)
      }
    }
    else {
      par_usr <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(par_usr))
    }
    distribs <- lapply(x$modelling$distributions, function(x) {
      if (is.data.frame(x)) {
        x[, 2] = resol * (x[, 2]%/%resol) + resol/2
        x[, 3] = resol * (x[, 3]%/%resol) + resol/2
        return(stats::aggregate(. ~ longitude + latitude + 
                                  taxonid, data = x, mean, na.action = NULL))
      }
      else {
        return(NA)
      }
    })
    veg_space <- do.call(rbind, distribs)[, c("longitude", 
                                              "latitude")]
    print("check2")
    veg_space <- plyr::count(veg_space)
    veg_space <- veg_space[!is.na(veg_space[, 1]), ]
    veg_space[, 3] <- base::log10(veg_space[, 3])
    veg_space <- raster::rasterFromXYZ(veg_space, crs = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))
    m3 <- rep(c(1:(2 * length(climate))), times = rep(c(1, 
                                                        3), times = length(climate)))
    m1 <- m3 + max(m3)
    m2 <- max(m1) + 1:(4 * length(climate))
    y3 <- y1 * 0.05
    y3 <- min(y0/2, y3)
    y3 <- max(y3, y0/3)
    graphics::layout(cbind(m1, m2, m3), width = c(x1, x2, 
                                                  x1), height = rep(c(y0, (y1 - y3)/2, y3, (y1 - y3)/2), 
                                                                    times = length(climate)))
    graphics::par(ps = 8 * 3/2)
    zlab = c(0, ceiling(max(raster::values(veg_space), na.rm = TRUE)))
    clab = c()
    i <- 0
    while (i <= max(zlab)) {
      clab <- c(clab, c(1, 2, 5) * 10^i)
      i <- i + 1
    }
    clab <- c(clab[log10(clab) <= max(raster::values(veg_space), 
                                      na.rm = TRUE)], clab[log10(clab) > max(raster::values(veg_space), 
                                                                             na.rm = TRUE)][1])
    print("check3")
    zlab[2] <- log10(clab[length(clab)])
    site_xy <- NA
    if (is.na(x$misc$site_info$long) | is.na(x$misc$site_info$lat)) 
      add_modern <- FALSE
    if (add_modern) {
      site_xy <- c(x$misc$site_info$long, x$misc$site_info$lat)
    }
    graphics::par(mar = c(0, 0, 0, 0), ps = 8 * 3/2)
    plot_map_eqearth(veg_space, ext, zlim = zlab, brks.pos = log10(clab), 
                     brks.lab = clab, col = viridis::plasma(20), title = "Number of unique species occurences", 
                     site_xy = site_xy, dim = c(x1 * width/sum(c(x1, 
                                                                 x2, x1)), height/length(climate)))
    ll <- do.call(rbind, distribs)
    ll <- ll[!is.na(ll[, 1]), ]
    ll_unique <- unique(ll[, colnames(ll) != "taxonid"])
    climate_space <- x$modelling$climate_space
    climate_space[, 1] = resol * (climate_space[, 1]%/%resol) + 
      resol/2
    climate_space[, 2] = resol * (climate_space[, 2]%/%resol) + 
      resol/2
    climate_space = stats::aggregate(. ~ longitude + latitude, 
                                     data = climate_space, mean)
    cs_colour <- rep(NA, nrow(climate_space))
    for (i in 1:length(cs_colour)) {
      w <- which(ll_unique[, 1] == climate_space[i, 1])
      cs_colour[i] <- ifelse(climate_space[i, 2] %in% 
                               ll_unique[w, 2], "black", "grey70")
    }
    if (length(climate) > 1) {
      oo <- order(cs_colour, decreasing = TRUE)
      for (clim in 1:(length(climate) - 1)) {
        miny <- min(climate_space[, climate[clim + 1]], 
                    na.rm = TRUE)
        maxy <- max(climate_space[, climate[clim + 1]], 
                    na.rm = TRUE)
        minx <- min(climate_space[, climate[clim]], 
                    na.rm = TRUE)
        maxx <- max(climate_space[, climate[clim]], 
                    na.rm = TRUE)
        dX <- maxx - minx
        dY <- maxy - miny
        minx <- minx - 0.03 * dX
        maxx <- maxx + 0.03 * dX
        miny <- miny - 0.03 * dY
        maxy <- maxy + 0.03 * dY
        xlim <- c(minx, maxx) + c(-0.005, 0.15) * dX
        ylim <- c(miny, maxy) + c(-0.1, 0.05^2) * dY
        graphics::par(mar = c(0, 0, 0, 0), ps = 8 * 
                        3/2)
        plot(NA, NA, type = "n", xlab = "", ylab = "", 
             main = "", axes = FALSE, frame = FALSE, xlim = xlim, 
             ylim = c(0, 1), xaxs = "i", yaxs = "i")
        graphics::text(mean(c(minx, maxx)), 0.25, paste(climate[clim], 
                                                        "(x-axis) vs.", climate[clim + 1], "(y-axis)"), 
                       adj = c(0.5, 0), cex = 1, font = 1)
        graphics::par(mar = c(0, 0.2, 0, 0.2))
        plot(NA, NA, type = "n", xaxs = "i", yaxs = "i", 
             axes = FALSE, frame = FALSE, xlim = xlim, 
             ylim = ylim)
        {
          for (yval in graphics::axTicks(4)) {
            if (yval >= miny & yval <= maxy) {
              graphics::segments(maxx, yval, minx, yval, 
                                 col = "grey90", lwd = 0.5)
              graphics::segments(maxx, yval, maxx - 
                                   diff(xlim) * 0.012, yval, lwd = 0.5)
              graphics::text(maxx + diff(xlim) * 0.015, 
                             yval, yval, cex = 6/8, adj = c(0, 0.4))
            }
          }
          for (xval in graphics::axTicks(1)) {
            if (xval >= minx & xval <= maxx) {
              graphics::segments(xval, miny, xval, maxy, 
                                 col = "grey90", lwd = 0.5)
              graphics::segments(xval, miny, xval, miny + 
                                   diff(ylim) * 0.012, lwd = 0.5)
              graphics::text(xval, miny - diff(ylim) * 
                               0.015, xval, cex = 6/8, adj = c(0.5, 
                                                               1))
            }
          }
          print("check4")
          graphics::rect(minx, miny, maxx, maxy, lwd = 0.5)
          graphics::points(climate_space[oo, climate[clim]], 
                           climate_space[oo, climate[clim + 1]], col = cs_colour[oo], 
                           pch = 20, cex = 0.5)
          if (add_modern) {
            if (is.numeric(x$misc$site_info$climate[, 
                                                    climate[clim]]) & is.numeric(x$misc$site_info$climate[, 
                                                                                                          climate[clim + 1]])) {
              graphics::points(x$misc$site_info$climate[, 
                                                        climate[clim]], x$misc$site_info$climate[, 
                                                                                                 climate[clim + 1]], pch = 23, cex = 2, 
                               lwd = 2, col = "white", bg = "red")
            }
          }
        }
      }
    }
    for (clim in climate) {
      brks <- c(x$modelling$ccs[[clim]]$k1, max(x$modelling$ccs[[clim]]$k1) + 
                  diff(x$modelling$ccs[[clim]]$k1[1:2]))
      R1 <- raster::rasterFromXYZ(cbind(climate_space[, 
                                                      1:2], climate_space[, clim]), crs = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))
      graphics::par(mar = c(0, 0, 0, 0), ps = 8 * 3/2)
      plot_map_eqearth(R1, ext, zlim = range(brks), col = viridis::viridis(length(brks) - 
                                                                             1), brks.pos = brks, brks.lab = brks, title = accClimateVariables(clim)[3], 
                       site_xy = site_xy, dim = c(x1 * width/sum(c(x1, 
                                                                   x2, x1)), height/length(climate)))
    }
    for (clim in climate) {
      h1 <- graphics::hist(climate_space[, clim], breaks = c(x$modelling$ccs[[clim]]$k1, 
                                                             max(x$modelling$ccs[[clim]]$k1) + diff(x$modelling$ccs[[clim]]$k1[1:2])), 
                           plot = FALSE)
      h2 <- graphics::hist(ll[, clim], breaks = c(x$modelling$ccs[[clim]]$k1, 
                                                  max(x$modelling$ccs[[clim]]$k1) + diff(x$modelling$ccs[[clim]]$k1[1:2])), 
                           plot = FALSE)
      xval <- range(h1$breaks)
      ext_factor_x <- max(graphics::strwidth(paste0("     ", 
                                                    c(h1$counts, h2$counts)), cex = 6/8, units = "inches")) + 
        graphics::strheight("Number of occurrences", 
                            cex = 6/8, units = "inches")
      w <- x2 - 2 * ext_factor_x
      ratio <- diff(xval)/w
      print("check5")
      xval <- xval + c(-1, 1) * (ext_factor_x * ratio)
      graphics::par(mar = c(0, 0, 0, 0), ps = 8 * 3/2)
      plot(NA, NA, type = "n", xlab = "", ylab = "", main = "", 
           axes = FALSE, frame = FALSE, xlim = c(0, 1), 
           ylim = c(0, 1), xaxs = "i", yaxs = "i")
      {
        s1 <- graphics::strwidth("Observed", cex = 1, 
                                 font = 2)
        s2 <- graphics::strwidth(" vs. ", cex = 1, font = 1)
        s3 <- graphics::strwidth("Sampled", cex = 1, 
                                 font = 2)
        graphics::text(0.5 - (s1 + s2 + s3)/2, 0.5, 
                       "Observed", cex = 1, font = 2, col = "grey70", 
                       adj = c(0, 0))
        graphics::text(0.5 - (s1 + s2 + s3)/2 + s1, 
                       0.5, "  vs.  ", cex = 1, font = 1, col = "black", 
                       adj = c(0, 0))
        graphics::text(0.5 - (s1 + s2 + s3)/2 + s1 + 
                         s2, 0.5, " Sampled", cex = 1, font = 2, col = "black", 
                       adj = c(0, 0))
        graphics::text(0.5, 0.3, paste(accClimateVariables(clim)[3], 
                                       " [", clim, "]", sep = ""), cex = 1, font = 1, 
                       col = "black", adj = c(0.5, 1))
      }
      graphics::par(mar = c(0, 0, 0.1, 0))
      opar <- graphics::par(lwd = 0.5)
      plot(NA, NA, type = "n", xlim = xval, ylim = c(0, 
                                                     max(h1$counts)), axes = FALSE, main = "", xaxs = "i", 
           yaxs = "i")
      {
        for (yval in graphics::axTicks(2)) {
          if (yval + graphics::strheight(yval, cex = 6/8)/2 < 
              max(h1$counts)) {
            graphics::text(h1$breaks[1] - diff(xval) * 
                             0.015, yval, yval, cex = 6/8, adj = c(1, 
                                                                   ifelse(yval == 0, 0, 0.4)))
          }
          graphics::segments(h1$breaks[1], yval, max(h1$breaks), 
                             yval, col = ifelse(yval == 0, "black", "grey90"))
          graphics::segments(h1$breaks[1], yval, h1$breaks[1] + 
                               diff(xval) * 0.012, yval)
        }
        graphics::segments(h1$breaks[1], 0, h1$breaks[1], 
                           max(h1$counts))
        plot(h1, add = TRUE, col = "grey70")
        graphics::text(xval[1], max(h1$counts)/2, "Number of occurrences", 
                       cex = 6/8, adj = c(0.5, 1), srt = 90)
      }
      graphics::par(opar)
      graphics::par(mar = c(0, 0, 0, 0))
      plot(NA, NA, type = "n", xlab = "", ylab = "", main = "", 
           axes = FALSE, frame = FALSE, xlim = xval, ylim = c(0, 
                                                              1), xaxs = "i", yaxs = "i")
      {
        d1 <- -9999999
        print("check6")
        if (add_modern) {
          if (is.numeric(x$misc$site_info$climate[, 
                                                  clim])) {
            graphics::points(x$misc$site_info$climate[, 
                                                      clim], 0.83, pch = 24, col = NA, bg = "red", 
                             cex = 0.9, lwd = 1.5)
            graphics::points(x$misc$site_info$climate[, 
                                                      clim], 0.17, pch = 25, col = NA, bg = "red", 
                             cex = 0.9, lwd = 1.5)
          }
        }
        for (i in 1:length(h1$breaks)) {
          d2 <- h1$breaks[i] - graphics::strwidth(paste0(" ", 
                                                         h1$breaks[i], " "), cex = 6/8, units = "user")/2
          if (d2 - d1 >= 0) {
            d1 <- h1$breaks[i] + graphics::strwidth(paste0(" ", 
                                                           h1$breaks[i], " "), cex = 6/8, units = "user")/2
            graphics::text(h1$breaks[i], 0.5, h1$breaks[i], 
                           cex = 6/8, adj = c(0.5, 0.5))
            graphics::segments(h1$breaks[i], 1, h1$breaks[i], 
                               0.9, lwd = 0.5)
            graphics::segments(h1$breaks[i], 0, h1$breaks[i], 
                               0.1, lwd = 0.5)
          }
        }
      }
      graphics::par(mar = c(0.1, 0, 0, 0))
      opar <- graphics::par(lwd = 0.5)
      plot(NA, NA, type = "n", xlim = xval, ylim = c(max(h2$counts), 
                                                     0), axes = FALSE, main = "", xaxs = "i", yaxs = "i")
      {
        M <- max(h2$breaks)
        for (yval in graphics::axTicks(4)) {
          if (yval - graphics::strheight(yval, cex = 6/8)/2 < 
              max(h2$counts)) {
            graphics::text(M + diff(xval) * 0.015, yval, 
                           yval, cex = 6/8, adj = c(0, ifelse(yval == 
                                                                0, 1, 0.6)))
          }
          graphics::segments(M, yval, M - diff(xval) * 
                               0.012, yval)
          graphics::segments(h2$breaks[1], yval, M, 
                             yval, col = ifelse(yval == 0, "black", "grey90"))
        }
        plot(h2, add = TRUE, col = "black", border = "grey70")
        graphics::segments(M, 0, M, max(h2$counts))
        graphics::segments(h1$breaks[1], 0, max(h1$breaks), 
                           0)
        graphics::par(opar)
      }
      graphics::text(xval[2], max(h2$counts)/2, "Number of occurrences", 
                     cex = 6/8, adj = c(0.5, 1), srt = -90)
    }
    if (save) {
      grDevices::dev.off()
    }
  }
  else {
    cat("This function only works with a crestObj.\n\n")
  }
  invisible()
}

