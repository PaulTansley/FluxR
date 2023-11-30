
#' @param site_name # site name as chacter string
#' @export ds_meld
ds_meld <- function(site_name){

  require(vroom)
  require(tidyr)
  require(readr)

  wd <- getwd()
  setwd("..")
  wd1 <- getwd()
  setwd(wd)


  mydir <- paste0(wd1, "/FluxR/")
  if (!dir.exists(mydir)) {
    dir.create(mydir)
  }
  mydir <- paste0(mydir, site_name)
  if (!dir.exists(mydir)) {
    dir.create(mydir)}

  mydir <- paste0(mydir, "/Datasheet")
  if (!dir.exists(mydir)) {
    dir.create(mydir)}

  ls <- list.files(pattern = "ds", full.names = T, recursive = T)

  ds <- vroom(ls, col_types = c(Rep = "c"))

  ds_path <<- paste0(mydir, "/", site_name, "_ds_melded.csv")

  write_csv(ds, ds_path)
}
#' @param directory
#' @param site_name # site name as chacter string
#' @param hr # hour as numeric
#' @param dy # day as numeric
#' @export lgr_meld
lgr_meld <- function(directory, site_name, hr = 0, dy = 0) {
  require(vroom)
  require(readr)
  require(lubridate)

  message("Removing empty files")

  wd <- getwd()
  setwd("..")
  wd1 <- getwd()
  setwd(wd)

  mydir <- paste0(wd1, "/FluxR/")
  if (!dir.exists(mydir)) {
    dir.create(mydir)
  }
  mydir <- paste0(mydir, site_name)
  if (!dir.exists(mydir)) {
    dir.create(mydir)}

  mydir <- paste0(mydir, "/Data")
  if (!dir.exists(mydir)) {
    dir.create(mydir)}

  l<- list.files(directory, pattern = ".txt", full.names = T, recursive = T)
  for (i in l) {
    if (file.size(i) == 0) {unlink(i)}
  }

  message("Melding lgr data...")
  l <- list.files(directory, full.names = T)

  l <-
    list.files(
      directory,
      pattern = ".txt",
      full.names = T,
      recursive = T
    )
  b <-  suppressMessages(vroom(l, skip = 1, delim = ","))
  l <- list.files(directory, pattern = ".txt")
  b$Time <- dmy_hms(b$Time)

  lgr_path <<- paste0(mydir, "/",  site_name,  "_melded.csv")

  write_csv(b, lgr_path)

  message("Files Melded.")
}
#' @param x # directory path
#' @param ch4 # is ch4 boolean
#' @export check
check <- function(x, ch4) {

  if (ch4 == T) {
    flux_edit_1 <- flux_edit_ch4 %>%
      filter(id == i)

    linear <- lm(mass_ch4_umol ~ t,
                 data = flux_edit_1, na.action = NULL)
    lin_sum <- summary(linear)
    slope <- lin_sum$coefficients[2]
    r2 <- lin_sum$r.squared
    p <- lin_sum$coefficients[8]
  }
  else{
    flux_edit_1 <- flux_edit_co2 %>%
      filter(id == i)
    linear <- lm(mass_co2_umol ~ t,
                 data = flux_edit_1, na.action = na.exclude)
    lin_sum <- summary(linear)
    slope <- lin_sum$coefficients[2]
    r2 <- lin_sum$r.squared
    p <- lin_sum$coefficients[8]

  }
  conf1 <- 3
  assign("conf1", conf1, .GlobalEnv)
  if (p >= 0.5 | between(r2, -0.7, 0.7)) {
    print(paste0(i, ": Check"))
    if (ch4 == T) {
      plot(
        flux_edit_1$t,
        flux_edit_1$mass_ch4_umol,
        xlab = "",
        main = i,
        sub = paste0("\nslope = ", slope, "\n", "r2 = ", r2,
                     "\n", "p = ", p)
      )
    }
    else{
      plot(
        flux_edit_1$t,
        flux_edit_1$mass_co2_umol,
        xlab = "",
        main = i,
        sub = paste0("\nslope = ", slope, "\n", "r2 = ", r2,
                     "\n", "p = ", p)
      )
    }
    conf_1 <- readline(prompt = "Choose options:
    1:Keep
    2:Discard")
    assign("conf1", conf1, .GlobalEnv)
  }
  else{
    fl1 <- flux_edit_1
    conf_1 <<- 3
    assign("conf1", conf1, .GlobalEnv)
    print(paste0(i, ": Plot within statistical boundaries"))
  }
  if (conf_1 == 1) {
    fl1 <- flux_edit_1
  }
  else if (conf_1 == 2) {
    f <- flux_edit_1
    f[5:58] <- NA
    fl1 <- f
  }

  flux_edit_2 <<- rbind(flux_edit_2, fl1)
  assign("flux_edit_2",  flux_edit_2, .GlobalEnv)
  if (ch4 == T) {
    assign("flux_edit_ch4_1", flux_edit_2, .GlobalEnv)
  }
  else{
    assign("flux_edit_co2_1", flux_edit_2, .GlobalEnv)
  }

}
#' @param lgr_data # lgr data path
#' @param start_times # Directory path to start times
#' @param site_name # site name as chacter string
#' @param chamber # boolean
#' @param chamb_volume # numeric
#' @param chamb_area # numeric
#' @export lgr_flux
lgr_flux <- function(lgr_data,
                     start_times,
                     site_name,
                     chamber = T,
                     chamb_volume,
                     chamb_area) {

  require(lubridate)
  require(ggpubr)
  require(nlme)
  require(qdapRegex)
  require(magrittr)
  require(dplyr)
  require(crayon)
  require(tibble)


  grph <- function(x, ch4) {
    require(crayon)
    require(dplyr)

    ch4 <<- ch4

    #check linear regression coefficients
    linear <- lm(mass_ch4_umol ~ t,
                 data = lgr1)

    lin_sum <- summary(linear)

    slope <- lin_sum$coefficients[2]
    r2 <- lin_sum$r.squared
    p <- lin_sum$coefficients[8]

    #plot and identify
    if (ch4 == T) {
      win.graph()
      plot(
        lgr1$t,
        lgr1$mass_ch4_umol,
        xlab = "",
        main = lgr1$id[1],
        sub = paste0("\nslope = ", slope, "\n", "r2 = ", r2,
                     "\n", "p = ", p)
      )
      a <<-
        identify(
          lgr1$t,
          lgr1$mass_ch4_umol,
          atpen = T,
          n = 2,
          offset = 15
        )
      invisible(dev.off())
      confirm()
    }

    else{
      win.graph()
      plot(
        lgr1$t,
        lgr1$mass_co2_umol,
        xlab = "",
        main = lgr1$id[1],
        sub = paste0("\nslope = ", slope, "\n", "r2 = ", r2,
                     "\n", "p = ", p)
      )
      a <<-
        identify(
          lgr1$t,
          lgr1$mass_co2_umol,
          atpen = T,
          n = 2,
          offset = 15
        )
      invisible(dev.off())
      confirm()
    }
  }
  confirm <- function() {
    require(magrittr)
    require(dplyr)
    require(crayon)

    retry <- "a"
    keepog <- "a"

    conf <<- readline(prompt = "Choose options:
    1:Keep Changes
    2:Keep Original
    3:Retry
    4:Discard Flux\n")

    if (conf == 1) {
      fl <- lgr1 %>%
        filter(t > min(a) & t < max(a))
      assign("conf", conf, .GlobalEnv)
    }
    else if (conf == 2) {
      fl <- lgr1
      assign("conf", conf, .GlobalEnv)
    }

    else if(conf ==4){
      f <- lgr1
      f[5:58] <- NA
      fl <- f
      assign("conf", conf, .GlobalEnv)
    }
    else{
      fl <- lgr1
      conf <- 3
      assign("conf", conf, .GlobalEnv)
    }
    flux_edit <<- rbind(flux_edit, fl)


    if (ch4 == T) {
      assign("flux_edit_ch4", flux_edit, .GlobalEnv)
    }
    else{
      assign("flux_edit_co2", flux_edit, .GlobalEnv)

    }

  }
  check <- function(x, ch4) {

    if (ch4 == T) {
      flux_edit_1 <- flux_edit_ch4 %>%
        filter(id == i)

      linear <- lm(mass_ch4_umol ~ t,
                   data = flux_edit_1, na.action = NULL)
      lin_sum <- summary(linear)
      slope <- lin_sum$coefficients[2]
      r2 <- lin_sum$r.squared
      p <- lin_sum$coefficients[8]
    }
    else{
      flux_edit_1 <- flux_edit_co2 %>%
        filter(id == i)
      linear <- lm(mass_co2_umol ~ t,
                   data = flux_edit_1, na.action = na.exclude)
      lin_sum <- summary(linear)
      slope <- lin_sum$coefficients[2]
      r2 <- lin_sum$r.squared
      p <- lin_sum$coefficients[8]

    }
    conf1 <- 3
    assign("conf1", conf1, .GlobalEnv)
    if (p >= 0.5 | between(r2, -0.7, 0.7)) {
      print(paste0(i, ": Check"))
      if (ch4 == T) {
        plot(
          flux_edit_1$t,
          flux_edit_1$mass_ch4_umol,
          xlab = "",
          main = i,
          sub = paste0("\nslope = ", slope, "\n", "r2 = ", r2,
                       "\n", "p = ", p)
        )
      }
      else{
        plot(
          flux_edit_1$t,
          flux_edit_1$mass_co2_umol,
          xlab = "",
          main = i,
          sub = paste0("\nslope = ", slope, "\n", "r2 = ", r2,
                       "\n", "p = ", p)
        )
      }
      conf_1 <- readline(prompt = "Choose options:
    1:Keep
    2:Discard")
      assign("conf1", conf1, .GlobalEnv)
    }
    else{
      fl1 <- flux_edit_1
      conf_1 <<- 3
      assign("conf1", conf1, .GlobalEnv)
      print(paste0(i, ": Plot within statistical boundaries"))
    }
    if (conf_1 == 1) {
      fl1 <- flux_edit_1
    }
    else if (conf_1 == 2) {
      f <- flux_edit_1
      f[5:58] <- NA
      fl1 <- f
    }

    flux_edit_2 <<- rbind(flux_edit_2, fl1)
    assign("flux_edit_2",  flux_edit_2, .GlobalEnv)
    if (ch4 == T) {
      assign("flux_edit_ch4_1", flux_edit_2, .GlobalEnv)
    }
    else{
      assign("flux_edit_co2_1", flux_edit_2, .GlobalEnv)
    }

  }




  #`Read in data, presumes if .txt then a column needs to be skipped as data untreated,
  #if .csv presumes data has been generated through the lgr_meld function and doesn't skip a column
  if (grepl(".txt", lgr_data)) {
    ghg <-
      suppressMessages(read_csv(lgr_data, skip = 1, col_types = cols()))
  } else{
    ghg <- suppressMessages(read_csv(lgr_data, col_types = cols()))
  }

  b <-  separate(ghg,
                 col = Time,
                 into = c("date", "time"),
                 sep = " ")


  dates <- gsub("/", "-", unique(b$date))


  #Create directories -- must be a better way to deal with recurisive directory creation...

  mydir <- "Fluxr/"
  if (!dir.exists(mydir)) {
    dir.create(mydir)
  }
  mydir <- paste0("Fluxr/", site_name)
  if (!dir.exists(mydir)) {
    dir.create(mydir)
  }
  for (i in dates) {
    mydir <- paste0("Fluxr/", site_name, "/", i)
    if (!dir.exists(mydir)) {
      dir.create(mydir)
    }
  }

  if(start_times == "FULL DATA.csv"){
  time1 <- suppressMessages(read_csv(start_times)[-1,])
  time1$date <- gsub("/", "-", time1$date)
  }else{
    time1 <- suppressMessages(read_csv(start_times)) %>%
      rowid_to_column("Sample_ID")
    time1 <- time1[1:10]
    names(time1) <- c("Sample_ID", "site", "plot", "date", "wtd", "rep_nb",
                      "start_temp", "starttime_lgr", "end_temp", "PAR")

  }


  #Define save file names


  #########Start time loop
  for (i in dates) {

    if(exists("flux_edit")){
    rm(list = c("flux_edit", "flux_edit_2", "flux_edit_ch4",
                "flux_edit_ch4_1", "flux_edit_co2", "flux_edit_co2_1", "a",
                "ch4", "conf", "conf1", "lsf"))}


    mydir <- "Fluxr/"

    ff <- i

    csv_lgr <-
      paste0(mydir,
             site_name,
             "/",
             ff,
             "/",
             site_name,
             "_",
             ff,
             "_extracted.csv")
    csv_ch4_flux <-
      paste0(mydir,
             site_name,
             "/",
             ff,
             "/",
             site_name,
             "_",
             ff,
             "_ch4_flux.csv")

    csv_co2_flux <-
      paste0(mydir,
             site_name,
             "/",
             ff,
             "/",
             site_name,
             "_",
             ff,
             "_co2_flux.csv")

    grph_nm <-
      paste0(mydir,
             "/",
             site_name,
             "/",
             ff,
             "/",
             site_name,
             "_",
             ff)
    grph_f <-
      paste0(mydir,
             "/",
             site_name,
             "/",
             ff,
             "/",
             site_name,
             "_",
             ff,
             "_flux.png")


    #Load Data
    if (grepl(".txt", lgr_data)) {
      ghg <- read_csv(lgr_data, skip = 1, col_types = cols())
    } else{
      ghg <- read_csv(lgr_data, col_types = cols())
    }


    b <-  separate(ghg,
                   col = Time,
                   into = c("date", "time"),
                   sep = " ")

    cat(red(paste0(
      "\nCalculating for ", site_name, " ", i, "..."
    )))

    if(start_times == "FULL DATA.csv"){
    time <- time1 %>%
      filter(site == site_name & date == i)}
    else{time <- time1 %>%
      filter(date == i)}



    #set up time mapping from start times file
    times <- time %>%
      filter(site == site_name) %>%
      select(sample_ID, date, starttime_lgr, endtime_lgr, site) %>%
      rename(id = sample_ID)

    starts <- hms(paste0(times$starttime_lgr, ":00"))
    t1 <- hms(paste0(times$endtime_lgr, ":00"))

    mynames <- times$id

    #Seperate date time column
    ghg1 <- ghg %>%
      mutate(date.time = Time) %>%
      separate(
        col = Time,
        into = c("date", "time"),
        sep = " ",
        remove = F
      ) %>%
      mutate(date1 = date)

    ghg1$date1 <- gsub("/", "-", ghg1$date1)
    ghg1 <- filter(ghg1, date1 == i)
    date <- ghg1$date[1]

    #remove X from names
    names(ghg1) <- gsub("[][]", replacement = "", x = names(ghg1))

    #hms time
    ghg1$time <- hms(ghg1$time)

    #select time periods needed by group
    message("\nMapping Times...")
    ghg_out <-
      map2(starts, t1, ~ filter(ghg1, time >= .x & time <= .y))

    ghg_out <- set_names(ghg_out, mynames)

    ghg <- map_df(ghg_out, ~ as.data.frame(.x), .id = "id")



    veg_sml <- c(0.002097341, 0.010029)
    veg_lrg <- c(0.022166568, 0.064692)
    ditch <- c(0.032951432, 0.098980)



    ghg_renew <- data.frame()
    for (f in unique(ghg$id)) {
      idchmb <- filter(time1, sample_ID == f)

      chamber_type <- idchmb$chamber[1]

      ghg_f <- ghg %>%
        filter(id == f)

      if (chamber == T) {
        if (chamber_type == "veg_sml") {
          chmbr <- veg_sml
          ghg_f$chamber_volume <- chmbr[1]
          ghg_f$chamber_area  <- chmbr[2]
          chamb_volume <- chmbr[1]
          chamb_area  <- chmbr[2]
        }
        else if (chamber_type == "veg_lrg") {
          chmbr <- veg_lrg
          ghg_f$chamber_volume <- chmbr[1]
          ghg_f$chamber_area  <- chmbr[2]
          chamb_volume <- chmbr[1]
          chamb_area  <- chmbr[2]
        }
        else if (chamber_type == "ditch") {
          chmbr <- ditch
          ghg_f$chamber_volume <- chmbr[1]
          ghg_f$chamber_area  <- chmbr[2]
          chamb_volume <- chmbr[1]
          chamb_area  <- chmbr[2]
        }
      }

      else{
        ghg_f$chamber_volume <- chamb_volume
      ghg_f$chamber_area <- chamb_area
      }


      ghg_renew <- rbind(ghg_renew, ghg_f)
    }


    #Calculate and graph concentrations and flux

    #Standard pressure(kpa)
    sp <- 101.325
    #Standard temp(k)
    st <- 273.15
    # vol 1 mol/m3 ideal gas at stp
    gas <- 0.022414

    #Convert and calculate c to k, torr to kpa, ppm to umol, and time in seconds

    lgr <- ghg_renew

    lgr <- lgr %>%
      group_by(id) %>%
      mutate(
        GasP_kpa = GasP_torr / 7.5006157818041,
        Amb.t.k = AmbT_C + 273.15,
        time = hms(time),
        date.time = dmy_hms(date.time),
        vol_co2 = (CO2d_ppm / 1000000) * chamber_volume,
        vol_co2_stp = vol_co2 * (sp / Amb.t.k) * (st / sp),
        mol_co2 = (vol_co2_stp / gas),
        mass_co2_g = mol_co2 * 44.0095,
        mass_co2_umol = mol_co2 * 1000000,
        vol_ch4 = (CH4d_ppm / 1000000) * chamber_volume,
        vol_ch4_stp = vol_ch4 * (sp / Amb.t.k) * (st / sp),
        mol_ch4 = (vol_ch4_stp / gas),
        mass_ch4_g = mol_ch4 * 44.0095,
        mass_ch4_umol = mol_ch4 * 1000000,
        t = as.numeric(time - lag(time)),
        t = replace_na(t, 0),
        t = cumsum(t)
      )

    lgr$site <- site_name

    #Edit and check fluxes

    flux_edit <<- data.frame()
    for (i in  unique(lgr$id)) {
      lgr1 <- lgr %>%
        filter(id == i)
      repeat {
        grph(ch4 = T)
        if (!conf == "3") {
          break
        }
      }
    }

    flux_edit_2 <<-  data.frame()

    for (i in unique(flux_edit_ch4$id)) {
      f <- flux_edit_ch4 %>%
        filter(id == i)
      if (anyNA(f$CH4_ppm))
      {
        next
      }
      tryCatch({
        check(i, ch4 = T)
      },
      error = function(e) {
        print(paste0(i, ": Flux stats are too poor, Flux will be added as NA"))
      })
    }

    chk_ids <- unique(flux_edit_ch4$id)
    chk_ids_2 <- unique(flux_edit_ch4_1$id)
    dif <- setdiff(chk_ids, chk_ids_2)
    if (!is_empty(dif)) {
      for (h in dif) {
        add <- flux_edit_ch4 %>%
          filter(id == h)
        add[5:58] <- NA
        flux_edit_ch4_1 <- rbind(flux_edit_ch4_1, add)
        assign("flux_edit_ch4_1", flux_edit_ch4_1, .GlobalEnv)
      }
    }


    flux_edit <<- data.frame()
    for (i in  unique(lgr$id)) {
      lgr1 <- lgr %>%
        filter(id == i)
      repeat {
        grph(ch4 = F)
        if (!conf == "3") {
          break
        }
      }
    }



    flux_edit_2 <<-  data.frame()
    for (i in unique(flux_edit_co2$id)) {
      f <- flux_edit_co2 %>%
        filter(id == i)
      if (anyNA(f$CO2_ppm))
      {
        next
      }
      tryCatch({
        check(i, ch4 = F)
      },
      error = function(e) {
        print(paste0(i, ": Flux stats are too poor, Flux will be added as NA"))
      })
    }

    chk_ids <- unique(flux_edit_co2$id)
    chk_ids_2 <- unique(flux_edit_co2_1$id)
    dif <- setdiff(chk_ids, chk_ids_2)
    dif
    if (!is_empty(dif)) {
      for (h in dif) {
        add <- flux_edit_co2 %>%
          filter(id == h)
        add[5:58] <- NA
        flux_edit_co2_1 <- rbind(flux_edit_co2_1, add)
        assign("flux_edit_co2_1", flux_edit_co2_1, .GlobalEnv)
      }
    }

    #calculate regression coefficient and flux
    names_co2 <-
      c("id",
        "co2_intercept",
        "co2_slope",
        "co2_r2",
        "co2_flux")
    names_ch4 <-
      c("id",
        "ch4_intercept",
        "ch4_slope",
        "ch4_r2",
        "ch4_flux")

    # Function to Extract coefs and rsquare from lmlist into df
    sumfun <- function(x)
    {
      aux <- function(x)
        c(coef(x), summary(x)$r.squared)
      t(sapply(x, aux))
    }

    #Generate Orginal
    co2_coefs <- lmList(mass_co2_g ~ t |
                          id,
                        data = flux_edit_co2_1,
                        na.action = na.omit)

    co2_coefs <- as.data.frame(sumfun(co2_coefs)) %>%
      mutate(co2_flux =  t / chamb_area) %>%
      rownames_to_column(var =  "id")

    names(co2_coefs) <- names_co2

    ch4_coefs <- lmList(mass_ch4_g ~ t |
                          id,
                        data = flux_edit_ch4_1,
                        na.action = na.omit)

    ch4_coefs <- as.data.frame(sumfun(ch4_coefs)) %>%
      mutate(ch4_flux = t / chamb_area) %>%
      rownames_to_column(var =  "id")
    names(ch4_coefs) <- names_ch4

    flux_ch4 <- ch4_coefs
    flux_ch4$date <- date
    flux_ch4$site <- site_name
    flux_ch4 <- flux_ch4 %>%
      relocate(date, site, id)

    flux_co2 <- co2_coefs
    flux_co2$date <- date
    flux_co2$site <- site_name
    flux_co2 <- flux_co2 %>%
      relocate(date, site, id)

    write_csv(flux_ch4, csv_ch4_flux)
    write_csv(flux_co2, csv_co2_flux)

    #Graph concentrations
    message("\nMapping CO2 Regressions...")
    co2_conc <-  flux_edit_co2_1 %>%  ggscatter(
      x = "t",
      y = "mass_co2_g",
      add = "reg.line",
      # Add regressin line
      add.params = list(color = "blue", fill = "lightgray"),
      # Customize reg. line
      conf.int = TRUE,
      # Add confidence interval
      cor.coef = TRUE,
      # Add correlation coefficient
      cor.coeff.args = list(method = "pearson", label.sep = "\n"),
      scales = "free",
      facet.by = "id"
    ) +
      ggtitle("CO2 Concentration Edited")
    suppressMessages(ggsave(
      co2_conc,
      width = 45,
      height = 45,
      units = "cm",
      filename = paste0(grph_nm,"_grph_co2_edited.png"
    )))

    co2_conc_org <-  lgr %>%  ggscatter(
      x = "t",
      y = "mass_co2_g",
      add = "reg.line",
      # Add regressin line
      add.params = list(color = "blue", fill = "lightgray"),
      # Customize reg. line
      conf.int = TRUE,
      # Add confidence interval
      cor.coef = TRUE,
      # Add correlation coefficient
      cor.coeff.args = list(method = "pearson", label.sep = "\n"),
      scales = "free",
      facet.by = "id"
    ) +
      ggtitle("CO2 Concentration Original")
    suppressMessages(ggsave(
      co2_conc_org,
      width = 45,
      height = 45,
      units = "cm",
      filename = paste0(grph_nm,"_grph_co2_original.png")
    ))
    message("
Mapping CH4 Regressions...")
    ch4_conc<- flux_edit_ch4_1 %>% ggscatter(
      x = "t",
      y = "mass_ch4_g",
      add = "reg.line",
      # Add regressin line
      add.params = list(color = "blue", fill = "lightgray"),
      # Customize reg. line
      conf.int = TRUE,
      # Add confidence interval
      cor.coef = TRUE,
      # Add correlation coefficient. see ?stat_cor
      cor.coeff.args = list(method = "pearson", label.sep = "\n"),
      scales = "free",
      facet.by = "id"
    ) +
      ggtitle("CH4 Concentration")
    suppressMessages(ggsave(
      ch4_conc,
      width = 45,
      height = 45,
      units = "cm",
      filename = paste0(grph_nm, "_grph_ch4_edited.png")
    ))

    ch4_conc_org <- lgr %>% ggscatter(
      x = "t",
      y = "mass_ch4_g",
      add = "reg.line",
      # Add regressin line
      add.params = list(color = "blue", fill = "lightgray"),
      # Customize reg. line
      conf.int = TRUE,
      # Add confidence interval
      cor.coef = TRUE,
      # Add correlation coefficient. see ?stat_cor
      cor.coeff.args = list(method = "pearson", label.sep = "\n"),
      scales = "free",
      facet.by = "id"
    ) +
      ggtitle("CH4 Concentration")
    suppressMessages(ggsave(
      ch4_conc_org,
      width = 45,
      height = 45,
      units = "cm",
      filename = paste0(grph_nm, "_grph_ch4_original.png")
    ))


    cat(red(bold(paste0("\nFluxes for ", ff, " complete.\n"))))
    Sys.sleep(2)
    if(exists("flux_edit")){
    rm(list = c("flux_edit", "flux_edit_2", "flux_edit_ch4",
               "flux_edit_ch4_1", "flux_edit_co2", "flux_edit_co2_1", "a",
               "ch4", "conf", "conf1", "lsf"))}

  }
  cat(red(bold("\nAnalysis Complete\n")))

  if(exists("flux_edit")){
    rm(list = c("flux_edit", "flux_edit_2", "flux_edit_ch4",
                "flux_edit_ch4_1", "flux_edit_co2", "flux_edit_co2_1", "a",
                "ch4", "conf", "conf1", "lsf"))}

}
#' @param lgr_data # lgr data path
#' @param start_times # Directory path to start times
#' @param site_name # site name as chacter string
#' @param chamber # boolean
#' @param chamb_volume # numeric
#' @param chamb_area # numeric
#' @param flux_duration # numeric
#' @export lgr_flux2
lgr_flux2 <- function(lgr_data,
                     start_times,
                     site_name,
                     chamber_p = T,
                     chamb_volume,
                     chamb_area,
                     flux_duration = 3) {

  require(lubridate)
  require(ggpubr)
  require(nlme)
  require(qdapRegex)
  require(magrittr)
  require(dplyr)
  require(crayon)
  require(stringr)
  require(purrr)
  require(tidyr)
  require(readr)
  require(tibble)
  cat(red("\n Initialising Functions and Reading Data \n"))

  grph <- function(x, ch4) {
    require(crayon)
    require(dplyr)

    ch4 <<- ch4

    #check linear regression coefficients
    linear <- lm(mass_ch4_umol ~ t,
                 data = lgr1)

    lin_sum <- summary(linear)

    slope <- lin_sum$coefficients[2]
    r2 <- lin_sum$r.squared
    p <- lin_sum$coefficients[8]

    #plot and identify
    if (ch4 == T) {
      win.graph()
      plot(
        lgr1$t,
        lgr1$mass_ch4_umol,
        xlab = "",
        main = lgr1$id[1],
        sub = paste0("\nslope = ", slope, "\n", "r2 = ", r2,
                     "\n", "p = ", p)
      )
      a <-
        identify(
          lgr1$t,
          lgr1$mass_ch4_umol,
          atpen = T,
          n = 2,
          offset = 15
        )
      invisible(dev.off())
      assign("a", a, envir =  .GlobalEnv)
      confirm()

    }

    else{
      win.graph()
      plot(
        lgr1$t,
        lgr1$mass_co2_umol,
        xlab = "",
        main = lgr1$id[1],
        sub = paste0("\nslope = ", slope, "\n", "r2 = ", r2,
                     "\n", "p = ", p)
      )
      a <-
        identify(
          lgr1$t,
          lgr1$mass_co2_umol,
          atpen = T,
          n = 2,
          offset = 15
        )
      invisible(dev.off())
      assign("a", a, envir =  .GlobalEnv)
      confirm()

    }
  }
  confirm <- function() {
    require(magrittr)
    require(dplyr)
    require(crayon)

    retry <- "a"
    keepog <- "a"

    conf <<- readline(prompt = "Choose options:
    1:Keep Changes
    2:Keep Original
    3:Retry
    4:Discard Flux\n")

    if (conf == 1) {
      fl <- lgr1 %>%
        filter(t > min(a) & t < max(a))
      assign("conf", conf, .GlobalEnv)
    }
    else if (conf == 2) {
      fl <- lgr1
      assign("conf", conf, .GlobalEnv)
    }

    else if(conf ==4){
      f <- lgr1
      f[5:58] <- NA
      fl <- f
      assign("conf", conf, .GlobalEnv)
    }
    else{
      fl <- lgr1
      conf <- 3
      assign("conf", conf, .GlobalEnv)
    }
    flux_edit <<- rbind(flux_edit, fl)


    if (ch4 == T) {
      assign("flux_edit_ch4", flux_edit, .GlobalEnv)
      write_csv(flux_edit_ch4, csv_ch4_flux, append = T, col_names = T)

    }
    else{
      assign("flux_edit_co2", flux_edit, .GlobalEnv)
      write_csv(flux_edit_co2, csv_co2_flux, append = T, col_names = T)

    }

  }
  check <- function(x, ch4) {

    if (ch4 == T) {
      flux_edit_1 <- flux_edit_ch4 %>%
        filter(id == i)

      linear <- lm(mass_ch4_umol ~ t,
                   data = flux_edit_1, na.action = NULL)
      lin_sum <- summary(linear)
      slope <- lin_sum$coefficients[2]
      r2 <- lin_sum$r.squared
      p <- lin_sum$coefficients[8]
    }
    else{
      flux_edit_1 <- flux_edit_co2 %>%
        filter(id == i)
      linear <- lm(mass_co2_umol ~ t,
                   data = flux_edit_1, na.action = na.exclude)
      lin_sum <- summary(linear)
      slope <- lin_sum$coefficients[2]
      r2 <- lin_sum$r.squared
      p <- lin_sum$coefficients[8]

    }
    conf1 <- 3
    assign("conf1", conf1, .GlobalEnv)
    if (p >= 0.5 | between(r2, -0.7, 0.7)) {
      print(paste0(i, ": Check"))
      if (ch4 == T) {
        plot(
          flux_edit_1$t,
          flux_edit_1$mass_ch4_umol,
          xlab = "",
          main = i,
          sub = paste0("\nslope = ", slope, "\n", "r2 = ", r2,
                       "\n", "p = ", p)
        )
      }
      else{
        plot(
          flux_edit_1$t,
          flux_edit_1$mass_co2_umol,
          xlab = "",
          main = i,
          sub = paste0("\nslope = ", slope, "\n", "r2 = ", r2,
                       "\n", "p = ", p)
        )
      }
      conf_1 <- readline(prompt = "Choose options:
    1:Keep
    2:Discard")
      assign("conf1", conf1, .GlobalEnv)
    }
    else{
      fl1 <- flux_edit_1
      conf_1 <<- 3
      assign("conf1", conf1, .GlobalEnv)
      print(paste0(i, ": Plot within statistical boundaries"))
    }
    if (conf_1 == 1) {
      fl1 <- flux_edit_1
    }
    else if (conf_1 == 2) {
      f <- flux_edit_1
      f[5:58] <- NA
      fl1 <- f
    }

    flux_edit_2 <<- rbind(flux_edit_2, fl1)
    assign("flux_edit_2",  flux_edit_2, .GlobalEnv)
    if (ch4 == T) {
      assign("flux_edit_ch4_1", flux_edit_2, .GlobalEnv)
    }
    else{
      assign("flux_edit_co2_1", flux_edit_2, .GlobalEnv)
    }

  }




  #`Read in data, presumes if .txt then a column needs to be skipped as data untreated,
  #if .csv presumes data has been generated through the lgr_meld function and doesn't skip a column
  if (grepl(".txt", lgr_data)) {
    ghg <-
      suppressMessages(read_csv(lgr_data, skip = 1, col_types = cols()))
  } else{
    ghg <- suppressMessages(read_csv(lgr_data, col_types = cols()))
  }

  b <-  separate(ghg,
                 col = Time,
                 into = c("date", "time"),
                 sep = " ")

  if(start_times == "FULL DATA.csv"){
    time1 <- suppressMessages(read_csv(start_times)[-1,])
    time1$date <- gsub("/", "-", time1$date)
  } else{
    time1 <- suppressMessages(read_csv(start_times, col_types = c(Rep = "c"))) %>%
      drop_na(`Start time`) %>%
      rowid_to_column("sample_ID")


    names(time1) <- c("sample_ID", "site", "plot", "date", "wtd", "rep_nb",
                      "starttime_lgr", "start_temp", "end_temp", "PAR",
                      "Soil temp", "Soil moisture","EC","Humidity",
                      "chamber")

  }

  time1$date <- parse_date_time(time1$date, c("dmy", "ymd"))

  dates <- gsub("/", "-", unique(time1$date))


  #Create directories -- must be a better way to deal with recurisive directory creation...


  wd <- getwd()
  setwd("..")
  wd1 <- getwd()
  setwd(wd)


  mydir <- paste0(wd1, "/FluxR/")
  if (!dir.exists(mydir)) {
    dir.create(mydir)}

  mydir <- paste0(mydir, site_name)
  if (!dir.exists(mydir)) {
    dir.create(mydir)
  }

  mydir1 <- paste0(mydir, "/Processed")
  if(!dir.exists(mydir1)){
    dir.create(mydir1)
  }

  for (i in dates) {
    mydir <- paste0(mydir1, "/", i)
    if (!dir.exists(mydir)) {
      dir.create(mydir)
    }
  }



  #Define save file names


  #########Start time loop
  for (i in dates) {

    if(exists("flux_edit")){
      rm(list = c("flux_edit", "flux_edit_2", "flux_edit_ch4",
                  "flux_edit_ch4_1", "flux_edit_co2", "flux_edit_co2_1", "a",
                  "ch4", "conf", "conf1", "lsf"))}


     mydir <- paste0(wd1, "/FluxR/")

     ff <- i

     csv_lgr <- paste0(mydir1, "/", ff, "/", site_name, "_", ff, "_extracted.csv")

     csv_ch4_flux <- paste0(mydir1, "/", ff, "/", site_name, "_",  ff, "_ch4_flux.csv")

     csv_co2_flux <- paste0(mydir1, "/", ff, "/", site_name,"_",  ff,"_co2_flux.csv")

     grph_nm <- paste0(mydir1, "/", ff, "/",site_name, "_", ff)

     grph_f <- paste0(mydir1,"/", ff, "/",site_name, "_", ff, "_flux.png")


    #Load Data
    if (grepl(".txt", lgr_data)) {
      ghg <- read_csv(lgr_data, skip = 1, col_types = cols())
    } else{
      ghg <- read_csv(lgr_data, col_types = cols())
    }

       cat(red(paste0(
      "\nCalculating for ", site_name, " ", i, "..."
    )))


    i <- parse_date_time(i, c("dmy", "ymd"))

     time <- time1 %>%
        filter(site == site_name & date == i)



    #set up time mapping from start times file
      if(start_times == "FULL DATA.csv"){
    times <- time %>%
      filter(site == site_name) %>%
      select(sample_ID, date, starttime_lgr, endtime_lgr, site) %>%
      rename(id = sample_ID)
      }else{
        times <- time %>%
          filter(site == site_name) %>%           #rename(id = sample_ID) %>%
          mutate(id = paste0(plot,"-", rep_nb, "-", date))


        g <- parse_date_time(times$starttime_lgr, c("HM", "HMS"))
      times$starttime_lgr <- hms(str_sub(g, 11))
      times <- times %>%
        unite("time", date:starttime_lgr, remove = F)
      times$starttime_lgr <- hms(times$starttime_lgr)
      times$endtime_lgr <- times$starttime_lgr + hms(paste0("00:", flux_duration, ":00"))
      }


      starts <- hms(times$starttime_lgr)
      t1 <- hms(paste0(times$endtime_lgr))

    mynames <- times$id

    #Seperate date time column
    ghg1 <- ghg %>%
      mutate(date.time = Time) %>%
      separate(
        col = Time,
        into = c("date", "time"),
        sep = " ",
        remove = F
      ) %>%
      mutate(date1 = date)

    ghg1$date1 <- gsub("/", "-", ghg1$date1)
    ghg1$date <- ymd(ghg1$date)
    ghg1 <- filter(ghg1, date == i)
    date <- ghg1$date[1]

    #remove X from names
    names(ghg1) <- gsub("[][]", replacement = "", x = names(ghg1))

    #hms time
    ghg1$time <- hms(ghg1$time)

    #select time periods needed by group
    message("\nMapping Times...")



    ghg_out <- map2(starts, t1, ~ filter(ghg1, time >= .x & time <= .y))

    ghg_out <- set_names(ghg_out, mynames)

    ghg_out <<- ghg_out

    ghg <- map_df(ghg_out, ~ as.data.frame(.x), .id = "id")


    #Container information, volume, area

    chmb <- tibble(chamber = c("veg_sml", "veg_lrg", "ditch", "bottle", "box"),
                   chamber_volume = c(0.002097341,0.022166568,0.032951432,0.016,0.036),
                   chamber_area = c(0.010029, 0.064692, 0.098980, 0.0707,0.09))

    names(time1) <- gsub("dipwell", "sample_ID", names(time1))
    names(time1) <- gsub("Chamber Type", "chamber", names(time1))


    ghg_renew <- data.frame()


    for (f in unique(ghg$id)) {


      ghg_f <- ghg %>%
        filter(id == f)

      if (chamber_p == T) {

        tm1 <- times %>%
          select(id, chamber) %>%
          left_join(chmb, by = "chamber")

        ghg_f <- ghg_f %>%
          left_join(tm1, by = "id")

        }else{
        ghg_f$chamber_volume <- chamb_volume
        ghg_f$chamber_area <- chamb_area
      }


      ghg_renew <- rbind(ghg_renew, ghg_f)
    }


    #Calculate and graph concentrations and flux

    #Standard pressure(kpa)
    sp <- 101.325
    #Standard temp(k)
    st <- 273.15
    # vol 1 mol/m3 ideal gas at stp
    gas <- 0.022414

    #Convert and calculate c to k, torr to kpa, ppm to umol, and time in seconds

    lgr <- ghg_renew

    lgr <- lgr %>%
      group_by(id) %>%
      mutate(
        GasP_kpa = GasP_torr / 7.5006157818041,
        Amb.t.k = AmbT_C + 273.15,
        time = hms(time),
        date.time = ymd_hms(date.time),
        vol_co2 = (CO2d_ppm / 1000000) * chamber_volume,
        vol_co2_stp = vol_co2 * (sp / Amb.t.k) * (st / sp),
        mol_co2 = (vol_co2_stp / gas),
        mass_co2_g = mol_co2 * 44.0095,
        mass_co2_umol = mol_co2 * 1000000,
        vol_ch4 = (CH4d_ppm / 1000000) * chamber_volume,
        vol_ch4_stp = vol_ch4 * (sp / Amb.t.k) * (st / sp),
        mol_ch4 = (vol_ch4_stp / gas),
        mass_ch4_g = mol_ch4 * 44.0095,
        mass_ch4_umol = mol_ch4 * 1000000,
        t = as.numeric(time - lag(time)),
        t = replace_na(t, 0),
        t = cumsum(t)
      )

    lgr$site <- site_name



    #Edit and check fluxes

    flux_edit <<- data.frame()
    for (i in  unique(lgr$id)) {
      lgr1 <- lgr %>%
        filter(id == i)
      repeat {
        grph(ch4 = T)
        if (!conf == "3") {
          break
        }
      }
    }

    flux_edit_2 <<-  data.frame()

    for (i in unique(flux_edit_ch4$id)) {
      f <- flux_edit_ch4 %>%
        filter(id == i)
      if (anyNA(f$CH4_ppm))
      {
        next
      }
      tryCatch({
        check(i, ch4 = T)
      },
      error = function(e) {
        print(paste0(i, ": Flux stats are too poor, Flux will be added as NA"))
      })
    }

    chk_ids <- unique(flux_edit_ch4$id)
    chk_ids_2 <- unique(flux_edit_ch4_1$id)
    dif <- setdiff(chk_ids, chk_ids_2)
    if (!is_empty(dif)) {
      for (h in dif) {
        add <- flux_edit_ch4 %>%
          filter(id == h)
        add[5:58] <- NA
        flux_edit_ch4_1 <- rbind(flux_edit_ch4_1, add)
        assign("flux_edit_ch4_1", flux_edit_ch4_1, .GlobalEnv)
      }
    }


    flux_edit <<- data.frame()
    for (i in  unique(lgr$id)) {
      lgr1 <- lgr %>%
        filter(id == i)
      repeat {
        grph(ch4 = F)
        if (!conf == "3") {
          break
        }
      }
    }



    flux_edit_2 <<-  data.frame()
    for (i in unique(flux_edit_co2$id)) {
      f <- flux_edit_co2 %>%
        filter(id == i)
      if (anyNA(f$CO2_ppm))
      {
        next
      }
      tryCatch({
        check(i, ch4 = F)
      },
      error = function(e) {
        print(paste0(i, ": Flux stats are too poor, Flux will be added as NA"))
      })
    }

    chk_ids <- unique(flux_edit_co2$id)
    chk_ids_2 <- unique(flux_edit_co2_1$id)
    dif <- setdiff(chk_ids, chk_ids_2)
    dif
    if (!is_empty(dif)) {
      for (h in dif) {
        add <- flux_edit_co2 %>%
          filter(id == h)
        add[5:58] <- NA
        flux_edit_co2_1 <- rbind(flux_edit_co2_1, add)
        assign("flux_edit_co2_1", flux_edit_co2_1, .GlobalEnv)
      }
    }

    #calculate regression coefficient and flux
    names_co2 <-
      c("id",
        "co2_intercept",
        "co2_slope",
        "co2_r2",
        "co2_flux")
    names_ch4 <-
      c("id",
        "ch4_intercept",
        "ch4_slope",
        "ch4_r2",
        "ch4_flux")

    # Function to Extract coefs and rsquare from lmlist into df
    sumfun <- function(x)
    {
      aux <- function(x)
        c(coef(x), summary(x)$r.squared)
      t(sapply(x, aux))
    }


    #Generate Orginal
    co2_coefs <- lmList(mass_co2_g ~ t |
                          id,
                        data = flux_edit_co2_1,
                        na.action = na.omit)


    co2_coefs <- as.data.frame(coef(co2_coefs)) %>%
      rownames_to_column("id") %>%
      left_join(tm1, by = "id") %>%
      mutate(co2_flux =  t / chamber_area)

    names(co2_coefs) <- names_co2

    ch4_coefs <- lmList(mass_ch4_g ~ t |
                          id,
                        data = flux_edit_ch4_1,
                        na.action = na.omit)

    ch4_coefs <- as.data.frame(coef(ch4_coefs)) %>%
      rownames_to_column("id") %>%
      left_join(tm1, by = "id") %>%
      mutate(co2_flux =  t / chamber_area)


    flux_ch4 <- ch4_coefs
    flux_ch4$date <- date
    flux_ch4$site <- site_name
    flux_ch4 <- flux_ch4 %>%
      relocate(date, site, id)

    flux_co2 <- co2_coefs
    flux_co2$date <- date
    flux_co2$site <- site_name
    flux_co2 <- flux_co2 %>%
      relocate(date, site, id)

    write_csv(flux_ch4, csv_ch4_flux)
    write_csv(flux_co2, csv_co2_flux)

    #Graph concentrations
    message("\nMapping CO2 Regressions...")
    co2_conc <-  flux_edit_co2_1 %>%  ggscatter(
      x = "t",
      y = "mass_co2_g",
      add = "reg.line",
      # Add regressin line
      add.params = list(color = "blue", fill = "lightgray"),
      # Customize reg. line
      conf.int = TRUE,
      # Add confidence interval
      cor.coef = TRUE,
      # Add correlation coefficient
      cor.coeff.args = list(method = "pearson", label.sep = "\n"),
      scales = "free",
      facet.by = "id"
    ) +
      ggtitle("CO2 Concentration Edited")
    suppressMessages(ggsave(
      co2_conc,
      width = 45,
      height = 45,
      units = "cm",
      filename = paste0(grph_nm,"_grph_co2_edited.png"
      )))

    co2_conc_org <-  lgr %>%  ggscatter(
      x = "t",
      y = "mass_co2_g",
      add = "reg.line",
      # Add regressin line
      add.params = list(color = "blue", fill = "lightgray"),
      # Customize reg. line
      conf.int = TRUE,
      # Add confidence interval
      cor.coef = TRUE,
      # Add correlation coefficient
      cor.coeff.args = list(method = "pearson", label.sep = "\n"),
      scales = "free",
      facet.by = "id"
    ) +
      ggtitle("CO2 Concentration Original")
    suppressMessages(ggsave(
      co2_conc_org,
      width = 45,
      height = 45,
      units = "cm",
      filename = paste0(grph_nm,"_grph_co2_original.png")
    ))
    message("
Mapping CH4 Regressions...")
    ch4_conc<- flux_edit_ch4_1 %>% ggscatter(
      x = "t",
      y = "mass_ch4_g",
      add = "reg.line",
      # Add regressin line
      add.params = list(color = "blue", fill = "lightgray"),
      # Customize reg. line
      conf.int = TRUE,
      # Add confidence interval
      cor.coef = TRUE,
      # Add correlation coefficient. see ?stat_cor
      cor.coeff.args = list(method = "pearson", label.sep = "\n"),
      scales = "free",
      facet.by = "id"
    ) +
      ggtitle("CH4 Concentration")
    suppressMessages(ggsave(
      ch4_conc,
      width = 45,
      height = 45,
      units = "cm",
      filename = paste0(grph_nm, "_grph_ch4_edited.png")
    ))

    ch4_conc_org <- lgr %>% ggscatter(
      x = "t",
      y = "mass_ch4_g",
      add = "reg.line",
      # Add regressin line
      add.params = list(color = "blue", fill = "lightgray"),
      # Customize reg. line
      conf.int = TRUE,
      # Add confidence interval
      cor.coef = TRUE,
      # Add correlation coefficient. see ?stat_cor
      cor.coeff.args = list(method = "pearson", label.sep = "\n"),
      scales = "free",
      facet.by = "id"
    ) +
      ggtitle("CH4 Concentration")
    suppressMessages(ggsave(
      ch4_conc_org,
      width = 45,
      height = 45,
      units = "cm",
      filename = paste0(grph_nm, "_grph_ch4_original.png")
    ))


    cat(red(bold(paste0("\nFluxes for ", ff, " complete.\n"))))
    Sys.sleep(2)

    suppressWarnings(if(exists("flux_edit")){
    rm(list = c("flux_edit", "flux_edit_2", "flux_edit_ch4",
                "flux_edit_ch4_1", "flux_edit_co2", "flux_edit_co2_1", "a",
                "ch4", "conf", "conf1", "lsf"))})

  }
  cat(red(bold("\nAnalysis Complete\n")))

  if(exists("flux_edit")){
    rm(list = c("flux_edit", "flux_edit_2", "flux_edit_ch4",
                "flux_edit_ch4_1", "flux_edit_co2", "flux_edit_co2_1", "a",
                "ch4", "conf", "conf1", "lsf"))}

}
