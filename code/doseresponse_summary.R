require(tidyverse)
# devtools::install_github("AndersenLab/easyXpress")
require(easyXpress)
require(drc)
require(ddpcr)
require(RColorBrewer)
require(cowplot)
require(ggbeeswarm)
require(ggrepel)
require(ggnewscale)
require(magrittr)
require(lmtest)
require(sandwich)
require(sommer)
library(ggh4x)
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))
load("data/FileS1.RData")
tx.classes <- data.table::fread(file = "data/tx.classes.csv") %>%
  dplyr::rename(Toxicant = drug)
strain.data <- data.table::fread(file = "data/CelegansStrainData.tsv")

#############
# FUNCTIONS #
#############
strain_colors       <- c("blue",   "orange","#5A0C13","#C51B29",  "#a37000","#627264","#67697C","purple")
names(strain_colors) <- c("CB4856","N2", "ECA36",  "ECA396"  ,"CB4855", "RC301",   "MY16", "XZ1516")
concatenate.assays <- function(x){
  assay <- get(x)
  assay.df <- assay$summarized_processed %>%
    dplyr::mutate(Metadata_Plate = gsub(x = Metadata_Plate, 
                                        pattern = "p", 
                                        replacement = ""),
                  well_censor_reason = as.character(well_censor_reason),
                  notes = as.character(notes))
  
  nested.doses <- assay$summarized_processed %>%
    dplyr::select(drug, concentration_um) %>%
    dplyr::distinct() %>%
    dplyr::group_by(drug) %>%
    tidyr::nest()
  dose.rank.list <- list()
  for(i in 1:length(nested.doses$drug)){
    dose.rank.list[[i]] <- nested.doses$data[[i]] %>%
      dplyr::mutate(dose.rank = 1:nrow(nested.doses$data[[i]]),
                    drug = nested.doses$drug[[i]])
  }
  dose.ranks.df <- dose.rank.list %>%
    Reduce(rbind,.)
  assay.df %>%
    dplyr::full_join(.,dose.ranks.df)
  
  
}
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = na.rm, 
                         ...)
  H <- 1.5 * stats::IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
collapseDoses <- function(nested.raw.data){
  if(length(unique(nested.raw.data$concentration_um)) == 1){
    return(nested.raw.data)
  } else {
    nested.raw.data$concentration_um <- unique(nested.raw.data$concentration_um)[1]
    return(nested.raw.data)
  }
}
cleanData <- function(raw.data){
  
  nested.raw.data <- raw.data %>%
    dplyr::group_by(dose.rank) %>%
    tidyr::nest()
  
  raw.data <- purrr::map(nested.raw.data$data, collapseDoses) %>%
     Reduce(rbind,.)
  
  if(length(levels(as.factor(raw.data$Metadata_Experiment))) > 1){
    
    # Remove Low n and Censored Wells
    n.drug.df <- raw.data %>%
      dplyr::filter(n < 30,
                    n > 5,
                    is.na(well_censor))
    
    # Outlier Removal
    outlier.removed.drug.df <- n.drug.df %>%
      dplyr::group_by(concentration_um, strain, drug) %>%
      dplyr::mutate(median_wormlength_um = remove_outliers(median_wormlength_um)) %>%
      dplyr::filter(!is.na(median_wormlength_um))
    
    # Complete dose representation
    complete.doses <- outlier.removed.drug.df %>%
      dplyr::ungroup() %>%
      dplyr::select(Metadata_Experiment, concentration_um) %>%
      dplyr::distinct() %>%
      dplyr::group_by(concentration_um) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::filter(n > length(unique(outlier.removed.drug.df$Metadata_Experiment))*0.8)
    
    outlier.removed.drug.df  %<>%
      dplyr::filter(concentration_um %in% complete.doses$concentration_um)
    
    # Regression
    median.reg <- lm(data = outlier.removed.drug.df, formula = median_wormlength_um ~ Metadata_Experiment + bleach)
    mean.reg <- lm(data = outlier.removed.drug.df, formula = mean_wormlength_um ~ Metadata_Experiment + bleach)
    outlier.removed.drug.df <- cbind(outlier.removed.drug.df, median.reg$residuals, mean.reg$residuals)
    colnames(outlier.removed.drug.df) <- c(colnames(outlier.removed.drug.df)[-c(length(colnames(outlier.removed.drug.df))-1,length(colnames(outlier.removed.drug.df)))],
                                           "median_wormlength_um_reg","mean_wormlength_um_reg")
    
    
    # Control Delta
    if(!0 %in% outlier.removed.drug.df$concentration_um){
      full_df_proc <-  outlier.removed.drug.df
    } else {
      control_values <- outlier.removed.drug.df %>%
        dplyr::filter(concentration_um == 0) %>% # filter to control wells
        dplyr::group_by(strain, bleach, Metadata_Experiment, drug) %>%
        dplyr::mutate(control_pheno = mean(median_wormlength_um_reg)) %>% # get mean control value for each trait and strain
        dplyr::distinct(control_pheno, bleach, strain, Metadata_Experiment) # make it neat
      
      full_df_proc <-  outlier.removed.drug.df %>%
        dplyr::ungroup() %>%
        dplyr::left_join(., control_values) %>% # join control values for each trait
        dplyr::mutate(median_wormlength_um_reg = median_wormlength_um_reg - control_pheno)
    }
    
    n.raw <- raw.data %>%
      dplyr::group_by(Metadata_Experiment, bleach, drug, strain) %>%
      dplyr::summarise(raw = n())
    
    n.cleaned <- n.drug.df %>%
      dplyr::group_by(Metadata_Experiment, bleach, drug, strain) %>%
      dplyr::summarise(cleaned = n())
    
    n.outlier.removed <- outlier.removed.drug.df %>%
      dplyr::group_by(Metadata_Experiment, bleach, drug, strain) %>%
      dplyr::summarise(outlier.removed = n())
    
    filtering.summary <- n.raw %>%
      dplyr::full_join(., n.cleaned) %>%
      dplyr::full_join(., n.outlier.removed)
    filtering.summary[is.na(filtering.summary)] <- 0
    
    drug.filtering.summary <- filtering.summary %>%
      dplyr::mutate(pct.retained = (outlier.removed/raw)*100,
                    pct.outlier = ((cleaned-outlier.removed)/raw)*100,
                    pct.cleaned = ((raw-cleaned)/raw)*100,
                    check = pct.retained + pct.outlier + pct.cleaned) %>%
      tidyr::pivot_longer(cols = -c(Metadata_Experiment, bleach, drug, strain), names_to = "stat")
    
  } else {
    
    # Remove Low n and Censored Wells
    n.drug.df <- raw.data %>%
      dplyr::filter(n < 30,
                    n > 5,
                    is.na(well_censor))
    
    
    # Outlier Removal
    outlier.removed.drug.df <- n.drug.df %>%
      dplyr::group_by(concentration_um, strain, drug) %>%
      # dplyr::group_by(concentration_um) %>%
      dplyr::mutate(median_wormlength_um = remove_outliers(median_wormlength_um)) %>%
      dplyr::filter(!is.na(median_wormlength_um))
    
    complete.doses <- outlier.removed.drug.df %>%
      dplyr::ungroup() %>%
      dplyr::select(Metadata_Experiment, concentration_um) %>%
      dplyr::distinct() %>%
      dplyr::group_by(concentration_um) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::filter(n > length(unique(outlier.removed.drug.df$Metadata_Experiment))*0.8)
    
    outlier.removed.drug.df  %<>%
      dplyr::filter(concentration_um %in% complete.doses$concentration_um)
    
    # Regression
    median.reg <- lm(data = outlier.removed.drug.df, formula = median_wormlength_um ~ bleach)
    mean.reg <- lm(data = outlier.removed.drug.df, formula = mean_wormlength_um ~ bleach)
    outlier.removed.drug.df <- cbind(outlier.removed.drug.df, median.reg$residuals, mean.reg$residuals)
    colnames(outlier.removed.drug.df) <- c(colnames(outlier.removed.drug.df)[-c(length(colnames(outlier.removed.drug.df))-1,length(colnames(outlier.removed.drug.df)))],
                                           "median_wormlength_um_reg","mean_wormlength_um_reg")
    
    # Control Delta
    if(!0 %in% outlier.removed.drug.df$concentration_um){
      full_df_proc <-  outlier.removed.drug.df
    } else {
      control_values <- outlier.removed.drug.df %>%
        dplyr::filter(concentration_um == 0) %>% # filter to control wells
        dplyr::group_by(strain, bleach, Metadata_Experiment, drug) %>%
        dplyr::mutate(control_pheno = mean(median_wormlength_um_reg)) %>% # get mean control value for each trait and strain
        dplyr::distinct(control_pheno, bleach, strain, Metadata_Experiment) # make it neat
      
      full_df_proc <-  outlier.removed.drug.df %>%
        dplyr::ungroup() %>%
        dplyr::left_join(., control_values) %>% # join control values for each trait
        dplyr::mutate(median_wormlength_um_reg = median_wormlength_um_reg - control_pheno)
    }
    
    n.raw <- raw.data %>%
      dplyr::group_by(Metadata_Experiment, bleach, drug, strain) %>%
      dplyr::summarise(raw = n())
    
    n.cleaned <- n.drug.df %>%
      dplyr::group_by(Metadata_Experiment, bleach, drug, strain) %>%
      dplyr::summarise(cleaned = n())
    
    n.outlier.removed <- outlier.removed.drug.df %>%
      dplyr::group_by(Metadata_Experiment, bleach, drug, strain) %>%
      dplyr::summarise(outlier.removed = n())
    
    filtering.summary <- n.raw %>%
      dplyr::full_join(., n.cleaned) %>%
      dplyr::full_join(., n.outlier.removed)
    filtering.summary[is.na(filtering.summary)] <- 0
    
    drug.filtering.summary <- filtering.summary %>%
      dplyr::mutate(pct.retained = (outlier.removed/raw)*100,
                    pct.outlier = ((cleaned-outlier.removed)/raw)*100,
                    pct.cleaned = ((raw-cleaned)/raw)*100,
                    check = pct.retained + pct.outlier + pct.cleaned) %>%
      tidyr::pivot_longer(cols = -c(Metadata_Experiment, bleach, drug, strain), names_to = "stat")
  }
  
  full_df_proc %<>%
    dplyr::ungroup() %>%
    dplyr::group_by(concentration_um, strain, drug) %>%
    dplyr::mutate(median_wormlength_um_reg = remove_outliers(median_wormlength_um_reg))
  
  return(list(full_df_proc,drug.filtering.summary))
}
# dose-response model inference
safe_LL4_strain <- purrr::safely(.f = ~ drc::drm(data = ., 
                                                 formula = value ~ concentration_um, strain,
                                                 pmodels=list(~strain-1,  ~1, ~1, ~strain-1),
                                                 fct = LL.4(fixed=c(NA, -600, NA, NA))),
                                 otherwise = "Unable to optimize model")

safe_LL4_no_strain <- purrr::safely(.f = ~ drc::drm(data = ., 
                                                    formula = value ~ concentration_um,
                                                    pmodels=list(~1,  ~1, ~1, ~1),
                                                    fct = LL.4(fixed=c(NA, -600, NA, NA))),
                                    otherwise = "Unable to optimize model")
safe_ED <- function(fit){
  if(is.character(fit$result)){print("Unable to optimize model")} 
  else {
    drc::ED(object = fit$result, c(10, 50, 90), interval = "delta", display = FALSE)
  }
}
safe_slope_strain <- function(fit){
  if(is.character(fit$result)){print("Unable to extract slope")}
  else {
    slopes<- data.frame(fit$result$coefficients[1:8]) %>%
      dplyr::rename(slope = `fit.result.coefficients.1.8.`) %>%
      dplyr::mutate(strain = rownames(.),
                    strain = gsub(strain, pattern = "b:strain", replacement = ""))
    rownames(slopes) <- NULL
    return(slopes)
  }
}
safe_lower_asym <- function(fit){
  if(is.character(fit$result)){print("Unable to extract slope")}
  else {
    
    asyms <- data.frame(fit$result$coefficients[10:17]) %>%
      dplyr::rename(lower.asymp = `fit.result.coefficients.10.17.`) %>%
      dplyr::mutate(strain = rownames(.),
                    strain = gsub(strain, pattern = "d:strain", replacement = ""))
    rownames(asyms) <- NULL
    return(asyms)
  }
}
safe_predict <- function(data, fit){
  if(is.character(fit$result)){print("Unable to optimize model")}
  else{
    expand.grid(conc = exp(seq(log(max(data$concentration_um)),
                               log(min(data$concentration_um[data$concentration_um!=min(data$concentration_um)])),
                               length = 100))) %>%
      dplyr::bind_cols(., tibble::as.tibble(predict(fit$result, newdata=., interval="confidence")))
  }
  
}
fit_DRC_models <- function(data, traits = "median_wormlength_um_reg", ...){
  # extract concentration variable name
  conc_var_name <- str_subset(names(data), pattern = "concentration_um")
  
  
  ##### old code #####
  # complete.doses <- data %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select(strain, drug, concentration_um) %>%
  #   dplyr::distinct() %>%
  #   dplyr::group_by(concentration_um) %>%
  #   dplyr::summarise(n()) %>%
  #   dplyr::filter(`n()` == 8)
  # 
  
  # # gather trait data, nest by grouping variables and trait, fit drc, find ED10, ED50, ED90, and make predictions.
  
  # WORKING
  # LL4_fit_data <- data %>%
  #   dplyr::select(Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain, 
  #                 contains("wormlength"), concentration_um := !!glue::glue("{conc_var_name}")) %>%
  #   tidyr::gather(trait, value, -drug, -strain, -concentration_um, -Metadata_Plate, -Metadata_Well, -Metadata_Experiment, -bleach) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select(drug, strain, trait, value, concentration_um) %>% # reorder for drc::drm
  #   dplyr::filter(trait %in% traits) %>% # filter to traits provided in arguments
  #   # dplyr::group_by(drug, strain, trait, Metadata_Experiment) %>% 
  #   dplyr::group_by(drug, strain, trait) %>% # should eventually be "..." passed from function + trait
  #   tidyr::nest() %>% 
  #   mutate(fit = purrr::map(data, safe_LL4),
  #          slope = purrr::map(fit, safe_slope),
  #          lower_asym = purrr::map(fit, safe_lower_asym),
  #          ED = purrr::map(fit, safe_ED),
  #          predictions = purrr::map2(data, fit, safe_predict),
  #          model = "LL.4")
  ##### 
  
  
  LL4_fit_data <- data %>%
    dplyr::select(Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain, 
                  contains("wormlength"), concentration_um := !!glue::glue("{conc_var_name}")) %>%
    tidyr::gather(trait, value, -drug, -strain, -concentration_um, -Metadata_Plate, -Metadata_Well, -Metadata_Experiment, -bleach) %>%
    dplyr::ungroup() %>%
    dplyr::select(drug, strain, trait, value, concentration_um) %>% # reorder for drc::drm
    dplyr::filter(trait %in% traits) %>% # filter to traits provided in arguments
    # dplyr::group_by(drug, strain, trait, Metadata_Experiment) %>% 
    dplyr::group_by(drug, trait) %>% # should eventually be "..." passed from function + trait
    tidyr::nest() %>% 
    mutate(fit = purrr::map(data, safe_LL4_strain),
           slope = purrr::map(fit, safe_slope_strain),
           lower_asym = purrr::map(fit, safe_lower_asym),
           ED = purrr::map(fit, safe_ED),
           model = "LL.4")
  return(LL4_fit_data)
}
logticks <- function(datavar,type) {
  minimum <- 1/10^abs(floor(log10(min(datavar, na.rm=TRUE))))
  maximum <- 1*10^abs(floor(log10(max(datavar, na.rm=TRUE)))+1)
  multiple <- floor(log10(maximum/minimum))
  
  yourtickvector <- c()
  
  if (type=="breaks") {
    
    yourtickvector <- c(minimum)
    
    for (x in seq(0,multiple)) {
      
      andadd <- seq(minimum*10^x,minimum*10^(x+1),minimum*10^x)[-1]
      
      yourtickvector <- c(yourtickvector,andadd)
      
    }
  } else if (type=="labels") {
    
    for (x in seq(0,multiple)) {
      
      andadd <- c(minimum*10^x,rep("",8))
      
      yourtickvector <- c(yourtickvector,andadd)
      
    }
    
    yourtickvector <- c(yourtickvector,minimum*10^multiple)
    
  }
  
  return(yourtickvector)
  
}
dose.check.plot <- function(dat){
  strain_colors       <- c("blue",   "orange","#5A0C13","#C51B29","#a37000","#627264","#67697C","purple")
  names(strain_colors) <- c("CB4856","N2", "ECA36",  "ECA396"  ,"CB4855", "RC301",   "MY16", "XZ1516")
  dose.check <- ggplot(dat) + 
    theme_bw(base_size = 8) + 
    # geom_jitter(mapping = aes(x = strain, y = median_wormlength_um_reg, colour = Metadata_Experiment),size = 1) + 
    geom_quasirandom(mapping = aes(x = strain, y = median_wormlength_um_reg, colour = Metadata_Experiment), size = 1) + 
    geom_boxplot(mapping = aes(x = strain, y = median_wormlength_um_reg), 
                 outlier.shape = NA, alpha = 0.5) + 
    facet_wrap(.~round(concentration_um,2), scales = "free") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "top") + 
    labs(fill = "Strain",
         y = "Relative Worm Length (um)",
         x = "Strain") + 
    theme(plot.title = element_text(hjust = 0.5), legend.box="vertical") + 
    # ylim(c(-810,200)) +
    scale_colour_brewer(palette = "Dark2") + 
    ggtitle(unique(dat$drug))
  
  ggsave(dose.check, filename = paste0("output/", paste(paste(strsplit(unique(dat$drug),split = " ")[[1]], collapse = "_"), "dose.check.plot", sep = "_"),".png"), 
         width = 10, height = 6)
}
data.retention <- function(data){
  
  A <- data %>%
    dplyr::filter(stat %in% c("pct.retained","pct.outlier","pct.cleaned")) %>%
    dplyr::mutate(bleach = paste("Bleach", bleach, sep = " ")) %>%
    ggplot(., mapping = aes(x = strain, y = value, fill = stat)) + 
    theme_bw() +
    geom_col() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "top",
          panel.grid = element_blank()) +
    scale_fill_brewer(palette = "Set1", name = "Cleaning Category") + 
    facet_wrap(.~Metadata_Experiment+bleach) + 
    ggtitle("Data Retention Cleaning: Before") + 
    labs(x = "Strain", y = "Experimental Retention (% of Expected Wells)")
  
  data.retention.df <- data %>%
    tidyr::pivot_wider(names_from = stat, values_from = value) %>%
    dplyr::filter(pct.cleaned < 35) %>%
    dplyr::group_by(Metadata_Experiment, bleach) %>%
    dplyr::summarise(n.strains = n()) %>%
    dplyr::filter(n.strains > 4)
  
  
  B <-  data.retention.df %>%
    dplyr::left_join(.,data) %>%
    dplyr::filter(stat %in% c("pct.retained","pct.outlier","pct.cleaned")) %>%
    # tidyr::pivot_longer(cols = -c(Metadata_Experiment, bleach), names_to = "stat") %>%
    ggplot(., mapping = aes(x = strain, y = value, fill = stat)) + 
    theme_bw() +
    geom_col() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "top",
          panel.grid = element_blank()) +
    scale_fill_brewer(palette = "Set1", name = "Cleaning Category") + 
    facet_wrap(.~Metadata_Experiment+bleach) + 
    ggtitle("Data Retention Cleaning: After") + 
    labs(x = "Strain", y = "")
  
  AB <- cowplot::plot_grid(A + theme(legend.position = "none"), 
                           B + theme(legend.position = "none"))
  
  title <- ggdraw() + draw_label(unique(data$drug), fontface = "bold")
  data.retention.plot <- cowplot::plot_grid(title, AB, ncol = 1, 
                     rel_heights = c(0.1, 1))
  ggsave(data.retention.plot + theme(plot.background = element_rect(fill = "white")), filename = paste0("output/", paste(paste(strsplit(unique(data$drug),split = " ")[[1]], collapse = "_"), "data.retention.plot", sep = "_"),".png"), 
         width = 10, height = 6)
  
  return(data.retention.df[,1:2])
  
}
se.DRC <- function(x, dodge.width = 0.15, point.size = 0.4){
  
  strain_colors       <- c("blue",   "orange","#5A0C13","#C51B29",  "#a37000","#627264","#67697C","purple")
  names(strain_colors) <- c("CB4856","N2", "ECA36",  "ECA396"  ,"CB4855", "RC301",   "MY16", "XZ1516")
  
  norm <- x %>%
    dplyr::group_by(strain) %>%
    dplyr::filter(concentration_um == 0) %>%
    dplyr::summarize(mean_ctrl = mean(mean_wormlength_um_reg)) %>%
    dplyr::ungroup() %>%
    dplyr::select(strain, mean_ctrl) %>% ##added trait
    dplyr::full_join(x, by = c("strain")) %>% ## "trait" instead of "condition"
    dplyr::distinct() %>%
    dplyr::mutate(norm_pheno = mean_wormlength_um_reg - mean_ctrl) 
  
  dose_avg <- norm %>%
    dplyr::group_by(strain, concentration_um) %>%
    dplyr::summarise(avg = mean(norm_pheno), st_dev = sd(norm_pheno))

  
  
  strains <- dose_avg %>%
    dplyr::mutate(log.concentration = round(log(concentration_um), 2))
  
  # conc.range <- range(strains[which(strains$log.concentration != -Inf),]$concentration_um)[2] - range(strains[which(strains$log.concentration != -Inf),]$concentration_um)[1] 
  # dodge.width = conc.range/650
  
  SE.DRC <- strains[which(strains$log.concentration != -Inf),] %>%
    ggplot(.,mapping = aes(x=concentration_um, y=avg)) +
    theme_bw(base_size = 18) + 
    geom_line(mapping = aes(group = strain, colour = strain), 
              position = position_dodge2(width = dodge.width), size = 0.15) +
    geom_pointrange(mapping = aes(ymin=avg-st_dev, 
                                  ymax = avg+st_dev,
                                  colour = strain), 
                    position = position_dodge2(width = dodge.width),
                    size = point.size) + 
    
    scale_colour_manual(values = strain_colors, name = "Strain") +
    scale_y_continuous("Normalized Worm Length (µm)", limits = c(-700,100), breaks = seq(-700, 0, 100)) + 
    scale_x_log10() +
    theme(panel.grid = element_blank(),
          legend.position = c(0.15,0.25),
          legend.background = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text = element_text(colour = "black")) + 
    labs(x="Concentration (μM)") + 
    ggtitle(unique(x$drug))
  ggsave(SE.DRC, filename = paste0("output/", paste(paste(strsplit(unique(x$drug),split = " ")[[1]], collapse = "_"), "SE.DRC.plot", sep = "_"),".png"), 
         width = 6, height = 6)
  print(SE.DRC)
  
}
# data = tx.nested.dose.responses$data[[11]]
# toxicant = tx.nested.dose.responses$drug[[11]]
dose.response.summary <- function(data, toxicant){
  
  # Data cleaning
  raw.data <- data %>%
    dplyr::mutate(drug = toxicant)
  cleaned.data <- suppressMessages(cleanData(raw.data))
  
  
  retained.bleaches <- data.retention(cleaned.data[[2]])
  dose.check.plot(dat = cleaned.data[[1]])
  
  cleaned.filtered.data <- retained.bleaches %>%
    dplyr::left_join(.,cleaned.data[[1]]) %>%
    dplyr::group_by(strain, concentration_um) %>%
    dplyr::count() %>%
    dplyr::arrange(n) %>%
    dplyr::filter(n > 10) %>%
    dplyr::select(-n) %>%
    dplyr::left_join(.,cleaned.data[[1]])

  se.DRC(cleaned.filtered.data)
  
  drug.drm <- suppressMessages(fit_DRC_models(data = cleaned.filtered.data,
                                              traits = "median_wormlength_um_reg"))
  
  
  # Population-wide LOAELs
  population.wide.LOAEL <- aov(data = cleaned.filtered.data, formula = median_wormlength_um_reg ~ as.factor(concentration_um))
  LOAEL.tukey <- TukeyHSD(population.wide.LOAEL)
  
  LOAEL.tests <- rownames(LOAEL.tukey$`as.factor(concentration_um)`) %>% 
    data.frame() %>%
    tidyr::separate(data = ., col = `.`, c("dose1","dose2"), sep = "-") %>%
    dplyr::filter(dose2 == 0) %>%
    rownames() %>% as.numeric()
  
  adj.Ps <- LOAEL.tukey$`as.factor(concentration_um)`[LOAEL.tests,] %>%
    data.frame()
  
  doses.adj.Ps <- rownames(LOAEL.tukey$`as.factor(concentration_um)`) %>% 
    data.frame() %>%
    tidyr::separate(data = ., col = `.`, c("dose1","dose2"), sep = "-") %>%
    dplyr::filter(dose2 == 0) %>%
    cbind(.,adj.Ps)
  rownames(doses.adj.Ps) <- NULL
  
  
  LOAEL <-  doses.adj.Ps %>%
    dplyr::filter(p.adj < 0.05 & diff < 0) %>%
    dplyr::mutate(tx = toxicant) %>%
    dplyr::rename(LOAEL = dose1) %>%
    dplyr::select(-dose2)
  rownames(LOAEL) <- NULL
  overall.tx.LOAEL <- LOAEL[which(LOAEL$LOAEL == min(LOAEL$LOAEL)),]
  
  
  # Strain-specific LOAELs
  LOAEL.strain.nested <- cleaned.filtered.data %>%
    dplyr::group_by(strain) %>%
    tidyr::nest()
  strain.LOAELs <- list()
  for(i in 1:length(LOAEL.strain.nested$strain)){
    strain.LOAEL <- aov(data = LOAEL.strain.nested$data[[i]], formula = median_wormlength_um_reg ~ as.factor(concentration_um))
    strain.LOAEL.tukey <- TukeyHSD(strain.LOAEL)
    
    strain.LOAEL.tests <- rownames(strain.LOAEL.tukey$`as.factor(concentration_um)`) %>% 
      data.frame() %>%
      tidyr::separate(data = ., col = `.`, c("dose1","dose2"), sep = "-") %>%
      dplyr::filter(dose2 == 0) %>%
      rownames() %>% as.numeric()
    
    strain.adj.Ps <- strain.LOAEL.tukey$`as.factor(concentration_um)`[strain.LOAEL.tests,] %>%
      data.frame() %>%
      dplyr::select(p.adj)
    
    strain.doses.adj.Ps <- rownames(strain.LOAEL.tukey$`as.factor(concentration_um)`) %>% 
      data.frame() %>%
      tidyr::separate(data = ., col = `.`, c("dose1","dose2"), sep = "-") %>%
      dplyr::filter(dose2 == 0) %>%
      cbind(.,strain.adj.Ps)
    rownames(strain.doses.adj.Ps) <- NULL
    
    
    strain.LOAEL.df <-  strain.doses.adj.Ps %>%
      dplyr::filter(p.adj < 0.05) %>%
      dplyr::mutate(tx = toxicant) %>%
      dplyr::rename(LOAEL = dose1) %>%
      dplyr::select(-dose2) %>%
      dplyr::mutate(strain = LOAEL.strain.nested$strain[[i]])
    rownames(strain.LOAEL.df) <- NULL
    strain.LOAELs[[i]] <- strain.LOAEL.df
  }
  
  strain.tx.LOAEL <- strain.LOAELs %>%
    Reduce(rbind,.)
  
  
  # dose.response.parameters <- list()
  # for(j in 1:length(drug.drm$data)){
  #   rownames(drug.drm$ED[[j]]) <- c("EC10","EC50","EC90")
  #   ECs <- data.frame(drug.drm$ED[[j]])
  #   ECs$metric <- rownames(ECs)
  #   rownames(ECs) <- NULL
  #   
  #   slope <- data.frame(drug.drm$slope[[j]]) %>%
  #     dplyr::rename(Estimate = `drug.drm.slope..j..`) %>%
  #     dplyr::mutate(metric = "Slope")
  #   
  #   lower_asymp <- data.frame(drug.drm$lower_asym[[j]]) %>%
  #     dplyr::rename(Estimate = `drug.drm.lower_asym..j..`) %>%
  #     dplyr::mutate(metric = "Lower Asymptote")
  #   
  #   high_dose <- max(unique(drug.drm$data[[j]]$concentration_um))
  #   
  #   EC.df <- ECs %>%
  #     dplyr::mutate(metric = gsub(metric, pattern = "e:1:", replacement = "EC"),
  #                   above.max.conc = if_else(Estimate > high_dose, true = TRUE, false = FALSE)) %>%
  #     dplyr::full_join(.,slope) %>%
  #     dplyr::full_join(.,lower_asymp) %>%
  #     dplyr::mutate(strain = drug.drm$strain[[j]],
  #                   drug = drug.drm$drug[[j]])
  #   
  #   dose.response.parameters[[j]] <- EC.df
  #   
  # }
  # dose.response.parameters.df <- Reduce(rbind, dose.response.parameters)
  
  
  # Dose-response parameters
  ED.table.df <- data.frame(drug.drm$ED) 
  ED.table.df$params <- rownames(ED.table.df)
  rownames(ED.table.df) <- NULL
  ED.df <- ED.table.df %>%
    tidyr::separate(params, c("x","strain","metric"),sep = ":") %>%
    dplyr::select(-x) %>%
    dplyr::mutate(metric = paste0("EC",metric)) %>%
    dplyr::mutate(drug = toxicant)
  
  
  
  # Dose-response statistics
  strain.model.fit <- safe_LL4_strain(drug.drm$data[[1]])
  no.strain.model.fit <- safe_LL4_no_strain(drug.drm$data[[1]])
  summarized.model.fit <- summary(strain.model.fit$result)
  strain.coefs <- data.frame(summarized.model.fit$coefficients)
  strain.coefs$params <- rownames(strain.coefs)
  rownames(strain.coefs) <- NULL
  slope.df <- strain.coefs %>%
    tidyr::separate(params, c("parameter","strain"), sep = ":") %>%
    dplyr::mutate(strain = gsub(strain, pattern = "strain", replacement = "")) %>%
    dplyr::filter(parameter == "b") %>%
    dplyr::rename(metric = parameter) %>%
    dplyr::mutate(drug = toxicant) %>%
    dplyr::select(Estimate, Std..Error, metric, strain, drug, p.value)
  
  intercept <- strain.coefs %>%
    tidyr::separate(params, c("parameter","strain"), sep = ":") %>%
    # dplyr::mutate(strain = gsub(strain, pattern = "strain", replacement = "")) %>%
    dplyr::filter(parameter == "c") %>%
    dplyr::select(Estimate) %>%
    as.numeric()
  
  asymp.df <- strain.coefs %>%
    tidyr::separate(params, c("parameter","strain"), sep = ":") %>%
    dplyr::mutate(strain = gsub(strain, pattern = "strain", replacement = "")) %>%
    dplyr::filter(parameter == "d") %>%
    dplyr::rename(metric = parameter) %>%
    dplyr::mutate(drug = toxicant,
                  int = intercept,
                  lower.asymp = Estimate + int) %>%
    dplyr::select(Estimate, Std..Error, metric, lower.asymp, strain, drug, p.value)

  
  dose.response.parameters.df <- ED.df %>%
    dplyr::full_join(.,slope.df) %>%
    dplyr::full_join(.,asymp.df)
  
  
  drc.anova <- anova(no.strain.model.fit$result, strain.model.fit$result)
  
  
  EDcomps.df <- suppressMessages(data.frame(drc::EDcomp(strain.model.fit$result, c(10,10))))
  EDcomps.df$comps <- rownames(EDcomps.df)
  rownames(EDcomps.df) <- NULL
  ed.comps <- EDcomps.df %>%
    dplyr::mutate(drug = toxicant)
  
  
  Slope.comps.df <- suppressMessages(data.frame(drc::compParm(strain.model.fit$result, "b")))
  Slope.comps.df$comps <- rownames(Slope.comps.df)
  rownames(Slope.comps.df) <- NULL
  slope.comps <- Slope.comps.df %>%
    dplyr::mutate(drug = toxicant)
  
  
  
  return(list(overall.tx.LOAEL, strain.tx.LOAEL, 
              dose.response.parameters.df, 
              summarized.model.fit, drc.anova, ed.comps, slope.comps,
              cleaned.filtered.data))

}
dose.response.plots.only <- function(data, toxicant, SE.drc.dodge.width = 0.09, SE.drc.point.size = 0.4){
  raw.data <- data %>%
    dplyr::mutate(drug = toxicant)
  cleaned.data <- suppressMessages(cleanData(raw.data))
  
  
  retained.bleaches <- data.retention(cleaned.data[[2]])
  dose.check.plot(dat = cleaned.data[[1]])
  
  cleaned.filtered.data <- retained.bleaches %>%
    dplyr::left_join(.,cleaned.data[[1]]) %>%
    dplyr::group_by(strain, concentration_um) %>%
    dplyr::count() %>%
    dplyr::arrange(n) %>%
    dplyr::filter(n > 10) %>%
    dplyr::select(-n) %>%
    dplyr::left_join(.,cleaned.data[[1]])
  
  se.DRC(cleaned.filtered.data, dodge.width = SE.drc.dodge.width, point.size = SE.drc.point.size)
}
make.reps <- function(x){
  x$replicate <- seq(1,nrow(x))
  x %>%
    dplyr::rename(value = median_wormlength_um_reg)
}
genos <- data.table::fread("data/Genotype_Matrix.tsv") %>%
  as.data.frame()
H2 <- function(d){
  strain.fit <- lme4::lmer(data = d, formula = value ~ 1 + (1|strain))
  variances <- lme4::VarCorr(x = strain.fit)
  A <- as.data.frame(variances)
  Vg <- A$vcov[1]
  Ve <- A$vcov[2]
  H2 <- Vg/(Vg+Ve)
  return(H2)
}
h2 <- function(d, geno_matrix){
  pheno_strains <- unique(d$strain)
  A <- sommer::A.mat(t(geno_matrix[,colnames(geno_matrix) %in% pheno_strains]))
  
  df_y <- d %>%
    dplyr::arrange(strain) %>%
    dplyr::select(strain, value) %>%
    dplyr::mutate(strain = as.character(strain))
  h2_res <- sommer::mmer(value~1,
                         random=~vs(strain,Gu=A), 
                         data=df_y)
  h2 <- vpredict(h2_res, h2 ~ (V1) / ( V1+V2))[[1]][1]
  h2_SE <- vpredict(h2_res, h2 ~ (V1) / ( V1+V2))[[2]][1]
  
  h2_df <- data.frame(h2) %>%
    dplyr::mutate(h2.upper = h2 + h2_SE,
                  h2.lower = h2 - h2_SE)
  return(h2_df)
}
H2.bootstrapping.calc <- function(d, nreps = 100, boot = T){
  
  if(boot == T){
    # Broad-sense Heritability
    H2.point <- H2(d = d)
    h2.point <- h2(d = d, geno_matrix = genos)
    H2.boots <- list()
    for(i in 1:nreps) {
      if(i %% 10 == 0){
        print(paste0((i/nreps)*100,"%"))
      }
      #################################
      # Bootstrap within strain ##
      #################################
      nested <- d %>%
        dplyr::group_by(strain) %>%
        tidyr::nest()
      boot.strain <- list()
      for(j in 1:length(nested$strain)){
        boot.strain[[j]] <- nested$data[[j]][sample(seq(1:nrow(nested$data[[j]])),replace = T),] %>%
          dplyr::mutate(strain = nested$strain[[j]])
      }
      boot <- boot.strain %>%
        Reduce(rbind,.)
      
      ##################################
      ## Bootstrap agnostic of strain ##
      ##################################
      # boot <- d[sample(seq(1:nrow(d)), replace = T),]
      
      check <- boot %>%
        dplyr::group_by(strain) %>%
        dplyr::summarise(n())
      if(1 %in% check$`n()`){
        print("Only 1 Strain Sampled in Bootstrap - Skipping")
        next
      }
      # Broad-Sense Heritability
      H2.est <- H2(d = boot)
      H2.boots[i] <- H2.est
    }
    
    H2.boots.vec <- unlist(H2.boots)
    H2.quantiles <- c(quantile(H2.boots.vec, probs = seq(0,1,0.05)))
    H2.CI <- data.frame(H2.point, 
                        as.numeric(H2.quantiles[2]), 
                        as.numeric(H2.quantiles[21])) %>%
      `colnames<-`(c("H2.Point.Estimate","H2.5.perc","H2.95.perc"))
    
    
    
    return(list(H2.CI,H2.boots.vec,h2.point))
    
  } else {
    
    H2.point <- H2(d = d)
    # h2.point <- h2(d = d)
    H2.CI <- data.frame(H2.point, 
                        NA, 
                        NA) %>%
      `colnames<-`(c("H2.Point.Estimate","H2.5.perc","H2.95.perc"))
    return(H2.CI)
  }
  
}
drugHeritability <- function(rep.data, concentration_um){
  
  # Broad-sense Heritability
  no.reps <- rep.data %>%
    dplyr::group_by(strain) %>%
    dplyr::summarise(n = n())
  
  
  if(1 %in% no.reps$n){
    H2.est <- data.frame("NA", "NA", "NA") %>%
      dplyr::mutate(drug = unique(rep.data$drug),
                    concentration = concentration_um) %>%
      `colnames<-`(c("H2", "lower", "upper", "drug", "concentration_um"))
    
    h2.est <- data.frame("NA", "NA", "NA") %>%
      dplyr::mutate(drug = unique(rep.data$drug),
                    concentration = concentration_um) %>%
      `colnames<-`(c("h2", "upper", "lower", "drug", "concentration_um"))
    
    list(H2.est, h2.est)
    
  } else {

    herits <- suppressMessages(H2.bootstrapping.calc(d = rep.data, boot = T))
    
    H2.est <- as.data.frame(herits[[1]]) %>%
      dplyr::mutate(drug = unique(rep.data$drug),
                    concentration = concentration_um) %>%
      `colnames<-`(c("H2", "lower", "upper", "drug", "concentration_um"))
    
    H2.replicates <- as.data.frame(herits[[2]]) %>%
      dplyr::mutate(drug = unique(rep.data$drug),
                    concentration = concentration_um) %>%
      `colnames<-`(c("H2.rep", "drug", "concentration_um"))
    
    h2.est <- as.data.frame(herits[[3]]) %>%
      dplyr::mutate(drug = unique(rep.data$drug),
                    concentration = concentration_um) %>%
      `colnames<-`(c("h2", "upper", "lower", "drug", "concentration_um"))
    
    list(H2.est,H2.replicates,h2.est)
  }
}
heritability.calculation <- function(data, toxicant){
  # Data cleaning
  raw.data <- data %>%
    dplyr::mutate(drug = toxicant)
  cleaned.data <- suppressMessages(cleanData(raw.data))
  retained.bleaches <- data.retention(cleaned.data[[2]])
  cleaned.filtered.data <- retained.bleaches %>%
    dplyr::left_join(.,cleaned.data[[1]])
  
  dose.nested <- cleaned.filtered.data %>%
    dplyr::group_by(concentration_um) %>%
    tidyr::nest() 
  
  dose.nested$data <- purrr::map(dose.nested$data, make.reps)
  dose.data.nested.reps <- dose.nested %>%
    tidyr::unnest(cols = c(data)) %>%
    dplyr::group_by(concentration_um) %>%
    tidyr::nest()
  

  herits <- purrr::map2(dose.data.nested.reps$data, 
                        dose.data.nested.reps$concentration_um,
                        drugHeritability)
  
}
plot_HH <- function(herit.data){
  ggplot(herit.data, mapping = aes(x = H2, y = h2, colour = log(concentration_um))) + 
    theme_bw(base_size = 11) +
    geom_abline(slope = 1, alpha = 0.5, linetype = 3, size = 0.25) +
    # geom_smooth(method = "lm", se = F, colour = "black", linetype = 1, size = 0.5) +
    # stat_smooth(method="lm_right", fullrange=TRUE, colour = "black", se = F, size = 0.5, linetype = 1,alpha = 0.5) + 
    geom_errorbar(aes(ymin = h2.lower, ymax = h2.upper), size = 0.2) +
    geom_errorbarh(aes(xmin = H2.lower, xmax = H2.upper), size = 0.2) + 
    facet_grid(. ~ drug) + 
    scale_colour_gradient2(low = "darkblue", mid = "violet", high = "darkred", midpoint = 2, 
                           name = expression(log[10](Concentration (µM)))) + 
    theme(panel.grid = element_blank(),
          axis.title.y = element_text(size = 8.5, angle = 0, vjust = 0.5),
          axis.title.x = element_text(size = 8.5),
          axis.text = element_text(size = 7), 
          strip.text = element_text(size = 7),
          legend.position="bottom",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7)) + 
    scale_x_continuous(breaks = seq(0,1,0.5), limits = c(0,1)) + 
    scale_y_continuous(breaks = seq(0,1,0.5), limits = c(0,1)) + 
    labs(x = bquote(italic(H^2)),
         y = bquote(italic(h^2)))
}
gather.herit.ranges <- function(all.herits){
  print(all.herits[[1]][[3]]$drug)
  dose.herits <- list()
  for(i in 1:length(all.herits)){
    if(length(all.herits[[i]]) == 3){
      H2.estimate <- all.herits[[i]][[1]] %>%
        dplyr::rename(H2.upper = upper, H2.lower = lower)
      
      h2.estimate <- all.herits[[i]][[3]] %>%
        dplyr::rename(h2.upper = upper, h2.lower = lower)
      
      dose.herits[[i]] <- H2.estimate %>%
        dplyr::full_join(h2.estimate,.) %>%
        suppressMessages() %>%
        dplyr::select(drug, concentration_um, everything())
    } else {
      H2.estimate <- all.herits[[i]][[1]] %>%
        dplyr::rename(H2.upper = upper, H2.lower = lower)
      
      h2.estimate <- all.herits[[i]][[2]] %>%
        dplyr::rename(h2.upper = upper, h2.lower = lower)
      
      dose.herits[[i]] <- H2.estimate %>%
        dplyr::full_join(h2.estimate,.) %>%
        suppressMessages() %>%
        dplyr::select(drug, concentration_um, everything())
    }
    
  }
  dose.herits %>% Reduce(rbind,.)
}
h2.comps.to.individual.strains <- function(EC10, this.strain){
  strain.EC10.comp <- EC10 %>%
    dplyr::filter(drug %in% h2.comp.table$Toxicant) %>%
    dplyr::select(drug, Estimate) %>%
    dplyr::rename(Toxicant = drug) %>%
    dplyr::left_join(., h2.comp.table %>%
                       dplyr::select(Toxicant, `Top Heritable Dose`))
  strain.h2.comp.r2 <- summary(lm(data = strain.EC10.comp, formula = log10(`Top Heritable Dose`) ~ log10(Estimate)))
  lm.params <- list(r2 = format(strain.h2.comp.r2$r.squared, digits = 3),
                    pval = format(as.numeric(pf(strain.h2.comp.r2$fstatistic[1],
                                                strain.h2.comp.r2$fstatistic[2],
                                                strain.h2.comp.r2$fstatistic[3],
                                                lower.tail=FALSE)), digits=3)) %>% 
    data.frame() %>%
    dplyr::mutate(strain = this.strain)
  
}








MDHD.data.list <- as.list(ls(pattern = "MDHD"))
dose.responses.df <- purrr::map(MDHD.data.list, concatenate.assays) %>%
  Reduce(rbind,.) %>%
  dplyr::filter(!drug %in% c("Control","Bacteria OD","Mancozeb_Dilute")) %>%  # extra conditions not pertaining to current experiment
  dplyr::mutate(drug = str_to_sentence(drug)) %>%
  dplyr::mutate(drug = if_else(drug == "2,4-d", true = "2,4-D", false = drug),
                drug = if_else(drug %in% c("Cadmium","Cadmium dichloride"), true = "Cadmium chloride", false = drug),
                drug = if_else(drug %in% c("Copper","Copper(ii) chloride"), true = "Copper(II) chloride", false = drug),
                drug = if_else(drug == "Lead(ii) nitrate", true = "Lead(II) nitrate", false = drug),
                drug = if_else(drug %in% c("Methyl_mercury","Methylmercury dichloride"), true = "Methylmercury chloride", false = drug),
                drug = if_else(drug %in% c("Nickel","Nickel dichloride"), true = "Nickel chloride", false = drug),
                drug = if_else(drug == "Silver", true = "Silver nitrate", false = drug),
                drug = if_else(drug %in% c("Zinc","Zinc dichloride"), true = "Zinc chloride", false = drug),
                drug = if_else(drug %in% c("Manganese dichloride"), true = "Manganese(II) chloride", false = drug))
dose.responses.df <- dose.responses.df[-which(dose.responses.df$drug == "Chlorpyrifos" & dose.responses.df$Metadata_Experiment == "toxin17A"),] # inconsistent response in assay
dose.responses.df <- dose.responses.df[-which(dose.responses.df$drug == "Chlorfenapyr" & dose.responses.df$concentration_um >= 1.2),] # inconsistent response in assay
dose.responses.df[which(dose.responses.df$strain == "PD1074"),]$strain <- "N2"
dose.responses.df[which(dose.responses.df$strain == "ECA248"),]$strain <- "CB4855"

#############################
### Supplemental Figure 1 ###
#############################
control.wells.variance <- dose.responses.df %>%
  dplyr::filter(concentration_um == 0) %>%
  dplyr::group_by(Metadata_Experiment, bleach) %>%
  dplyr::summarise(mean.n = mean(n),
                   sd.n = sd(n),
                   cv.n = (sd.n/mean.n),
                   mean.median_wormlength_um = mean(median_wormlength_um),
                   sd.median_wormlength_um = sd(median_wormlength_um),
                   cv.worm.length = sd(median_wormlength_um)/mean(median_wormlength_um))

cv.n.cutoff = 0.6
annotated.control.wells.variance <- control.wells.variance %>%
  dplyr::mutate(censor = if_else(condition = cv.n > cv.n.cutoff, true = "CENSOR", false = "NO CENSOR")) %>% 
  dplyr::select(Metadata_Experiment, bleach, censor, cv.n, cv.worm.length) %>%
  dplyr::rename(censor.cv.n = cv.n, censor.cvWL = cv.worm.length)
cv.well.n.block.figure <- annotated.control.wells.variance %>%
  ggplot(., mapping = aes(x = Metadata_Experiment, y = censor.cv.n, colour = censor)) +
  theme_bw(base_size = 9) +
  geom_hline(yintercept = cv.n.cutoff, linetype = 3) + 
  geom_jitter() + 
  scale_color_manual(values = c("red","black"), name = "Data Censored") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "top",
        panel.grid = element_blank()) + 
  labs(y = "Well Titer Coefficient of Variation",
       x = "Assay")
cv.n.cv.WL.figure <- annotated.control.wells.variance %>%
  ggplot(., mapping = aes(y = censor.cv.n, x = censor.cvWL, color = as.factor(censor))) + 
  theme_bw(base_size = 9) + 
  geom_hline(yintercept = cv.n.cutoff, linetype = 3) + 
  geom_point() + 
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank()) + 
  scale_color_manual(values = c("red","black"), name = "Data Censored") + 
  labs(y = "Well Titer Coefficient of Variation",
       x = "Median Worm Length Coefficient of Variation")
cv.cleaning.supp.fig <- cowplot::plot_grid(cv.well.n.block.figure + theme(legend.position = "none"),
                   cv.n.cv.WL.figure, align = "h", labels = "AUTO")
ggsave(cv.cleaning.supp.fig + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "manuscript_figures/supp.fig.2.png", width = 7.5, height = 3.5)


# save(dose.responses.df, file = "data/all.dose.response.data.Rdata")

tx.nested.dose.responses <- dose.responses.df %>%
  dplyr::full_join(., annotated.control.wells.variance) %>%
  dplyr::filter(censor == "NO CENSOR") %>%
  dplyr::group_by(drug) %>%
  tidyr::nest()
dose.response.summaries <- purrr::map2(tx.nested.dose.responses$data,
                                       tx.nested.dose.responses$drug,
                                       dose.response.summary)
tr.dose.response.summaries <- purrr::transpose(dose.response.summaries)

dose.response.parameter.summaries <- purrr::map(.x = tr.dose.response.summaries,
                                                .f = function(x){Reduce(rbind,x)})
# 
# 
# test.DR <- dose.response.plots.only(tx.nested.dose.responses$data[[3]],
#                                     tx.nested.dose.responses$drug[[3]], 0.15)


###########
## LOAEL ##
###########

overall.LOAEL.summary.df <- dose.response.parameter.summaries[[1]] %>%
  # Reduce(rbind,.) %>%
  dplyr::mutate(LOAEL = as.numeric(LOAEL))
strain.LOAEL.summary.df <- dose.response.parameter.summaries[[2]] %>%
  # Reduce(rbind,.) %>%
  dplyr::select(tx, LOAEL, strain) %>%
  dplyr::group_by(tx, strain) %>% 
  dplyr::summarise(min(LOAEL)) %>%
  dplyr::mutate(`min(LOAEL)` = as.numeric(`min(LOAEL)`)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(values_from = `min(LOAEL)`, names_from = strain) %>%
  dplyr::rename(Toxicant = tx)
LOAEL.summary <- overall.LOAEL.summary.df %>%
  dplyr::select(tx, LOAEL) %>%
  dplyr::rename(Toxicant = tx,
                `Population-wide LOAEL` = LOAEL) %>%
  dplyr::full_join(.,strain.LOAEL.summary.df)

LOAEL.table <- LOAEL.summary %>%
  dplyr::full_join(.,tx.classes) %>%
  dplyr::arrange(big_class) %>%
  dplyr::select(big_class, everything(), -diluent, -class) %>%
  dplyr::rename(`Toxicant Class` = big_class)
write.csv(LOAEL.table, file = "manuscript_tables/supp.table.2.csv", row.names = F)




##############
## FIGURE 2 ##
##############

max.doses <- dose.responses.df %>%
  dplyr::select(drug, concentration_um) %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(max.dose = max(concentration_um))
 
all.DR.params <- dose.response.parameter.summaries[[3]]

EC10 <- dose.response.parameter.summaries[[3]] %>%
  dplyr::filter(metric == "EC10") %>%
  dplyr::mutate(null.model = is.na(Std..Error)) %>%
  dplyr::left_join(.,tx.classes %>%
                     dplyr::rename(drug = Toxicant)) %>%
  dplyr::left_join(.,max.doses) %>%
  dplyr::mutate(above.max.conc = Estimate >= max.dose) %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " "),
                flag = if_else(null.model == TRUE | Lower < 0, 
                               true = TRUE, false = FALSE))

  
above.max.concs <- EC10  %>%
  dplyr::group_by(drug, above.max.conc) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(above.max.conc == TRUE)


flagged <- EC10  %>%
  dplyr::group_by(drug, null.model) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(null.model == TRUE)
  
EC10.unflagged <- EC10 %>%
  dplyr::filter(null.model == FALSE,
                above.max.conc == FALSE)
EC10.filtered <- EC10.unflagged %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n == 8) %>%
  dplyr::left_join(.,EC10.unflagged) %>%
  dplyr::select(-n)

#####
## EC50
EC50 <- dose.response.parameter.summaries[[3]] %>%
  Reduce(rbind,.) %>%
  dplyr::filter(metric == "EC50") %>%
  dplyr::mutate(null.model = is.na(Std..Error)) %>%
  dplyr::left_join(.,tx.classes %>%
                     dplyr::rename(drug = Toxicant)) %>%
  dplyr::left_join(.,max.doses) %>%
  dplyr::mutate(above.max.conc = Estimate >= max.dose) %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " "),
                flag = if_else(null.model == TRUE | above.max.conc == TRUE | Lower < 0, true = TRUE, false = FALSE))

flagged <- EC50 %>%
  dplyr::group_by(drug, flag) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = flag, values_from = n) %>%
  dplyr::filter(!is.na(`TRUE`)) %>%
  dplyr::left_join(.,EC50) %>%
  dplyr::ungroup() %>%
  dplyr::select(drug, strain, null.model, above.max.conc, flag) %>%
  data.frame()

above.max.concs <- EC10 %>%
  dplyr::group_by(drug, above.max.conc) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = above.max.conc, values_from = n) %>%
  dplyr::filter(!is.na(`TRUE`))

EC50.unflagged <- EC50 %>%
  dplyr::filter(flag == FALSE)
EC50.filtered <- EC50.unflagged %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n == 8) %>%
  dplyr::left_join(.,EC50.unflagged) %>%
  dplyr::select(-n)


## EC90
EC90 <- dose.response.parameter.summaries[[3]] %>%
  Reduce(rbind,.) %>%
  dplyr::filter(metric == "EC90") %>%
  dplyr::mutate(null.model = is.na(Std..Error)) %>%
  dplyr::left_join(.,tx.classes %>%
                     dplyr::rename(drug = Toxicant)) %>%
  dplyr::left_join(.,max.doses) %>%
  dplyr::mutate(above.max.conc = Estimate >= max.dose) %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " "),
                flag = if_else(null.model == TRUE | above.max.conc == TRUE | Lower < 0, true = TRUE, false = FALSE))

flagged <- EC90 %>%
  dplyr::group_by(drug, flag) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = flag, values_from = n) %>%
  dplyr::filter(!is.na(`TRUE`)) %>%
  dplyr::left_join(.,EC90) %>%
  dplyr::ungroup() %>%
  dplyr::select(drug, strain, null.model, above.max.conc, flag) %>%
  data.frame()

above.max.concs <- EC90 %>%
  dplyr::group_by(drug, above.max.conc) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = above.max.conc, values_from = n) %>%
  dplyr::filter(!is.na(`TRUE`))

EC90.unflagged <- EC90 %>%
  dplyr::filter(flag == FALSE)
EC90.filtered <- EC90.unflagged %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n == 8) %>%
  dplyr::left_join(.,EC90.unflagged) %>%
  dplyr::select(-n)


reduced.EC10.filtered <- EC10.filtered %>%
  dplyr::select(drug, strain, metric, Estimate, Std..Error)
reduced.EC50.filtered <- EC50.filtered %>%
  dplyr::select(drug, strain, metric, Estimate, Std..Error)
reduced.EC90.filtered <- EC90.filtered %>%
  dplyr::select(drug, strain, metric, Estimate, Std..Error)

reduced.endpoint <- reduced.EC10.filtered %>%
  dplyr::full_join(., reduced.EC50.filtered) %>%
  dplyr::full_join(., reduced.EC90.filtered) %>%
  dplyr::arrange(metric, drug) %>%
  dplyr::rename(estimate = Estimate,
                SE = Std..Error)
write.csv(reduced.endpoint, "data/20220329_widmayer_endpoints.csv", row.names = F)
#####

#############################
### Supplemental Figure 3 ###
#############################\
proc_img_dir <- c("~/Documents/projects/toxin_dose_responses/data/supplemental.image.files")
dose.responses.df %>%
  dplyr::filter(drug %in% c("Arsenic trioxide")) %>%
  dplyr::distinct(Metadata_Experiment, Metadata_Plate, drug) %>% data.frame() 
proc_dose_data <- MDHD.toxin19A$processed_data
raw_dose_data <- MDHD.toxin19A$raw_data
plot_raw <- easyXpress::viewDose(raw_dose_data %>%
                                   dplyr::filter(Metadata_Plate == "p09"), 
                                 strain_name = "PD1074", 
                                 drug_name = "Arsenic trioxide", 
                                 proc_img_dir = proc_img_dir) 
ggsave(plot_raw, filename = "manuscript_figures/raw.arsenic.png", width = 20, height = 4)

#############################
### Supplemental Figure 4 ###
#############################
dose.responses.df %>%
  dplyr::filter(drug %in% c("Deltamethrin")) %>%
  dplyr::distinct(Metadata_Experiment, Metadata_Plate, drug) %>%
  dplyr::select(Metadata_Plate) %>% c()
proc_dose_data <- MDHD.toxin27A$processed_data
raw_dose_data <- MDHD.toxin27A$raw_data


PD1074_deltamethrin <- easyXpress::viewDose(proc_dose_data, 
                                  strain_name = "PD1074", 
                                  drug_name = "Deltamethrin", 
                                  proc_img_dir = proc_img_dir)
ggsave(PD1074_deltamethrin, filename = "manuscript_figures/PD1074_deltametrhin.png", width = 20, height = 4)


dose.responses.df %>%
  dplyr::filter(drug %in% c("Malathion")) %>%
  dplyr::distinct(Metadata_Experiment, Metadata_Plate, drug) %>% data.frame()
proc_dose_data <- MDHD.toxin23A$processed_data
raw_dose_data <- MDHD.toxin23A$raw_data
RC301_Malathion <- easyXpress::viewDose(proc_dose_data, 
                                        strain_name = "RC301", 
                                        drug_name = "Malathion", 
                                        proc_img_dir = proc_img_dir)
ggsave(RC301_Malathion, filename = "manuscript_figures/RC301_malathion.png", width = 20, height = 4)





#################################
### Supplemental Figures 5-29 ###
#################################
supp.DR.dodge.widths <- c(rep(0.15, 12),0.04, 0.08,rep(0.04, 2),0.09, 0.04, 0.02,rep(0.10, 2),0.04,0.02,rep(0.04,2))
for(i in 1:length(tx.nested.dose.responses$data)){
  supp.DR <- dose.response.plots.only(tx.nested.dose.responses$data[[i]],
                                      tx.nested.dose.responses$drug[[i]],
                                      supp.DR.dodge.widths[[i]])
  ggsave(supp.DR, filename = paste0("manuscript_figures/supp.fig.", i + 4,".png"), width = 6, height = 6)
}


################
### Figure 2 ###
################
DR.theme <- ggplot2::theme(axis.text = element_text(size = 7),
                           axis.title = element_text(size = 7), 
                           title = element_text(size = 7))

nickelDR <- dose.response.plots.only(tx.nested.dose.responses$data[[6]],
                                     tx.nested.dose.responses$drug[[6]],
                                     supp.DR.dodge.widths[[6]], 
                                     SE.drc.point.size = 0.05) + 
  DR.theme + 
  ggplot2::theme(legend.position = "none")


carbarylDR <- dose.response.plots.only(tx.nested.dose.responses$data[[15]],
                                     tx.nested.dose.responses$drug[[15]],
                                     supp.DR.dodge.widths[[15]], 
                                     SE.drc.point.size = 0.05) + 
  DR.theme + 
  ggplot2::theme(legend.position = "none")

x24DDR <- dose.response.plots.only(tx.nested.dose.responses$data[[19]],
                                   tx.nested.dose.responses$drug[[19]],
                                   supp.DR.dodge.widths[[19]], 
                                   SE.drc.point.size = 0.05) + 
  DR.theme + 
  ggplot2::theme(legend.position = "none")

pyraclostrobinDR <- dose.response.plots.only(tx.nested.dose.responses$data[[11]],
                                             tx.nested.dose.responses$drug[[11]],
                                   supp.DR.dodge.widths[[11]], 
                                   SE.drc.point.size = 0.05) + 
  DR.theme + 
  ggplot2::theme(legend.position = "none")


for.legend <- dose.response.plots.only(tx.nested.dose.responses$data[[6]],
                                       tx.nested.dose.responses$drug[[6]],
                                       supp.DR.dodge.widths[[6]], 
                                       SE.drc.point.size = 0.05)
strain.legend <- cowplot::get_legend(plot = for.legend + 
                                       ggplot2::theme(legend.text = element_text(size = 7),
                                                      legend.title = element_text(size = 9),
                                                      legend.position = "bottom",
                                                      legend.background = element_rect(fill = "white", linetype = "blank")))

DR.overview.fig <- cowplot::plot_grid(x24DDR, carbarylDR, nickelDR, pyraclostrobinDR, 
                   ncol = 2, 
                   nrow = 2, labels = "AUTO")
DR.overview.fig2 <- cowplot::plot_grid(DR.overview.fig, strain.legend, nrow = 2, rel_heights = c(1,0.15))

ggsave(DR.overview.fig2 + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "manuscript_figures/fig.2.png", width = 5, height = 5)





################
### Figure 3 ###
################
complete.EC10.plot <- EC10.filtered %>%
  dplyr::filter(!drug %in% c("Malathion","Deltamethrin")) %>%
  ggplot(., mapping = aes(y = drug, x = Estimate, xmin = Lower, xmax = Upper,
                          color = strain)) + 
  theme_bw(base_size = 11) + 
  geom_pointrange(position = position_quasirandom(groupOnX = FALSE),
                  size = 0.25) +
  scale_x_log10() + 
  scale_color_manual(values = strain_colors, name = "Strain") + 
  facet_grid(big_class~., scales = "free", space = "free") + 
  theme(strip.text.y = element_text(angle = 0, size = 9),
        axis.text = element_text(color = "black", size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top") + 
  labs(x = "EC10 Estimate (µM)")

complete.EC10.plot


n.EC.comp.tests <- dose.response.parameter.summaries[[6]] %>%
  dplyr::filter(drug %in% EC10.filtered$drug) %>%
  dplyr::filter(!drug %in% c("Deltamethrin","Malathion")) %>%  
  tidyr::separate(comps, c("strains","fracs"), sep = ":") %>%
  tidyr::separate(strains, c("strain1","strain2"), sep = "/") %>%
  dplyr::filter(strain2 == "N2" | strain1 == "N2") %>%
  nrow()
BF <- 0.05/n.EC.comp.tests
options(scipen = 999999)


EC10.relative.potency <- dose.response.parameter.summaries[[6]] %>%
  tidyr::separate(comps, c("strains","fracs"), sep = ":") %>%
  tidyr::separate(strains, c("strain1","strain2"), sep = "/") %>%
  dplyr::filter(drug %in% unique(EC10.filtered$drug),
                strain2 == "N2" | strain1 == "N2") %>%
  dplyr::mutate(N2.normed.estimate = if_else(strain1 == "N2", true = 1/Estimate, false = Estimate)) %>%
  dplyr::mutate(N2.normed.diff = N2.normed.estimate-1,
                focal.strain = if_else(strain1 == "N2", true = strain2, false = strain1)) %>%
  dplyr::left_join(.,tx.classes %>%
                     dplyr::rename(drug = Toxicant)) %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " ")) %>%
  dplyr::mutate(start = 1,
                sig = if_else(condition = p.value < BF, true = "SIG", false = "NONSIG"))

EC10.relative.potency.figure <- EC10.relative.potency %>%
  dplyr::filter(!drug %in% c("Deltamethrin","Malathion")) %>%
  ggplot(., mapping = aes(y = drug, x = N2.normed.diff, 
                          xmin = N2.normed.diff-Std..Error, 
                          xmax = N2.normed.diff+Std..Error,
                          color = focal.strain,
                          alpha = sig)) + 
  theme_bw(base_size = 11) +
  geom_vline(xintercept = 0, linetype = 3, colour = "orange") + 
  geom_pointrange(position = position_dodge(width = 0.2),
                  size = 0.25) + 
  scale_color_manual(values = strain_colors[c(1,3:8)], name = "Strain") +
  scale_alpha_manual(values = c(0.15,1), guide = "none") +
  facet_grid(big_class~., scales = "free", space = "free") + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 9),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0,size = 9),
        legend.position = "top") + 
  # xlim(c(-2,2)) +
  labs(x = "Resistance Compared to Reference (N2)")
EC10.relative.potency.figure

fig.2.legend <- cowplot::get_legend(complete.EC10.plot)
EC10.fig.2 <- cowplot::plot_grid(complete.EC10.plot +  
                     theme(strip.background = element_blank(),strip.text.y = element_blank(),
                           plot.margin = unit(c(0.5,0.5,0.25,0.5), "cm"),
                           legend.position = "none"),
                   EC10.relative.potency.figure + 
                     theme(plot.margin = unit(c(0.5,0.5,0.25,0.5), "cm"),
                           legend.position = "none"), 
                   rel_widths = c(0.9,1), labels = "AUTO")
complete.EC10.fig.2.plot <- cowplot::plot_grid(EC10.fig.2, fig.2.legend, rel_heights = c(20,2), ncol = 1)
complete.EC10.fig.2.plot
ggsave(complete.EC10.fig.2.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), filename = "manuscript_figures/fig.3.png", width = 7.5, height = 5)



## EC10 Stats by Toxicant Class ##
EC10.aov.nested <- EC10.filtered %>%
  dplyr::filter(!drug %in% c("Deltamethrin","Malathion")) %>%
  dplyr::select(drug, Estimate, strain, big_class) %>%
  dplyr::group_by(big_class) %>%
  tidyr::nest()

EC10.aov <- function(estimate, class){
  if(length(unique(estimate$drug)) > 1){
    big.class.aov <- aov(estimate, formula = Estimate ~ strain + drug)
    aov.df <- data.frame(summary(big.class.aov)[[1]][[5]][1],
                         summary(big.class.aov)[[1]][[5]][2])
    colnames(aov.df) <- c("strain.p","tox.p")
    aov.df$class <- class
    return(list(aov.df, TukeyHSD(big.class.aov)))
  } else {
    return("Only one toxicant!")
  }
  
}
EC10.aov.list <- purrr::map2(EC10.aov.nested$data,
            EC10.aov.nested$big_class,
            EC10.aov)
EC10.aov.list.tr <- EC10.aov.list %>%
  purrr::keep(., is.list) %>%
  purrr::transpose()

# Overall anova results for EC10s
Reduce(rbind,EC10.aov.list.tr[[1]])

# Fungicide Tukey's HSD Results
EC10.aov.list.tr[[2]][[4]][[2]] %>%
  data.frame() %>%
  dplyr::filter(p.adj < 0.05)
# Herbicide Tukey's HSD Results
EC10.aov.list.tr[[2]][[1]][[2]] %>%
  data.frame() %>%
  dplyr::filter(p.adj < 0.05)

# Insecticide Tukey's HSD Results
EC10.aov.list.tr[[2]][[2]][[2]] %>%
  data.frame() %>%
  dplyr::filter(p.adj > 0.05)

# Metal Tukey's HSD Results
EC10.aov.list.tr[[2]][[3]][[2]] %>%
  data.frame() %>%
  dplyr::filter(p.adj < 0.05)





## EC10 Relative Potency Tests ##
EC10.relative.potency.supp.table <- dose.response.parameter.summaries[[6]] %>%
  dplyr::filter(!drug %in% c("Deltamethrin","Malathion")) %>%
  tidyr::separate(comps, c("strains","fracs"), sep = ":") %>%
  tidyr::separate(strains, c("strain1","strain2"), sep = "/") %>%
  dplyr::mutate(`p < 0.05` = if_else(condition = p.value < 0.05, true = "T", false = "F"),
                p.value = round(p.value,6),
                p.value = if_else(p.value < 0.000001, true =  "<< 0.00001", false = as.character(p.value))) %>%
  dplyr::rename(`Relative Potency Estimate` = Estimate,
                SE = Std..Error,
                `Fraction Comparison` = fracs,
                Toxicant = drug) %>%
  dplyr::select(Toxicant, strain1, strain2, everything()) %>%
  data.frame()
write.csv(EC10.relative.potency.supp.table, "manuscript_tables/supp.table.3.csv", row.names = F)

sig.EC.comps <- dose.response.parameter.summaries[[6]] %>%
  dplyr::filter(drug %in% EC10.filtered$drug) %>%
  dplyr::filter(!drug %in% c("Deltamethrin","Malathion")) %>%
  tidyr::separate(comps, c("strains","fracs"), sep = ":") %>%
  tidyr::separate(strains, c("strain1","strain2"), sep = "/") %>%
  dplyr::filter(strain2 == "N2" | strain1 == "N2") %>%
  dplyr::mutate(N2.normed.estimate = if_else(strain1 == "N2", true = 1/Estimate, false = Estimate)) %>%
  dplyr::mutate(N2.normed.diff = N2.normed.estimate-1,
                focal.strain = if_else(strain1 == "N2", true = strain2, false = strain1)) %>%
  dplyr::left_join(.,tx.classes %>%
                     dplyr::rename(drug = Toxicant)) %>%
  dplyr::filter(p.value < BF)

nrow(sig.EC.comps)

sig.EC.comps %>%
  dplyr::group_by(drug) %>%
  dplyr::count()

strain.resistance <- dose.response.parameter.summaries[[6]] %>%
  dplyr::filter(drug %in% EC10.filtered$drug) %>%
  dplyr::filter(!drug %in% c("Deltamethrin","Malathion")) %>%
  tidyr::separate(comps, c("strains","fracs"), sep = ":") %>%
  tidyr::separate(strains, c("strain1","strain2"), sep = "/") %>%
  dplyr::filter(strain2 == "N2" | strain1 == "N2") %>%
  dplyr::mutate(N2.normed.estimate = if_else(strain1 == "N2", true = 1/Estimate, false = Estimate)) %>%
  dplyr::mutate(N2.normed.diff = N2.normed.estimate-1,
                focal.strain = if_else(strain1 == "N2", true = strain2, false = strain1)) %>%
  dplyr::left_join(.,tx.classes %>%
                     dplyr::rename(drug = Toxicant)) %>%
  dplyr::mutate(sens = if_else(N2.normed.diff > 0, "resistant", "sensitive")) %>%
  dplyr::group_by(focal.strain, sens) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = sens, values_from = n) %>%
  dplyr::mutate(ratio = resistant/sensitive)


strain.resistance.sig <- sig.EC.comps %>%
  dplyr::mutate(N2.normed.estimate = if_else(strain1 == "N2", true = 1/Estimate, false = Estimate)) %>%
  dplyr::mutate(N2.normed.diff = N2.normed.estimate-1,
                focal.strain = if_else(strain1 == "N2", true = strain2, false = strain1)) %>%
  dplyr::left_join(.,tx.classes %>%
                     dplyr::rename(drug = Toxicant)) %>%
  dplyr::mutate(sens = if_else(N2.normed.diff > 0, "resistant", "sensitive")) %>%
  dplyr::group_by(focal.strain, sens) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = sens, values_from = n)
strain.resistance.sig[is.na(strain.resistance.sig)] <- 0
strain.resistance.sig %<>%
  dplyr::mutate(ratio = sensitive/resistant)

total.resistant <- sum(strain.resistance.sig$resistant)
strain.resistance.sig %>%
  dplyr::mutate(pct.of.resistant.strains = resistant/total.resistant)
total.sensitive <- sum(strain.resistance.sig$sensitive)
strain.resistance.sig %>%
  dplyr::mutate(pct.of.sensitive.strains = sensitive/total.sensitive)


resistance <- strain.resistance.sig %>%
  dplyr::select(focal.strain, resistant) %>%
  dplyr::mutate(not.resistant = 21-resistant)

strain.resistance.fisher.reps <- list()
for(i in 1:1000){
  rep <- fisher.test(resistance[,2:3], simulate.p.value = T)
  strain.resistance.fisher.reps[[i]] <- rep$p.value
}
c(mean(unlist(strain.resistance.fisher.reps)),sd(unlist(strain.resistance.fisher.reps)))


sensitive <- strain.resistance.sig %>%
  dplyr::select(focal.strain, sensitive) %>%
  dplyr::mutate(not.sensitive = 21-sensitive)
strain.sensitive.fisher.reps <- list()
for(i in 1:1000){
  rep <- fisher.test(sensitive[,2:3], simulate.p.value = T)
  strain.sensitive.fisher.reps[[i]] <- rep$p.value
}
c(mean(unlist(strain.sensitive.fisher.reps)),sd(unlist(strain.sensitive.fisher.reps)))



relative.potency.table <- strain.resistance %>%
  dplyr::mutate(ratio = signif(ratio,2)) %>%
  dplyr::rename(Strain = focal.strain,
                Resistant = resistant,
                Sensitive = sensitive,
                `Resistant:Sensitive` = ratio)
write.csv(relative.potency.table, "manuscript_tables/table.1.csv", row.names = F)






################
### Figure 4 ###
################

slopes.filtered <- EC10.filtered %>%
  dplyr::select(strain, drug) %>%
  dplyr::left_join(.,all.DR.params) %>%
  dplyr::left_join(.,tx.classes %>%
                     dplyr::rename(drug = Toxicant)) %>%
  dplyr::filter(metric == "b") %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " "))


n.slope.comp.tests <- dose.response.parameter.summaries[[7]] %>%
  dplyr::filter(!drug %in% c("Deltamethrin","Malathion")) %>%
  tidyr::separate(comps, c("strain1","strain2"), sep = "/") %>%
  dplyr::mutate(strain1 = gsub(strain1, pattern = "strain", replacement = "")) %>%
  dplyr::mutate(strain2 = gsub(strain2, pattern = "strain", replacement = "")) %>% 
  dplyr::filter(drug %in% unique(EC10.filtered$drug),
                strain2 == "N2" | strain1 == "N2") %>%
  dplyr::mutate(N2.normed.estimate = if_else(strain1 == "N2", true = 1/Estimate, false = Estimate)) %>% 
  nrow()
BF <- 0.05/n.slope.comp.tests
options(scipen = 999999)

relative.slope <- dose.response.parameter.summaries[[7]] %>%
  tidyr::separate(comps, c("strain1","strain2"), sep = "/") %>%
  dplyr::mutate(strain1 = gsub(strain1, pattern = "strain", replacement = "")) %>%
  dplyr::mutate(strain2 = gsub(strain2, pattern = "strain", replacement = "")) %>% 
  dplyr::filter(drug %in% unique(EC10.filtered$drug),
                strain2 == "N2" | strain1 == "N2") %>%
  dplyr::mutate(N2.normed.estimate = if_else(strain1 == "N2", true = 1/Estimate, false = Estimate)) %>%
  dplyr::mutate(N2.normed.diff = N2.normed.estimate-1,
                focal.strain = if_else(strain1 == "N2", true = strain2, false = strain1)) %>%
  dplyr::left_join(.,tx.classes %>%
                     dplyr::rename(drug = Toxicant)) %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " ")) %>%
  dplyr::mutate(start = 1,
                sig = if_else(condition = p.value < BF, true = "SIG", false = "NONSIG"))




complete.slope.plot <- slopes.filtered %>%
  dplyr::filter(!drug %in% c("Deltamethrin","Malathion")) %>%
  ggplot(., mapping = aes(y = drug, x = Estimate, xmin = Estimate-Std..Error, xmax = Estimate+Std..Error,
                          color = strain)) + 
  theme_bw(base_size = 11) + 
  geom_pointrange(position = position_quasirandom(groupOnX = FALSE),
                  size = 0.25) +
  scale_color_manual(values = strain_colors, name = "Strain") + 
  facet_grid(big_class~., scales = "free", space = "free") + 
  theme(strip.text.y = element_text(angle = 0, size = 9),
        axis.text = element_text(color = "black", size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top") + 
  labs(x = "Estimated Slope")
complete.slope.plot

relative.slope.figure <- relative.slope %>%
  dplyr::filter(!drug %in% c("Deltamethrin","Malathion")) %>%
  ggplot(., mapping = aes(y = drug, x = N2.normed.diff, 
                          xmin = N2.normed.diff-Std..Error, 
                          xmax = N2.normed.diff+Std..Error,
                          color = focal.strain,
                          alpha = sig)) + 
  theme_bw(base_size = 11) +
  geom_vline(xintercept = 0, linetype = 3, colour = "orange") + 
  geom_pointrange(position = position_dodge(width = 0.2),
                  size = 0.25) + 
  scale_color_manual(values = strain_colors[c(1,3:8)], name = "Strain") +
  scale_alpha_manual(values = c(0.15,1), guide = "none") +
  facet_grid(big_class~., scales = "free", space = "free") + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 9),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0,size = 9),
        legend.position = "top") + 
  labs(x = "Relative Slope Compared to Reference (N2)")
relative.slope.figure
fig.3.legend <- cowplot::get_legend(complete.slope.plot)
slope.fig.3 <- cowplot::plot_grid(complete.slope.plot +  
                                   theme(strip.background = element_blank(),strip.text.y = element_blank(),
                                         plot.margin = unit(c(0.5,0.5,0.25,0.5), "cm"),
                                         legend.position = "none"),
                                 relative.slope.figure + 
                                   theme(plot.margin = unit(c(0.5,0.5,0.25,0.5), "cm"),
                                         legend.position = "none"), 
                                 rel_widths = c(0.9,1), labels = "AUTO")
complete.slope.fig.3.plot <- cowplot::plot_grid(slope.fig.3, fig.3.legend, rel_heights = c(20,2), ncol = 1)
complete.slope.fig.3.plot
ggsave(complete.slope.fig.3.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), filename = "manuscript_figures/fig.4.png", width = 7.5, height = 5)


slope.aov.nested <- slopes.filtered %>%
  dplyr::filter(!drug %in% c("Deltamethrin","Malathion")) %>%
  dplyr::select(drug, Estimate, strain, big_class) %>%
  dplyr::group_by(big_class) %>%
  tidyr::nest()

slope.aov <- function(estimate, class){
  if(length(unique(estimate$drug)) > 1){
    big.class.aov <- aov(estimate, formula = Estimate ~ strain + drug)
    aov.df <- data.frame(summary(big.class.aov)[[1]][[5]][1],
                         summary(big.class.aov)[[1]][[5]][2])
    colnames(aov.df) <- c("strain.p","tox.p")
    aov.df$class <- class
    return(list(aov.df, TukeyHSD(big.class.aov)))
  } else {
    return("Only one toxicant!")
  }
  
}
slope.aov.list <- purrr::map2(slope.aov.nested$data,
                             slope.aov.nested$big_class,
                             slope.aov)
slope.aov.list.tr <- slope.aov.list %>%
  purrr::keep(., is.list) %>%
  purrr::transpose()


# Overall anova results for slope\
Reduce(rbind,slope.aov.list.tr[[1]])


# Fungicide Tukey's HSD Results
slope.aov.list.tr[[2]][[4]][[2]] %>%
  data.frame() %>%
  dplyr::filter(p.adj < 0.05)

# Herbicide Tukey's HSD Results
slope.aov.list.tr[[2]][[1]][[2]] %>%
  data.frame() %>%
  dplyr::filter(p.adj < 0.05)

# Insecticide Tukey's HSD Results
slope.aov.list.tr[[2]][[2]][[2]] %>%
  data.frame() %>%
  dplyr::filter(p.adj < 0.05)

# Metal Tukey's HSD Results
slope.aov.list.tr[[2]][[3]][[2]] %>%
  data.frame() %>%
  dplyr::filter(p.adj < 0.05)




## Slope Relative Potency Tests ##
dose.response.parameter.summaries[[7]] %>%
  dplyr::filter(!drug %in% c("Deltamethrin","Malathion")) %>%
  tidyr::separate(comps, c("strain1","strain2"), sep = "/") %>%
  dplyr::mutate(strain1 = gsub(strain1, pattern = "strain", replacement = "")) %>%
  dplyr::mutate(strain2 = gsub(strain2, pattern = "strain", replacement = "")) %>% 
  dplyr::filter(drug %in% unique(EC10.filtered$drug),
                strain2 == "N2" | strain1 == "N2")

sig.EC.comps <- dose.response.parameter.summaries[[7]] %>%
  dplyr::filter(!drug %in% c("Deltamethrin","Malathion")) %>%
  tidyr::separate(comps, c("strain1","strain2"), sep = "/") %>%
  dplyr::mutate(strain1 = gsub(strain1, pattern = "strain", replacement = "")) %>%
  dplyr::mutate(strain2 = gsub(strain2, pattern = "strain", replacement = "")) %>% 
  dplyr::filter(drug %in% unique(EC10.filtered$drug),
                strain2 == "N2" | strain1 == "N2") %>%
  dplyr::mutate(N2.normed.estimate = if_else(strain1 == "N2", true = 1/Estimate, false = Estimate)) %>%
  dplyr::mutate(N2.normed.diff = N2.normed.estimate-1,
                focal.strain = if_else(strain1 == "N2", true = strain2, false = strain1)) %>%
  dplyr::left_join(.,tx.classes %>%
                     dplyr::rename(drug = Toxicant)) %>%
  dplyr::filter(p.value < BF)

nrow(sig.EC.comps)


sig.slope.comps <- relative.slope %>%
  dplyr::filter(sig == "SIG")

nrow(sig.slope.comps)
sig.slope.comps %>%
  dplyr::group_by(big_class) %>%
  dplyr::count()

sig.slope.comps %>%
  dplyr::group_by(focal.strain) %>%
  dplyr::count()



slope.comp.sig <- sig.slope.comps %>%
  dplyr::mutate(sens = if_else(N2.normed.diff > 0, "steep", "shallow")) %>%
  dplyr::group_by(focal.strain, sens) %>%
  dplyr::count() %>%
  tidyr::pivot_wider(names_from = sens, values_from = n)


resistance <- slope.comp.sig %>%
  dplyr::select(focal.strain, shallow) %>%
  dplyr::mutate(not.shallow = 21-shallow)

strain.resistance.fisher.reps <- list()
for(i in 1:1000){
  rep <- fisher.test(resistance[,2:3], simulate.p.value = T)
  strain.resistance.fisher.reps[[i]] <- rep$p.value
}
c(mean(unlist(strain.resistance.fisher.reps)),sd(unlist(strain.resistance.fisher.reps)))


sensitive <- slope.comp.sig %>%
  dplyr::select(focal.strain, steep) %>%
  dplyr::mutate(not.steep = 21-steep)
strain.sensitive.fisher.reps <- list()
for(i in 1:1000){
  rep <- fisher.test(sensitive[,2:3], simulate.p.value = T)
  strain.sensitive.fisher.reps[[i]] <- rep$p.value
}
c(mean(unlist(strain.sensitive.fisher.reps)),sd(unlist(strain.sensitive.fisher.reps)))





################
### Figure 5 ###
################
all.heritabilities <- purrr::map2(tx.nested.dose.responses$data,
                                  tx.nested.dose.responses$drug,
                                  heritability.calculation)
summarized.heritabilities <- purrr::map(all.heritabilities, gather.herit.ranges)
summarized.heritabilities.df <- Reduce(rbind,summarized.heritabilities) %>%
  dplyr::filter(h2 != "NA",
                drug %in% unique(EC10.filtered$drug),
                !drug %in% c("Deltamethrin","Malathion")) %>% 
  dplyr::mutate(h2 = as.numeric(h2),
                H2 = as.numeric(H2),
                h2.upper= as.numeric(h2.upper),
                h2.lower= as.numeric(h2.lower),
                H2.lower= as.numeric(H2.lower),
                H2.upper= as.numeric(H2.upper)) %>%
  dplyr::left_join(.,tx.classes %>%
                     dplyr::rename(drug = Toxicant)) %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " "))



flame.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Flame Retardant") %>%
  plot_HH(.)

herb.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Herbicide" ) %>%
  plot_HH(.)

metal.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Metal") %>%
  plot_HH(.)

fung.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Fungicide") %>%
  plot_HH(.)

insect.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Insecticide") %>%
  plot_HH(.)

fig.4.legend <- cowplot::get_legend(insect.HH)
A <- cowplot::plot_grid(flame.HH + theme(legend.position = "none"), 
                        herb.HH + theme(legend.position = "none"), 
                        labels = c("A","B"), rel_widths = c(0.45, 1))
B <- cowplot::plot_grid(fung.HH + 
                          theme(legend.position = "none"),
                        labels = c("C"))
C <- cowplot::plot_grid(insect.HH + 
                          theme(legend.position = "none") + 
                          facet_wrap(facets = drug ~., nrow = 2), 
                        metal.HH + 
                          theme(legend.position = "none") +
                          facet_wrap(facets = drug ~., nrow = 3), 
                        labels = c("D","E"),
                        ncol = 1, rel_heights = c(0.65,1))
heritability.summary.plot <- cowplot::plot_grid(A, B, C, fig.4.legend, ncol = 1, rel_heights = c(0.3, 
                                                                                                 0.28, 
                                                                                                 1, 
                                                                                                 0.1))
ggsave(heritability.summary.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "manuscript_figures/fig.5.png", width = 5, height = 8)



################
### Figure 6 ###
################

nested.sum.herits <- summarized.heritabilities.df %>%
  dplyr::group_by(drug) %>%
  tidyr::nest()
ranked.herits <- purrr::map2(nested.sum.herits$data,
            nested.sum.herits$drug, 
            function(x,y){
              x %>%
                dplyr::mutate(rank = seq(1:nrow(x))) %>%
                dplyr::mutate(drug = y)}) %>%
  Reduce(rbind,.)
max.H2 <- ranked.herits %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(H2 = max(H2)) %>%
  dplyr::left_join(., ranked.herits) %>%
  dplyr::select(drug, concentration_um, H2, rank) %>%
  dplyr::arrange(rank) %>%
  data.frame()
max.H2
max.h2 <- ranked.herits %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(h2 = max(h2)) %>%
  dplyr::left_join(., ranked.herits) %>%
  dplyr::select(drug, concentration_um, h2, rank) %>%
  dplyr::arrange(rank) %>%
  data.frame()
max.h2

max.h2 %>%
  dplyr::group_by(rank) %>%
  dplyr::summarise(n = n())

avg.herits <- ranked.herits %>%
  dplyr::filter(concentration_um != 0) %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(mean(H2), 
                   sd(H2),
                   mean(h2),
                   sd(h2))
avg.herits %>%
  dplyr::filter(`mean(H2)` %in% range(avg.herits$`mean(H2)`))

avg.herits %>%
  dplyr::filter(`mean(h2)` %in% range(avg.herits$`mean(h2)`))


  
EC.10.diff.H2 <- EC10.filtered %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(mean.EC10 = mean(Estimate)) %>%
  dplyr::full_join(ranked.herits,.) %>%
  dplyr::mutate(herit.EC10.diff = abs(concentration_um - mean.EC10))

herit.at.EC10.closest.dose <- EC.10.diff.H2 %>%
  dplyr::select(drug, concentration_um, rank, mean.EC10, herit.EC10.diff) %>%
  dplyr::filter(!is.na(mean.EC10)) %>%
  dplyr::group_by(drug) %>% 
  dplyr::summarise(herit.EC10.diff = min(herit.EC10.diff)) %>%
  dplyr::left_join(.,EC.10.diff.H2) %>%
  dplyr::select(drug, concentration_um, H2, h2) %>%
  dplyr::left_join(., max.H2 %>%
                     dplyr::rename(max.H2 = H2,
                                   max.H2.conc = concentration_um) %>%
                     dplyr::select(-rank)) %>%
  dplyr::left_join(., max.h2 %>%
                     dplyr::rename(max.h2 = h2,
                                   max.h2.conc = concentration_um))

h2.comp.table <- herit.at.EC10.closest.dose %>%
  dplyr::filter(drug %in% unique(EC10.filtered$drug),
                !drug %in% c("Deltamethrin","Malathion")) %>%
  dplyr::select(-max.h2.conc) %>%
  dplyr::mutate(`EC10 Dose - Top Heritable Dose` =  concentration_um - max.H2.conc) %>%
  dplyr::select(drug, concentration_um,max.H2.conc, `EC10 Dose - Top Heritable Dose`, H2, max.H2, h2, max.h2, rank) %>%
  dplyr::rename(Toxicant = drug,
                `Closest Dosage to EC10` = concentration_um,
                `Top Heritable Dose` = max.H2.conc,
                `H2 EC10` = H2,
                `Top H2` = max.H2,
                `h2 EC10` = h2,
                `Top h2` = max.h2,
                `Top Heritable Dose Rank` = rank)
colnames(h2.comp.table)

strain.nested.EC10 <- EC10 %>%
  dplyr::filter(drug %in% unique(EC10.filtered$drug),
                !drug %in% c("Deltamethrin","Malathion")) %>%
  dplyr::select(drug, strain, Estimate) %>%
  dplyr::group_by(strain) %>%
  tidyr::nest()



# 
# purrr::map2(strain.nested.EC10$data,
#             strain.nested.EC10$strain,
#             h2.comps.to.individual.strains)

h2.comp.r2 <- summary(lm(data = h2.comp.table, formula = log10(`Top Heritable Dose`) ~ log10(`Closest Dosage to EC10`)))
l <- list(r2 = format(h2.comp.r2$r.squared, digits = 3),
          pval = format(as.numeric(pf(h2.comp.r2$fstatistic[1],
                           h2.comp.r2$fstatistic[2],
                           h2.comp.r2$fstatistic[3],
                           lower.tail=FALSE)), digits=3)
)

eq <- substitute(italic(r)^2 == r2*","~~italic(p) == pval,l)
eqstr <- as.character(as.expression(eq))
h2.comp.figure <- ggplot() +
  theme_bw(base_size = 11) + 
  geom_smooth(h2.comp.table, mapping = aes(x = log10(`Closest Dosage to EC10`), 
                                           y = log10(`Top Heritable Dose`)),
              method = "lm", colour = "darkred") + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_jitter(h2.comp.table, mapping = aes(x = log10(`Closest Dosage to EC10`), 
                                           y = log10(`Top Heritable Dose`)),
              width = 0.06, height = 0.06) + # some values are overplotted on log scale
  
  annotate(geom = "text", x = 1.5, y = 3, label = eqstr, parse = TRUE, colour = "darkred", size = 3) + 
  theme(panel.grid = element_blank()) + 
  labs(x = expression(log[10](`Closest Dosage to EC10 `(µM))),
       y = expression(log[10](`Most Heritable Dosage `(µM))))
ggsave(h2.comp.figure + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "manuscript_figures/fig.6.png", width = 3.5, height = 3.5)


pivoted.top.herits <- max.H2 %>%
  tidyr::pivot_longer(cols = H2, names_to = "measurement") %>%
  dplyr::full_join(., max.h2 %>%
                     tidyr::pivot_longer(cols = h2, names_to = "measurement"))

pivoted.herits.at.EC10 <- herit.at.EC10.closest.dose %>%
  dplyr::rename(H2.EC10 = H2, 
                h2.EC10 = h2) %>%
  tidyr::pivot_longer(cols = c(H2.EC10,h2.EC10), names_to = "measurement")
  
top.v.EC10.herit <- pivoted.top.herits %>%
  dplyr::full_join(., pivoted.herits.at.EC10)

top.v.EC10.herit %>%
  dplyr::arrange(drug) %>%
  dplyr::filter(measurement %in% c("h2", "h2.EC10")) %>%
  dplyr::select(-rank) %>%
  dplyr::group_by(drug) 

