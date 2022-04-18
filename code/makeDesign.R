
#Load necessary packages
library(tidyverse)
library(easyXpress)
##############################################
###            Generate Thumbs             ### 
##############################################
easyXpress::tidyProject("/Volumes/ECA_Image/sam/20210701-toxin28A/", rm_other = TRUE)


###############################################
###            Assign Variables             ### 
###############################################
n_plates <- 84 # how many plates in assay
assay <- "48h" # 48h or 96h
food_od <- 10 # typically 5 - 15 (measured in OD)
food_t <- "15hHB101_20210223" # growth hour, e. coli strain, date of prep
rep_plates <- 4 # number of replicate plates in assay
bleaches <- 3 # number of bleaches performed in assay
levels <- 12 # number of doses per drug
drugs <- c("Propoxur",
           "Chlorothalonil",
           "Mancozeb",
           "Silver nitrate",
           "2,4-D",
           "Aldicarb",
           "Bacteria OD") # in order administered to plates
diluent_concentration_perc <- 1 # 1% dilutent should be enters as 1 not 0.01
drug_prep_method <- "serial; Rscript"
drug_sheet <- read.csv("~/Documents/projects/NemaDose/drug_stock_sheet.csv") # read in drug sheet

###############################################
###            Make Plates and Wells        ### 
###############################################
Metadata_Plate = stringr::str_pad(rep(1:n_plates, each = 96), width = 2, side = "left", pad = 0) # NOTE: p prefix not added yet
Metadata_Well = tibble(row = rep(LETTERS[1:8], each = 12)) %>%
  dplyr::mutate(col = stringr::str_pad(rep(1:12, 8), width = 2, side = "left", pad = 0)) %>%
  dplyr::mutate(full = paste0(row, col)) %>%
  dplyr::pull(full) %>%
  rep(., times = n_plates)

###############################################
###            Make drug doses              ### 
###############################################
clean_drug_sheet <- drug_sheet %>%
  dplyr::select(Drug, High_Dose_uM, drug_stock_prep_date, diluent, DilutionFactors)

drug_doses <- tibble(drug = rep(drugs, each = 96*rep_plates*bleaches)) %>%
  dplyr::left_join(clean_drug_sheet, by = c("drug" = "Drug")) %>%
  dplyr::mutate(plate = rep(1:n_plates, each = 96),
                row = str_extract(Metadata_Well, pattern = "[A-Z]"),
                col = rep(1:12, times = 8*n_plates),
                bleach = rep(1:3, each = 96*rep_plates, length.out = n_plates*96),
                concentration_um = case_when(col == 1 ~ 0,
                                             col == levels ~ as.double(High_Dose_uM),
                                             col != 1 | col != levels ~ as.double(High_Dose_uM / DilutionFactors^(levels - col))),
                concentration_um = round(concentration_um, digits = 2))

###############################################
###            Assign strains               ### 
###############################################
strains = dplyr::tibble(Metadata_Plate, 	Metadata_Well,	assay_type = assay,	food_conc_OD = food_od,	food_type = food_t, diluent_concentration_perc, drug_prep_method) %>%
  dplyr::mutate(plate = rep(1:n_plates, each = 96),
                row = str_extract(Metadata_Well, pattern = "[A-Z]"),
                col = rep(1:12, times = 8*n_plates),
                bleach = rep(1:3, each = 96*rep_plates, length.out = n_plates*96),
                strain = case_when(bleach == 1 & row == "A" ~ "RC301", bleach == 1 & row == "B" ~ "ECA396", bleach == 1 & row == "C" ~ "PD1074", bleach == 1 & row == "D" ~ "ECA248", bleach == 1 & row == "E" ~ "XZ1516", bleach == 1 & row == "F" ~ "CB4856", bleach == 1 & row == "G" ~ "ECA36", bleach == 1 & row == "H" ~ "MY16",
                                   bleach == 2 & row == "A" ~ "RC301", bleach == 2 & row == "B" ~ "MY16", bleach == 2 & row == "C" ~ "CB4856", bleach == 2 & row == "D" ~ "ECA396", bleach == 2 & row == "E" ~ "ECA36", bleach == 2 & row == "F" ~ "XZ1516", bleach == 2 & row == "G" ~ "PD1074", bleach == 2 & row == "H" ~ "ECA248",
                                   bleach == 3 & row == "A" ~ "MY16", bleach == 3 & row == "B" ~ "RC301", bleach == 3 & row == "C" ~ "PD1074", bleach == 3 & row == "D" ~ "CB4856", bleach == 3 & row == "E" ~ "ECA36", bleach == 3 & row == "F" ~ "XZ1516", bleach == 3 & row == "G" ~ "ECA396", bleach == 3 & row == "H" ~ "ECA248",),
                drug = rep(drugs, each = 96*rep_plates*bleaches),
                well_censor = NA_integer_,
                well_censor_reason = NA_character_,
                notes = NA_character_)

###############################################
###            Finalize design              ### 
############################################### 
design1 <- strains %>%
  left_join(drug_doses) %>%
  dplyr::mutate(Metadata_Plate = paste0("p", Metadata_Plate)) %>%
  dplyr::select(Metadata_Plate,
                Metadata_Well,
                assay_type,
                food_conc_OD,
                food_type,
                diluent,
                diluent_concentration_perc,
                drug,
                concentration_um,
                drug_stock_prep_date,
                drug_prep_method,
                strain,
                bleach,
                well_censor,
                well_censor_reason,
                notes)

###############################################
###            Manual design  edits         ### 
###############################################
design2 <- design1 %>%
  dplyr::mutate(notes = if_else(strain == "ECA36" & bleach == 2, true = "males possible, check images", false = notes))

###############################################
###            Export design                ### 
###############################################
rio::export(design2, file = "~/Dropbox/HTA_ImageXpress/Projects/20210701_toxin28A/design/20210701_toxin28A_design.csv")
