#'Takes a raw dataset for emocemp teleconsults and clean it
#'@description The raw data need to be exported from qualtrics EMOCEMP Tele in .csv with the numeric values option checked
#'@import dplyr
#'@import stringr
#'@importFrom lubridate dmy
#'@importFrom utils read.csv
#'@param datafile a EMOCEMP teleconsult .csv file exported from qualtrics with numeric data option
#'@return a data.frame with cleaned data which may be used to merge with other visits
#'@export
#'@author R.C.S

clean_tele2 <- function(datafile){
  
  emocemp_messy <- read.csv(
    datafile,
    skip = 1,
    encoding = "UTF-8",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  emocemptele <- emocemp_messy[-c(1), -c(1:17)]
  
  colnames(emocemptele) <- c(
    "consentimento","acompanhante","nome_acompanhante_1",
    "nome_acompanhante_2","visita_tele",
    "instrumento_utilizado","id_centro","id_paciente",
    "data_televisita","novo_surto_tele",
    "novas_sind_clinica","data_nob",
    "data_nou","data_mp",
    "data_mt","data_adem",
    "data_rombo","data_outras",
    "lcr_cel","lcr_dif",
    "lcr_prot","lcr_boc",
    "lcr_igg",
    "tto_fa","tto_fa_mpiv",
    "tto_fa_co","tto_fa_igiv",
    "tto_fa_plex","obsv_surtos",
    "tabagismo","peso_kg",
    "altura_m","dmd",
    "dmd_ifb1a","dmd_ifb1a22",
    "dmd_ifb1a44","dmd_ifb1b",
    "dmd_glat","dmd_terif",
    "dmd_dmf","dmd_fingo",
    "dmd_nataliz","dmd_alent","dmd_ritux",
    "dmd_aza","dmd_outros",
    "ea","ea_outros",
    "medica_continuo","fan",
    "fr","vhs","ssa",
    "ssb","vit_d","vit_b12",
    "labs_outros","obs","data_encefalo",
    "data_orbita","data_medula",
    "DICOM","brain_mr","brain_realce_mr",
    "brain_mr_outros", "orbit_mr",
    "orbit_realce_mr","orbit_mr_outros",
    "spinal_mr","spinal_mr_realce",
    "spinal_mr_outros","mr_outros",
    "deambul","sf_pira",
    "sf_cereb","sf_sens",
    "sf_tronco","sf_visual",
    "sf_vesicalintest","sf_cerebral",
    "cv_visual","cv_vesicalintest",
    "edss","criterio_diag_tele",
    "criterio_outros_tele", "obsv_geral"
  )
  
  message("Dumming multiple choice variables")
  
  clean <- function(var, database, n, names){
    df <- dummie_mcv(
      names = names,
      var = var,
      n = n,
      dataf = database
    )
    database <- database[, names(database) != var]
    database <- cbind(database, df)
    return(database)
  }
  
  # -------------------------------------------------
  # Clinical
  # -------------------------------------------------
  message("Clinical")
  
  name_clinical <- c(
    "neurite_b",
    "neurite_u",
    "mielite_p",
    "mielite_t",
    "adem",
    "romboencefalite",
    "outras",
    "nenhumas_das_sindromes"
  )
  
  emocemptele <- clean("novas_sind_clinica", emocemptele, 8, name_clinical)
  
  # -------------------------------------------------
  # Classificação de surto
  # -------------------------------------------------
  message("Classifying relapse")
  
  emocemptele$novo_surto_tele <- dplyr::case_when(
    emocemptele$novo_surto_tele %in% c("1", 1, "sim", "Sim", "SIM") ~ "sim",
    emocemptele$novo_surto_tele %in% c("2", 2, "nao", "não", "Nao", "Não", "NAO", "NÃO") ~ "nao",
    TRUE ~ NA_character_
  )
  
  # Garantir que as dummies clínicas estejam numéricas/lógicas
  for (nm in name_clinical) {
    emocemptele[[nm]] <- suppressWarnings(as.numeric(as.character(emocemptele[[nm]])))
    emocemptele[[nm]][is.na(emocemptele[[nm]])] <- 0
  }
  
  emocemptele$n_sindromes_surto <- rowSums(emocemptele[, c(
    "neurite_b",
    "neurite_u",
    "mielite_p",
    "mielite_t",
    "adem",
    "romboencefalite",
    "outras"
  )], na.rm = TRUE)
  
  emocemptele$class_surto_tele <- dplyr::case_when(
    emocemptele$novo_surto_tele == "nao" ~ "sem_surto",
    
    emocemptele$novo_surto_tele == "sim" & emocemptele$n_sindromes_surto == 0 ~ "surto_sem_classificacao",
    
    emocemptele$novo_surto_tele == "sim" & emocemptele$n_sindromes_surto > 1 ~ "multifocal_ou_multissindromico",
    
    emocemptele$novo_surto_tele == "sim" & emocemptele$neurite_b == 1 ~ "neurite_optica_bilateral",
    emocemptele$novo_surto_tele == "sim" & emocemptele$neurite_u == 1 ~ "neurite_optica_unilateral",
    emocemptele$novo_surto_tele == "sim" & emocemptele$mielite_p == 1 ~ "mielite_parcial",
    emocemptele$novo_surto_tele == "sim" & emocemptele$mielite_t == 1 ~ "mielite_transversa",
    emocemptele$novo_surto_tele == "sim" & emocemptele$adem == 1 ~ "adem",
    emocemptele$novo_surto_tele == "sim" & emocemptele$romboencefalite == 1 ~ "romboencefalite",
    emocemptele$novo_surto_tele == "sim" & emocemptele$outras == 1 ~ "outras",
    
    TRUE ~ NA_character_
  )
  
  emocemptele$inconsistencia_surto <- dplyr::case_when(
    emocemptele$novo_surto_tele == "nao" & emocemptele$n_sindromes_surto > 0 ~ TRUE,
    emocemptele$novo_surto_tele == "sim" & emocemptele$nenhumas_das_sindromes == 1 ~ TRUE,
    TRUE ~ FALSE
  )
  
  # -------------------------------------------------
  # MRI
  # -------------------------------------------------
  message("Image")
  
  name_brain_mr <- c(
    "perp_cc","apenas_bem_deli","5ou_mais_T2",
    "2ou_mais_peri","1_tronco","black_holes",
    "difusas_bilat","uma_ou_mais_justa","uma_ou_mais_peri",
    "1_infra","outros_realce",
    "nenhuma_acima_brain_mr",
    "brain_mr_nr",
    "outras_brain_mr","realce_anel_incompleto"
  )
  
  df_alt <- split_mcv("brain_mr", 15, emocemptele)
  df_alt[df_alt == 17] <- 15
  df_alt <- mutate_mcv(name_brain_mr, df_alt)
  emocemptele <- emocemptele[, names(emocemptele) != "brain_mr"]
  emocemptele <- cbind(emocemptele, df_alt)
  
  # Spinal
  name_spinal_mr <- c(
    "spinal_mr_mt","spinal_mr_mcm","spinal_mr_mp",
    "spinal_mr_<3","spinal_mr_>3","spinal_mr_cervical",
    "spinal_mr_dorsal","spinal_mr_lomb","spinal_mr_nenhuma_ac",
    "spinal_mr_nr","spinal_mr_realce",
    "spinal_mr_bright_spot",
    "spinal_mr_outras"
  )
  
  emocemptele <- clean("spinal_mr", emocemptele, 13, name_spinal_mr)
  
  # Orbit
  name_orbit_mr <- c(
    "hiper_t2_flair_retro","hiper_t2_flair_qui","hiper_t2_flair_med",
    "hiper_t2_flair_ext","orbit_realce","orbit_mr_nenhuma_ac",
    "orbit_mr_nr","orbit_mr_outras"
  )
  
  emocemptele <- clean("orbit_mr", emocemptele, 8, name_orbit_mr)
  
  # -------------------------------------------------
  # Treatment
  # -------------------------------------------------
  message("Treatment")
  
  name_dmd <- c(
    "ifn_b1a_IM","ifn_b1a_22","ifn_b1a_44","ifn_b1b",
    "ac_glat","teriflunomida","dmf","fingo","nataliz",
    "alentuz","rituximab","azat","outros_dmd","nenhum_dmd"
  )
  
  emocemptele <- clean("dmd", emocemptele, 14, name_dmd)
  
  # -------------------------------------------------
  # Factors
  # -------------------------------------------------
  message("Cleaning factors")
  
  emocemptele$criterio_diag_tele <- as.factor(emocemptele$criterio_diag_tele)
  
  vals_criterio <- sort(unique(stats::na.omit(as.character(emocemptele$criterio_diag_tele))))
  
  if (identical(vals_criterio, c("1", "2", "3", "4", "5"))) {
    emocemptele$criterio_diag_tele <- factor(
      emocemptele$criterio_diag_tele,
      levels = c("1", "2", "3", "4", "5"),
      labels = c("em", "nmosd", "isolated_event", "outras", "cis")
    )
  } else if (identical(vals_criterio, c("1", "3", "4", "5"))) {
    emocemptele$criterio_diag_tele <- factor(
      emocemptele$criterio_diag_tele,
      levels = c("1", "3", "4", "5"),
      labels = c("em", "isolated_event", "outras", "cis")
    )
  }
  
  emocemptele$novo_surto_tele <- factor(
    emocemptele$novo_surto_tele,
    levels = c("nao", "sim")
  )
  
  emocemptele$class_surto_tele <- as.factor(emocemptele$class_surto_tele)
  
  # -------------------------------------------------
  # Dates
  # -------------------------------------------------
  message("Cleaning dates")
  
  dates <- c(
    "data_televisita","data_nob",
    "data_nou","data_mp",
    "data_mt","data_adem",
    "data_rombo","data_outras",
    "data_encefalo",
    "data_orbita","data_medula"
  )
  
  for (i in dates) {
    emocemptele[[i]] <- suppressWarnings(lubridate::dmy(emocemptele[[i]]))
  }
  
  # -------------------------------------------------
  # Numerics
  # -------------------------------------------------
  message("Cleaning numerics")
  
  # helper numérico sem readr
  extrai_num <- function(x) {
    x <- stringr::str_replace_all(x, ",", ".")
    x <- stringr::str_extract(x, "\\d+(?:\\.\\d+)?")
    suppressWarnings(as.numeric(x))
  }
  
  # EDSS
  emocemptele$edss <- suppressWarnings(as.numeric(as.character(emocemptele$edss))) / 2
  emocemptele$edss[emocemptele$edss == 0.5] <- 0
  
  # Peso
  emocemptele$peso_kg <- extrai_num(emocemptele$peso_kg)
  
  # Altura
  altura_original <- emocemptele$altura_m
  emocemptele$altura_m <- extrai_num(emocemptele$altura_m)
  
  altura_i <- stringr::str_detect(altura_original, "^\\d{2}([,.]\\d+)?$")
  altura_i[is.na(altura_i)] <- FALSE
  emocemptele$altura_m[altura_i & !is.na(emocemptele$altura_m)] <-
    emocemptele$altura_m[altura_i & !is.na(emocemptele$altura_m)] / 100
  
  emocemptele$altura_m[!is.na(emocemptele$altura_m) & emocemptele$altura_m > 3] <-
    emocemptele$altura_m[!is.na(emocemptele$altura_m) & emocemptele$altura_m > 3] / 100
  
  # -------------------------------------------------
  # LCR
  # -------------------------------------------------
  message("Cleaning LCR")
  
  emocemptele$lcr_boc <- trimws(emocemptele$lcr_boc)
  emocemptele$lcr_boc[emocemptele$lcr_boc == ""] <- NA
  
  emocemptele$lcr_cel <- stringr::str_extract(
    emocemptele$lcr_cel,
    "\\d{1,3}(?:[\\.,]\\d{1,3})?"
  )
  emocemptele$lcr_cel <- suppressWarnings(
    as.numeric(stringr::str_replace_all(emocemptele$lcr_cel, ",", "."))
  )
  
  emocemptele$lcr_prot <- stringr::str_extract(
    emocemptele$lcr_prot,
    "\\d{1,3}(?:[\\.,]\\d{1,3})?"
  )
  emocemptele$lcr_prot <- suppressWarnings(
    as.numeric(stringr::str_replace_all(emocemptele$lcr_prot, ",", "."))
  )
  
  return(emocemptele)
}
