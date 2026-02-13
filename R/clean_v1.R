#' Takes the raw .csv file exported from qualtrics and clean it
#' @description The raw data need to be exported from qualtrics EMOCEMP V1 in .csv with the numeric values option checked
#' @param antibodies an autoantibodies .csv file to be merged with the visit data
#' @param exclude If TRUE, exclude patients excluded from the study in the final cleaned data frame
#' @param datafile a V1 .csv file exported from qualtrics with numeric data option
#' @return a data.frame with cleaned data which may be used to merge with other visits
#' @export
#' @import dplyr
#' @import stringr
#' @importFrom lubridate dmy
#' @importFrom utils read.csv
#' @importFrom zscorer addWGSR
#' @author R.C.S

clean_v1 <- function (datafile, antibodies = NULL, exclude = TRUE){
  emocemp_messy <- read.csv(datafile,
                            skip = 1,
                            encoding = "UTF-8") # Load Data
  emocemp <- emocemp_messy[-c(1),-c(1:17)] # Exclude irrelevant rows and collumns from table

  print("starting script")

  # Renaming collumns
  colnames(emocemp) <- c("data_visita","id_centro","id_paciente",
                        "registro_centro", "tcle", "inclusao",
                        "exclusao", "nascimento", "sex",
                        "etnic", "peso", "altura",
                        "imc", "data_onset", "clinica_onset",
                        "clinica_onset_outros","clinica_obs", "historia_familiar",
                        "historia_familiar_outros",
                        "infeccao_2meses", "infeccao_especificar", "data_vacina",
                        "nome_vacina", "tabagismo", "comorbidade",
                        "medicacao", "sf_piramidal", "sf_cerebelar",
                        "sf_sensitivo", "sf_tronco", "sf_visual",
                        "sf_vesical_intestinal", "sf_cerebral", "conver_visual",
                        "conver_vesical_intestinal", "deambulacao", "edss",
                        "data_brain_rm", "data_spinal_rm", "data_orbit_rm",
                        "brain_mr", "brain_mr_realce", "brain_outras",
                        "orbit_rm", "orbit_realce", "spinal_rm",
                        "spinal_realce", "spinal_outros", "rm_outros",
                        "fan", "fr", "vhs", "ssa", "ssb", "vitd",
                        "b12", "ebv", "cmv",
                        "soro_outros", "lcr_cel", "lcr_dif",
                        "lcr_prot", "lcr_boc", "lcr_igg",
                        "tto_fa","tto_fa_ivmp","tto_fa_co","tto_fa_igg",
                        "tto_fa_plex","dmd","dmd_ifb_im","dmd_ifb_sc22",
                        "dmd_ifb_sc44","dmd_ifb","dmd_glatir","dmd_teriflu",
                        "dmd_fuma","dmd_fingo","dmd_nat","dmd_alen",
                        "dmd_ritux","dmd_azt","dmd_outros","obsv"

                        )

  print("Dumming multiple choice variables")
  # Clinical dumming
  name_c <- c("neurite_b",
              "neurite_u",
              "mielite_parcial",
              "mielite_transversa",
              "adem",
              "romboencefalite",
              "outros_clinic")

  df_clinica <- split_mcv("clinica_onset",7,emocemp)
  df_clinica <- mutate_mcv(name_c, df_clinica)

  # Acute treatment dumming
  name_ta <- c("ivmp",
               "corticode_vo",
               "igiv",
               "plex",
               "nenhum_ta")

  df_tta <- split_mcv("tto_fa",5,emocemp)
  df_tta <- mutate_mcv(name_ta, df_tta)

  # Brain MRI dumming
  name_brain_mri <- c("uma_ou_mais_perpendicular",
                      "apenas_bem_delimitadas",
                      "cinco_ou_mais_t2",
                      "duas_ou_mais_peri",
                      "uma_ou_mais_tronco",
                      "black_holes",
                      "difusas_bilat",
                      "uma_ou_mais_justa",
                      "uma_ou_mais_peri",
                      "uma_ou_mais_infra",
                      "realce_anel_inc",
                      "outro_realce",
                      "nenhum_brm",
                      "naorealizada_brm",
                      "outras_brm")

  df_brm <- split_mcv("brain_mr",15,emocemp)
  df_brm <- mutate_mcv(name_brain_mri, df_brm)

  # Spinal MRI dumming
  name_spinal_mri <- c("mielite_trans_mr",
                       "mielite_centromedular",
                       "mielite_perif",
                       "mielite_menos3",
                       "letm",
                       "mielite_cerv",
                       "mielite_dorsal",
                       "mielite_lombo",
                       "nenhuma_spinalmr",
                       "nao_realizada_srm",
                       "realce_medular",
                       "bright_spot",
                       "outras_srm")

  df_srm <- split_mcv("spinal_rm",13,emocemp)
  df_srm <- mutate_mcv(name_spinal_mri, df_srm)

  # Orbital MRI dumming
  name_orbital_mri <- c("hipersinal_retrob",
                        "hipersinal_quiasma",
                        "hipersinal_media",
                        "hipersinal_extensa",
                        "realce_orbit",
                        "nenhuma_orbitrm",
                        "nao_realizada_orbitrm",
                        "outras_orbitrm")

  df_orm <- split_mcv("orbit_rm",8,emocemp)
  df_orm <- mutate_mcv(name_orbital_mri, df_orm)

  # Merge dfs
  dummed <- cbind(df_orm,df_srm,df_brm,df_tta,df_clinica)
  emocemp <- cbind(emocemp, dummed)
  print("Done")
  print("cleaning factors")
  # Cleaning factors
  # -------------------------------------------
  # Demographics ------------------------------
  # -------------------------------------------
  # Etnic
  emocemp$etnic <- as.character(emocemp$etnic)
  emocemp$etnic <- as.numeric(emocemp$etnic)
  emocemp$etnic <- as.factor(emocemp$etnic)

  # Infecção últimos 2 meses
  emocemp$infeccao_2meses <- as.character(emocemp$infeccao_2meses)
  emocemp$infeccao_2meses <- as.numeric(emocemp$infeccao_2meses)
  emocemp$infeccao_2meses <- as.factor(emocemp$infeccao_2meses)
  levels(emocemp$infeccao_2meses) <- c("sim", "nao")

  # Sex
  emocemp$sex <- as.character(emocemp$sex)
  emocemp$sex <- as.numeric(emocemp$sex)
  emocemp$sex <- as.factor(emocemp$sex)
  levels(emocemp$sex) <- c("m", "f")
  print("Done")


  print("Cleaning numerics")

  # Altura
  emocemp$altura <- as.character(emocemp$altura)
  altura_cm <- str_detect(emocemp$altura, "^\\d{2}")
  emocemp$altura <- str_replace(emocemp$altura, ",", ".")
  emocemp$altura <- as.numeric(emocemp$altura)
  emocemp$altura[altura_cm] <- emocemp$altura[altura_cm] / 100 # Transform in meters

  # Peso
  emocemp$peso <- as.character(emocemp$peso)
  emocemp$peso <- as.numeric(str_replace(emocemp$peso, ",",","))

  # IMC
  emocemp$imc <- as.character(emocemp$imc)
  emocemp$imc <- as.numeric(str_replace(emocemp$imc, ",","."))

  # Tabagismo
  emocemp$tabagismo <- as.character(emocemp$tabagismo)
  emocemp$tabagismo <- as.numeric(emocemp$tabagismo)
  emocemp$tabagismo <- factor(emocemp$tabagismo)
  levels(emocemp$tabagismo) <- c("ativo","passivo","ausente")

  # -----------------------------------------------------
  # Clinical analysis -----------------------------------
  # -----------------------------------------------------

  # EDSS ------------------------------------------------
  # Change to numeric and divided by 2 (it is doubled from the exportation)
  emocemp$edss <- as.character(emocemp$edss)
  emocemp$edss <- as.numeric(emocemp$edss)
  emocemp$edss <- emocemp$edss / 2
  #Create index to verify 0.5 EDSS
  emocemp <- emocemp %>%
    mutate(edss_i = if_else(edss < 1, TRUE, FALSE))
  # Change 0.5 for zeros
  emocemp$edss[emocemp$edss_i] <- 0
  # SF
  fs <- c( "sf_piramidal", "sf_cerebelar",
           "sf_sensitivo", "sf_tronco", "sf_visual",
           "sf_vesical_intestinal", "sf_cerebral")

  for (i in 1:length(fs)) {
    emocemp[,fs[i]] <- as.numeric(emocemp[,fs[i]])
    emocemp[,fs[i]] <- emocemp[,fs[i]] - 1
  }

  print("Done")
  # Dates ------------------------------------------------

  print("Cleaning dates")
  # Nascimento
  emocemp$nascimento <- dmy(emocemp$nascimento)
  # Visit
  emocemp$data_visita <- dmy(emocemp$data_visita)
  # Onset
  emocemp$data_onset <- dmy(emocemp$data_onset)
  # Visit age
  emocemp <- emocemp %>%
    mutate(idade_visita1 = as.integer(
      (emocemp$data_visita - emocemp$nascimento)/365 ) )
  # Onset age
  emocemp <- emocemp %>%
    mutate(idade_onset = as.integer(
      (emocemp$data_onset - emocemp$nascimento)/365 ) )
  # Disease duration at visit (years)
  emocemp <- emocemp %>%
    mutate(disease_duration_y = as.integer(
      (emocemp$idade_visita1 - emocemp$idade_onset)
    ))
  # Disease duration at visit (months)
  emocemp <- emocemp %>%
    mutate(disease_duration_m = (emocemp$data_visita - emocemp$data_onset) / 30 )

  # Disease duration at visit (days)
  emocemp <- emocemp %>%
    mutate(disease_duration_d = emocemp$data_visita - emocemp$data_onset)

  # Vaccine date
  emocemp$data_vacina <- dmy(emocemp$data_vacina)

  # Replace negative values with NAs
  wrong_dd <- emocemp$disease_duration_m < 0
  emocemp$disease_duration_m[wrong_dd] <- NA
  wrong_ageonset <- emocemp$idade_onset < 0
  emocemp$idade_onset[wrong_ageonset] <- NA
  # Create age factor > | < 10
  emocemp <- emocemp %>%
    mutate(idade_10 = ifelse(idade_onset <= 10, "m", "M") )
  print("Done")
  print("Cleaning strings")

  # LCR
  lcr_messy <- emocemp %>%
    select(lcr_cel, lcr_dif, lcr_prot, lcr_boc)

  lcr_messy <- lcr_messy %>%
    rename (
      cel = lcr_cel ,
      dif = lcr_dif ,
      prot = lcr_prot ,
      boc = lcr_boc
    )

  lcr_messy$boc <- as.character(lcr_messy$boc)
  lcr_messy$cel <- as.character(lcr_messy$cel)
  lcr_messy$prot <- as.character(lcr_messy$prot)
  lcr_messy$dif <- as.character(lcr_messy$dif)


  lcr_messy <- lcr_messy %>%
    mutate(across(c(boc, cel, prot),
                ~ na_if(str_trim(.), "")))


  correct_prot <- str_extract(lcr_messy$prot, "^\\d{1,3}[\\. ,]?\\d{0,3}$")
  correct_prot <- str_replace(correct_prot, ",", ".")
  correct_prot <- as.numeric(correct_prot)
  emocemp$lcr_prot <- correct_prot

  correct_cel <- str_extract(lcr_messy$cel, "^\\d{1,4}[\\. ,]?\\d{0,3}")
  correct_cel <- str_replace(correct_cel, "," , ".")
  correct_cel <- as.numeric(correct_cel)
  emocemp$lcr_cel <- correct_cel

  lcr_messy$boc <- tolower(lcr_messy$boc)


  # Create new IMC column with calculated from weight and height
  emocemp <- emocemp %>%
    mutate(imc_c = peso / altura^2)


  # Z-Score and nutritional status
  emocemp <- emocemp %>%
    mutate(altura_cm = altura * 100) %>%
    mutate(idade_dias = idade_visita1 * 365.25) %>%
    mutate(sex_coded = ifelse(sex == "m", 1, 2))

   emocemp <- zscorer::addWGSR(emocemp, sex="sex_coded",
                                   firstPart="peso",
                                   secondPart="altura_cm",
                                   thirdPart = "idade_dias",
                                   index="bfa")
   emocemp <- emocemp %>%
     mutate(obese =ifelse(bfaz >= 2,1,0)) %>%
     mutate(sobrepe = ifelse(bfaz >= 1,1,0))

   # Join vaccination and infection data
   emocemp <- emocemp %>%
     mutate(vacina_m = as.integer ((data_onset - data_vacina) / 12 ) )

   # NA Wrong dates
   wrong_dtvacina <- emocemp$vacina_m < 0
   emocemp$vacina_m[wrong_dtvacina] <- NA
   # Create inf_ou_vac
   emocemp <- emocemp %>%
     mutate(vacina_recente = ifelse(vacina_m <= 2, "1", "0") ) %>%
     mutate (inf_ou_vac_rec = ifelse(vacina_m==1 | infeccao_2meses == "sim", "1", "0" ))

   emocemp$inf_ou_vac_rec <- as.factor (emocemp$inf_ou_vac_rec)

  print("Done")
  print("Writing final data frame")
  col_2_exclude <- c("edss_i",
                     "clinica_onset",
                     "tto_fa",
                     "brain_mr",
                     "orbit_rm",
                     "spinal_rm"
                     )

  cleaned <- emocemp[, ! names(emocemp) %in% col_2_exclude, drop = F]

  # Merge Antibodies
  if(!is.null(antibodies)){
  print("Merging Autoantibodies")
    autoantibody <- read.csv2(antibodies)
    autoantibody <- rename(autoantibody,
                           id_paciente = id ,
                           mog = mog.positivo)
    cleaned <- merge(emocemp, autoantibody, by = "id_paciente", all.x = TRUE)}

  # Exclude patients
  if(exclude){
  cleaned <- subset (cleaned, !(id_paciente %in% excluded))
  }

  print("Script completed")
  return(cleaned)

}
