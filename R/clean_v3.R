#'Takes the raw dataset from qualtrics and clean it
#'@description The raw data need to be exported from qualtrics EMOCEMP V3 in .csv with the numeric values option checked
#' @import dplyr
#' @import stringr
#' @importFrom lubridate dmy
#' @importFrom utils read.csv
#'@param datafile a V3 .csv file exported from qualtrics with numeric data option
#'@return a data.frame with cleaned data which may be used to merge with other visits
#'@export
#'@author R.C.S

clean_v3 <- function(datafile){
        # Read input
        emocemp_messy <- read.csv(datafile,
                                  skip = 1,
                                  encoding = "UTF-8",
                                  stringsAsFactors = FALSE)


        # Exclude irrelevant rows and columns from table
        emocempv3 <- emocemp_messy[-c(1),-c(1:17)]
        emocempv3 <- emocempv3[,-79]
# -------------------------------------------------------------------------------- #
        # Rename
         colnames(emocempv3) <- c("data_visita",
                                 "id_centro",
                                 "id_paciente",
                                 "peso",
                                 "altura",
                                 "imc",
                                 "tabagismo",
                                 "sf_piramidal",
                                 "sf_cerebelar",
                                 "sf_sensitivo",
                                 "sf_tronco",
                                 "sf_visual",
                                 "sf_vesical_intestinal",
                                 "sf_cerebral",
                                 "conv_sf_visual",
                                 "conv_sf_vesical_intestinal",
                                 "deambula",
                                 "edss",
                                 "data_cerebro_rm_v3",
                                 "data_orbitas_rm_v3",
                                 "data_medula_rm_v3",
                                 "rm_encefalo",
                                 "rm_encefalo_realce",
                                 "rm_encefalo_coment",
                                 "rm_orbita",
                                 "rm_orbita_realce_text",
                                 "rm_medula",
                                 "rm_medula_realce_text",
                                 "rm_medula_coment",
                                 "rm_outros_texto",
                                 "medicacoes_cont",
                                 "criterio_diag_v3",
                                 "criterio_diag_outros_v3_text",
                                 "fan_titulo_padrao",
                                 "fr",
                                 "vhs",
                                 "anti-ssa",
                                 "anti-ssb",
                                 "vit_d",
                                 "vit_b12",
                                 "outros_labs",
                                 "dmd_v3",
                                 "ifn_b1a_im_dose",
                                 "ifn_b1a_sc22_dose",
                                 "ifn_b1a_sc44_dose",
                                 "ifn_b1b_dose",
                                 "acetato_glati_dose",
                                 "teriflunomida_dose",
                                 "fumarato_dose",
                                 "fingolimod_dose",
                                 "natalizumba_dose",
                                 "alentuzumab_dose",
                                 "rituximab_dose",
                                 "azatioprina_dose",
                                 "outros_dose",
                                 "efeitos_adversos",
                                 "efeitos_adversos_texto",
                                 "observacoes",
                                 "novo_surto_v3",
                                 "novas_sindromes_clinicas",
                                 "nob_data",
                                 "nou_data",
                                 "mielite_p_data",
                                 "mielite_t_data",
                                 "adem_data",
                                 "romboencefalite_data",
                                 "outras_data",
                                 "lcr_cel",
                                 "lcr_dif",
                                 "lcr_prot",
                                 "lcr_boc",
                                 "lcr_igg_i",
                                 "tto_fa_v3",
                                 "mpiv_dose",
                                 "corticoide_vo_dose",
                                 "igiv_dose",
                                 "plex_dose")
# -------------------------------------------------------------------------------- #
        # Multiple choice variables
        print("Dumming multiple choice variables")
        # Implement the cleaning function
        clean <- function(var,database,n,names){
                df <- dummie_mcv(names = names,
                                 var = var,
                                 n = n,
                                 dataf = database)
                database <- database[, names(database) != var]
                database <- cbind(database, df)
                return(database)
        }
# -------------------------------------------------------------------------------- #
        print("Clinical Data")
         # Clinical syndromes
        name_clinical <- c("neurite_b",
                           "neurite_u",
                           "mielite_p",
                           "mielite_t",
                           "adem" ,
                           "romboencefalite",
                           "outras",
                           "nenhumas_das_sindromes")

        emocempv3 <- clean("novas_sindromes_clinicas",emocempv3,8,name_clinical)
# -------------------------------------------------------------------------------- #
        print("Brain Data")
        # Brain imaging
        name_brm <- c("1_lesao_perpendicular_cc_brm",
                      "apenas_bem_delimitadas_brm",
                      "5_ou_mais_t2_brm",
                      "2_ou_mais_peri_brm",
                      "1_ou_mais_tronco_brm" ,
                      "black_holes_brm",
                      "lesoes_difusas_bilat_brm",
                      "1_ou_mais_justa_brm",
                      "1_ou_mais_peri_brm",
                      "1_ou_mais_infra_brm",
                      "realce_anel_inc_brm",
                      "outros_realces_brm",
                      "nenhuma_acima_brm",
                      "nao_realizada_brm",
                      "outras_brm")



        emocempv3 <- clean("rm_encefalo",emocempv3,15,name_brm)

        print("Done")
# -------------------------------------------------------------------------------- #
        print("Spinal Data")
        # Spinal imaging
        name_srm <- c("mielite_t_srm",
                      "mielite_centromed_srm",
                      "mielite_perif_srm",
                      "mielite_menos_3seg_srm",
                      "mielite_mais_3seg_srm" ,
                      "mielite_cerv_srm",
                      "mielite_dors_srm",
                      "mielite_lombo_srm",
                      "nenhuma_das_acima_srm",
                      "nao_realizada_srm",
                      "realce_srm",
                      "brigth_spot_srm",
                      "outras_srm")

        emocempv3 <- clean("rm_medula",emocempv3,13,name_srm)
        print("Done")
# -------------------------------------------------------------------------------- #
         # Orbit imaging
        print("Orbital Data")
        name_orm <- c("t2_hiper_retrobulbar_orm",
                      "t2_hiper_quiasma_orm",
                      "t2_hiper_medportion_orm",
                      "t2_hiper_mais_2tercos_orm",
                      "realce_orm" ,
                      "nenhuma_das_acima_orm",
                      "nao_realizada_orm",
                      "outras_orm")

        emocempv3 <- clean("rm_orbita",emocempv3,8,name_orm)
        print("Done")
# -------------------------------------------------------------------------------- #
        print("DMD Data")
        # DMD dumming
        names_dmd <- c("ifnb_1a_im",
                       "ifnb_1a_22sc",
                       "ifnb_1a_44sc",
                       "ifnb_1b",
                       "glatiramer",
                       "teriflunomida",
                       "fumarato_dimetil",
                       "fingolimod",
                       "natalizumab",
                       "alentuzumab",
                       "rituximab",
                       "azatioprina",
                       "outros_dmds",
                       "nenhum_dmd")

        emocempv3 <- clean("dmd_v3", emocempv3, 14, names_dmd)
        print("Done")
# -------------------------------------------------------------------------------- #
        print("Acute treatment Data")
        # Acute treatment dumming
        names_tta <- c("mpiv",
                       "corticoide_vo",
                       "igiv",
                       "plex",
                       "nenhum_tto_a")

        emocempv3 <- clean("tto_fa_v3",emocempv3,5,names_tta)
        print("Done")
# -------------------------------------------------------------------------------- #
        print("Adverse Events Data")
        # Adverse events dumming
        names_ae <- c("dor_local_ae",
                      "lipodistrofia_ae",
                      "flu_like_ae",
                      "infec__ae",
                      "lemp__ae",
                      "leuco_linfopenia_ae",
                      "desordem_autoimune_ae",
                      "infusao_react_ae",
                      "rubor_facial_ae",
                      "gastrointestinal_ae",
                      "elev_transaminases_ae",
                      "outros")

         emocempv3 <- clean("efeitos_adversos", emocempv3, 12, names_ae)
         print("Done")
# -------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------- #

        # Cleaning dates
        print("Cleaning dates")
        dates <- c("data_visita","romboencefalite_data",
                   "data_cerebro_rm_v3","data_medula_rm_v3",
                   "data_orbitas_rm_v3","nob_data",
                   "nou_data","adem_data","outras_data","mielite_p_data",
                   "mielite_t_data")

        for (i in dates) {
                emocempv3[[i]] <- suppressWarnings(dmy(emocempv3[[i]]))
        }


        print("Done")
# -------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------- #

        # Cleaning factors
        print("Cleaning factors")

        # Smoking status
        emocempv3$tabagismo <- as.factor(as.numeric(emocempv3$tabagismo))
        levels(emocempv3$tabagismo) <- c("ativo","passivo","ausente")

# -------------------------------------------------------------------------------- #
        # Diagnostic Criteria
        emocempv3$criterio_diag_v3 <- factor(as.numeric(emocempv3$criterio_diag_v3))
        levels(emocempv3$criterio_diag_v3) <- c("em","nmosd","isolated_event",
                                                "outras","cis")

# -------------------------------------------------------------------------------- #
        # Relapses
        emocempv3$novo_surto_v3 <- as.numeric(emocempv3$novo_surto_v3)
        emocempv3$novo_surto_v3 <- as.factor(emocempv3$novo_surto_v3)
        levels(emocempv3$novo_surto_v3) <- c("sim","nao")

        print("Done")
# -------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------- #

        # Cleaning numeric variables
        print("Cleaning numeric variables")

        # EDSS
        emocempv3$edss <- as.numeric(emocempv3$edss) / 2
        emocempv3$edss[emocempv3$edss == 0.5] <- 0
# -------------------------------------------------------------------------------- #

        # Weight
        emocempv3$peso <- as.numeric(str_replace(emocempv3$peso, ",","."))
# -------------------------------------------------------------------------------- #
        # Height
        altura_i <- str_detect(emocempv3$altura, "^\\d{2}")
        emocempv3$altura <- as.numeric(str_replace(emocempv3$altura,",","."))
        altura_i[is.na(altura_i)] <- FALSE
        emocempv3$altura[altura_i] <- emocempv3$altura[altura_i] / 100
# -------------------------------------------------------------------------------- #

        # IMC
        emocempv3$imc <- as.numeric(str_replace(emocempv3$imc,",","."))
# -------------------------------------------------------------------------------- #


        # Cleaning LCR strings

        # Input NAs
          emocempv3$lcr_boc   <- stringr::str_squish(emocempv2$lcr_boc)
          emocempv3$lcr_cel   <- stringr::str_squish(emocempv2$lcr_cel)
          emocempv3$lcr_prot  <- stringr::str_squish(emocempv2$lcr_prot)
          emocempv3$lcr_igg_i <- stringr::str_squish(emocempv2$lcr_igg_i)
          
          emocempv3$lcr_boc[emocempv2$lcr_boc == ""] <- NA
          emocempv3$lcr_cel[emocempv2$lcr_cel == ""] <- NA
          emocempv3$lcr_prot[emocempv2$lcr_prot == ""] <- NA
          emocempv3$lcr_igg_i[emocempv2$lcr_igg_i == ""] <- NA
  

        # Clean strings
        # Protein
        emocempv3$lcr_prot <- str_extract(emocempv3$lcr_prot,
                                          "^\\d{1,3}[\\. ,]?\\d{0,3}$")

        # Cel
        emocempv3$lcr_cel <- str_extract(emocempv3$lcr_cel,
                                         "^\\d{1,3}[\\. ,]?\\d{0,3}")
        emocempv3$lcr_cel <- as.numeric(str_replace(emocempv3$lcr_cel,
                                                    ",","."))

        return(emocempv3)

        }
