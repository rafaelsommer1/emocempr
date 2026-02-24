#'Takes the raw dataset from qualtrics and clean it
#'@description The raw data need to be exported from qualtrics EMOCEMP V2 in .csv with the numeric values option checked
#' @import dplyr
#' @import stringr
#' @importFrom lubridate dmy
#' @importFrom utils read.csv
#'@param datafile a V2 .csv file exported from qualtrics with numeric data option
#'@return a data.frame with cleaned data which may be used to merge with other visits
#'@export
#'@author R.C.S

clean_v2 <- function (datafile) {
        print("Starting script")

        # Load Data
        emocemp_messy <- read.csv(datafile,
                                  skip = 1,
                                  encoding = "UTF-8",
                                  stringsAsFactors = FALSE)


        # Exclude irrelevant rows and collumns from table
        emocempv2 <- emocemp_messy[-c(1),-c(1:17)]
        emocempv2 <- emocempv2[,-79]

        # Rename columns
        colnames(emocempv2) <- c("id_centro",
                                 "id_paciente",
                                 "data_visita",
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
                                 "data_cerebro_rm_v2",
                                 "data_orbitas_rm_v2",
                                 "data_medula_rm_v2",
                                 "rm_encefalo",
                                 "rm_encefalo_realce",
                                 "rm_encefalo_coment",
                                 "rm_orbita",
                                 "rm_orbita_realce_text",
                                 "rm_orbita_coment",
                                 "rm_medula",
                                 "rm_medula_realce_text",
                                 "rm_medula_coment",
                                 "rm_outros_texto",
                                 "medicacoes_cont",
                                 "criterio_diag_v2",
                                 "criterio_diag_outros_v2_text",
                                 "fan_titulo_padrao",
                                 "fr",
                                 "vhs",
                                 "anti-ssa",
                                 "anti-ssb",
                                 "vit_d",
                                 "vit_b12",
                                 "outros_labs",
                                 "dmd_v2",
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
                                 "novo_surto_v2",
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
                                 "tto_fa_v2",
                                 "mpiv_dose",
                                 "corticoide_vo_dose",
                                 "igiv_dose",
                                 "plex_dose")
        rm(emocemp_messy)

        #Start cleaning
        print("Dumming multiple choice variables")

        # Clinical syndromes
        name_clinical <- c("neurite_b",
                           "neurite_u",
                           "mielite_p",
                           "mielite_t",
                           "adem" ,
                           "romboencefalite",
                           "outras",
                           "nenhumas_das_sindromes")

        clinical_df <- split_mcv("novas_sindromes_clinicas",8,emocempv2)
        clinical_df <- mutate_mcv(name_clinical,clinical_df)
        # Include and exclude
        emocempv2 <- emocempv2[, names(emocempv2) != "novas_sindromes_clinicas"]
        emocempv2 <- cbind(emocempv2, clinical_df)


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
                           "outros_realce_brm",
                           "nenhuma_das_alt_brm",
                           "nao_realizada_brm",
                           "outras_brm",
                           "realce_anel_inc_brm")

        brm_df <- split_mcv("rm_encefalo",15,emocempv2)
        # Correct for qualtrics inconsistency
        brm_df[brm_df == 17] <- 15
        # Continue
        brm_df <- mutate_mcv(name_brm,brm_df)
        # Include and exclude
        emocempv2 <- emocempv2[, names(emocempv2) != "rm_encefalo"]
        emocempv2 <- cbind(emocempv2, brm_df)

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

        srm_df <- split_mcv("rm_medula",13,emocempv2)
        srm_df <- mutate_mcv(name_srm,srm_df)
        # Include and exclude
        emocempv2 <- emocempv2[, names(emocempv2) != "rm_medula"]
        emocempv2 <- cbind(emocempv2, srm_df)

        # Orbit imaging
        # Spinal imaging
        name_orm <- c("t2_hiper_retrobulbar_orm",
                      "t2_hiper_quiasma_orm",
                      "t2_hiper_medportion_orm",
                      "t2_hiper_mais_2tercos_orm",
                      "realce_orm" ,
                      "nenhuma_das_acima_orm",
                      "nao_realizada_orm",
                      "outras_orm")

        orm_df <- split_mcv("rm_orbita",8,emocempv2)
        orm_df <- mutate_mcv(name_orm,orm_df)
        # Include and exclude
        emocempv2 <- emocempv2[, names(emocempv2) != "rm_orbita"]
        emocempv2 <- cbind(emocempv2, orm_df)


        emocempv2$criterio_diag_v2 <- factor(emocempv2$criterio_diag_v2)
        if (length(levels(emocempv2$criterio_diag_v2)) == 4){
                levels(emocempv2$criterio_diag_v2) <- c("em","nmosd",
                                                        "isolated_event",
                                                        "outras")
                } else if (length(levels(emocempv2$criterio_diag_v2)) == 5) {
                        levels(emocempv2$criterio_diag_v2) <- c("em","nmosd",
                                                                "isolated_event",
                                                                "outras","cis")}
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

        dmd_df <- split_mcv("dmd_v2",14,emocempv2)
        dmd_df <- mutate_mcv(names_dmd,dmd_df)
        # Include and exclude
        emocempv2 <- emocempv2[, names(emocempv2) != "dmd_v2"]
        emocempv2 <- cbind(emocempv2, dmd_df)

        # Acute treatment dumming
        names_tta <- c("mpiv",
                       "corticoide_vo",
                       "igiv",
                       "plex",
                       "nenhum_tto_a")

        tta_df <- split_mcv("tto_fa_v2",5,emocempv2)
        tta_df <- mutate_mcv(names_tta,tta_df)
        # Include and exclude
        emocempv2 <- emocempv2[, names(emocempv2) != "tto_fa_v2"]
        emocempv2 <- cbind(emocempv2, tta_df)


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

        ae_df <- split_mcv("efeitos_adversos",12,emocempv2)
        ae_df <- mutate_mcv(names_ae,ae_df)
        # Include and exclude
        emocempv2 <- emocempv2[, names(emocempv2) != "efeitos_adversos"]
        emocempv2 <- cbind(emocempv2, ae_df)

        print("Ok!")

        # Cleaning dates
        print("Cleaning dates")
        dates <- c("data_visita","romboencefalite_data",
                   "data_cerebro_rm_v2","data_medula_rm_v2",
                   "data_orbitas_rm_v2","nob_data",
                   "nou_data","adem_data","outras_data","mielite_p_data",
                   "mielite_t_data")

         for (i in dates) {
                emocempv2[[i]] <- suppressWarnings(dmy(emocempv2[[i]]))
        }


        print("Ok!")

        # Cleaning factors
        print("Cleaning factors")
        # Tabagismo
        emocempv2$tabagismo <- as.factor(as.numeric(emocempv2$tabagismo))
        levels(emocempv2$tabagismo) <- c("ativo","passivo","ausente")


        # Cleaning numerics
        print("Cleaning numerics")

        emocempv2$edss <- as.numeric(emocempv2$edss) / 2
        emocempv2$edss[emocempv2$edss == 0.5] <- 0
        # Weigth
        emocempv2$peso <- as.numeric(str_replace(emocempv2$peso, ",","."))
        # Heigth
        altura_i <- str_detect(emocempv2$altura, "^\\d{2}")
        emocempv2$altura <- as.numeric(str_replace(emocempv2$altura,",","."))
        altura_i[is.na(altura_i)] <- FALSE
        emocempv2$altura[altura_i] <- emocempv2$altura[altura_i] / 100
        # IMC
        emocempv2$imc <- as.numeric(str_replace(emocempv2$imc,",","."))
        # Cleaning LCR strings

                # Input NAs
                  emocempv2$lcr_boc   <- stringr::str_squish(emocempv2$lcr_boc)
                  emocempv2$lcr_cel   <- stringr::str_squish(emocempv2$lcr_cel)
                  emocempv2$lcr_prot  <- stringr::str_squish(emocempv2$lcr_prot)
                  emocempv2$lcr_igg_i <- stringr::str_squish(emocempv2$lcr_igg_i)
          
                  emocempv2$lcr_boc[emocempv2$lcr_boc == ""] <- NA
                  emocempv2$lcr_cel[emocempv2$lcr_cel == ""] <- NA
                  emocempv2$lcr_prot[emocempv2$lcr_prot == ""] <- NA
                  emocempv2$lcr_igg_i[emocempv2$lcr_igg_i == ""] <- NA
  
                # Clean strings
                # Protein
                emocempv2$lcr_prot <- str_extract(emocempv2$lcr_prot,
                                                  "^\\d{1,3}[\\. ,]?\\d{0,3}$")

                # Cel
                emocempv2$lcr_cel <- str_extract(emocempv2$lcr_cel,
                                                  "^\\d{1,3}[\\. ,]?\\d{0,3}")
                emocempv2$lcr_cel <- as.numeric(str_replace(emocempv2$lcr_cel,
                                                 ",","."))
                # BOC

        # Binary variable
        emocempv2$novo_surto_v2 <- as.numeric(emocempv2$novo_surto_v2)
        emocempv2$novo_surto_v2 <- as.factor(emocempv2$novo_surto_v2)
        levels(emocempv2$novo_surto_v2) <- c("sim","nao")

        return(emocempv2)

        }



















