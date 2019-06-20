testthat::context("Test format_data output (approx. 40 seconds)")

################################################################################
## Generate testing data ##

# Load pmartR test data #
isobaric_object    <- pmartRdata::isobaric_object
lipid_object       <- pmartRdata::lipid_object
metab_object       <- pmartRdata::metab_object
pep_object         <- pmartRdata::pep_object
pro_object         <- pmartRdata::pro_object
techrep_pep_object <- pmartRdata::techrep_pep_object

# Set appropriate group designations #
isobaric_object    <- pmartR::group_designation(isobaric_object, "Group")
lipid_object       <- pmartR::group_designation(lipid_object, "Condition")
metab_object       <- pmartR::group_designation(metab_object, "Condition")
pep_object         <- pmartR::group_designation(pep_object, "Condition")
pro_object         <- pmartR::group_designation(pro_object, "Condition")
techrep_pep_object <- pmartR::group_designation(techrep_pep_object,
                                                c("FACTOR", "DILUTION"))

# Log transform untransformed data #
# isobaric_object2   <- pmartR::edata_transform(isobaric_object, "log2") ################## Error
lipid_object2      <- pmartR::edata_transform(lipid_object, "log2")
metab_object2      <- pmartR::edata_transform(metab_object, "log2")
pep_object2        <- pmartR::edata_transform(pep_object, "log2")

###
# Error in .as.isobaricpepData(e_data, f_data, e_meta, edata_cname, fdata_cname,  : 
# 'refpool_notation= Yes  is not in every experiment. See Details and Examples for more information.
###


### map test for different protein quants? anova/gtest/combined?

# Protein quantification from pepData with mappings #

qpro1 <- pmartR::protein_quant(isobaric_object, "rrollup")
# qpro2 <- pmartR::protein_quant(isobaric_object2, "rrollup")  ##################### Error from above
qpro3 <- pmartR::protein_quant(pep_object, "rrollup")
qpro4 <- pmartR::protein_quant(pep_object2, "rrollup")

## Note: qpor1 and qpro 3 are in "abundance" as is isobaric and pep objects, thus no stats or stats filters run

# Stats Filters #

# Not run on non-transformed data

# isobaric_object2 <- pmartR::applyFilt(pmartR::imdanova_filter(isobaric_object2),    ########
#                                       isobaric_object2,
#                                       min_nonmiss_anova = 2,
#                                       min_nonmiss_gtest = 3)
lipid_object2 <- pmartR::applyFilt(pmartR::imdanova_filter(lipid_object2),
                                   lipid_object2,
                                   min_nonmiss_anova = 2,
                                   min_nonmiss_gtest = 3)
metab_object2 <- pmartR::applyFilt(pmartR::imdanova_filter(metab_object2),
                                  metab_object2,
                                  min_nonmiss_anova = 2,
                                  min_nonmiss_gtest = 3)
pep_object2 <- pmartR::applyFilt(pmartR::imdanova_filter(pep_object2),
                                pep_object2,
                                min_nonmiss_anova = 2,
                                min_nonmiss_gtest = 3)
pro_object <- pmartR::applyFilt(pmartR::imdanova_filter(pro_object),
                                pro_object,
                                min_nonmiss_anova = 2,
                                min_nonmiss_gtest = 3)
techrep_pep_object <- pmartR::applyFilt(
  pmartR::imdanova_filter(techrep_pep_object),
  techrep_pep_object,
  min_nonmiss_anova = 2,
  min_nonmiss_gtest = 3)

# qpro2 <- pmartR::applyFilt(pmartR::imdanova_filter(qpro2),     ############
#                            qpro2,
#                            min_nonmiss_anova = 2,
#                            min_nonmiss_gtest = 3)
qpro4 <- pmartR::applyFilt(pmartR::imdanova_filter(qpro4),
                           qpro4,
                           min_nonmiss_anova = 2,
                           min_nonmiss_gtest = 3)


# Generate stats (requires log-transformed data) #
# isobaric_stats     <- pmartR::imd_anova(isobaric_object2,     ##################### Error from above
#                                         test_method = "combined")

# lipid_stats        <- pmartR::imd_anova(lipid_object2, 
#                                         test_method = "combined")   Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found
metab_stats        <- pmartR::imd_anova(metab_object2, 
                                        test_method = "combined")
# pep_stats          <- pmartR::imd_anova(pep_object2, 
#                                         test_method = "combined") ##  Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found
# pro_stats          <- pmartR::imd_anova(pro_object,
#                                         test_method = "combined")  ##  Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found
# techrep_pep_stats  <- pmartR::imd_anova(techrep_pep_object,
#                                         test_method = "combined")  ## Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found

# qpro2_stats        <- pmartR::imd_anova(qpro2, test_method = "combined")   ###################
# qpro4_stats        <- pmartR::imd_anova(qpro4, test_method = "combined") ## Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found



################################################################################

# Manipulated input #

badstat <- metab_stats
badstat$Full_results <-  badstat$Full_results[1:3,]  # Missing rows
badstat2 <- metab_stats
badstat2$Full_results <-  NULL  # Missing df
badstat3 <- metab_stats
badstat3$Flags <-  NULL  # Missing df
badstat3$P_values <- NULL  # Missing df
badstat4  <- metab_stats
badstat4$Flags$Metabolite <- NULL  # Missing e_data cname
badstat5 <- metab_stats
badstat5$Full_results[3] <- NULL  # Missing expected columns

baddat <- metab_object2
baddat$e_data <- NULL
baddat2 <- metab_object2
baddat2$f_data$SampleID <- NULL
baddat3 <- metab_object2
baddat3$e_data <- baddat3$e_data[1:3]

# Expected function input error throwing #

testthat::test_that("Correct format_data() error throwing for arguments omicsData, omicsStats", {
  testthat::expect_error(pmartR::format_data(omicsData = 23, omicsStats = 34))
  testthat::expect_error(pmartR::format_data(omicsStats = 34))
  testthat::expect_error(pmartR::format_data(omicsData = "23")) 
  testthat::expect_error(pmartR::format_data(omicsData = list()))
  testthat::expect_error(pmartR::format_data(omicsData = list(), 
                                             omicsStats = list()))
  testthat::expect_error(pmartR::format_data(omicsData = NULL, 
                                             omicsStats = NULL))
  testthat::expect_error(pmartR::format_data(omicsData = c(), 
                                             omicsStats = c()))
  testthat::expect_error(
    pmartR::format_data(omicsData = pro_object, omicsStats = pro_object))
  testthat::expect_error(
    pmartR::format_data(omicsData = metab_stats, omicsStats = metab_stats))
  testthat::expect_warning(
    pmartR::format_data(omicsData = NULL, omicsStats = pro_object))
  testthat::expect_error(
    pmartR::format_data(omicsData = metab_stats, omicsStats = pro_object))
  testthat::expect_error(
    pmartR::format_data(omicsData = pro_object, omicsStats = metab_stats))
  
  # Manipulated data
  testthat::expect_error(pmartR::format_data(badstat))
  testthat::expect_error(pmartR::format_data(badstat2))
  testthat::expect_error(pmartR::format_data(badstat3))
  testthat::expect_error(pmartR::format_data(badstat4))
  testthat::expect_error(pmartR::format_data(badstat5))
  testthat::expect_error(pmartR::format_data(badstat, metab_object2))
  testthat::expect_error(pmartR::format_data(badstat2, metab_object2))
  testthat::expect_error(pmartR::format_data(badstat3, metab_object2))
  testthat::expect_error(pmartR::format_data(badstat4, metab_object2))
  testthat::expect_error(pmartR::format_data(badstat5, metab_object2))
  testthat::expect_error(pmartR::format_data(baddat))
  testthat::expect_error(pmartR::format_data(baddat2))
  testthat::expect_error(pmartR::format_data(baddat3))
  testthat::expect_error(pmartR::format_data(baddat, metab_stats))
  testthat::expect_error(pmartR::format_data(baddat2, metab_stats))
  testthat::expect_error(pmartR::format_data(baddat3, metab_stats))
  
})


# Invalid until errors are resolved
testthat::test_that("Subfunction recursive_format correctly throws errors", {
  
  testthat::expect_error(format_data(omicsData = list(pep_object2, qpro4), omicsStats = list(pep_stats)))
  testthat::expect_error(format_data(omicsData = list(pep_object2, qpro4), omicsStats = list(pep_stats, pep_stats)))
  testthat::expect_error(format_data(omicsData = list(pep_object2, qpro4), omicsStats = list(pep_stats, pep_stats, pep_stats)))
  testthat::expect_error(format_data(omicsData = list(pep_object2, qpro4), omicsStats = list(metab_stats)))
  testthat::expect_error(format_data(omicsData = c(pep_object2, qpro4), omicsStats = list(pep_stats, pep_stats)))
  
  testthat::expect_error(format_data(omicsData = list(pep_object2, pep_object2)))
  testthat::expect_error(format_data(omicsData = list(pep_object2, pep_object2, pep_object2)))
  
  testthat::expect_error(format_data(omicsStats = list(pep_stats, pep_stats)))
  testthat::expect_error(format_data(omicsStats = list(pep_object2, pep_object2)))
})

################################################################################

## Generate function outputs from test data ##

get_edata_cname(isobaric_object)
get_emeta_cname(isobaric_object)
get_fdata_cname(isobaric_object)


# Data only #
format_isoobject   <- pmartR::format_data(isobaric_object)
# format_isoobject2  <- pmartR::format_data(isobaric_object2) ############
format_lipobject   <- pmartR::format_data(lipid_object)
format_lipobject2  <- pmartR::format_data(lipid_object2)
format_metobject   <- pmartR::format_data(metab_object)
format_metobject2  <- pmartR::format_data(metab_object2)
format_pepobject   <- pmartR::format_data(pep_object)
format_pepobject2  <- pmartR::format_data(pep_object2)
format_proobject   <- pmartR::format_data(pro_object)
format_tecobject   <- pmartR::format_data(techrep_pep_object)

# format_qpro1       <- pmartR::format_data(qpro1)              ################
# format_qpro2       <- pmartR::format_data(qpro2)               #################
# format_qpro3       <- pmartR::format_data(qpro3)                ################
format_qpro4       <- pmartR::format_data(qpro4)

# Stats only #
# format_isostats  <- pmartR::format_data(isobaric_stats) #############
# format_lipstats    <- pmartR::format_data(lipid_stats) #############
format_metstats    <- pmartR::format_data(metab_stats)
# format_pepstats    <- pmartR::format_data(pep_stats) #############
# format_prostats    <- pmartR::format_data(pro_stats) #############
# format_tecstats    <- pmartR::format_data(techrep_pep_stats) #############

# format_qpro2stats  <- pmartR::format_data(qpro2_stats)         ################
# format_qpro4stats  <- pmartR::format_data(qpro4_stats)     #############

# Both data and stats #
## Only log transformed data? ##  ADDD for matching log transforms

# format_isoboth    <- pmartR::format_data(isobaric_object, isobaric_stats) ############
# format_isoboth2   <- pmartR::format_data(isobaric_object2, isobaric_stats) ###########
# format_lipboth    <- pmartR::format_data(lipid_object, lipid_stats) ##### Error: different cnames (techrep)
# format_lipboth2   <- pmartR::format_data(lipid_object2, lipid_stats) ## Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found
# format_metboth    <- pmartR::format_data(metab_object, metab_stats) ##### Error: different cnames (techrep)
format_metboth2   <- pmartR::format_data(metab_object2, metab_stats)
# format_pepboth    <- pmartR::format_data(pep_object, pep_stats)        ##### Error: different cnames (techrep)
# format_pepboth2   <- pmartR::format_data(pep_object2, pep_stats) ## Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found
# format_proboth    <- pmartR::format_data(pro_object, pro_stats) ## Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found
# format_tecboth    <- pmartR::format_data(techrep_pep_object, techrep_pep_stats) ## Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found

# format_qpro2both <- pmartR::format_data(qpro2_stats, qpro2)   ###
# format_qpro4both <- pmartR::format_data(qpro4_stats, qpro4)   ###

# List based on pepData and quantified proData #
# Data
# forlist_qpro1dat  <- pmartR::format_data(omicsData = list(isobaric_object,    ################
#                                                           qpro1))          
# forlist_qpro2dat  <- pmartR::format_data(omicsData = list(isobaric_object2, ##############
#                                                           qpro2))
# forlist_qpro3dat  <- pmartR::format_data(omicsData = list(pep_object,       ################
#                                                           qpro3))
forlist_qpro4dat  <- pmartR::format_data(omicsData = list(pep_object2,
                                                          qpro4))

# Stats                                                                                                                    qpro4))
# forlist_qpro1stat <- pmartR::format_data(omicsStats = list(qpro1_stats,   ############
#                                                           isobaric_stats))
# forlist_qpro2stat <- pmartR::format_data(omicsStats = list(qpro2_stats,  ############
#                                                           isobaric_stats))
# forlist_qpro3stat  <- pmartR::format_data(omicsStats = list(qpro3_stats,   ################
#                                                           pep_stats))
# forlist_qpro4stat  <- pmartR::format_data(omicsStats = list(qpro4_stats, #######
#                                                           pep_stats))

#Both     qpro4))

# forlist_qpro1both <- pmartR::format_data(omicsStats = list(qpro1_stats,     #########
#                                                            isobaric_stats),
#                                          omicsData = list(qpro1, 
#                                                           isobaric_object))
# forlist_qpro2both <- pmartR::format_data(omicsStats = list(qpro2_stats,    ########
#                                                            isobaric_stats),
#                                          omicsData = list(qpro2, 
#                                                           isobaric_object2))

# forlist_qpro3both  <- pmartR::format_data(omicsStats = list(qpro3_stats,   ################ 
#                                                             pep_stats),
#                                           omicsData = list(qpro3,
#                                                            pep_object))
# forlist_qpro4both  <- pmartR::format_data(omicsStats = list(qpro4_stats,   #######
#                                                             pep_stats),
#                                           omicsData = list(qpro4,
#                                                            pep_object2))


################################################################################

Valid_format_dat <- list(
  format_isoobject, format_lipobject, 
  format_lipobject2, format_metobject, 
  format_metobject2, format_pepobject, 
  format_pepobject2, format_proobject, 
  format_tecobject, format_qpro4
)

Vfd_obs <- list(
  isobaric_object, lipid_object,
  lipid_object2, metab_object,
  metab_object2, pep_object,
  pep_object2, pro_object,
  techrep_pep_object, qpro4)

##

Valid_format_stat <- list(
  format_metstats 
  # format_pepstats, format_prostats, 
  # format_tecstats, format_qpro4stats
)

Vfs_obs <- list(
  metab_stats
)

##

Valid_format_both <- list(
  format_metboth2
  # format_pepboth2, format_proboth,
  # format_tecboth
)

Vfb_obs1 <- list(
  metab_object2
)

Vfb_obs2 <- list(
  metab_stats
)

#####

Valid_forlist_dat <- list(
  forlist_qpro4dat
)

Vfld_list <- list(
  list(pep_object2, qpro4)
)

##

Valid_forlist_stat <- list(
  # forlist_qpro4stat
)

Vfls_list <- list(
  # list(pro4_stat, pep_stat)
)

##

Valid_forlist_both <- list(
  # forlist_qpro4both
)

Vflb_list1 <- list(
  # list(pro4, pep_object2)
)

Vflb_list2 <- list(
  # list(pro4_stat, pep_stat)
)


################################################################################

## Test format_data outputs ##
#### Dimensions and correct columns in output ####

testthat::test_that("Correct dataframe population and dimensions", {
                      
  purrr::map2(Valid_format_dat, Vfd_obs, function(formDat, parOb){
    
    # Data frame population #
    testthat::expect_length(formDat, 3)
    testthat::expect_null(formDat$comp_stats)
    testthat::expect_null(formDat$summary_stats)
    testthat::expect_false(is.null(formDat$data_value))
    testthat::expect_gt(nrow(formDat$data_value), 0)
    
    # Data value correct columns and number of rows #
    form_cols <- colnames(formDat$data_values)
    testthat::expect_match(toString(form_cols), get_edata_cname(parOb))
    testthat::expect_match(toString(form_cols), get_fdata_cname(parOb))
    testthat::expect_true(all(colnames(parOb$f_data) %in% form_cols))
    form_cols <- form_cols[!(form_cols %in% colnames(parOb$f_data))]
    testthat::expect_match(toString(form_cols), "abundance")
    testthat::expect_match(toString(form_cols), "Group")
    
    if(!is.null(parOb$e_meta)){
      testthat::expect_match(toString(form_cols), get_emeta_cname(parOb))
      testthat::expect_true(all(colnames(parOb$e_meta) %in% form_cols))
      mult_row <- nrow(unique(
        parOb$e_meta[c(get_emeta_cname(parOb), get_edata_cname(parOb))]))
    } else {
      mult_row <- length(parOb$e_data[[get_edata_cname(parOb)]])
    }
    
    testthat::expect_identical(mult_row * 
                                 length(parOb$f_data[[get_fdata_cname(parOb)]]),
                               nrow(formDat$data_values))
    
  })
  
##
  
  purrr::map2(Valid_format_stat, Vfs_obs, function(formStat, parOb){
    
    # Data frame population #
    testthat::expect_length(formStat, 3)
    testthat::expect_null(formStat$data_value)
    testthat::expect_false(is.null(formStat$summary_stats))
    testthat::expect_false(is.null(formStat$comp_stats))
    testthat::expect_gt(nrow(formStat$summary_stats), 0)
    testthat::expect_gt(nrow(formStat$comp_stats), 0)
    
    # Comps and summary stats correct columns and number of rows #
    
    form_cols_summ <- colnames(formStat$summary_stats)
    #expect_cols_summ <- c(get_edata_cname(parOb), "Group", "Count", "Mean")
    expect_cols_summ <- c(attr(parOb, "cnames")$edata_cname, "Group", "Count", "Mean")
    # expect_rows_summ <- length(parOb$Full_results[[get_edata_cname(parOb)]]) *
    #   length(unique(attr(parOb, "group_DF")$Group))
    expect_rows_summ <- length(parOb$Full_results[[attr(parOb, "cnames")$edata_cname]]) *
      length(unique(attr(parOb, "group_DF")$Group))
    testthat::expect_identical(form_cols_summ, expect_cols_summ)
    testthat::expect_identical(nrow(formStat$summary_stats), expect_rows_summ)
    
    form_cols_comp <- colnames(formStat$comp_stats)
    test <- attr(parOb, "statistical_test")
    if (test == "gtest"){
      # expect_cols_comp <- c(get_edata_cname(parOb), "Comparison",
      #                       "P_value_G", "Fold_change", "Flag")
      expect_cols_comp <- c(attr(parOb, "cnames")$edata_cname, "Comparison",
                            "P_value_G", "Fold_change", "Flag")
    } else if (test == "anova"){
      # expect_cols_comp <- c(get_edata_cname(parOb), "Comparison",
      #                       "P_value_T", "Fold_change", "Flag")
      expect_cols_comp <- c(attr(parOb, "cnames")$edata_cname, "Comparison",
                            "P_value_T", "Fold_change", "Flag")
    } else {
      # expect_cols_comp <- c(get_edata_cname(parOb), "Comparison", 
      #                       "P_value_G", "P_value_T", 
      #                       "Fold_change", "Flag")
      expect_cols_comp <- c(attr(parOb, "cnames")$edata_cname, "Comparison", 
                            "P_value_G", "P_value_T", 
                            "Fold_change", "Flag")
    }
    # expect_rows_comp <- length(parOb$Full_Results[get_edata_cname(parOb)]) *
    #   length(unique(get_group_info(par_Ob)$Group))
    expect_rows_comp <- length(
      parOb$Full_results[[attr(parOb, "cnames")$edata_cname]]) *
      length(attr(parOb, "comparisons"))
    
    testthat::expect_identical(form_cols_comp, expect_cols_comp)
    testthat::expect_identical(nrow(formStat$comp_stats), expect_rows_comp)
    
  })

##  
  
  purrr::pmap(list(Valid_format_both, 
                   Vfb_obs1, 
                   Vfb_obs2), 
              function(formBoth, parDat, parStat){
                
    testthat::expect_length(formBoth, 3)
    testthat::expect_false(is.null(formBoth$data_value))
    testthat::expect_false(is.null(formBoth$summary_stats))
    testthat::expect_false(is.null(formBoth$comp_stats))
    testthat::expect_gt(nrow(formBoth$data_value), 0)
    testthat::expect_gt(nrow(formBoth$summary_stats), 0)
    testthat::expect_gt(nrow(formBoth$comp_stats), 0)
    
    # Data value correct columns and number of rows #
    form_cols <- colnames(formBoth$data_values)
    testthat::expect_match(toString(form_cols), get_edata_cname(parDat))
    testthat::expect_match(toString(form_cols), get_fdata_cname(parDat))
    testthat::expect_true(all(colnames(parDat$f_data) %in% form_cols))
    form_cols <- form_cols[!(form_cols %in% colnames(parDat$f_data))]
    testthat::expect_match(toString(form_cols), "abundance")
    testthat::expect_match(toString(form_cols), "Group")
    
    if(!is.null(parDat$e_meta)){
      testthat::expect_match(toString(form_cols), get_emeta_cname(parDat))
      testthat::expect_true(all(colnames(parDat$e_meta) %in% form_cols))
      mult_row <- nrow(unique(
        parDat$e_meta[c(get_emeta_cname(parDat), get_edata_cname(parDat))]))
    } else {
      mult_row <- length(parDat$e_data[[get_edata_cname(parDat)]])
    }
    
    testthat::expect_identical(mult_row * 
                                 length(parDat$f_data[[get_fdata_cname(parDat)]]),
                               nrow(formBoth$data_values))
    
    # Comps and summary stats correct columns and number of rows #
    
    form_cols_summ <- colnames(formBoth$summary_stats)
    #expect_cols_summ <- c(get_edata_cname(parStat), "Group", "Count", "Mean")
    expect_cols_summ <- c(attr(parStat, "cnames")$edata_cname, "Group", "Count", "Mean")
    # expect_rows_summ <- length(parStat$Full_results[[get_edata_cname(parStat)]]) *
    #   length(unique(attr(parStat, "group_DF")$Group))
    expect_rows_summ <- length(parStat$Full_results[[attr(parStat, "cnames")$edata_cname]]) *
      length(unique(attr(parStat, "group_DF")$Group))
    testthat::expect_identical(form_cols_summ, expect_cols_summ)
    testthat::expect_identical(nrow(formBoth$summary_stats), expect_rows_summ)
    
    form_cols_comp <- colnames(formBoth$comp_stats)
    test <- attr(parStat, "statistical_test")
    if (test == "gtest"){
      # expect_cols_comp <- c(get_edata_cname(parStat), "Comparison",
      #                       "P_value_G", "Fold_change", "Flag")
      expect_cols_comp <- c(attr(parStat, "cnames")$edata_cname, "Comparison",
                            "P_value_G", "Fold_change", "Flag")
    } else if (test == "anova"){
      # expect_cols_comp <- c(get_edata_cname(parStat), "Comparison",
      #                       "P_value_T", "Fold_change", "Flag")
      expect_cols_comp <- c(attr(parStat, "cnames")$edata_cname, "Comparison",
                            "P_value_T", "Fold_change", "Flag")
    } else {
      # expect_cols_comp <- c(get_edata_cname(parStat), "Comparison", 
      #                       "P_value_G", "P_value_T", 
      #                       "Fold_change", "Flag")
      expect_cols_comp <- c(attr(parStat, "cnames")$edata_cname, "Comparison", 
                            "P_value_G", "P_value_T", 
                            "Fold_change", "Flag")
    }
    # expect_rows_comp <- length(parStat$Full_Results[get_edata_cname(parStat)]) *
    #   length(unique(get_group_info(par_Ob)$Group))
    expect_rows_comp <- length(
      parStat$Full_results[[attr(parStat, "cnames")$edata_cname]]) *
      length(attr(parStat, "comparisons"))
    
    testthat::expect_identical(form_cols_comp, expect_cols_comp)
    testthat::expect_identical(nrow(formBoth$comp_stats), expect_rows_comp)
    
  })
  
##  
  
  purrr::map2(Valid_forlist_dat, Vfld_list, function(index_form, index_ob){
    testthat::expect_length(index_form, 2)
    purrr::map2(index_form, index_ob, function(formDat, parOb){
      # Data frame population #
      testthat::expect_length(formDat, 3)
      testthat::expect_null(formDat$comp_stats)
      testthat::expect_null(formDat$summary_stats)
      testthat::expect_false(is.null(formDat$data_value))
      testthat::expect_gt(nrow(formDat$data_value), 0)
      
      # Data value correct columns and number of rows #
      form_cols <- colnames(formDat$data_values)
      testthat::expect_match(toString(form_cols), get_edata_cname(parOb))
      testthat::expect_match(toString(form_cols), get_fdata_cname(parOb))
      testthat::expect_true(all(colnames(parOb$f_data) %in% form_cols))
      form_cols <- form_cols[!(form_cols %in% colnames(parOb$f_data))]
      testthat::expect_match(toString(form_cols), "abundance")
      testthat::expect_match(toString(form_cols), "Group")
      
      if(!is.null(parOb$e_meta)){
        testthat::expect_match(toString(form_cols), get_emeta_cname(parOb))
        testthat::expect_true(all(colnames(parOb$e_meta) %in% form_cols))
        mult_row <- nrow(unique(
          parOb$e_meta[c(get_emeta_cname(parOb), get_edata_cname(parOb))]))
      } else {
        mult_row <- length(parOb$e_data[[get_edata_cname(parOb)]])
      }
      
      testthat::expect_identical(mult_row * 
                                   length(parOb$f_data[[get_fdata_cname(parOb)]]),
                                 nrow(formDat$data_values))
      
    })
  })
  
  ## Update as other changes happen before errors are fixed
  # purrr::map(Valid_forlist_stat, function(index){
  #   testthat::expect_length(index, 2)
  #   purrr::map(index, function(formStat){
  #     testthat::expect_length(formStat, 3)
  #     testthat::expect_null(formStat$data_value)
  #     testthat::expect_false(is.null(formStat$summary_stats))
  #     testthat::expect_false(is.null(formStat$comp_stats))
  #   })
  # })
  # 
  # purrr::map(Valid_forlist_both, function(index){
  #   testthat::expect_length(index, 2)
  #   purrr::map(index, function(formBoth){
  #     testthat::expect_length(formBoth, 3)
  #     testthat::expect_false(is.null(formBoth$data_value))
  #     testthat::expect_false(is.null(formBoth$summary_stats))
  #     testthat::expect_false(is.null(formBoth$comp_stats))
  #   })
  # })
  
})

testthat::test_that("Format list slices are equal to non-lists", {
  testthat::expect_equal(forlist_qpro4dat[[1]], format_pepobject2)
  testthat::expect_equal(forlist_qpro4dat[[2]], format_qpro4)
})

#### Testing data specific expectations #####

## Random stats?


################################################################################
################################################################################

testthat::context("Test format_plot output")

## Generate nested plot structures ##
plot_isoobject <- format_plot(format_isoobject)
plot_lipobject <- format_plot(format_lipobject)
plot_lipobject2 <- format_plot(format_lipobject2)
plot_metobject <- format_plot(format_metobject)
plot_metobject2 <- format_plot(format_metobject2)
plot_pepobject <- format_plot(format_pepobject)
plot_pepobject2 <- format_plot(format_pepobject2)
plot_proobject <- format_plot(format_proobject)
plot_tecobject <- format_plot(format_tecobject)
plot_qpro4 <- format_plot(format_qpro4)

plot_metstats <- format_plot(format_metstats)

plot_metboth2 <- format_plot(format_metboth2)

# Plot object lists #
Valid_plot_dat <- list(
  plot_isoobject, plot_lipobject, 
  plot_lipobject2, plot_metobject, 
  plot_metobject2, plot_pepobject,
  plot_pepobject2, plot_proobject,  
  plot_tecobject, plot_qpro4)

Valid_plot_stat <- list(
  plot_metstats
)

##### Test error throwing #####

## "Bad" input data ##
badform <- format_pepobject
class(badform) <- "blah"  # incorrect class 
badform2 <- format_pepobject
badform2$data_values <- NULL  # Empty dataframes
badform3 <- format_pepobject
badform3$data_values <- NA    # Empty dataframes
badform4 <- format_pepobject
badform4$data_values <- "rr"  # not a dataframe
badform5 <- format_metstats
attr(badform5, "statistical_test") <- "blah"  # Bad test entry 
badform6 <- format_metstats
attr(badform6, "statistical_test") <- NULL  # Bad test entry 
badform7 <- format_metstats
badform7$comp_stats <- NULL  # missing stats dataframe
badform8 <- format_metstats
badform8$summary_stats <- NULL  # missing stats dataframe

## Error throwing ##
testthat::test_that("Correct error/warning throwing", {
  # Bad input #
  testthat::expect_error(format_plot(badform))
  testthat::expect_error(format_plot(badform2))
  testthat::expect_error(format_plot(badform3))
  testthat::expect_error(format_plot(badform4))
  testthat::expect_error(format_plot(badform5))
  testthat::expect_error(format_plot(badform6))
  testthat::expect_error(format_plot(badform7))
  testthat::expect_error(format_plot(badform8))
  
  # p val error #
  testthat::expect_error(format_plot(format_pepobject, p_val = NULL))
  testthat::expect_error(format_plot(format_pepobject, p_val = "ghjk"))
  testthat::expect_error(format_plot(format_pepobject, p_val = c(0.3, 1)))
  
  # Comps error #
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_y_range = 4, 
                                     comps_y_limits = "free"))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_y_limits = 4))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_y_range = 4, 
                                     comps_y_max = 5,
                                     comps_y_min = 1))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_y_limits = "free", 
                                     comps_y_max = 5,
                                     comps_y_min = 1))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_y_limits = c("free", "free")))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_y_range = c("free", "free")))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_y_range = c(1, 2)))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_y_range = -5))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_y_range = 0))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_y_max = "0"))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_y_max = c(1,3)))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_y_min = "0"))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_y_min = c(1,3)))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_plot_type = c(1,3)))
  testthat::expect_error(format_plot(format_pepobject, 
                                     comps_plot_type = c("box", 
                                                         "boxpoint", 
                                                         "raster")))
  testthat::expect_warning(format_plot(format_metstats, 
                                      comps_panel_y_axis = "P_value_G",
                                      comps_panel_x_axis = "P_value_G"), 
                           "identical")
  testthat::expect_error(format_plot(format_metstats, 
                                     comps_panel_x_axis = "Peptide",
                                     panel_variable = "Peptide"))
  testthat::expect_error(format_plot(format_metstats, 
                                     comps_panel_y_axis = "Peptide",
                                     panel_variable = "Peptide"))
  testthat::expect_error(format_plot(format_metstats, 
                                     comps_color_variable = "Peptide",
                                     panel_variable = "Peptide"))
  testthat::expect_error(format_plot(format_metstats, 
                                     comps_color_variable = "blue"))
  testthat::expect_error(format_plot(format_metstats, 
                                     panel_variable = "blue"))
  testthat::expect_error(format_plot(format_metstats, 
                                     comps_panel_y_axis = "blue"))
  testthat::expect_error(format_plot(format_metstats, 
                                     comps_panel_x_axis = "blue"))
  
  
  # Value error #
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_y_range = 4, 
                                     value_y_limits = "free"))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_y_limits = 4))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_y_range = 4, 
                                     value_y_max = 5,
                                     value_y_min = 1))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_y_limits = "free", 
                                     value_y_max = 5,
                                     value_y_min = 1))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_y_limits = c("free", "free")))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_y_range = c("free", "free")))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_y_range = c(1, 2)))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_y_range = -5))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_y_range = 0))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_y_max = "0"))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_y_max = c(1,3)))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_y_min = "0"))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_y_min = c(1,3)))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_plot_type = c(1,3)))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_plot_type = c("box", 
                                                         "boxpoint", 
                                                         "raster")))
  testthat::expect_warning(format_plot(format_pepobject, 
                                       value_panel_y_axis = "abundance",
                                       value_panel_x_axis = "abundance"), 
                           "identical")
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_panel_x_axis = "Peptide",
                                     panel_variable = "Peptide"))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_panel_y_axis = "Peptide",
                                     panel_variable = "Peptide"))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_color_variable = "Peptide",
                                     panel_variable = "Peptide"))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_color_variable = "blue"))
  testthat::expect_error(format_plot(format_pepobject, 
                                     panel_variable = "blue"))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_panel_y_axis = "blue"))
  testthat::expect_error(format_plot(format_pepobject, 
                                     value_panel_x_axis = "blue"))
})


## Test y-value messages

testthat::test_that("Correct messages", {
  testthat::expect_message(
    format_plot(format_pepobject2),
    "No specified value y-axis parameters. Axis y-limits will be scaled per plot, as per y_limits = 'free'.")
  testthat::expect_message(
    format_plot(format_pepobject2, value_y_max = 3),
    "No range or limits specified; Axis y-limits will be scaled per plot with a maximum of y_max.")
  testthat::expect_message(
    format_plot(format_pepobject2, value_y_min = 3), 
    "No range or limits specified; Axis y-limits will be scaled per plot with a minimum of y_min.")
  testthat::expect_message(
    format_plot(format_pepobject2, value_y_limits = "free"), 
    "Specified value y-limit: 'free'. Axis y-limits will be scaled per plot.")
  testthat::expect_message(
    format_plot(format_pepobject2, value_y_limits = "fixed"), 
    "Specified value y-limit: 'fixed'. Axis y-limits will fixed for all plots based on maximum and minimum y-values.")
  testthat::expect_message(
    format_plot(format_pepobject2, value_y_limits = "fixed",  value_y_max = 3), 
   "Specified value y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a maximum of y_max. Specified y_max")
  testthat::expect_message(
    format_plot(format_pepobject2, value_y_limits = "fixed",  value_y_min = 3), 
    "Specified value y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a minimum of y_min. Specified y_min")
  testthat::expect_message(
    format_plot(format_pepobject2, value_y_range = 10), 
    " units, split over the median.")
  testthat::expect_message(
    format_plot(format_pepobject2, value_y_range = 10,  value_y_min = 3), 
    " units from the y_min. Specified y_min: ")
  testthat::expect_message(
    format_plot(format_pepobject2, value_y_range = 10,  value_y_max = 3), 
    " units from the y_max. Specified y_max: ")

  
  testthat::expect_message(
    format_plot(format_metstats),
    "No specified comparison y-axis parameters. Axis y-limits will be scaled per plot, as per y_limits = 'free'.")
  testthat::expect_message(
    format_plot(format_metstats, comps_y_max = 3),
    "No range or limits specified; Axis y-limits will be scaled per plot with a maximum of y_max.")
  testthat::expect_message(
    format_plot(format_metstats, comps_y_min = 3), 
    "No range or limits specified; Axis y-limits will be scaled per plot with a minimum of y_min.")
  testthat::expect_message(
    format_plot(format_metstats, comps_y_limits = "free"), 
    "Specified comparison y-limit: 'free'. Axis y-limits will be scaled per plot.")
  testthat::expect_message(
    format_plot(format_metstats, comps_y_limits = "fixed"), 
    "Specified comparison y-limit: 'fixed'. Axis y-limits will fixed for all plots based on maximum and minimum y-values.")
  testthat::expect_message(
    format_plot(format_metstats, comps_y_limits = "fixed",  comps_y_max = 3), 
    "Specified comparison y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a maximum of y_max. Specified y_max")
  testthat::expect_message(
    format_plot(format_metstats, comps_y_limits = "fixed",  comps_y_min = 3), 
    "Specified comparison y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a minimum of y_min. Specified y_min")
  testthat::expect_message(
    format_plot(format_metstats, comps_y_range = 10), 
    " units, split over the median.")
  testthat::expect_message(
    format_plot(format_metstats, comps_y_range = 10,  comps_y_min = 3), 
    " units from the y_min. Specified y_min: ")
  testthat::expect_message(
    format_plot(format_metstats, comps_y_range = 10,  comps_y_max = 3), 
    " units from the y_max. Specified y_max: ")
})

##### Test dimensions #####
testthat::test_that("Correct nested table dimensions", {
  purrr::map2(c(Valid_plot_dat, Valid_plot_stat), 
              c(Valid_format_dat, Valid_format_stat), 
              function(plotDat, parForm){
                if(is.null(parForm$data_values)){
                  panelVar <- parForm$comp_stats[[colnames(plotDat)[1]]]
                  } else {
                    panelVar <- parForm$data_values[[colnames(plotDat)[1]]]
                    }
                expect_cols <- c(colnames(plotDat)[1], "data", "panel")
                # testthat::expect_identical(nrow(plotDat), length(panelVar))
                testthat::expect_equal(nrow(plotDat), 10)     ######### Subsetted version
                testthat::expect_identical(colnames(plotDat), expect_cols)
                testthat::expect_true(inherits(plotDat$data, "list"))
                testthat::expect_true(inherits(plotDat$panel, c(
                  "trelliscope_panels",
                  "list")))
                testthat::expect_true(all(lapply(plotDat$data, is.data.frame)))
                testthat::expect_true(all(
                  map(plotDat$panel, function(plot) inherits(
                    plot, c("gg", "plotly")))))
  })
})


##### Test random stats? #####

#Consistant numerics of over trnsformations

################################################################################
################################################################################

testthat::context("Test data_cogs output")

################################################################################
################################################################################

testthat::context("Test set_increment and set_ylimits output")

#### Testing data: y-values #### 
ylist1 <- data.frame(NA, "NA", 8, c(3,2), NA, NA)
ylist2 <- c(NA, "r"," r", "r", NA, NA)
ylist3 <- c(NA, "NA"," NA", "NA", NA, NA)
ylist4 <- c(0)
ylist5 <- data.frame(NA, NA, NA, NA, NA, NA)
ylist6 <- c(NULL)
ylist7 <- c(NA, NA, NA, NA, NA, NA)
ylist8 <- c(4)
ylist9 <- c(-4, -4, NA, NA, NA, NA)
ylist10 <- c(-4, -5, NA, NA, NA, NA)
ylist11 <- c(4, NA, NA, 6, NA, 5)
ylist12 <- data.frame(NA, NA, NA, 3, NA, NA)
ylist13 <- data.frame(NA, NA, NA, c(3,2), NA, NA)
ylist14 <- data.frame(NA, NA, -8, c(-9,-10), NA, NA)

#### Test set_increment #### 
testthat::test_that("Subfunction set_increment correctly processes", {
  testthat::expect_error(pmartR::set_increment(ylist1))
  testthat::expect_error(pmartR::set_increment(ylist2))
  testthat::expect_error(pmartR::set_increment(ylist3))
  testthat::expect_error(pmartR::set_increment(ylist4, include_zero = "z"))
  testthat::expect_error(pmartR::set_increment(ylist4, include_zero = 1))
  testthat::expect_identical(pmartR::set_increment(ylist4), 0)
  testthat::expect_identical(pmartR::set_increment(ylist5), 0)
  testthat::expect_identical(pmartR::set_increment(ylist6), 0)
  testthat::expect_identical(pmartR::set_increment(ylist7), 0)
  testthat::expect_identical(pmartR::set_increment(ylist8), 4/20)
  testthat::expect_identical(pmartR::set_increment(
    ylist8, include_zero = FALSE), 4/20)
  testthat::expect_identical(pmartR::set_increment(ylist9), 4/20)
  testthat::expect_identical(pmartR::set_increment(
    ylist9, include_zero = FALSE), 4/20)
  testthat::expect_identical(pmartR::set_increment(ylist10), 5/20)
  testthat::expect_identical(pmartR::set_increment(
    ylist10, include_zero = FALSE), 1/20)
  testthat::expect_identical(pmartR::set_increment(ylist11), 6/20)
  testthat::expect_identical(pmartR::set_increment(
    ylist11, include_zero = FALSE), 2/20)
  testthat::expect_identical(pmartR::set_increment(ylist12), 3/20)
  testthat::expect_identical(pmartR::set_increment(
    ylist12, include_zero = FALSE), 3/20)
  testthat::expect_identical(pmartR::set_increment(ylist13), 3/20)
  testthat::expect_identical(pmartR::set_increment(
    ylist13, include_zero = FALSE), 1/20)
  testthat::expect_identical(pmartR::set_increment(ylist14), 10/20)
  testthat::expect_identical(pmartR::set_increment(
    ylist14, include_zero = FALSE), 2/20)
})

#### Test set_ylimits #### 
testthat::test_that("Subfunction set_ylimits correctly processes", {
  
  temp_inc <- 1    # Testing increment value
  temp_ymax <- 20  # Testing y_max
  temp_ymin <- -16 # Testing y_min
  temp_yrange <- 6 # Testing y_range
  
  #### Errors ####
  testthat::expect_error(pmartR::set_ylimits(ylist1, temp_inc))
  testthat::expect_error(pmartR::set_ylimits(ylist2, temp_inc))
  testthat::expect_error(pmartR::set_ylimits(ylist3, temp_inc))
  testthat::expect_error(pmartR::set_ylimits(ylist4, "z"))
  testthat::expect_error(pmartR::set_ylimits(ylist4, NA))
  testthat::expect_error(pmartR::set_ylimits(ylist4, c(1,2)))
  testthat::expect_error(pmartR::set_ylimits(ylist4, temp_inc, include_zero = "z"))
  testthat::expect_error(pmartR::set_ylimits(ylist4, temp_inc, include_zero = 1))
  testthat::expect_error(pmartR::set_ylimits(ylist4, temp_inc, y_max = "z"))
  testthat::expect_error(pmartR::set_ylimits(ylist4, temp_inc, y_min = "z"))
  testthat::expect_error(pmartR::set_ylimits(ylist4, temp_inc, y_range = "z"))
  testthat::expect_error(pmartR::set_ylimits(ylist4, temp_inc, 
                                             y_range = temp_yrange,
                                             y_max = temp_ymax,
                                             y_min = temp_ymin))
  
  #### Value match w/o limits ####
  testthat::expect_identical(pmartR::set_ylimits(ylist4, temp_inc), c(-5, 5))
  testthat::expect_identical(pmartR::set_ylimits(ylist5, temp_inc), c(-5, 5))
  testthat::expect_identical(pmartR::set_ylimits(ylist6, temp_inc), c(-5, 5))
  testthat::expect_identical(pmartR::set_ylimits(ylist7, temp_inc), c(-5, 5))
  testthat::expect_identical(pmartR::set_ylimits(ylist8, temp_inc), c(-1, 9))
  testthat::expect_identical(pmartR::set_ylimits(
    ylist8, temp_inc, include_zero = FALSE), c(-1, 9))
  testthat::expect_identical(pmartR::set_ylimits(ylist9, temp_inc), c(-9, 1))
  testthat::expect_identical(pmartR::set_ylimits(
    ylist9, temp_inc, include_zero = FALSE), c(-9, 1))
  testthat::expect_identical(pmartR::set_ylimits(ylist10, temp_inc), c(-10, 1))
  testthat::expect_identical(pmartR::set_ylimits(
    ylist10, temp_inc, include_zero = FALSE), c(-10, 1))
  testthat::expect_identical(pmartR::set_ylimits(ylist11, temp_inc), c(-1, 11))
  testthat::expect_identical(pmartR::set_ylimits(
    ylist11, temp_inc, include_zero = FALSE), c(-1, 11))
  testthat::expect_identical(pmartR::set_ylimits(ylist12, temp_inc), c(-2, 8))
  testthat::expect_identical(pmartR::set_ylimits(
    ylist12, temp_inc, include_zero = FALSE), c(-2, 8))
  testthat::expect_identical(pmartR::set_ylimits(ylist13, temp_inc), c(-3, 8))
  testthat::expect_identical(pmartR::set_ylimits(
    ylist13, temp_inc, include_zero = FALSE), c(-3, 8))
  testthat::expect_identical(pmartR::set_ylimits(ylist14, temp_inc), c(-15, 5))
  testthat::expect_identical(pmartR::set_ylimits(
    ylist14, temp_inc, include_zero = FALSE), c(-15, -3))
  
  
  #### Value match w/ limits ####
  testthat::expect_identical(pmartR::set_ylimits(
    ylist4, temp_inc, y_max = temp_ymax), c(-5, 20))
  testthat::expect_identical(pmartR::set_ylimits(
    ylist4, temp_inc, y_min = temp_ymin), c(-16, 5))
  testthat::expect_identical(pmartR::set_ylimits(
    ylist4, temp_inc, y_min = temp_ymin), c(-16, 5))
  testthat::expect_identical(pmartR::set_ylimits(
    ylist4, temp_inc, y_min = temp_ymin, y_max = temp_ymax), c(-16, 20))
  testthat::expect_identical(pmartR::set_ylimits(
    ylist4, temp_inc, y_range = temp_yrange, y_max = temp_ymax), c(14, 20))
  testthat::expect_identical(pmartR::set_ylimits(
    ylist4, temp_inc, y_min = temp_ymin, y_range = temp_yrange), c(-16, -10))
  
})

################################################################################
################################################################################

testthat::context("Test main TrelliVis function output")

##  covr::file_coverage("./R/as.omicsPlotter.R", "./tests/testthat/test-as.R")  42.27%