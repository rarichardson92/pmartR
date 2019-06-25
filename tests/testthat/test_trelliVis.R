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
# isobaric_object   <- pmartR::edata_transform(isobaric_object, "log2") ################## Error
lipid_object      <- pmartR::edata_transform(lipid_object, "log2")
metab_object      <- pmartR::edata_transform(metab_object, "log2")
pep_object        <- pmartR::edata_transform(pep_object, "log2")

###
# Error in .as.isobaricpepData(e_data, f_data, e_meta, edata_cname, fdata_cname,  : 
# 'refpool_notation= Yes  is not in every experiment. See Details and Examples for more information.
###

# Protein quantification from pepData with mappings #

# qpro1 <- pmartR::protein_quant(isobaric_object, "rrollup")  ##################### Error from above
qpro2 <- pmartR::protein_quant(pep_object, "rrollup")

# Stats Filters #

# Not run on non-transformed data

# isobaric_object <- pmartR::applyFilt(pmartR::imdanova_filter(isobaric_object),    ########
#                                       isobaric_object,
#                                       min_nonmiss_anova = 2,
#                                       min_nonmiss_gtest = 3)
lipid_object <- pmartR::applyFilt(pmartR::imdanova_filter(lipid_object),
                                   lipid_object,
                                   min_nonmiss_anova = 2,
                                   min_nonmiss_gtest = 3)
metab_object <- pmartR::applyFilt(pmartR::imdanova_filter(metab_object),
                                  metab_object,
                                  min_nonmiss_anova = 2,
                                  min_nonmiss_gtest = 3)
pep_object <- pmartR::applyFilt(pmartR::imdanova_filter(pep_object),
                                pep_object,
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

# qpro1 <- pmartR::applyFilt(pmartR::imdanova_filter(qpro1),     ############
#                            qpro1,
#                            min_nonmiss_anova = 2,
#                            min_nonmiss_gtest = 3)
qpro2 <- pmartR::applyFilt(pmartR::imdanova_filter(qpro2),
                           qpro2,
                           min_nonmiss_anova = 2,
                           min_nonmiss_gtest = 3)


# Generate stats (requires log-transformed data) #
# isobaric_stats     <- pmartR::imd_anova(isobaric_object,     ##################### Error from above
#                                         test_method = "combined")

# lipid_stats        <- pmartR::imd_anova(lipid_object, 
#                                         test_method = "combined")   Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found
metab_stats        <- pmartR::imd_anova(metab_object, 
                                        test_method = "combined")
# pep_stats          <- pmartR::imd_anova(pep_object, 
#                                         test_method = "combined") ##  Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found
# pro_stats          <- pmartR::imd_anova(pro_object,
#                                         test_method = "combined")  ##  Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found
# techrep_pep_stats  <- pmartR::imd_anova(techrep_pep_object,
#                                         test_method = "combined")  ## Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found

# qpro1_stats        <- pmartR::imd_anova(qpro1, test_method = "combined")   ###################
# qpro2_stats        <- pmartR::imd_anova(qpro2, test_method = "combined") ## Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found



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

baddat <- metab_object
baddat$e_data <- NULL
baddat2 <- metab_object
baddat2$f_data$SampleID <- NULL
baddat3 <- metab_object
baddat3$e_data <- baddat3$e_data[1:3]

# Expected function input error throwing - format_data to validate_input #

testthat::test_that("Correct format_data() error throwing for arguments omicsData, omicsStats", {
  testthat::expect_error(
    pmartR::format_data(omicsData = 23, omicsStats = 34), 
    "must be of class")
  testthat::expect_error(
    pmartR::format_data(omicsStats = 34), 
    "must be of the class")
  testthat::expect_error(
    pmartR::format_data(omicsData = "23"),  
    "must be of class") 
  testthat::expect_error(
    pmartR::format_data(omicsData = list()), 
    "Empty list")
  testthat::expect_error(
    pmartR::format_data(omicsData = list(), omicsStats = list()), 
    "Empty list")
  testthat::expect_error(
    pmartR::format_data(omicsData = NULL, omicsStats = NULL), 
    "requires at least")
  testthat::expect_error(
    pmartR::format_data(omicsData = c(), omicsStats = c()), 
    "requires at least")
  testthat::expect_error(
    pmartR::format_data(omicsData = pro_object, omicsStats = pro_object), 
    "class 'statRes'")
  testthat::expect_error(
    pmartR::format_data(omicsData = metab_stats, omicsStats = metab_stats), 
    "omicsData must be of class")
  testthat::expect_warning(
    pmartR::format_data(omicsData = NULL, omicsStats = pro_object), 
    "Input reordered")
  testthat::expect_error(
    pmartR::format_data(omicsData = metab_stats, omicsStats = pro_object), 
    "Non-matching cname attributes")
  testthat::expect_error(
    pmartR::format_data(omicsData = pro_object, omicsStats = metab_stats), 
    "Non-matching cname attributes")
  
  # Manipulated data
  testthat::expect_error(
    pmartR::format_data(badstat), 
    "Mismatched rows")
  testthat::expect_error(
    pmartR::format_data(badstat2), 
    "should contain only and all of the following dataframes:")
  testthat::expect_error(
    pmartR::format_data(badstat3),  
    "should contain only and all of the following dataframes:")
  testthat::expect_error(
    pmartR::format_data(badstat4), 
    "column must be present in all omicsStats dataframes")
  testthat::expect_error(
    pmartR::format_data(badstat5), 
    "Number of columns in omicsStats dataframes is different than expected")
  testthat::expect_error(
    pmartR::format_data(badstat, metab_object), 
    "Mismatched rows")
  testthat::expect_error(
    pmartR::format_data(badstat2, metab_object), 
    "Requires compatible identifiers")
  testthat::expect_error(
    pmartR::format_data(badstat3, metab_object), 
    "should contain only and all of the following dataframes:")
  testthat::expect_error(
    pmartR::format_data(badstat4, metab_object), 
    "column must be present in all omicsStats dataframes")
  testthat::expect_error(
    pmartR::format_data(badstat5, metab_object), 
    "Number of columns in omicsStats dataframes is different than expected")
  testthat::expect_error(
    pmartR::format_data(baddat), 
    "Omicsdata requires both e_data and f_data.")
  testthat::expect_error(
    pmartR::format_data(baddat2), 
    "column must be present in omicsData f_data.")
  testthat::expect_error(
    pmartR::format_data(baddat3), 
    "column in f_data does not match column names in e_data" )
  testthat::expect_error(
    pmartR::format_data(baddat, metab_stats), 
    "Biomolecules in omicsStats do not match biomolecules in omicsData.")
  testthat::expect_error(
    pmartR::format_data(baddat2, metab_stats), 
    "SampleID column does not match between omicsData and omicsStats")
  testthat::expect_error(
    pmartR::format_data(baddat3, metab_stats), 
    "SampleID column in f_data does not match column names in e_data")
})

rm(baddat, baddat2, baddat3, badstat, badstat2, badstat3, badstat4, badstat5)

# Invalid until errors are resolved; 
# format_data -> recursive_format -> format_data (checks list input) -> format_data -> format_data (checks individual omics objects)
testthat::test_that("Subfunction recursive_format correctly throws errors", {
  
  # testthat::expect_error(
  #   format_data(omicsData = list(pep_object, qpro2), omicsStats = list(pep_stats)), 
  #   "Biomolecules in omicsStats do not match biomolecules in omicsData.")
  # testthat::expect_error(
  #   format_data(omicsData = list(pep_object, qpro2), omicsStats = list(pep_stats, pep_stats)), 
  #   "Lists in omicsData and omicsStats have mismatched cname attributes.")
  # testthat::expect_error(
  #   format_data(omicsData = list(pep_object, qpro2), omicsStats = list(pep_stats, pep_stats, pep_stats)),
  #   "List length does not match;")
  testthat::expect_error(
    format_data(omicsData = list(pep_object, qpro2), omicsStats = list(metab_stats)),
    "Non-matching cname attributes in omicsStats and omicsData.")
  # testthat::expect_error(
  #   format_data(omicsData = c(pep_object, qpro2), omicsStats = list(pep_stats, pep_stats)),
  #   "List length does not match;")
  # 
  testthat::expect_error(
    format_data(omicsData = list(pep_object, pep_object)), 
    "Only one pepData")
  testthat::expect_error(
    format_data(omicsData = list(pep_object, pep_object, pep_object)), 
    "List length != 2")
  
  # testthat::expect_error(
  #   format_data(omicsStats = list(pep_stats, pep_stats)), 
  #   "Only one stats object derived from pepData")
  testthat::expect_error(
    format_data(omicsStats = list(pep_object, pep_object)),
    "Only one stats object derived from pepData")
})

################################################################################

## Generate function outputs from test data ##

# Data only #
format_isoobject   <- pmartR::format_data(isobaric_object)
format_lipobject   <- pmartR::format_data(lipid_object)
format_metobject   <- pmartR::format_data(metab_object)
format_pepobject   <- pmartR::format_data(pep_object)
format_proobject   <- pmartR::format_data(pro_object)
format_tecobject   <- pmartR::format_data(techrep_pep_object)

# format_qpro1       <- pmartR::format_data(qpro1)              ################
format_qpro2       <- pmartR::format_data(qpro2)

# Stats only #
# format_isostats  <- pmartR::format_data(isobaric_stats) #############
# format_lipstats    <- pmartR::format_data(lipid_stats) #############
format_metstats    <- pmartR::format_data(metab_stats)
# format_pepstats    <- pmartR::format_data(pep_stats) #############
# format_prostats    <- pmartR::format_data(pro_stats) #############
# format_tecstats    <- pmartR::format_data(techrep_pep_stats) #############

# format_qpro1stats  <- pmartR::format_data(qpro1_stats)         ################
# format_qpro1stats  <- pmartR::format_data(qpro2_stats)     #############

# Both data and stats #
## Only log transformed data? ##  ADDD for matching log transforms

# format_isoboth    <- pmartR::format_data(isobaric_object, isobaric_stats) ############
# format_lipboth    <- pmartR::format_data(lipid_object, lipid_stats) ##### Error: different cnames (techrep)
format_metboth   <- pmartR::format_data(metab_object, metab_stats)
# format_pepboth    <- pmartR::format_data(pep_object, pep_stats)        ##### Error: different cnames (techrep)
# format_proboth    <- pmartR::format_data(pro_object, pro_stats) ## Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found
# format_tecboth    <- pmartR::format_data(techrep_pep_object, techrep_pep_stats) ## Error in Peptide %in% as.character(to_fix) : object 'Peptide' not found

# format_qpro1both <- pmartR::format_data(qpro1_stats, qpro1)   ###
# format_qpro2both <- pmartR::format_data(qpro2_stats, qpro2)   ###

# List based on pepData and quantified proData #
# Data
# forlist_qpro1dat  <- pmartR::format_data(omicsData = list(isobaric_object,    ################
#                                                           qpro1))          
forlist_qpro2dat  <- pmartR::format_data(omicsData = list(pep_object, qpro2))

# Stats                                                                       
# forlist_qpro1stat <- pmartR::format_data(omicsStats = list(qpro1_stats,   ############
#                                                           isobaric_stats))
# forlist_qpro2stat <- pmartR::format_data(omicsStats = list(qpro2_stats,  ############
#                                                           pep_stats))


# Both #

# forlist_qpro1both <- pmartR::format_data(omicsStats = list(qpro1_stats,     #########
#                                                            isobaric_stats),
#                                          omicsData = list(qpro1, 
#                                                           isobaric_object))
# forlist_qpro2both <- pmartR::format_data(omicsStats = list(qpro2_stats,    ########
#                                                            pep_stats),
#                                          omicsData = list(qpro2, 
#                                                           pep_object))


################################################################################

Valid_format_dat <- list(
  format_isoobject, format_lipobject, 
  format_metobject, format_pepobject, 
  format_proobject, format_tecobject, 
  format_qpro2
)

Vfd_obs <- list(
  isobaric_object, lipid_object, 
  metab_object, pep_object,
  pro_object, techrep_pep_object, 
  qpro2)

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
  format_metboth
)

Vfb_obs1 <- list(
  metab_object
)

Vfb_obs2 <- list(
  metab_stats
)

#####

Valid_forlist_dat <- list(
  forlist_qpro2dat
)

Vfld_list <- list(
  list(pep_object, qpro2)
)

##

Valid_forlist_stat <- list(
  # forlist_qpro2stat
)

Vfls_list <- list(
  # list(qpro2_stat, pep_stat)
)

##

Valid_forlist_both <- list(
  # forlist_qpro2both
)

Vflb_list1 <- list(
  # list(qpro2, pep_object)
)

Vflb_list2 <- list(
  # list(qpro2_stat, pep_stat)
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
  testthat::expect_equal(forlist_qpro2dat[[1]], format_pepobject)
  testthat::expect_equal(forlist_qpro2dat[[2]], format_qpro2)
})

#### Testing data specific expectations #####

## Random stats?


################################################################################
################################################################################

testthat::context("Test format_plot output")

## Generate nested plot structures ##
plot_isoobject <- pmartR::format_plot(format_isoobject)
plot_lipobject <- pmartR::format_plot(format_lipobject)
plot_metobject <- pmartR::format_plot(format_metobject)
plot_pepobject <- pmartR::format_plot(format_pepobject)
plot_proobject <- pmartR::format_plot(format_proobject)
plot_tecobject <- pmartR::format_plot(format_tecobject)
# plot_qpro1   <- pmartR::format_plot(format_qpro1)
plot_qpro2     <- pmartR::format_plot(format_qpro2)

plot_metstats  <- pmartR::format_plot(format_metstats)$panel[[1]]

plot_metboth   <- pmartR::format_plot(format_metboth)

# Plot object lists #
Valid_plot_dat <- list(
  plot_isoobject, plot_lipobject, 
  plot_metobject, plot_pepobject, 
  plot_proobject, plot_tecobject, 
  plot_qpro2)

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
  testthat::expect_error(
    pmartR::format_plot(badform),
    "trellData must be of the class 'trellData'")
  testthat::expect_error(
    pmartR::format_plot(badform2),
    "No data values or comparison statistics in trellData to plot")
  testthat::expect_error(
    pmartR::format_plot(badform3))
  testthat::expect_error(pmartR::format_plot(badform4))
  testthat::expect_error(pmartR::format_plot(badform5))
  testthat::expect_error(pmartR::format_plot(badform6))
  testthat::expect_error(pmartR::format_plot(badform7))
  testthat::expect_error(pmartR::format_plot(badform8))
  
  rm(badform, badform2, badform3, badform4, badform5, badform6, badform7, badform8)
  
  # p val error #
  testthat::expect_error(pmartR::format_plot(format_pepobject, p_val = NULL))
  testthat::expect_error(pmartR::format_plot(format_pepobject, p_val = "ghjk"))
  testthat::expect_error(pmartR::format_plot(format_pepobject, p_val = c(0.3, 1)))
  
  # Comps error #
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_y_range = 4, 
                                     comps_y_limits = "free"))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_y_limits = 4))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_y_range = 4, 
                                     comps_y_max = 5,
                                     comps_y_min = 1))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_y_limits = "free", 
                                     comps_y_max = 5,
                                     comps_y_min = 1))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_y_limits = c("free", "free")))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_y_range = c("free", "free")))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_y_range = c(1, 2)))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_y_range = -5))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_y_range = 0))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_y_max = "0"))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_y_max = c(1,3)))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_y_min = "0"))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_y_min = c(1,3)))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_plot_type = c(1,3)))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     comps_plot_type = c("box", 
                                                         "boxpoint", 
                                                         "raster")))
  testthat::expect_warning(
    pmartR:::validate_format_plot_input(
      format_metstats, 
      comps_panel_y_axis = "P_value_G",
      comps_panel_x_axis = "P_value_G",
      comps_color_variable = "Comparison",
      panel_variable = "Metabolite"),
    "identical")
  testthat::expect_error(pmartR::format_plot(format_metstats, 
                                     comps_panel_x_axis = "Peptide",
                                     panel_variable = "Peptide"))
  testthat::expect_error(pmartR::format_plot(format_metstats, 
                                     comps_panel_y_axis = "Peptide",
                                     panel_variable = "Peptide"))
  testthat::expect_error(pmartR::format_plot(format_metstats, 
                                     comps_color_variable = "Peptide",
                                     panel_variable = "Peptide"))
  testthat::expect_error(pmartR::format_plot(format_metstats, 
                                     comps_color_variable = "blue"))
  testthat::expect_error(pmartR::format_plot(format_metstats, 
                                     panel_variable = "blue"))
  testthat::expect_error(pmartR::format_plot(format_metstats, 
                                     comps_panel_y_axis = "blue"))
  testthat::expect_error(pmartR::format_plot(format_metstats, 
                                     comps_panel_x_axis = "blue"))
  
  
  # Value error #
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_y_range = 4, 
                                     value_y_limits = "free"))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_y_limits = 4))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_y_range = 4, 
                                     value_y_max = 5,
                                     value_y_min = 1))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_y_limits = "free", 
                                     value_y_max = 5,
                                     value_y_min = 1))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_y_limits = c("free", "free")))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_y_range = c("free", "free")))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_y_range = c(1, 2)))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_y_range = -5))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_y_range = 0))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_y_max = "0"))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_y_max = c(1,3)))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_y_min = "0"))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_y_min = c(1,3)))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_plot_type = c(1,3)))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_plot_type = c("box", 
                                                         "boxpoint", 
                                                         "raster")))
  
  testthat::expect_warning(
    pmartR:::validate_format_plot_input(
      format_pepobject, 
      value_panel_y_axis = "log2_abundance",
      value_panel_x_axis = "log2_abundance",
      value_color_variable = "Group",
      panel_variable = "Mass_Tag_ID"),
    "identical")

  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_panel_x_axis = "Peptide",
                                     panel_variable = "Peptide"))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_panel_y_axis = "Peptide",
                                     panel_variable = "Peptide"))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_color_variable = "Peptide",
                                     panel_variable = "Peptide"))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_color_variable = "blue"))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     panel_variable = "blue"))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_panel_y_axis = "blue"))
  testthat::expect_error(pmartR::format_plot(format_pepobject, 
                                     value_panel_x_axis = "blue"))
})


## Test y-value messages

testthat::test_that("Correct messages in generate_plot_message", {
  testthat::expect_message(
    pmartR:::generate_plot_message(format_pepobject),
    "No specified value y-axis parameters. Axis y-limits will be scaled per plot, as per y_limits = 'free'.")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_pepobject, value_y_max = 3),
    "No range or limits specified; Axis y-limits will be scaled per plot with a maximum of y_max.")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_pepobject, value_y_min = 3), 
    "No range or limits specified; Axis y-limits will be scaled per plot with a minimum of y_min.")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_pepobject, value_y_limits = "free"), 
    "Specified value y-limit: 'free'. Axis y-limits will be scaled per plot.")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_pepobject, value_y_limits = "fixed"), 
    "Specified value y-limit: 'fixed'. Axis y-limits will fixed for all plots based on maximum and minimum y-values.")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_pepobject, value_y_limits = "fixed",  value_y_max = 3), 
   "Specified value y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a maximum of y_max. Specified y_max")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_pepobject, value_y_limits = "fixed",  value_y_min = 3), 
    "Specified value y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a minimum of y_min. Specified y_min")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_pepobject, value_y_range = 10), 
    " units, split over the median.")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_pepobject, value_y_range = 10,  value_y_min = 3), 
    " units from the y_min. Specified y_min: ")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_pepobject, value_y_range = 10,  value_y_max = 3), 
    " units from the y_max. Specified y_max: ")

  
  testthat::expect_message(
    pmartR:::generate_plot_message(format_metstats),
    "No specified comparison y-axis parameters. Axis y-limits will be scaled per plot, as per y_limits = 'free'.")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_metstats, comps_y_max = 3),
    "No range or limits specified; Axis y-limits will be scaled per plot with a maximum of y_max.")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_metstats, comps_y_min = 3), 
    "No range or limits specified; Axis y-limits will be scaled per plot with a minimum of y_min.")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_metstats, comps_y_limits = "free"), 
    "Specified comparison y-limit: 'free'. Axis y-limits will be scaled per plot.")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_metstats, comps_y_limits = "fixed"), 
    "Specified comparison y-limit: 'fixed'. Axis y-limits will fixed for all plots based on maximum and minimum y-values.")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_metstats, comps_y_limits = "fixed",  comps_y_max = 3), 
    "Specified comparison y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a maximum of y_max. Specified y_max")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_metstats, comps_y_limits = "fixed",  comps_y_min = 3), 
    "Specified comparison y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a minimum of y_min. Specified y_min")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_metstats, comps_y_range = 10), 
    " units, split over the median.")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_metstats, comps_y_range = 10,  comps_y_min = 3), 
    " units from the y_min. Specified y_min: ")
  testthat::expect_message(
    pmartR:::generate_plot_message(format_metstats, comps_y_range = 10,  comps_y_max = 3), 
    " units from the y_max. Specified y_max: ")
})

## Test expected dataframe generation

## Generate plotting mid-point test data ##
## Nesting is avoided due to large memory allotment (84 Mb to 108 Gb)

testnestdat <- purrr::map(Valid_format_dat, function(trellData){
  if(!is.null(trellData$comp_stats)){
    ## Generate nested data ##
    if ("Group_DF" %in% colnames(trellData$summary_stats)) {
      group_df_name <- "Group_DF"
    } else {
      group_df_name <- "Group"
    }
    plotter <- tidyr::separate(trellData$comp_stats, Comparison,
                               c("comp1", "comp2"), sep = "_vs_", 
                               remove = FALSE) %>%
      reshape2::melt(id.vars = names(trellData$comp_stats),
                     value.name = group_df_name)
    plotter <- suppressWarnings(dplyr::left_join(plotter, trellData$summary_stats))
    if(!is.null(trellData$data_values)){
      plotter <- suppressWarnings(dplyr::left_join(trellData$data_values, plotter))
      plotter[[attr(format_isoobject, "cname")$edata_cname]] <- NULL
      return(plotter)
    }
  } else {
    plotter <- trellData$data_values
    plotter[[attr(format_isoobject, "cname")$edata_cname]] <- NULL
    return(plotter)
  }
})
    
testnestdat[[1]]

testthat::test_that("Expected dataframes are generated before plotting", {
 
  purrr::map2(Valid_format_dat, testnestdat, function(trellData, nestdata){
    
    print("hang here 0")
    df <- pmartR:::add_plot_features(trellData, nestdata, 
                               value_panel_y_axis = grep("abundance", colnames(nestdata)),
                               panel_variable = attr(trellData, "cname")$edata_cname)
    
    print("hang here 1")
    # Expected columns #
    if (!is.null(trellData$comp_stats)){
      addcol <- c("bord", "text", "label")
      expect_col <- c(colnames(trellData$comp_stats), "comp1", "comp2", group_df_name, colnames(trellData$summary_stats))
      if (!is.null(trellData$data_values)){
        expect_col <- c(expect_col, colnames(trellData$data_values))
      }
    } else {
      print("hang here 2")
      addcol <- "text"
      expect_col <- colnames(trellData$data_values)
    }
    
    # Expected rows #
    print("hang here 3")
    print(colnames(df))
    print(c(expect_col, addcol))
    # testthat::expect_true(all(colnames(df) %in% c(expect_col, addcol)))
    #testthat::expect_true()
    
  })
})
    
  
  format_isoobject
  format_lipobject
  format_metboth
  format_metobject
  format_pepobject
  format_proobject

  format_isoobject
  
  testisoval <- tidyr::nest(format_isoobject$data_values, -!!attr(format_isoobject, "cname")$edata_cname)
  
  pmartR:::add_plot_features(format_isoobject, 
                             testiso$data[[1]], 
                             value_panel_y_axis = "abundance",
                             panel_variable = attr(format_isoobject, "cname")$edata_cname)
   
  ?pmartR:::add_plot_features
  
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
                # testthat::expect_identical(nrow(plotDat), length(panelVar)) ######### non-Subsetted version
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