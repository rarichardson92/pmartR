testthat::context("Test as.trellData output (approx. 1 minute)")
library(pmartR)
library(pmartRdata)
################################################################################
## Generate testing data ##

# Load pmartR test data #
isobaric_object    <- pmartRdata::isobaric_object
lipid_object       <- pmartRdata::lipid_object
metab_object       <- pmartRdata::metab_object
pep_object         <- pmartRdata::pep_object
pro_object         <- pmartRdata::pro_object
techrep_pep_object <- pmartRdata::techrep_pep_object

# Log transform untransformed data #
isobaric_object    <- pmartR::edata_transform(isobaric_object, "log2")
lipid_object       <- pmartR::edata_transform(lipid_object, "log10")
metab_object       <- pmartR::edata_transform(metab_object, "log")
pep_object         <- pmartR::edata_transform(pep_object, "log2")
techrep_pep_object <- pmartR::edata_transform(techrep_pep_object, "log10")

isobaric_object    <- pmartR::normalize_isobaric(isobaric_object, apply_norm = T)

# Set appropriate group designations #
isobaric_object    <- pmartR::group_designation(isobaric_object, "Group")
lipid_object       <- pmartR::group_designation(lipid_object, "Condition")
metab_object       <- pmartR::group_designation(metab_object, "Condition")
pep_object         <- pmartR::group_designation(pep_object, "Condition")
pro_object         <- pmartR::group_designation(pro_object, "Condition")
techrep_pep_object <- pmartR::group_designation(techrep_pep_object,
                                                c("FACTOR", "DILUTION"))


# Protein quantification from pepData with mappings #

qpro1 <- pmartR::protein_quant(isobaric_object, "rrollup")
qpro2 <- pmartR::protein_quant(pep_object, "rrollup")

# Stats Filters #
# Not run on non-transformed data
isobaric_object <- pmartR::applyFilt(pmartR::imdanova_filter(isobaric_object),
                                      isobaric_object,
                                      min_nonmiss_anova = 2,
                                      min_nonmiss_gtest = 3)
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

qpro1 <- pmartR::applyFilt(pmartR::imdanova_filter(qpro1),
                           qpro1,
                           min_nonmiss_anova = 2,
                           min_nonmiss_gtest = 3)
qpro2 <- pmartR::applyFilt(pmartR::imdanova_filter(qpro2),
                           qpro2,
                           min_nonmiss_anova = 2,
                           min_nonmiss_gtest = 3)

# Normalize #

isobaric_object    <- pmartR::normalize_global(isobaric_object, "all", "median", apply_norm = TRUE, backtransform = TRUE)
lipid_object       <- pmartR::normalize_global(lipid_object, "all", "median", apply_norm = TRUE, backtransform = TRUE)
metab_object       <- pmartR::normalize_global(metab_object, "all", "median", apply_norm = TRUE, backtransform = TRUE)
pep_object         <- pmartR::normalize_global(pep_object, "all", "median", apply_norm = TRUE, backtransform = TRUE)
pro_object         <- pmartR::normalize_global(pro_object, "all", "median", apply_norm = TRUE, backtransform = TRUE)
techrep_pep_object <- pmartR::normalize_global(techrep_pep_object, "all", "median", apply_norm = TRUE, backtransform = TRUE)
qpro1              <- pmartR::normalize_global(qpro1, "all", "median", apply_norm = TRUE, backtransform = TRUE)
qpro2              <- pmartR::normalize_global(qpro2, "all", "median", apply_norm = TRUE, backtransform = TRUE)

# Generate stats (requires log-transformed data) #
isobaric_stats     <- pmartR::imd_anova(isobaric_object,  
                                        test_method = "anova")
lipid_stats        <- pmartR::imd_anova(lipid_object,
                                        test_method = "gtest")
metab_stats        <- pmartR::imd_anova(metab_object, 
                                        test_method = "combined")
pep_stats          <- pmartR::imd_anova(pep_object,
                                        test_method = "anova")
pro_stats          <- pmartR::imd_anova(pro_object,
                                        test_method = "gtest")
techrep_pep_stats  <- pmartR::imd_anova(techrep_pep_object,
                                        test_method = "combined")

qpro1_stats        <- pmartR::imd_anova(qpro1, test_method = "combined")
qpro2_stats        <- pmartR::imd_anova(qpro2, test_method = "combined")



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
badstat6 <- metab_stats
attr(badstat6 , "group_DF") <- NULL
badstat7 <- metab_stats
attr(badstat7 , "statistical_test") <- NULL

baddat <- metab_object
baddat$e_data <- NULL
baddat2 <- metab_object
baddat2$f_data$SampleID <- NULL
baddat3 <- metab_object
baddat3$e_data <- baddat3$e_data[1:3]
baddat4 <- metab_object
attr(baddat4, "group_DF") <- NULL
baddat5 <- metab_object
baddat5$e_data <- NULL
baddat6 <- metab_object
baddat6$f_data[[pmartR:::get_fdata_cname(baddat6)]] <- NULL

# Expected function input error throwing - as.trellData to validate_input #

testthat::test_that("Correct as.trellData() error throwing for arguments omicsData, omicsStats", {
  testthat::expect_error(
    pmartR::as.trellData(omicsData = 23, omicsStats = 34), 
    "must be of class")
  testthat::expect_error(
    pmartR::as.trellData(omicsStats = 34), 
    "must be of the class")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = "23"),  
    "must be of class") 
  testthat::expect_error(
    pmartR::as.trellData(omicsData = list()), 
    "Empty list")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = list(), omicsStats = list()), 
    "Empty list")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = NULL, omicsStats = NULL), 
    "requires at least")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = c(), omicsStats = c()), 
    "requires at least")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = pro_object, omicsStats = pro_object), 
    "class 'statRes'")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = metab_stats, omicsStats = metab_stats), 
    "omicsData must be of class")
  testthat::expect_warning(
    pmartR::as.trellData(omicsData = NULL, omicsStats = pro_object), 
    "Input reordered")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = metab_stats, omicsStats = pro_object), 
    "Non-matching cname attributes")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = pro_object, omicsStats = metab_stats), 
    "Non-matching cname attributes")
  
  # Manipulated data
  testthat::expect_error(
    pmartR::as.trellData(badstat), 
    "Mismatched rows")
  testthat::expect_error(
    pmartR::as.trellData(badstat2), 
    "should contain only and all of the following dataframes:")
  testthat::expect_error(
    pmartR::as.trellData(badstat3),  
    "should contain only and all of the following dataframes:")
  testthat::expect_error(
    pmartR::as.trellData(badstat4), 
    "column must be present in all omicsStats dataframes")
  testthat::expect_error(
    pmartR::as.trellData(badstat5), 
    "Number of columns in omicsStats dataframes is different than expected")
  testthat::expect_error(
    pmartR::as.trellData(badstat, metab_object), 
    "Mismatched rows")
  testthat::expect_error(
    pmartR::as.trellData(badstat2, metab_object), 
    "Requires compatible identifiers")
  testthat::expect_error(
    pmartR::as.trellData(badstat3, metab_object), 
    "should contain only and all of the following dataframes:")
  testthat::expect_error(
    pmartR::as.trellData(badstat4, metab_object), 
    "column must be present in all omicsStats dataframes")
  testthat::expect_error(
    pmartR::as.trellData(badstat5, metab_object), 
    "Number of columns in omicsStats dataframes is different than expected")
  testthat::expect_error(
    pmartR::as.trellData(baddat), 
    "Omicsdata requires both e_data and f_data.")
  testthat::expect_error(
    pmartR::as.trellData(baddat2), 
    "column must be present in omicsData f_data.")
  testthat::expect_error(
    pmartR::as.trellData(baddat3), 
    "column in f_data does not match column names in e_data" )
  testthat::expect_error(
    pmartR::as.trellData(baddat, metab_stats), 
    "Biomolecules in omicsStats do not match biomolecules in omicsData.")
  testthat::expect_error(
    pmartR::as.trellData(baddat2, metab_stats), 
    "SampleID column does not match between omicsData and omicsStats")
  testthat::expect_error(
    pmartR::as.trellData(baddat3, metab_stats), 
    "SampleID column in f_data does not match column names in e_data")
  testthat::expect_error(
    pmartR::as.trellData(baddat4), 
    "group_designation()")
  testthat::expect_error(
    pmartR::as.trellData(baddat4, metab_stats), 
    "group_designation()")
  testthat::expect_error(
    pmartR::as.trellData(baddat5), 
    "Omicsdata requires both e_data and f_data")
  testthat::expect_error(
    pmartR::as.trellData(baddat5, metab_stats), 
    "Biomolecules in omicsStats do")
  testthat::expect_error(
    pmartR::as.trellData(baddat6), 
    "column must be")
  testthat::expect_error(
    pmartR::as.trellData(baddat6, metab_stats), 
    "column does not match")
  testthat::expect_error(
    pmartR::as.trellData(badstat6), 
    "group_designation()")
  testthat::expect_error(
    pmartR::as.trellData(metab_object, badstat6), 
    "group_designation()")
  testthat::expect_error(
    pmartR::as.trellData(badstat7), 
    "must be combined, anova, or gtest")
  testthat::expect_error(
    pmartR::as.trellData(metab_object, badstat7), 
    "must be combined, anova, or gtest")
})

# as.trellData -> recursive_format -> as.trellData (checks list input) -> as.trellData -> as.trellData (checks individual omics objects)
testthat::test_that("Subfunction recursive_format correctly throws errors", {
  
  testthat::expect_error(
    pmartR::as.trellData(omicsData = list(pep_object, qpro2), omicsStats = list(pep_stats)),
    "List length does not match;")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = list(pep_object, qpro2), omicsStats = list(pep_stats, pep_stats)),
    "Lists in omicsData and omicsStats have mismatched cname attributes.")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = list(pep_object, qpro2), omicsStats = list(pep_stats, pep_stats, pep_stats)),
    "List length does not match;")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = list(pep_object, qpro2), omicsStats = list(metab_stats, metab_stats)),
    "Lists in omicsData and omicsStats have mismatched cname attributes")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = c(pep_object, qpro2), omicsStats = list(pep_stats, pep_stats)),
    "List length does not match;")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = list(pep_object, pep_object)), 
    "Only one pepData")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = list(pep_object, pep_object, pep_object)), 
    "List length != 2")
  testthat::expect_error(
    pmartR::as.trellData(omicsStats = list(pep_stats, pep_stats)),
    "Only one stats object derived from pepData")
  testthat::expect_error(
    pmartR::as.trellData(omicsStats = list(pep_object, pep_object)),
    "Only one stats object derived from pepData")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = c(isobaric_object, qpro1)),
    "List/vector entry error")
  testthat::expect_error(
    pmartR::as.trellData(omicsData = c(isobaric_object, qpro1), c(isobaric_stats, qpro1_stats)),
    "List/vector entry error")
  testthat::expect_error(
    pmartR::as.trellData(omicsStats = c(isobaric_stats, qpro1_stats)),
    "List/vector entry error")
  testthat::expect_error(pmartR::as.trellData(list(pep_object, lipid_object)),
                         "Only pepData and proData are valid")
  testthat::expect_error(pmartR::as.trellData(list(pep_object, lipid_object), list(pep_stats, lipid_stats)),
                         "Only pepData and proData are valid")
  testthat::expect_error(pmartR::as.trellData(omicsStats = list(pep_stats, lipid_stats)),
                         "Only stats derived from ")
  testthat::expect_error(pmartR::as.trellData(omicsStats = list(pep_stats, pep_stats), omicsData = list(pep_object, pep_object)),
                         "Only one pepData, one ")
  testthat::expect_error(pmartR::as.trellData(omicsStats = list(pep_stats, pro_stats, pep_stats)),
                         "List length != 2")
  testthat::expect_error(pmartR::as.trellData(omicsStats = list(pep_stats, lipid_stats)),
                         "Only stats derived")
  testthat::expect_error(pmartR::as.trellData(omicsStats = list(pep_stats),  omicsData = list(pep_object, pro_object)),
                         "List length does not match")
})

################################################################################

## Generate function outputs from test data ##

# Data only #
format_isoobject   <- pmartR::as.trellData(isobaric_object)
format_lipobject   <- pmartR::as.trellData(list(lipid_object))
format_metobject   <- pmartR::as.trellData(metab_object)
format_pepobject   <- pmartR::as.trellData(pep_object)
format_proobject   <- pmartR::as.trellData(pro_object)
format_tecobject   <- pmartR::as.trellData(techrep_pep_object)

format_qpro1       <- pmartR::as.trellData(qpro1)
format_qpro2       <- pmartR::as.trellData(qpro2)

# Stats only #
format_isostats  <- pmartR::as.trellData(isobaric_stats)
format_lipstats    <- pmartR::as.trellData(lipid_stats)
format_metstats    <- pmartR::as.trellData(metab_stats)
format_pepstats    <- pmartR::as.trellData(pep_stats)
format_prostats    <- pmartR::as.trellData(omicsStats = list(pro_stats))
format_tecstats    <- pmartR::as.trellData(techrep_pep_stats)

format_qpro1stats  <- pmartR::as.trellData(qpro1_stats)
format_qpro2stats  <- pmartR::as.trellData(qpro2_stats)

# Both data and stats #

format_isoboth    <- pmartR::as.trellData(isobaric_object, isobaric_stats)
format_lipboth    <- pmartR::as.trellData(list(lipid_object), list(lipid_stats))
format_metboth   <- pmartR::as.trellData(metab_object, metab_stats)
format_pepboth    <- pmartR::as.trellData(list(pep_object), list(pep_stats))
format_proboth    <- pmartR::as.trellData(pro_object, pro_stats)
format_tecboth    <- pmartR::as.trellData(techrep_pep_object, techrep_pep_stats)

format_qpro1both <- pmartR::as.trellData(qpro1_stats, qpro1)
format_qpro2both <- pmartR::as.trellData(qpro2_stats, qpro2)

# List based on pepData and quantified proData #
# Data
forlist_qpro1dat  <- pmartR::as.trellData(omicsData = list(isobaric_object, qpro1))
forlist_qpro2dat  <- pmartR::as.trellData(omicsData = list(pep_object, qpro2))

# Stats                                                                       
forlist_qpro1stat <- pmartR::as.trellData(omicsStats = list(qpro1_stats, isobaric_stats))
forlist_qpro2stat <- pmartR::as.trellData(omicsStats = list(qpro2_stats, pep_stats))

# Both #

forlist_qpro1both <- pmartR::as.trellData(omicsStats = list(qpro1_stats,
                                                           isobaric_stats),
                                         omicsData = list(qpro1,
                                                          isobaric_object))
forlist_qpro2both <- pmartR::as.trellData(omicsStats = list(qpro2_stats,
                                                           pep_stats),
                                         omicsData = list(qpro2,
                                                          pep_object))

################################################################################

Valid_format_dat <- list(
  format_isoobject, 
  format_lipobject, 
  format_metobject, 
  format_pepobject, 
  format_proobject, 
  format_tecobject, 
  format_qpro1, 
  format_qpro2
)

Vfd_obs <- list(
  isobaric_object, 
  lipid_object, 
  metab_object, 
  pep_object,
  pro_object, 
  techrep_pep_object, 
  qpro1, 
  qpro2)

##

Valid_format_stat <- list(
  format_isostats,
  format_lipstats,
  format_metstats, 
  format_pepstats, 
  format_prostats,
  format_tecstats, 
  format_qpro1stats, 
  format_qpro2stats
)

Vfs_obs <- list(
  isobaric_stats,
  lipid_stats,
  metab_stats,
  pep_stats,
  pro_stats,
  techrep_pep_stats,
  qpro1_stats,
  qpro2_stats
)

##

Valid_format_both <- list(
  format_isoboth,
  format_lipboth,
  format_metboth,
  format_pepboth,
  format_proboth,
  format_tecboth,
  format_qpro1both,
  format_qpro2both
)

Vfb_obs1 <- Vfd_obs

Vfb_obs2 <- Vfs_obs

#####

Valid_forlist_dat <- list(
  forlist_qpro1dat,
  forlist_qpro2dat
)

Vfld_list <- list(
  list(isobaric_object, qpro1),
  list(pep_object, qpro2)
)

##

Valid_forlist_stat <- list(
  forlist_qpro1stat,
  forlist_qpro2stat
)

Vfls_list <- list(
  list(qpro1_stats, isobaric_stats),
  list(qpro2_stats, pep_stats)
)

##

Valid_forlist_both <- list(
  forlist_qpro1both,
  forlist_qpro2both
)

Vflb_list1 <- Vfld_list

Vflb_list2 <- Vfls_list



################################################################################

## Test as.trellData outputs ##
#### Dimensions and correct columns in output ####

testthat::test_that("Correct dataframe population and dimensions", {
  
  purrr::map2(Valid_format_dat, Vfd_obs, function(formDat, parOb){
    
    ## Print dimension check ## 
    testthat::expect_true(length(capture.output(print(formDat))) == 12)
    testthat::expect_true(length(strsplit(capture.output(print(formDat))[2], 
                                          "[[:blank:]]+")[[1]]) == 6)
    temp <- formDat
    temp$data_values <- temp$data_values[1:5,1:3]
    testthat::expect_true(length(capture.output(print(temp))) == 8)
    testthat::expect_true(length(strsplit(capture.output(print(temp))[2], 
                                          "[[:blank:]]+")[[1]]) == 4)
    temp2 <- formDat
    temp2$data_values <- cbind(temp2$data_values, temp2$data_values)
    testthat::expect_message(capture.output(print(temp2)), "first 5")
    
    # abundance name if un-log transformed
    testthat::expect_true(
      "abundance" %in% colnames(
        pmartR::as.trellData(pmartR::edata_transform(parOb, "abundance"))$data_values))
    
    # Data frame population #
    testthat::expect_length(formDat, 3)
    testthat::expect_null(formDat$comp_stats)
    testthat::expect_null(formDat$summary_stats)
    testthat::expect_false(is.null(formDat$data_value))
    testthat::expect_gt(nrow(formDat$data_value), 0)
    
    # Data value correct columns and number of rows #
    form_cols <- colnames(formDat$data_values)
    
    testthat::expect_match(toString(form_cols), pmartR:::get_edata_cname(parOb))
    testthat::expect_match(toString(form_cols), pmartR:::get_fdata_cname(parOb))
    testthat::expect_true(all(colnames(parOb$f_data) %in% form_cols))
    form_cols <- form_cols[!(form_cols %in% colnames(parOb$f_data))]
    testthat::expect_match(toString(form_cols), "abundance")
    testthat::expect_match(toString(form_cols), "Group")
    
    if(!is.null(parOb$e_meta)){
      testthat::expect_match(toString(form_cols), pmartR:::get_emeta_cname(parOb))
      testthat::expect_true(all(colnames(parOb$e_meta) %in% form_cols))
      mult_row <- nrow(unique(
        parOb$e_meta[c(pmartR:::get_emeta_cname(parOb), pmartR:::get_edata_cname(parOb))]))
    } else {
      mult_row <- length(parOb$e_data[[pmartR:::get_edata_cname(parOb)]])
    }
    
    testthat::expect_identical(mult_row * 
                                 length(parOb$f_data[[pmartR:::get_fdata_cname(parOb)]]),
                               nrow(formDat$data_values))
    return(NULL)
  })
  
##
  purrr::map2(Valid_format_stat, Vfs_obs, function(formStat, parOb){
    
    ## Print dimension check ## 
    testthat::expect_error(pmartR:::print.trellData("blah"), "trellData object")
    testthat::expect_true(length(capture.output(print(formStat))) == 24)
    testthat::expect_true(length(strsplit(capture.output(print(formStat))[2], 
                                          "[[:blank:]]+")[[1]]) == 5)
    testthat::expect_true(length(strsplit(capture.output(print(formStat))[14], 
                                          "[[:blank:]]+")[[1]]) == 6)
    
    temp <- formStat
    temp$summary_stats <- temp$summary_stats[1:5,1:3]
    temp$comp_stats <- temp$comp_stats[1:5,1:3]
    testthat::expect_true(length(capture.output(print(temp))) == 16)
    testthat::expect_true(length(strsplit(capture.output(print(temp))[2], 
                                          "[[:blank:]]+")[[1]]) == 4)
    testthat::expect_true(length(strsplit(capture.output(print(temp))[10], 
                                          "[[:blank:]]+")[[1]]) == 4)
    
    temp2 <- temp
    temp2$summary_stats <- cbind(temp2$summary_stats, temp2$summary_stats)
    testthat::expect_message(capture.output(print(temp2)), "only first 5")
    
    temp3 <- temp
    temp3$comp_stats <- cbind(temp3$comp_stats, temp3$comp_stats)
    testthat::expect_message(capture.output(print(temp3)), "only first 5")
    
    # Data frame population #
    testthat::expect_length(formStat, 3)
    testthat::expect_null(formStat$data_value)
    testthat::expect_false(is.null(formStat$summary_stats))
    testthat::expect_false(is.null(formStat$comp_stats))
    testthat::expect_gt(nrow(formStat$summary_stats), 0)
    testthat::expect_gt(nrow(formStat$comp_stats), 0)
    
    # Comps and summary stats correct columns and number of rows #
    
    form_cols_summ <- colnames(formStat$summary_stats)
    expect_cols_summ <- c(pmartR:::get_edata_cname(parOb), "Group", "Count", "Mean")
    expect_rows_summ <- length(as.character(parOb$Full_results[[pmartR:::get_edata_cname(parOb)]])) *
      length(unique(attr(parOb, "group_DF")$Group))
    testthat::expect_identical(form_cols_summ, expect_cols_summ)
    testthat::expect_identical(nrow(formStat$summary_stats), expect_rows_summ)
    
    form_cols_comp <- colnames(formStat$comp_stats)
    test <- attr(parOb, "statistical_test")
    if (test == "gtest"){
      expect_cols_comp <- c(attr(parOb, "cnames")$edata_cname, "Comparison",
                            "P_value_G", "Fold_change", "Flag")
    } else if (test == "anova"){
      expect_cols_comp <- c(attr(parOb, "cnames")$edata_cname, "Comparison",
                            "P_value_T", "Fold_change", "Flag")
    } else {
      expect_cols_comp <- c(attr(parOb, "cnames")$edata_cname, "Comparison", 
                            "P_value_G", "P_value_T", 
                            "Fold_change", "Flag")
    }
    expect_rows_comp <- length(
      parOb$Full_results[[attr(parOb, "cnames")$edata_cname]]) *
      length(attr(parOb, "comparisons"))
    
    testthat::expect_true(all(form_cols_comp %in% expect_cols_comp))
    testthat::expect_identical(nrow(formStat$comp_stats), expect_rows_comp)
    return(NULL)
  })

 ##

  purrr::pmap(list(Valid_format_both, Vfb_obs1, Vfb_obs2),
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
                testthat::expect_match(toString(form_cols), pmartR:::get_edata_cname(parDat))
                testthat::expect_match(toString(form_cols), pmartR:::get_fdata_cname(parDat))
                testthat::expect_true(all(colnames(parDat$f_data) %in% form_cols))
                form_cols <- form_cols[!(form_cols %in% colnames(parDat$f_data))]
                testthat::expect_match(toString(form_cols), "abundance")
                testthat::expect_match(toString(form_cols), "Group")
                
                if(!is.null(parDat$e_meta)){
                  testthat::expect_match(toString(form_cols), pmartR:::get_emeta_cname(parDat))
                  testthat::expect_true(all(colnames(parDat$e_meta) %in% form_cols))
                  mult_row <- nrow(unique(
                    parDat$e_meta[c(pmartR:::get_emeta_cname(parDat), pmartR:::get_edata_cname(parDat))]))
                } else {
                  mult_row <- length(parDat$e_data[[pmartR:::get_edata_cname(parDat)]])
                }
                
                testthat::expect_identical(mult_row * 
                                             length(parDat$f_data[[pmartR:::get_fdata_cname(parDat)]]),
                                           nrow(formBoth$data_values))
                
                # Comps and summary stats correct columns and number of rows #
                
                form_cols_summ <- colnames(formBoth$summary_stats)
                
                expect_cols_summ <- c(attr(parStat, "cnames")$edata_cname, "Group", "Group_DF", "Count", "Mean")  ## OOD
                if(!is.null(attr(parStat, "cnames")$emeta_cname)){
                  expect_cols_summ <- c(expect_cols_summ, colnames(parDat$e_meta))
                }
                
                if(!is.null(parDat$e_meta)){
                  expect_rows_summ <- sum(parDat$e_meta[[attr(parDat, "cnames")$edata_cname]] %in% 
                                   parStat$Full_results[[attr(parDat, "cnames")$edata_cname]]) *
                    length(unique(as.character(attr(parDat, "group_DF")$Group)))
                } else{
                  expect_rows_summ <- length(parStat$Full_results[[attr(parStat, "cnames")$edata_cname]]) *
                    length(unique(attr(parStat, "group_DF")$Group))
                }
                
                testthat::expect_true(all(form_cols_summ %in% expect_cols_summ))
                testthat::expect_identical(nrow(formBoth$summary_stats), expect_rows_summ)
                
                form_cols_comp <- colnames(formBoth$comp_stats)
                test <- attr(parStat, "statistical_test")
                if (test == "gtest"){
                  expect_cols_comp <- c(attr(parStat, "cnames")$edata_cname, "Comparison",
                                        "P_value_G", "Fold_change", "Flag")
                } else if (test == "anova"){
                  expect_cols_comp <- c(attr(parStat, "cnames")$edata_cname, "Comparison",
                                        "P_value_T", "Fold_change", "Flag")
                } else {
                  expect_cols_comp <- c(attr(parStat, "cnames")$edata_cname, "Comparison", 
                                        "P_value_G", "P_value_T", 
                                        "Fold_change", "Flag")
                }
                
                if(!is.null(attr(parStat, "cnames")$emeta_cname)){
                  expect_cols_comp <- c(expect_cols_comp, colnames(parDat$e_meta))
                }
                if(!is.null(parDat$e_meta)){
                  expect_rows_comp <- sum(
                    parDat$e_meta[[attr(parStat, "cnames")$edata_cname]] %in% 
                      parStat$Full_result[[attr(parStat, "cnames")$edata_cname]]) *
                    length(attr(parStat, "comparisons"))
                } else{
                  expect_rows_comp <- length(
                    parStat$Full_results[[attr(parStat, "cnames")$edata_cname]]) *
                    length(attr(parStat, "comparisons"))
                }
                testthat::expect_true(all(form_cols_comp %in% expect_cols_comp))
                testthat::expect_identical(nrow(formBoth$comp_stats), expect_rows_comp)
                return(NULL)
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
      testthat::expect_match(toString(form_cols), pmartR:::get_edata_cname(parOb))
      testthat::expect_match(toString(form_cols), pmartR:::get_fdata_cname(parOb))
      testthat::expect_true(all(colnames(parOb$f_data) %in% form_cols))
      form_cols <- form_cols[!(form_cols %in% colnames(parOb$f_data))]
      testthat::expect_match(toString(form_cols), "abundance")
      testthat::expect_match(toString(form_cols), "Group")
      
      if(!is.null(parOb$e_meta)){
        testthat::expect_match(toString(form_cols), pmartR:::get_emeta_cname(parOb))
        testthat::expect_true(all(colnames(parOb$e_meta) %in% form_cols))
        mult_row <- nrow(unique(
          parOb$e_meta[c(pmartR:::get_emeta_cname(parOb), pmartR:::get_edata_cname(parOb))]))
      } else {
        mult_row <- length(parOb$e_data[[pmartR:::get_edata_cname(parOb)]])
      }
      
      testthat::expect_identical(mult_row * 
                                   length(parOb$f_data[[pmartR:::get_fdata_cname(parOb)]]),
                                 nrow(formDat$data_values))
      return(NULL)
    })
  })
  
  # Update as other changes happen before errors are fixed
  purrr::map(Valid_forlist_stat, function(index){
    testthat::expect_length(index, 2)
    purrr::map(index, function(formStat){
      testthat::expect_length(formStat, 3)
      testthat::expect_null(formStat$data_value)
      testthat::expect_false(is.null(formStat$summary_stats))
      testthat::expect_false(is.null(formStat$comp_stats))
      return(NULL)
    })
    return(NULL)
  })

  purrr::map(Valid_forlist_both, function(index){
    testthat::expect_length(index, 2)
    purrr::map(index, function(formBoth){
      testthat::expect_length(formBoth, 3)
      testthat::expect_false(is.null(formBoth$data_value))
      testthat::expect_false(is.null(formBoth$summary_stats))
      testthat::expect_false(is.null(formBoth$comp_stats))
      return(NULL)
    })
    return(NULL)
  })
})

testthat::test_that("Format list slices are equal to non-lists", {
  testthat::expect_equal(forlist_qpro2dat[[1]], format_pepobject)
  testthat::expect_equal(forlist_qpro2dat[[2]], format_qpro2)
})

################################################################################
################################################################################

## Test expected dataframe generation of add_plot_features
testdat <- function(dat){
  df <- dat$data_values
  random <- floor(runif(1, min = 1, max = nrow(df)))
  pick <- df[attr(dat, "cname")$edata_cname][random,]
  subset <-  df[df[attr(dat, "cname")$edata_cname] == as.character(pick),]
  abundance <- grep("abundance", colnames(subset), value = TRUE)
  df2 <- pmartR:::add_plot_features(dat,
                                    subset,
                                    y_axis = abundance,
                                    panel_variable = attr(dat, "cname")$edata_cname)
  testthat::expect_true(nrow(subset) == nrow(df2))
  testthat::expect_true(ncol(subset) == ncol(df2) - 1)
  testthat::expect_true(all(colnames(df2) %in% c(colnames(subset), "text")))
}

teststat <- function(stat){
  df <- stat$comp_stats
  random <- floor(runif(1, min = 1, max = nrow(df)))
  pick <- df[attr(stat, "cname")$edata_cname][random,]
  subset <-  df[df[attr(stat, "cname")$edata_cname] == as.character(pick),]
  df2 <- pmartR:::add_plot_features(stat,
                                    subset,
                                    y_axis = "Fold_change",
                                    panel_variable = attr(stat, "cname")$edata_cname)
  testthat::expect_true(nrow(subset) == nrow(df2))
  testthat::expect_true(ncol(subset) == ncol(df2) - 3)
  testthat::expect_true(all(colnames(df2) %in% c(colnames(subset), "text", "labels", "bord")))
}

testthat::test_that("Correct add_plot_features dimensions", {
  purrr::map(Valid_format_dat, testdat)
  purrr::map(Valid_format_stat, teststat)
  purrr::map(Valid_format_both, testdat)
  purrr::map(Valid_format_both, teststat)
})

##########################

# list_y_limits errors #
  testthat::test_that("Correct list_y_limits digestion", {
  # list_y_limits errors #
  testthat::expect_warning( 
    pmartR:::list_y_limits(plot_type = "abundance_global", y_limits = 'free'),
    "Scale option 'free' is not valid with global plots")
  testthat::expect_warning( 
    pmartR:::list_y_limits(plot_type = "abundance_heatmap", y_limits = 'free'),
    "y_limits are not supported ")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_global", y_limits = list('free')),
    "Invalid format. List inputs must be named")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_global", 
                           y_limits = list(abundance_global = list(scale = 'free', scale = 'free'))),
    "List inputs may not have duplicate names.")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_global", 
                           y_limits = list(abundance_global = list(truth = 'free'))),
    "List names for y_limits are restricted")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_global", 
                           y_limits = list(abundance_global = list(scale = 'free', scale = 'free'),
                                           abundance_boxplot = list(scale = 'free', scale = 'free'))),
    "List inputs may not have duplicate names.")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", 
                           y_limits = list(abundance_boxplot = list(truth = 'free'),
                                           missing_bar = list(scale = 'free'))),
    "List names for y_limits are restricted")
  
  testthat::expect_error(pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                                                y_limits = list(abundance_boxplot = list(scale = 'free'),
                                                                missing_bar = 10),
                                                "For mixed input of lists and non-lists"))
  
  testthat::expect_warning( 
    pmartR:::list_y_limits(plot_type = "foldchange_global", y_limits = 'free'),
    "Scale option 'free' is not valid with global plots")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                           y_limits = list(min = 4, max = 3)),
    "Set min must be less than set max")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                           y_limits = list(min = "lp", max = 3)),
    "Min must be a numeric")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                           y_limits = list(max = "lp")),
    "Max must be a numeric")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                           y_limits = list(range = -5)),
    "Range must be a numeric")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                           y_limits = list(scale = -5)),
    "Scale must be a character")
  testthat::expect_error(pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                                                y_limits = -5),
                         "Invalid input")
  
  testthat::expect_error( pmartR:::list_y_limits(plot_type = "abundance_boxplot", 
                                                 y_limits = list(blah = list(min = 1, max = 2))),
                          "Names in y_limits should match plot_type")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                           y_limits = list(bro = 4, max = 3)),
    "'min', 'max', 'range', and 'scale'")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", 
                           y_limits = list(min = 5, min = 6)),
    "List inputs may not have duplicate names.")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", 
                           y_limits = list(list(list(min = 3, max = 4)))),
    "List of lists of lists")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", 
                           y_limits = list(list(min = 4, max = 5))),
    "inputs must be named")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", 
                           y_limits = list(abundance_boxplot = list(min = 4, max = 5), 
                                           abundance_boxplot = list(min = 4, max = 5))),
    "duplicate names")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "blah", 
                           y_limits = "free"),
    "Must be a character string in the following")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = list(blah = "free")),
    "List names for y_limits are restricted")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = list(abundance_boxplot = list("free"))),
    "List inputs must be named")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = "blah"),
    "Invalid input")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = list(
      scale = 'free',
      range = 5,
      min = 5)),
    "are not supported")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = list(
      max = 8,
      range = 5,
      min = 5)),
    "are not supported")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = list(
      max = 8,
      scale = 'free',
      min = 5)),
    "are not supported")
})



################################################################################
################################################################################

testthat::context("Test format_plot output, (approx. 2 minute)") ## takes some time for pep/pro data (large), not used

## Generate nested plot structures ##

plot_isoobject <- pmartR:::format_plot(format_isoobject, panel_variable = "Protein", 
                                      plot_type = "abundance_heatmap", p_val = NULL,
                                      interactive = TRUE)
plot_lipobject <- pmartR:::format_plot(format_lipobject, plot_type = "abundance_global", 
                                      p_val = NULL, y_limits = c(3,4),
                                      interactive = TRUE)
plot_metobject <- pmartR:::format_plot(format_metobject, plot_type = "missing_bar",
                                      y_limits = list("missing_bar" = list(max = 0.5)),
                                      interactive = TRUE)
plot_qpro1     <- pmartR:::format_plot(format_qpro1, plot_type = "abundance_boxplot", 
                                      y_limits = "fixed", interactive = TRUE)

plot_isostats   <- pmartR:::format_plot(format_isostats, plot_type = "foldchange_bar", 
                                       y_limits = list(foldchange_bar = list(range = 10),
                                                       missing_bar = c(0.2, 0.5)),
                                       interactive = TRUE)
plot_lipstats   <- pmartR:::format_plot(format_lipstats, plot_type = "foldchange_global",
                                       y_limits = c(3,6), interactive = TRUE)
plot_metstats   <- pmartR:::format_plot(format_metstats, plot_type = "missing_bar", 
                                       y_limits = list("missing_bar" = list(max = 0.5)))
plot_qpro1stats <- pmartR:::format_plot(format_qpro1stats, plot_type = "foldchange_bar",
                                       y_limits = "fixed", interactive = TRUE, plot_text = TRUE)

plot_isoboth    <- pmartR:::format_plot(format_isoboth, panel_variable = "Protein",
                                       plot_type = "presence_heatmap", interactive = TRUE)
plot_isoboth2    <- pmartR:::format_plot(format_isoboth, panel_variable = "Protein",
                                       plot_type = "foldchange_heatmap", interactive = TRUE)
plot_lipboth    <- pmartR:::format_plot(format_lipboth, plot_type = "abundance_boxplot",
                                       y_limits = c(3,6), interactive = TRUE)
plot_metboth    <- pmartR:::format_plot(format_metboth, plot_type = "foldchange_bar",
                                       y_limits = c(3,6), interactive = TRUE)
plot_qpro1both  <- pmartR:::format_plot(format_qpro1both, custom_plot = "custom_fn")

# Plot object lists #
Valid_plot_dat <- list(
  plot_isoobject, 
  plot_lipobject, 
  plot_metobject, 
  plot_qpro1)

Valid_plot_stat <- list(
  plot_isostats, 
  plot_lipstats, 
  plot_metstats, 
  plot_qpro1stats
)

Valid_plot_both <- list(
  plot_isoboth, 
  plot_lipboth, 
  plot_metboth, 
  plot_qpro1both
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
    pmartR:::format_plot(format_isoobject, plot_type = "missing_bar",
                        y_limits = list("missing_bar" = list(max = 5))),
    "Minimum and maximum for missing_bar is only supported for proportions")
  testthat::expect_error(
    pmartR:::format_plot(format_isoobject, plot_type = "missing_bar",
                        y_limits = list("missing_bar" = list(min = 5))),
    "Minimum and maximum for missing_bar is only supported for proportions")
  testthat::expect_error(
    pmartR:::format_plot(format_isoobject, plot_type = NULL),
    "No plotting information found in plot_type or custom_plot.")
  testthat::expect_error(
    pmartR:::format_plot(format_isoobject, plot_type = "abundance_boxplot",
                        interactive = "blah"),
    "interactive must be a logical")
  testthat::expect_error(
    pmartR:::format_plot(format_isoobject, plot_type = "abundance_boxplot",
                        plot_text = "blah"),
    "must be a logical")
  testthat::expect_error(
    pmartR:::format_plot(format_isostats, plot_type = "foldchange_bar", panel_variable = "blah"),
    "Panel_variable is not valid")
  testthat::expect_error(
    pmartR:::format_plot(format_isoobject, custom_plot = "NULL"),
    "custom_plot function not found!")
  testthat::expect_error(
    pmartR:::format_plot(format_isoobject, custom_plot = c("NULL", "blah")),
    "length 1")
  testthat::expect_error(
    pmartR:::format_plot(format_isoobject, plot_type = "abundance_bxplot"),
    "plot_type specified is not supported")
  testthat::expect_error(
    pmartR:::format_plot(format_isoobject, plot_type = c("abundance_boxplot", "abundance_boxplot")),
    "Invalid plot_type input")
  testthat::expect_error(
    pmartR:::format_plot(format_isoobject, plot_type = "abundance_heatmap", panel_variable = "Peptide"),
    "Heatmaps require a panel_variable that is not the same as edata cname")
  testthat::expect_error(
    pmartR:::format_plot(format_isoobject, plot_type = "foldchange_heatmap", panel_variable = "Protein"),
    "Statistical comparisons")
  testthat::expect_error(
    pmartR:::format_plot(format_isostats, plot_type = "abundance_heatmap", panel_variable = "Count"),
    "Data values")
  testthat::expect_error(
    pmartR:::format_plot(format_isoobject, custom_plot = "mean"),
    "Generated custom plot")
  testthat::expect_warning(
    pmartR:::format_plot(format_lipobject, custom_plot = "custom_fn", y_limits = 'free'),
    "y-limits are not applied to custom plots.")
  testthat::expect_error(
    pmartR:::format_plot(badform, plot_type = "abundance_boxplot"),
    "trellData must be of the class 'trellData'")
  testthat::expect_error(
    pmartR:::format_plot(badform2, plot_type = "foldchange_bar"),
    "No data values or comparison statistics in trellData to plot")
  testthat::expect_error(pmartR:::format_plot(badform3, plot_type = "abundance_boxplot"),
                         "No data values or comparison statistics in trellData to plot")
  testthat::expect_error(pmartR:::format_plot(badform4, plot_type = "abundance_boxplot"),
                         " must be NULL or data.frames")
  testthat::expect_error(pmartR:::format_plot(badform5, plot_type = "abundance_boxplot"),
                         "combined, anova, or gtest")
  testthat::expect_error(pmartR:::format_plot(badform6, plot_type = "abundance_boxplot"),
                         "combined, anova, or gtest")
  testthat::expect_error(pmartR:::format_plot(badform7, plot_type = "abundance_boxplot"),
                         "No data values or comparison statistics in trellData to plot")
  testthat::expect_error(pmartR:::format_plot(badform8, plot_type = "abundance_boxplot"),
                         "Where comp_stats is generated, summary_stats are required and vice versa.")
  
  # p val error #
  testthat::expect_error(pmartR:::format_plot(format_pepobject, p_val = NA, plot_type = "abundance_boxplot"),
                         "p_val must be")
  testthat::expect_error(pmartR:::format_plot(format_pepobject, p_val = "ghjk", plot_type = "abundance_boxplot"),
                         "p_val must be")
  testthat::expect_error(pmartR:::format_plot(format_pepobject, p_val = c(0.3, 1), plot_type = "abundance_boxplot"),
                         "p_val must be")
})


##### Test dimensions #####

Valid_format_stat2 <- list(
  format_isostats,
  format_lipstats,
  format_metstats, 
  format_qpro1stats 
)

Valid_format_dat2 <- list(
  format_isoobject, 
  format_lipobject, 
  format_metobject, 
  format_qpro1
)

Valid_format_both2 <- list(
  format_isoboth,
  format_lipboth,
  format_metboth,
  format_qpro1both
)

plotDat <- Valid_plot_dat[[1]]
parForm <- Valid_format_dat2[[1]]

testthat::test_that("Correct format_plot nested table dimensions", {
  purrr::map2(c(Valid_plot_dat, Valid_plot_stat),
              c(Valid_format_dat2, Valid_format_stat2),
              function(plotDat, parForm){
                if(is.null(parForm$data_values)){
                  panelVar <- parForm$comp_stats[[colnames(plotDat)[1]]]
                  } else {
                    panelVar <- parForm$data_values[[colnames(plotDat)[1]]]
                    }
                expect_cols <- c(colnames(plotDat)[1], "panel")
                testthat::expect_identical(nrow(plotDat), length(unique(panelVar)))
                testthat::expect_true(all(as.character(unlist(plotDat[1])) %in% 
                                            as.character(unique(panelVar))))
                testthat::expect_identical(colnames(plotDat), expect_cols)
                testthat::expect_true(inherits(plotDat$panel, c(
                  "trelliscope_panels",
                  "list")))
                testthat::expect_true(all(
                  unlist(purrr::map(plotDat$panel, function(plot) inherits(
                    plot, c("gg", "plotly", "rbokeh"))))))
  })
})

# ################################################################################
# ################################################################################

testthat::context("Test data_cogs output")

testthat::test_that("Correct data_cogs error/warning throwing", {
  testthat::expect_error(
    pmartR:::data_cogs(format_lipobject, plot_lipobject), 
    "nested_plot must be a data.frame passed from format_plot")
  testthat::expect_error(
    pmartR:::data_cogs(plot_lipobject, format_lipobject, "blah"), 
    "p_val must be a numeric of length 1")
  testthat::expect_error(
    pmartR:::data_cogs(plot_lipobject, format_lipobject, try_URL = "blah"), 
    "try_URL must be a logical")
  testthat::expect_error(
    pmartR:::data_cogs(plot_lipobject, badform), 
    "trellData must be class trellData.")
  testthat::expect_error(
    pmartR:::data_cogs(plot_lipobject, format_lipobject, custom_cog = 5), 
    "custom_cog argument must")
})

## Generate cogs ##

cog_list <- purrr::map2(c(Valid_plot_dat, Valid_plot_stat, Valid_plot_both), 
            c(Valid_format_dat2, Valid_format_stat2, Valid_format_both2),
            function(plot, trelldata){
              pmartR:::data_cogs(plot, 
                                 trelldata, 
                                 try_URL = TRUE, 
                                 custom_cog = "custom_fn_cog")
              
            })

testthat::test_that("Correct data_cogs output", {
  purrr::map2(cog_list, c(Valid_plot_dat, Valid_plot_stat, Valid_plot_both), function(cog, plot){
    testthat::expect_true(nrow(cog) == nrow(plot))
    testthat::expect_true(ncol(cog) == 3)
    testthat::expect_true(colnames(cog)[2] == "panel")
    testthat::expect_true(colnames(cog)[3] == "cogs")
    testthat::expect_true(all(as.character(unlist(cog[1])) %in% as.character(unlist(plot[1]))))
    
    testthat::expect_true(inherits(cog$panel, c(
      "trelliscope_panels",
      "list")))
    testthat::expect_true(inherits(cog$cogs, c(
      "trelliscope_cogs",
      "list")))
    testthat::expect_true(all(
      unlist(purrr::map(cog$panel, function(plot) inherits(
        plot, c("gg", "plotly", "rbokeh"))))))
    testthat::expect_true(all(
      unlist(purrr::map(cog$cogs, function(cog) inherits(
        cog, c("data.frame"))))))
    return(NULL)
  })
})

################################################################################
################################################################################

testthat::context("Test list_y_limits output and error throwing")

testthat::test_that("Correct list_y_limits error/warning throwing", {
# list_y_limits errors #
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                           y_limits = list(min = 4, max = 3)),
    "Set min must be less than set max")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                           y_limits = list(min = "lp", max = 3)),
    "Min must be a numeric")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                           y_limits = list(max = "lp")),
    "Max must be a numeric")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                           y_limits = list(range = -5)),
    "Range must be a numeric")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                           y_limits = list(scale = -5)),
    "Scale must be a character")
  testthat::expect_error(pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                                                y_limits = -5),
                         "Invalid input")
  
  testthat::expect_error( pmartR:::list_y_limits(plot_type = "abundance_boxplot", 
                                                 y_limits = list(blah = list(min = 1, max = 2))),
                          "Names in y_limits should match plot_type")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "abundance_boxplot",
                           y_limits = list(bro = 4, max = 3)),
    "'min', 'max', 'range', and 'scale'")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", 
                           y_limits = list(min = 5, min = 6)),
    "List inputs may not have duplicate names.")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", 
                           y_limits = list(list(list(min = 3, max = 4)))),
    "List of lists of lists")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", 
                           y_limits = list(list(min = 4, max = 5))),
    "inputs must be named")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", 
                           y_limits = list(abundance_boxplot = list(min = 4, max = 5), 
                                           abundance_boxplot = list(min = 4, max = 5))),
    "duplicate names")
  testthat::expect_error( 
    pmartR:::list_y_limits(plot_type = "blah", 
                           y_limits = "free"),
    "Must be a character string in the following")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = list(blah = "free")),
    "List names for y_limits are restricted")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = list(abundance_boxplot = list("free"))),
    "List inputs must be named")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = "blah"),
    "Invalid input")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = list(
      scale = 'free',
      range = 5,
      min = 5)),
    "are not supported")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = list(
      max = 8,
      range = 5,
      min = 5)),
    "are not supported")
  testthat::expect_error(
    pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = list(
      max = 8,
      scale = 'free',
      min = 5)),
    "are not supported")
  
  
  testthat::expect_identical(
    names(pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = 'free')), "abundance_boxplot")
  testthat::expect_identical(
    names(pmartR:::list_y_limits(plot_type = "abundance_boxplot", y_limits = 'free')), "abundance_boxplot")
})

# ################################################################################
# ################################################################################
# 
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

# #### Test set_increment #### 
testthat::test_that("Subfunction set_increment correctly processes", {
  testthat::expect_error(pmartR:::set_increment(ylist1),
                         "numeric")
  testthat::expect_error(pmartR:::set_increment(ylist2),
                         "numeric")
  testthat::expect_error(pmartR:::set_increment(ylist3),
                         "numeric")
  testthat::expect_error(pmartR:::set_increment(ylist4, include_zero = "z"),
                         "logical")
  testthat::expect_error(pmartR:::set_increment(ylist4, include_zero = 1),
                         "logical")
  
  testthat::expect_identical(pmartR:::set_increment(ylist4), 0)
  testthat::expect_identical(pmartR:::set_increment(ylist5), 0)
  testthat::expect_identical(pmartR:::set_increment(ylist6), 0)
  testthat::expect_identical(pmartR:::set_increment(ylist7), 0)
  testthat::expect_identical(pmartR:::set_increment(ylist8), 4/20)
  testthat::expect_identical(pmartR:::set_increment(
    ylist8, include_zero = FALSE), 4/20)
  testthat::expect_identical(pmartR:::set_increment(ylist9), 4/20)
  testthat::expect_identical(pmartR:::set_increment(
    ylist9, include_zero = FALSE), 4/20)
  testthat::expect_identical(pmartR:::set_increment(ylist10), 5/20)
  testthat::expect_identical(pmartR:::set_increment(
    ylist10, include_zero = FALSE), 1/20)
  testthat::expect_identical(pmartR:::set_increment(ylist11), 6/20)
  testthat::expect_identical(pmartR:::set_increment(
    ylist11, include_zero = FALSE), 2/20)
  testthat::expect_identical(pmartR:::set_increment(ylist12), 3/20)
  testthat::expect_identical(pmartR:::set_increment(
    ylist12, include_zero = FALSE), 3/20)
  testthat::expect_identical(pmartR:::set_increment(ylist13), 3/20)
  testthat::expect_identical(pmartR:::set_increment(
    ylist13, include_zero = FALSE), 1/20)
  testthat::expect_identical(pmartR:::set_increment(ylist14), 10/20)
  testthat::expect_identical(pmartR:::set_increment(
    ylist14, include_zero = FALSE), 2/20)
})

#### Test set_ylimits ####
testthat::test_that("Subfunction set_ylimits correctly processes", {

  temp_inc <- 1    # Testing increment value
  temp_ymax <- 20  # Testing y_max
  temp_ymin <- -16 # Testing y_min
  temp_yrange <- 6 # Testing y_range

  #### Errors ####
  testthat::expect_error(pmartR:::set_ylimits(ylist1, temp_inc),
                         "numeric")
  testthat::expect_error(pmartR:::set_ylimits(ylist2, temp_inc),
                         "numeric")
  testthat::expect_error(pmartR:::set_ylimits(ylist3, temp_inc),
                         "numeric")
  testthat::expect_error(pmartR:::set_ylimits(ylist4, "z"),
                         "numeric")
  testthat::expect_error(pmartR:::set_ylimits(ylist4, NA),
                         "numeric")
  testthat::expect_error(pmartR:::set_ylimits(ylist4, c(1,2)),
                         "numeric")
  testthat::expect_error(pmartR:::set_ylimits(ylist4, temp_inc, include_zero = "z"),
                         "logical")
  testthat::expect_error(pmartR:::set_ylimits(ylist4, temp_inc, include_zero = 1),
                         "logical")
  testthat::expect_error(pmartR:::set_ylimits(ylist4, temp_inc, y_max = "z"),
                         "numeric")
  testthat::expect_error(pmartR:::set_ylimits(ylist4, temp_inc, y_min = "z"),
                         "numeric")
  testthat::expect_error(pmartR:::set_ylimits(ylist4, temp_inc, y_range = "z"),
                         "numeric")
  testthat::expect_error(pmartR:::set_ylimits(ylist4, temp_inc,
                                             y_range = temp_yrange,
                                             y_max = temp_ymax,
                                             y_min = temp_ymin),
                         "y_range must be NULL when y_max and y_min are assigned")

  #### Value match w/o limits ####
  testthat::expect_identical(pmartR:::set_ylimits(ylist4, temp_inc), c(-3, 3))
  testthat::expect_identical(pmartR:::set_ylimits(ylist5, temp_inc), c(-3, 3))
  testthat::expect_identical(pmartR:::set_ylimits(ylist6, temp_inc), c(-3, 3))
  testthat::expect_identical(pmartR:::set_ylimits(ylist7, temp_inc), c(-3, 3))
  testthat::expect_identical(pmartR:::set_ylimits(ylist8, temp_inc), c(-3, 7))
  testthat::expect_identical(pmartR:::set_ylimits(
    ylist8, temp_inc, include_zero = FALSE), c(1, 7))
  testthat::expect_identical(pmartR:::set_ylimits(ylist9, temp_inc), c(-7, 3))
  testthat::expect_identical(pmartR:::set_ylimits(
    ylist9, temp_inc, include_zero = FALSE), c(-7, -1))
  testthat::expect_identical(pmartR:::set_ylimits(ylist10, temp_inc), c(-8, 3))
  testthat::expect_identical(pmartR:::set_ylimits(
    ylist10, temp_inc, include_zero = FALSE), c(-8, -1))
  testthat::expect_identical(pmartR:::set_ylimits(ylist11, temp_inc), c(-3, 9))
  testthat::expect_identical(pmartR:::set_ylimits(
    ylist11, temp_inc, include_zero = FALSE), c(1, 9))
  testthat::expect_identical(pmartR:::set_ylimits(ylist12, temp_inc), c(-3, 6))
  testthat::expect_identical(pmartR:::set_ylimits(
    ylist12, temp_inc, include_zero = FALSE), c(0, 6))
  testthat::expect_identical(pmartR:::set_ylimits(ylist13, temp_inc), c(-1, 6))
  testthat::expect_identical(pmartR:::set_ylimits(
    ylist13, temp_inc, include_zero = FALSE), c(-1, 6))
  testthat::expect_identical(pmartR:::set_ylimits(ylist14, temp_inc), c(-13, 3))
  testthat::expect_identical(pmartR:::set_ylimits(
    ylist14, temp_inc, include_zero = FALSE), c(-13, -5))


  #### Value match w/ limits ####
  testthat::expect_identical(pmartR:::set_ylimits(
    ylist4, temp_inc, y_max = temp_ymax), c(-3, 20))
  testthat::expect_identical(pmartR:::set_ylimits(
    ylist4, temp_inc, y_min = temp_ymin), c(-16, 3))
  testthat::expect_identical(pmartR:::set_ylimits(
    ylist4, temp_inc, y_min = temp_ymin, y_max = temp_ymax), c(-16, 20))
  testthat::expect_identical(pmartR:::set_ylimits(
    ylist4, temp_inc, y_range = temp_yrange, y_max = temp_ymax), c(14, 20))
  testthat::expect_identical(pmartR:::set_ylimits(
    ylist4, temp_inc, y_min = temp_ymin, y_range = temp_yrange), c(-16, -10))

})

################################################################################
################################################################################

testthat::context("Test main TrelliVis function output (Takes )")

## Generate function outputs from test data ##

testthat::test_that("Validate error throwing", {
  testthat::expect_error(
  suppressWarnings(trelliVis(metab_object, metab_stats, custom_cog = "mean")), 
  "Validation failed")
  testthat::expect_error(
    suppressWarnings(trelliVis(list(isobaric_object, qpro1), list(isobaric_stats, qpro1_stats), custom_cog = "mean")), 
  "Validation failed")
  testthat::expect_error(
    suppressWarnings(trelliVis(metab_object, metab_stats, custom_plot = "mean")), 
    "Validation failed")
  testthat::expect_error(
    suppressWarnings(trelliVis(list(isobaric_object, qpro1), list(isobaric_stats, qpro1_stats), custom_plot = "mean")), 
    "Validation failed")
    
  testthat::expect_error(
    trelliVis(metab_object, metab_stats, panel_variable = "mean"), 
    "not present in input data")
  testthat::expect_error(
    trelliVis(list(isobaric_object, qpro1), list(isobaric_stats, qpro1_stats), panel_variable = c("mean", "mean")), 
    "not present in input data")
  testthat::expect_error(
    trelliVis(metab_object, metab_stats, state = "mean"), 
    "currently not supported")
  testthat::expect_error(
    trelliVis(omicsStats = NULL, omicsData = NULL, omicsFormat = NULL), 
    "least one of omicsData")
  testthat::expect_error(
    trelliVis(list(isobaric_object, qpro1), list(isobaric_stats, qpro1_stats), panel_variable = "Peptide"), 
    "specified for each index")
  testthat::expect_error(
    trelliVis(metab_object, metab_stats, try_URL = "blh"), 
    "try_URL must be a")

  x <- as.trellData(isobaric_object, isobaric_stats)
  attributes(x)$isobaric_info$norm_info$is_normalized <- FALSE
  attributes(x)$data_info$norm_info$is_normalized <- FALSE
  testthat::expect_error(
    trelliVis(omicsFormat = x), 
    "normalize_isobaric")
  x <- as.trellData(metab_object, metab_stats)
  attributes(x)$isobaric_info$norm_info$is_normalized <- FALSE
  attributes(x)$data_info$norm_info$is_normalized <- FALSE
  testthat::expect_error(
    trelliVis(omicsFormat = x), 
    "normalize_global")
  
  x <- as.trellData(list(isobaric_object, qpro1), list(isobaric_stats, qpro1_stats))
  attributes(x[[1]])$isobaric_info$norm_info$is_normalized <- FALSE
  attributes(x[[1]])$data_info$norm_info$is_normalized <- FALSE
  testthat::expect_error(
    trelliVis(omicsFormat = x), 
    "normalize_isobaric")
  x <- as.trellData(list(pep_object, qpro2), list(pep_stats, qpro2_stats))
  attributes(x[[2]])$isobaric_info$norm_info$is_normalized <- FALSE
  attributes(x[[2]])$data_info$norm_info$is_normalized <- FALSE
  testthat::expect_error(
    trelliVis(omicsFormat = x), 
    "normalize_global")
  
  x <- as.trellData(list(pep_object, qpro2), list(pep_stats, qpro2_stats))
  attributes(x[[1]])$isobaric_info$norm_info$is_normalized <- FALSE
  attributes(x[[1]])$data_info$norm_info$is_normalized <- FALSE
  testthat::expect_error(
    trelliVis(omicsFormat = x), 
    "normalize_global")
  x <- as.trellData(list(pep_object, qpro2), list(pep_stats, qpro2_stats))
  attributes(x[[2]])$isobaric_info$norm_info$is_normalized <- FALSE
  attributes(x[[2]])$data_info$norm_info$is_normalized <- FALSE
  testthat::expect_error(
    trelliVis(omicsFormat = x), 
    "normalize_global")
  
  # These warnings aren't being picked up for some reason??
  # testthat::expect_warning(testthat::expect_warning(
  #   x <- trelliVis(omicsData = metab_object, trelli_name = "la:la"), 
  #   "Non-word characters"))
  # testthat::expect_warning(
  #   x <- trelliVis(omicsData = metab_object, trelli_path_out = "la:la"), 
  #   "Non-word characters")
  
  testthat::expect_message(
    trelliVis(list(isobaric_object, qpro1), trelli_name = c("this", "that", "theother"), display = FALSE), 
    "Length of trelli_name")
})


## check all the options for correct output (For full checks, set display = TRUE)
iso_display    <- trelliVis(isobaric_object, isobaric_stats, custom_plot = "custom_fn", display = FALSE)
metab_display  <- trelliVis(metab_object, metab_stats, self_contained = TRUE, display = FALSE)
metab_display2  <- trelliVis(metab_stats, display = FALSE)
lipid_display  <- trelliVis(list(lipid_object), list(lipid_stats), display = FALSE)
pro_display    <- trelliVis(qpro1_stats, qpro1, display = FALSE)
linked_display <- trelliVis(list(isobaric_object, qpro1),
                            list(isobaric_stats, qpro1_stats),
                            custom_plot = "custom_fn", display = FALSE)
# linked_display2 <- trelliVis(list(isobaric_object, qpro1), trelli_name = c("this", "that", "the other"))
linked_display3 <- trelliVis(list(isobaric_object, qpro1),
                             custom_plot = "custom_fn",
                             trelli_name = c("Truth", "Happiness"),
                             self_contained = TRUE, display = FALSE)
linked_display4 <- trelliVis(list(isobaric_object, qpro1),
                             trelli_name = "blah",
                             plot_type = "abundance_boxplot",
                             custom_plot = "custom_fn", 
                             display = FALSE, 
                             panel_variable = c("Protein", "Protein"))

##

metab_display  <- trelliVis(metab_object, trelli_name = "test")
lipid_display  <- trelliVis(list(lipid_object), list(lipid_stats))

# devtools::test_coverage_file() => 
