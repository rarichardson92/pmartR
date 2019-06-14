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

################################################################################

## Generate function outputs from test data ##

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

#Both                                                                                                                   qpro4))

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
  format_metobject, format_metobject2,
  format_pepobject, format_pepobject2,
  format_proobject, format_tecobject,
  format_qpro4
)

Valid_format_stat <- list(
  format_metstats 
  # format_pepstats, format_prostats, 
  # format_tecstats, format_qpro4stats
)

Valid_format_both <- list(
  format_metboth2
  # format_pepboth2, format_proboth,
  # format_tecboth
)

Valid_forlist_dat <- list(
  forlist_qpro4dat
)

Valid_forlist_stat <- list(
  # forlist_qpro4stat
)

Valid_forlist_both <- list(
  # forlist_qpro4both
)

################################################################################

## Test format_data outputs ##
# Dimensions and correct columns in output #

testthat::test_that("Correct dataframe population format data", {
                      
  purrr::map(Valid_format_dat, function(formDat){
    testthat::expect_length(formDat, 3)
    testthat::expect_null(formDat$comp_stats)
    testthat::expect_null(formDat$summary_stats)
    testthat::expect_false(is.null(formDat$data_value))
    testthat::expect_gt(nrow(formDat$data_value), 0)
  })

  purrr::map(Valid_format_stat, function(formStat){
    testthat::expect_length(formStat, 3)
    testthat::expect_null(formStat$data_value)
    testthat::expect_false(is.null(formStat$summary_stats))
    testthat::expect_false(is.null(formStat$comp_stats))
    testthat::expect_gt(nrow(formStat$summary_stats), 0)
    testthat::expect_gt(nrow(formStat$comp_stats), 0)
  })

  purrr::map(Valid_format_both, function(formBoth){
    testthat::expect_length(formBoth, 3)
    testthat::expect_false(is.null(formBoth$data_value))
    testthat::expect_false(is.null(formBoth$summary_stats))
    testthat::expect_false(is.null(formBoth$comp_stats))
    testthat::expect_gt(nrow(formBoth$data_value), 0)
    testthat::expect_gt(nrow(formBoth$summary_stats), 0)
    testthat::expect_gt(nrow(formBoth$comp_stats), 0)
  })
  
  purrr::map(Valid_forlist_dat, function(index){
    testthat::expect_length(index, 2)
    purrr::map(index, function(formDat){
      testthat::expect_length(formDat, 3)
      testthat::expect_null(formDat$comp_stats)
      testthat::expect_null(formDat$summary_stats)
      testthat::expect_false(is.null(formDat$data_value))
      testthat::expect_gt(nrow(formDat$data_value), 0)
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

################################################################################

## Test format_data outputs ##
# Testing data specific expectations #

### Thus farrr with errors :(
# Is there a good way to estimate expected row number mathmatically?

# data_values:   cnames, abundance value, f_data colnames, group_DF group, e_meta colnames of associated omicsData
# comp stats:   e_data cname, "comparison", respective P values (T and G), "Fold_change", "Flag"
# summary stats : e_data cname, "Group" or "Group_DF", "Count", "Mean"

testthat::test_that("Correct columns and rows in format_data", {
  
  # Data-only generated ############
  testthat::expect_match(toString(colnames(format_isoobject$data_values)), 
                         "abundance")
  testthat::expect_match(toString(colnames(format_isoobject$data_values)), 
                         "Group_DF")
  testthat::expect_match(toString(colnames(format_isoobject$data_values)),
                        attr(format_isoobject, "cname")$edata_cname)
  testthat::expect_true(all(colnames(isobaric_object$f_data) %in%
                              colnames(format_isoobject$data_values)))
  testthat::expect_true(all(colnames(isobaric_object$e_meta) %in%
                              colnames(format_isoobject$data_values)))
  
  testthat::expect_match(toString(colnames(format_lipobject$data_values)), 
                         "abundance")
  testthat::expect_match(toString(colnames(format_lipobject$data_values)), 
                         "Group")
  testthat::expect_match(toString(colnames(format_lipobject$data_values)),
                         attr(format_lipobject, "cname")$edata_cname)
  testthat::expect_true(all(colnames(lipid_object$f_data) %in%
                              colnames(format_lipobject$data_values)))
  
  testthat::expect_match(toString(colnames(format_lipobject2$data_values)), 
                         "abundance")
  testthat::expect_match(toString(colnames(format_lipobject2$data_values)), 
                         "Group")
  testthat::expect_match(toString(colnames(format_lipobject2$data_values)),
                         attr(format_lipobject2, "cname")$edata_cname)
  testthat::expect_true(all(colnames(lipid_object2$f_data) %in%
                              colnames(format_lipobject2$data_values)))
  
  testthat::expect_match(toString(colnames(format_metobject$data_values)), 
                         "abundance")
  testthat::expect_match(toString(colnames(format_metobject$data_values)), 
                         "Group")
  testthat::expect_match(toString(colnames(format_metobject$data_values)),
                         attr(format_metobject, "cname")$edata_cname)
  testthat::expect_true(all(colnames(metab_object$f_data) %in%
                              colnames(format_metobject$data_values)))
  
  testthat::expect_match(toString(colnames(format_metobject2$data_values)), 
                         "abundance")
  testthat::expect_match(toString(colnames(format_metobject2$data_values)), 
                         "Group")
  testthat::expect_match(toString(colnames(format_metobject2$data_values)),
                         attr(format_metobject2, "cname")$edata_cname)
  testthat::expect_true(all(colnames(metab_object2$f_data) %in%
                              colnames(format_metobject2$data_values)))

  testthat::expect_match(toString(colnames(format_pepobject$data_values)), 
                         "abundance")
  testthat::expect_match(toString(colnames(format_pepobject$data_values)), 
                         "Group")
  testthat::expect_match(toString(colnames(format_pepobject$data_values)),
                         attr(format_pepobject, "cname")$edata_cname)
  testthat::expect_true(all(colnames(pep_object$f_data) %in%
                              colnames(format_pepobject$data_values)))
  testthat::expect_true(all(colnames(pep_object$e_meta) %in%
                              colnames(format_pepobject$data_values)))
  
  testthat::expect_match(toString(colnames(format_pepobject2$data_values)), 
                         "abundance")
  testthat::expect_match(toString(colnames(format_pepobject2$data_values)), 
                         "Group")
  testthat::expect_match(toString(colnames(format_pepobject2$data_values)),
                         attr(format_pepobject2, "cname")$edata_cname)
  testthat::expect_true(all(colnames(pep_object2$f_data) %in%
                              colnames(format_pepobject2$data_values)))
  testthat::expect_true(all(colnames(pep_object2$e_meta) %in%
                              colnames(format_pepobject2$data_values)))
  
  testthat::expect_match(toString(colnames(format_proobject$data_values)), 
                         "abundance")
  testthat::expect_match(toString(colnames(format_proobject$data_values)), 
                         "Group")
  testthat::expect_match(toString(colnames(format_proobject$data_values)),
                         attr(format_proobject, "cname")$edata_cname)
  testthat::expect_true(all(colnames(pro_object$f_data) %in%
                              colnames(format_proobject$data_values)))
  testthat::expect_true(all(colnames(pro_object$e_meta) %in%
                              colnames(format_proobject$data_values)))
  
  testthat::expect_match(toString(colnames(format_tecobject$data_values)), 
                         "abundance")
  testthat::expect_match(toString(colnames(format_tecobject$data_values)), 
                         "Group")
  testthat::expect_match(toString(colnames(format_tecobject$data_values)),
                         attr(format_tecobject, "cname")$edata_cname)
  testthat::expect_true(all(colnames(techrep_pep_object$f_data) %in%
                              colnames(format_tecobject$data_values)))
  testthat::expect_true(all(colnames(techrep_pep_object$e_meta) %in%
                              colnames(format_tecobject$data_values)))
  
  testthat::expect_match(toString(colnames(format_qpro4$data_values)), 
                         "abundance")
  testthat::expect_match(toString(colnames(format_qpro4$data_values)), 
                         "Group")
  testthat::expect_match(toString(colnames(format_qpro4$data_values)),
                         attr(format_qpro4, "cname")$edata_cname)
  testthat::expect_true(all(colnames(qpro4$f_data) %in%
                              colnames(format_qpro4$data_values)))
  testthat::expect_true(all(colnames(qpro4$e_meta) %in%
                              colnames(format_qpro4$data_values)))
  
  # Stats-only generated ##########
  
  testthat::expect_match(toString(colnames(format_metstats$comp_stats)), 
                         "Comparison")
  testthat::expect_match(toString(colnames(format_metstats$comp_stats)), 
                         "Fold_change")
  testthat::expect_match(toString(colnames(format_metstats$comp_stats)), 
                         "Flag")
  testthat::expect_match(toString(colnames(format_metstats$comp_stats)), 
                         "P_value_G")
  testthat::expect_match(toString(colnames(format_metstats$comp_stats)), 
                         "P_value_T")
  testthat::expect_match(toString(colnames(format_metstats$comp_stats)),
                         attr(metab_stats, "cname")$edata_cname)
  
  testthat::expect_match(toString(colnames(format_metstats$summary_stats)), 
                         "Count")
  testthat::expect_match(toString(colnames(format_metstats$summary_stats)), 
                         "Mean")
  testthat::expect_match(toString(colnames(format_metstats$summary_stats)), 
                         "Group")
  testthat::expect_match(toString(colnames(format_metstats$summary_stats)),
                         attr(metab_stats, "cname")$edata_cname)
  
  # Both Stats and Data generated ########
  
  testthat::expect_match(toString(colnames(format_metboth2$data_values)), 
                         "abundance")
  testthat::expect_match(toString(colnames(format_metboth2$data_values)), 
                         "Group")
  testthat::expect_match(toString(colnames(format_metboth2$data_values)),
                         attr(metab_object2, "cname")$edata_cname)
  testthat::expect_true(all(colnames(metab_object2$f_data) %in%
                              colnames(format_metboth2$data_values)))
  testthat::expect_true(all(colnames(metab_object2$e_meta) %in%
                              colnames(format_metboth2$data_values)))
  
  testthat::expect_match(toString(colnames(format_metboth2$comp_stats)), 
                         "Comparison")
  testthat::expect_match(toString(colnames(format_metboth2$comp_stats)), 
                         "Fold_change")
  testthat::expect_match(toString(colnames(format_metboth2$comp_stats)), 
                         "Flag")
  testthat::expect_match(toString(colnames(format_metboth2$comp_stats)), 
                         "P_value_G")
  testthat::expect_match(toString(colnames(format_metboth2$comp_stats)), 
                         "P_value_T")
  testthat::expect_match(toString(colnames(format_metboth2$comp_stats)),
                         attr(metab_stats, "cname")$edata_cname)
  
  testthat::expect_match(toString(colnames(format_metboth2$summary_stats)), 
                         "Count")
  testthat::expect_match(toString(colnames(format_metboth2$summary_stats)), 
                         "Mean")
  testthat::expect_match(toString(colnames(format_metboth2$summary_stats)), 
                         "Group")
  testthat::expect_match(toString(colnames(format_metboth2$summary_stats)),
                         attr(metab_stats, "cname")$edata_cname)
  
  # List of data generated ##########

  testthat::expect_match(toString(colnames(forlist_qpro4dat[[1]]$data_values)), 
                         "abundance")
  testthat::expect_match(toString(colnames(forlist_qpro4dat[[1]]$data_values)), 
                         "Group")
  testthat::expect_match(toString(colnames(forlist_qpro4dat[[1]]$data_values)),
                         attr(pep_object2, "cname")$edata_cname)
  testthat::expect_true(all(colnames(pep_object2$f_data) %in%
                              colnames(forlist_qpro4dat[[1]]$data_values)))
  testthat::expect_true(all(colnames(pep_object2$e_meta) %in%
                              colnames(forlist_qpro4dat[[1]]$data_values)))
  
  testthat::expect_match(toString(colnames(forlist_qpro4dat[[2]]$data_values)), 
                         "abundance")
  testthat::expect_match(toString(colnames(forlist_qpro4dat[[2]]$data_values)), 
                         "Group")
  testthat::expect_match(toString(colnames(forlist_qpro4dat[[2]]$data_values)),
                         attr(qpro4, "cname")$edata_cname)
  testthat::expect_true(all(colnames(qpro4$f_data) %in%
                              colnames(forlist_qpro4dat[[2]]$data_values)))
  testthat::expect_true(all(colnames(qpro4$e_meta) %in%
                              colnames(forlist_qpro4dat[[2]]$data_values)))
  
})

################################################################################
################################################################################

testthat::context("Test format_plot output")

################################################################################
################################################################################

testthat::context("Test data_cogs output")

################################################################################
################################################################################

testthat::context("Test recursive_format, set_increment, and set_ylimits output")

testthat::test_that("Subfunction recursive_format correctly processes", {
  
})

testthat::test_that("Subfunction set_increment correctly processes", {
  
})

testthat::test_that("Subfunction set_ylimits correctly processes", {
  
})

################################################################################
################################################################################

testthat::context("Test main function output")

