#' @name as.trellData
#' @rdname as.trellData
#' @title Convert Omics data and pairwise statistics to a plotting object
#'
#' @description Converts a ResObject and its respective OmicsData into an easily plottable object.
#'
#' @param omicsData A pmartR object of class pepData, lipidData, metabData, or proData
#' @param omicsStats A statistical results object produced by running \code{imd_anova} on omicsData.
#' @param ... further arguments
#'
#' @details Objects of class 'trellData' inherit attributes from omicsStats and omicsData. These attributes can be changed from their default value by manual specification. A list of these attributes as well as their default values are as follows:
#' \tabular{ll}{
#' data_scale \tab Scale of the data provided in \code{e_data}. Acceptable values are 'log2', 'log10', 'log', and 'abundance', which indicate data is log base 2, base 10, natural log transformed, and raw abundance, respectively. Default values is 'abundance'. \cr
#' \tab \cr
#' is_normalized \tab A logical argument, specifying whether the data has been normalized or not. Default value is FALSE. \cr
#' \tab \cr
#' norm_info \tab Default value is an empty list, which will be populated with a single named element \code{is_normalized = is_normalized}. When a normalization is applied to the data, this becomes populated with a list containing the normalization function, normalization subset and subset parameters, the location and scale parameters used to normalize the data, and the location and scale parameters used to backtransform the data (if applicable). \cr
#' \tab \cr
#' data_types \tab Character string describing the type of data (e.g.'Positive ion'). Default value is NULL. \cr
#' \tab\cr
#' check.names \tab Logical defaults to TRUE. Indicates whether 'check.names' attribute of returned omicsData object is TRUE or FALSE. \cr
#' }
#' Computed values included in the \code{data_info} attribute are as follows:
#' \tabular{ll}{
#' num_edata \tab The number of unique \code{edata_cname} entries.\cr
#' \tab \cr
#' num_miss_obs \tab The number of missing observations.\cr
#' \tab \cr
#' num_emeta \tab The number of unique \code{emeta_cname} entries. \cr
#' \tab \cr
#' prop_missing \tab The proportion of \code{e_data} values that are NA. \cr
#' \tab \cr
#' num_samps \tab The number of samples that make up the columns of \code{e_data}.\cr
#' \tab \cr
#' meta_info \tab A logical argument, specifying whether \code{e_meta} is provided.\cr
#' \tab \cr
#' }
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data("pro_edata")
#' data("pro_fdata")
#' myproData <- as.proData(e_data = pro_edata, f_data = pro_fdata, edata_cname = "Reference", fdata_cname = "SampleID", is_normalized=TRUE)
#'}
#'
#' @author Rachel Richardson
#'
#' @export
#' 
as.trellData <- function(...){
  .as.trellData(...)
}

.as.trellData <- suppressWarnings(function(omicsData = NULL, omicsStats = NULL){
  
  ## Checks and recursive for lists as input in as.trellData ##
  if ((class(omicsData) == "list") || (class(omicsStats) == "list")){
    return(recursive_format(omicsData = omicsData, omicsStats = omicsStats))
  }
  
  ## Switch Stats and Omics data as appropriate ##
  if (is.null(omicsStats) & inherits(omicsData, 'statRes')){
    omicsStats <- omicsData
    omicsData <- NULL
  }
  
  if (is.null(omicsData) && inherits(omicsStats, 
                                     c("proData", "pepData",
                                       "metabData", "lipidData"))){
    warning("Input reordered: omicsStats inherits omicsData attributes, input object reassigned to omicsData.")
    omicsData <- omicsStats
    omicsStats <- NULL
  }
  
  if (!is.null(omicsStats) && 
      !is.null(omicsData) && 
      inherits(omicsStats, c("proData", "pepData", 
                             "metabData", "lipidData")) &&
      inherits(omicsData, 'statRes')){
    temp <- omicsStats
    omicsStats <- omicsData
    omicsData <- temp
    rm(temp)
  }
  
  ## Validate input ##
  validate_omics_input(omicsStats = omicsStats, omicsData = omicsData)
  
  ## Moniker Variables ##
  uniqedatId <- attributes(omicsData)$cnames$edata_cname
  if(is.null(uniqedatId)){
    uniqedatId <- attributes(omicsStats)$cnames$edata_cname
  }
  sampID <- attributes(omicsData)$cnames$fdata_cname
  if(is.null(sampID)){
    sampID <- attributes(omicsStats)$cnames$fdata_cname
  }
  stats <- omicsStats$Full_results
 
  ## Formatting ##
  # --omicsData--#
  if (!is.null(omicsData)){
    ## Manipulate Dataframes ##
    # Original data values from omicsData #
    # 1) Melt to assign e_data values to each sample ID
    # 2) Combine with groups in f_data
    # 3) Save attributes to be used
    if (get_data_scale(omicsData) != "abundance"){
      value_name <- paste(get_data_scale(omicsData), "_abundance", sep = "")
    } else {
      value_name <- "abundance"
    }
    
    data_values <- suppressWarnings(
      reshape2::melt(omicsData$e_data, 
                     id.vars = uniqedatId,
                     value.name = value_name,
                     variable.name = sampID) %>%
        dplyr::left_join(omicsData$f_data, by = sampID))
    
    # Adjust for class differences when read in for group_DF and data values #
    joingroupDF <- attr(omicsData, "group_DF")
    fixclass <- data_values[colnames(data_values) %in% 
                              colnames(joingroupDF)]
    data_values[colnames(data_values) %in% 
                  colnames(joingroupDF)] <- apply(fixclass,
                                                  2,
                                                  as.character)
    
    # Adjust for "Group" in original data columns #
    if("Group" %in% colnames(fixclass)){
      joingroupDF$Group_DF <- joingroupDF$Group
      joingroupDF$Group <- joingroupDF[["Group.1"]]
      joingroupDF[["Group.1"]] <- NULL
    }
    
    # Join with group_DF
    data_values <- suppressWarnings(dplyr::left_join(data_values, joingroupDF, by = sampID))
    
    # Join with e_meta if present #
    if(!is.null(omicsData$e_meta)) {
      data_values <- suppressWarnings(dplyr::left_join(data_values, 
                               omicsData$e_meta, by = uniqedatId))
    }
  
    # Inherit attributes of omicsData
    inheritDatAtt <- attributes(omicsData)[grep("filters|^names", 
                                                names(attributes(omicsData)), 
                                                invert = T)]
    } else {
      data_values <- NULL
      inheritDatAtt <- NULL
  }
  
  if (!is.null(omicsStats)){
    # Comparison statistics by pairwise comparison from omicsStats #
    # 1) For each comparison, extract comparison relevant columns from omicsStats
    #  2) Bind extracted with unique IDs from omicsStats and a comparison column
    #  3) Rename columns by removing comparison designation
    # 4) Bind rows from each comparison into one dataframe 
    comp_stats <- suppressWarnings(
      dplyr::bind_rows(
        attr(omicsStats, "comparisons") %>%
          purrr::map(function(paircomp) {
            df <- cbind(
              stats[[uniqedatId]], 
              rep(paircomp, nrow(stats)),
              stats[stringr::str_detect(names(stats),
                                        pattern = paste(paircomp, 
                                                        "$", 
                                                        sep = ""))])
            trimname <- stringr::str_remove_all(colnames(df)[3:ncol(df)],
                                                paste("_", paircomp, sep = ""))
            colnames(df) <- c(uniqedatId, "Comparison", trimname)
            return(df)
            }
            )))
    
    # Comparison statistics by groups from omicsStats #
    # 1) Extract all rows w/o comparison information in the column names
    summary_stats <- stats[,!stringr::str_detect(
      names(stats),
      paste(attr(omicsStats, "comparisons"), collapse = '|')
      )]
    
    # Re-sort for any grouping variables #
    # 1) If there is a group_DF specified,
    # 2) Extract unique groups in group_DF, and for each group
    #  3) Bind extracted with unique IDs from omicsStats and a Group column
    #  4) Rename columns by removing group designation
    # 5) Bind rows from each group into one dataframe 
    if (!is.null(attr(omicsStats, "group_DF"))){
      groups <-  unique(attr(omicsStats, "group_DF")$Group)
      summary_stats <- suppressWarnings(
        dplyr::bind_rows(purrr::map(groups, function(group){
          df <- cbind(summary_stats[[uniqedatId]],
                      rep(group, nrow(summary_stats)),
                      summary_stats[,stringr::str_detect(
                        names(summary_stats), 
                        paste(group, "$", sep = ""))])
          trimname <- stringr::str_remove_all(colnames(df)[3:ncol(df)],
                                              paste("_", group, sep = ""))
          if(exists("fixclass") && "Group" %in% colnames(fixclass)){
            colnames(df) <- c(uniqedatId, "Group_DF", trimname)
          } else {
            colnames(df) <- c(uniqedatId, "Group", trimname)
          }
          return(df)
          }
        )))
    }
    
    # Join with e_meta if present #
    if(!is.null(omicsData$e_meta)) {
      comp_stats <- suppressWarnings(dplyr::left_join(comp_stats,
                              omicsData$e_meta, by = uniqedatId))
      summary_stats <- suppressWarnings(dplyr::left_join(summary_stats,
                              omicsData$e_meta, by = uniqedatId))
    }
    
    # Save attributes to be used #
    inheritStatAtt <- attributes(omicsStats)[grep("filters|^names", 
                                                  names(attributes(omicsStats)), 
                                                  invert = T)]
  } else {
    comp_stats <- NULL
    summary_stats <- NULL
    inheritStatAtt <- NULL
  }
  
  ## Store Results ##
  # Create object res from generated dataframes and inherited attributes
  # 1) Assign dataframes to res list
  # 2) Define inherited attributes from omicsData and omicsStats
  # 3) Assign attributes to res
  # 4) Assign class to res
  res <- list(data_values = data_values, 
              summary_stats = summary_stats,
              comp_stats = comp_stats)
  
  attributes(res) <- c(attributes(res), inheritDatAtt, inheritStatAtt)
  class(res) <- "trellData" 
  
  if (!is.null(omicsData) & !is.null(omicsStats)){
    attr(res, "parent_class") <- c(class(omicsStats), class(omicsData))
    attr(res, "data_class") <- class(omicsData)
  } else if (!is.null(omicsData)){
    attr(res, "parent_class") <- class(omicsData)
    attr(res, "data_class") <- class(omicsData)
  } else {
    attr(res, "parent_class") <- class(omicsStats)
    attr(res, "data_class") <- attr(omicsStats, "data_class")
  }

  return(res)
})


#' print.trellData
#' 
#' For printing an S3 object of type 'trellData':
#' 
#'@rdname print-trellData
#'@export
#'
print.trellData<- function(trellData){
  if(!inherits(trellData, "trellData")) stop("trellData object must be of the class 'trellData'")
  
  data_values <- as.data.frame(lapply(trellData$data_values, as.character), stringsAsFactors = FALSE, check.names = attr(trellData, "check.names"))
  summary_stats <- as.data.frame(lapply(trellData$summary_stats, as.character), stringsAsFactors = FALSE, check.names = attr(trellData, "check.names"))
  comp_stats <- as.data.frame(lapply(trellData$comp_stats, as.character), stringsAsFactors = FALSE, check.names = attr(trellData, "check.names"))
  data_values_ncols <- ncol(data_values)
  summary_stats_ncols <- ncol(summary_stats)
  comp_stats_ncols <- ncol(comp_stats)
  blank_row = rep("----", 5)
    
  ## Truncate visualization for many rows, columns ##
  if(nrow(data_values) >= 9){
    data_values_head = head(data_values, 4)[, 1:min(data_values_ncols, 5)]
    data_values_tail = tail(data_values, 4)[, 1:min(data_values_ncols, 5)]
    data_values = rbind(data_values_head, blank_row, data_values_tail)
  }else if (nrow(data_values) > 0) {
    data_values = data_values[, 1:min(data_values_ncols, 5)]
  }
  if(nrow(summary_stats) >= 9){
    summary_stats_head = head(summary_stats, 4)[, 1:min(summary_stats_ncols, 5)]
    summary_stats_tail = tail(summary_stats, 4)[, 1:min(summary_stats_ncols, 5)]
    summary_stats = rbind(summary_stats_head, blank_row, summary_stats_tail)
  }else if (nrow(summary_stats) > 0){
    summary_stats = summary_stats[, 1:min(summary_stats_ncols, 5)]
  }
  if(nrow(comp_stats) >= 9){
    comp_stats_head = head(comp_stats, 4)[, 1:min(comp_stats_ncols, 5)]
    comp_stats_tail = tail(comp_stats, 4)[, 1:min(comp_stats_ncols, 5)]
    comp_stats = rbind(comp_stats_head, blank_row, comp_stats_tail)
  }else if (nrow(comp_stats) > 0){
    comp_stats = comp_stats[, 1:min(comp_stats_ncols, 5)]
  }
  
  ## Print applicable results ##  
  if(nrow(data_values) > 0){
    if(data_values_ncols > 5) message("only first 5 columns are shown")
    cat("data_values\n")
    cat(capture.output(data_values), sep = "\n")
    cat("\n")
  }
  if(nrow(summary_stats) > 0){
    if(summary_stats_ncols > 5) message("only first 5 columns are shown")
    cat("summary_stats\n")
    cat(capture.output(summary_stats), sep = "\n")
    cat("\n")
  }
  if(nrow(comp_stats) > 0){
    if(comp_stats_ncols > 5) message("only first 5 columns are shown")
    cat("comp_stats\n")
    cat(capture.output(comp_stats), sep = "\n")
    cat("\n")
  }
}

#' @name recursive_format
#' @rdname recursive_format
#' @title Recursive call of as.trellData and associated checks for list inputs
#' 
#' @description Checks for validity of list inputs and handles different list combinations as input.
#'
#' @param yvalues y-values for plotting
#' @param testtype consideration for different statistical tests
#'
#' @author Rachel Richardson
#'

recursive_format <- function(...){
  .recursive_format(...)
}

.recursive_format <- function(omicsData = NULL, omicsStats = NULL){

  #--Both--#
  if (!is.null(omicsData) & !is.null(omicsStats)) {
    if(length(omicsData) != 1 & length(omicsStats) != 1) {

      validate_omics_input(omicsData = omicsData, omicsStats = omicsStats)
      classlist <- omicsData %>% purrr::map(function(omics) attr(omics, which = "class"))
      plotterlist <- purrr::map2(omicsData, omicsStats, as.trellData)
      
    } else {
      plotterlist <- as.trellData(omicsData[[1]], omicsStats[[1]])
      return(plotterlist)
    }
    
    #--omicsData--#
  } else if(!is.null(omicsData)){
    if (length(omicsData) != 1){
      
      validate_omics_input(omicsData = omicsData, omicsStats = omicsStats)
      classlist <- omicsData %>% purrr::map(function(omics) attr(omics, which = "class"))
      plotterlist <- purrr::map(omicsData, as.trellData)
      
    } else {
      plotterlist <- as.trellData(omicsData[[1]])
      return(plotterlist)
    }
    
    #--omisStats--#
  } else {
    if (length(omicsStats) > 1){
      
      validate_omics_input(omicsData = omicsData, omicsStats = omicsStats)
      classlist <- omicsStats %>% purrr::map(function(omics) attr(omics, which = "data_class"))
      plotterlist <- purrr::map(omicsStats, as.trellData)
      
    } else {
      plotterlist <- as.trellData(omicsStats[[1]])
      return(plotterlist)
    }
  }
  attr(plotterlist, "data_types") <- classlist
  return(plotterlist)
}

#' @name validate_omics_input
#' @rdname validate_omics_input
#' @title Validate inputs for omicsData and omicsStats in trelliVis processing
#' 
#' @description Checks for validity of omicsData and omicsStats inputs.
#'
#' @param omicsData A pmartR object of class pepData, lipidData, metabData, or proData
#' @param omicsStats A statistical results object produced by running \code{imd_anova} on omicsData.
#'
#' @author Rachel Richardson
#'

validate_omics_input <- function(...){
  .validate_omics_input(...)
}

.validate_omics_input <- function(omicsData = NULL, omicsStats = NULL){
  
  
  ## Checks and recursive for lists as input in as.trellData ##
  if ((class(omicsData) == "list") || (class(omicsStats) == "list")){
    
    # Check for empty lists #
    if(((length(omicsData) == 0) && (class(omicsData) == "list")) || 
       ((length(omicsStats) == 0) && (class(omicsStats) == "list"))) stop(
         "Empty lists are not supported in omicsData or omicsStats input."
       )
    
    #--Both--#
    if (!is.null(omicsData) & !is.null(omicsStats)) {
      if(length(omicsData) != 1 & length(omicsStats) != 1) {
        # Check that the list inputs are equal in length (1 and 1, or 2 and 2) #
        if (length(omicsData) != length(omicsStats)) stop(
          "List length does not match; lists for omicsData amd omicsStats should contain
          the same protien set(s) and peptide set(s).")
        
        # Check that omicsData input is a list() not a c() of omics data #
        if (any(unlist(purrr::map(c(omicsData), is.data.frame)))) stop(
          "List/vector entry error, dataframes in input list. Possible solutions: 
          1) use list() instead of c() to preverse integrity of omicsData 
          and omicsStats lists, 2) If using both omicsData and omicsStats, 
          both inputs must be lists.")
        
        # Check that omicsStats input is a list() not a c() #
        if (any(unlist(purrr::map(c(omicsStats), is.data.frame)))) stop(
          "List/vector entry error, dataframes in input list. Possible solutions: 
            1) use list() instead of c() to preverse integrity of omicsData 
            and omicsStats lists, 2) If using both omicsData and omicsStats, 
            both inputs must be lists.")
        
        #Generate cname lists #
        listData <- omicsData %>% purrr::map(function(omics) attr(omics, which = "cname"))
        listStats <- omicsStats %>% purrr::map(function(omics) attr(omics, which = "cname"))
        
        # Check if cnames match between lists #
        if (!identical(listData, listStats)) stop(
          "Lists in omicsData and omicsStats have mismatched cname attributes. 
          Order matters; correct the associated data and stats to have the same 
          index in each list.")
        
        # Generate class list and check for pro/pep classes #
        classlist <- omicsData %>% purrr::map(function(omics) attr(omics, which = "class"))
        classlgl <- classlist %>% 
          purrr::map_lgl(function(class) all(stringr::str_detect(class, 
                                                                 pattern = "pepData|proData")))
        if(!all(classlgl)) stop(
          "Only pepData and proData are valid for lists in omicsData. 
          For omicsStats, please use omics stats argument in this format: 
          omicsStats = list(omicsStats1, omicsStats2)")
        if(any(duplicated(classlist))) stop(
          "Only one pepData, one proData and their respective statistics supported in omicsData and omicsStats lists."
        )
        
      } 
      
      #--omicsData--#
    } else if(!is.null(omicsData)){
      if (length(omicsData) != 1){

        # Check that input is a list() not a c() of omics data #
        if (any(unlist(purrr::map(c(omicsData), is.data.frame)))) stop(
          "List/vector entry error, dataframes in input list. Possible solutions: 
      1) use list() instead of c() to preverse integrity of omicsData 
      and omicsStats lists, 2) If using both omicsData and omicsStats, 
        both inputs must be lists.")
        
        # Check that list is length 2 #
        if (length(omicsData) != 2) stop(
          "List length != 2; list for omicsData should contain one 
      protien set and one peptide set.")
        
        # Generate class list and check for pro/pep classes #
        classlist <- omicsData %>% purrr::map(function(omics) attr(omics, which = "class"))
        classlgl <- classlist %>% 
          purrr::map_lgl(function(class) all(stringr::str_detect(class, 
                                                                 pattern = "pepData|proData")))
        if(!all(classlgl)) stop(
          "Only pepData and proData are valid for lists in omicsData.  
          For omicsStats, please use omics stats argument in this format: 
          omicsStats = list(omicsStats1, omicsStats2)")
        if(any(duplicated(classlist))) stop(
          "Only one pepData and one proData supported in omicsData lists."
        )

      } 

      #--omisStats--#
    } else {
      if (length(omicsStats) > 1){
        
        # Check that input is a list() not a c() #
        if (any(unlist(purrr::map(c(omicsStats), is.data.frame)))) stop(
          "List/vector entry error, dataframes in input list. Possible solutions: 
        1) use list() instead of c() to preverse integrity of omicsData 
        and omicsStats lists, 2) If using both omicsData and omicsStats, 
          both inputs must be lists.")
        
        # Check that list is length 2 #
        if (length(omicsStats) != 2 ) stop(
          "List length != 2; list for omicsStats should contain one 
        protien set and one peptide set.")
        
        classlist <- omicsStats %>% purrr::map(function(omics) attr(omics, which = "data_class"))
        classlgl <- classlist %>% 
          purrr::map_lgl(function(class) all(stringr::str_detect(class, 
                                                                 pattern = "pepData|proData")))
        if(!all(classlgl)) stop(
          "Only stats derived from pepData and proData are valid for lists in omicsStats.")
        if(any(duplicated(classlist))) stop(
          "Only one stats object derived from pepData and one stats object derived from proData supported in omicsStats lists.")   
        
      }
    }
  } else {
    
    ## Initial Checks  for non-lists ##
    # Make sure at least one of omicsData or omicsStats is present #
    if(is.null(omicsStats) && is.null(omicsData)) stop(
      "as.trellData() requires at least one of the following: 
      omicsStats, omicsData")
    
    # Check that omicsData and omicsStats are the correct classes #
    if (!is.null(omicsData) && 
        !inherits(omicsData, c("proData", 
                               "pepData", 
                               "metabData", 
                               "lipidData"))) stop(
                                 "omicsData must be of class 'proData', 
                                 'pepData', 'metabData', or 'lipidData'")
    if(!is.null(omicsStats) && !inherits(omicsStats, "statRes")) stop(
      "omicsStats must be of the class 'statRes'")
    
    
    # Moniker Variables #
    uniqedatId <- attributes(omicsData)$cnames$edata_cname
    if(is.null(uniqedatId)){
      uniqedatId <- attributes(omicsStats)$cnames$edata_cname
    }
    sampID <- attributes(omicsData)$cnames$fdata_cname
    if(is.null(sampID)){
      sampID <- attributes(omicsStats)$cnames$fdata_cname
    }
    stats <- omicsStats$Full_results
    
    # --omicsStats and omicsData--#
    if(!is.null(omicsStats) & !is.null(omicsData)){
      # Check if omicsData and omicsStats have the same cname attributes. #
      if (!(identical(attr(omicsData, "cnames"), 
                      attr(omicsStats, "cnames")))) stop(
                        "Non-matching cname attributes in omicsStats and omicsData. 
                        Check that omicsStats is correctly derived from omicsData.")
      
      # Check that the biomolecule unique ID column exists in omicsData and omicsStats #
      if(!(uniqedatId %in% names(omicsStats$Full_results))) stop(
        paste("Column ", 
              uniqedatId,
              " not found in omicsStats. Requires compatible identifiers.", 
              sep = ""))
      
      # Check if stats biomolecules are (at least) a subset of omicsData biomolecules. #
      if(!all(
        omicsStats$Full_results[[uniqedatId]] %in% omicsData$e_data[[uniqedatId]])) stop(
          paste("Biomolecules in omicsStats do not match biomolecules in omicsData.", 
                sep = ""))
      
      # Check if group_designation has been run #
      if(is.null(attr(omicsStats, "group_DF")) | 
         is.null(attr(omicsData, "group_DF"))) stop(
           "Function group_designation() has not been run on data, 
           please run group_designation().")
      
      # Check if stats sampleIDs are (at least) a subset of omicsData f_data sample IDs. #
      if(!all(
        attr(omicsStats, "group_DF")[[sampID]] %in% omicsData$f_data[[sampID]])) stop(
          paste(sampID, "column does not match between omicsData and omicsStats."))
      
      # Check correct object length #
      if(is.null(omicsStats$Full_results) || 
         is.null(omicsStats$Flag) || 
         is.null(omicsStats$P_value)) stop(
           "OmicsStats should contain only and all of the following dataframes: Full_results, P_value, and Flags."
         )
      
      # Check correct lengths of dataframes #
      if(nrow(omicsStats$Full_results) != nrow(omicsStats$P_value) ||
         nrow(omicsStats$Full_results) != nrow(omicsStats$Flags) ||
         nrow(omicsStats$Flags) != nrow(omicsStats$P_value)) stop(
           "Mismatched rows between omicsStats dataframes Full_results, P_value, and Flags. Check integrity of omicsStats object."
         )
      
      # Check attribute statistical test is populated #
      if(!(attr(omicsStats, "statistical_test") %in% c("combined", "anova", "gtest"))) stop(
        "OmicsStats statistical_test attribute incorrect; must be combined, anova, or gtest."
      )
      
      # Check cname e_data is in all columns #
      if(!(uniqedatId %in% colnames(omicsStats$Flags) &&
           uniqedatId %in% colnames(omicsStats$Full_results) &&
           uniqedatId %in% colnames(omicsStats$P_values))) stop(
             paste(uniqedatId, "column must be present in all omicsStats dataframes.")
           )
      
      # Check correct column number #
      n_comps <- length(attr(omicsStats, "comparisons"))
      n_groups <- length(unique(attr(omicsStats, "group_DF")$Group))
      if (attr(omicsStats, "statistical_test") == "combined"){
        # identifier, count for each group, mean for each group, p-val g for each comp, 
        # p-val t for each comp, fold change for each comp, flag for each comp
        full_col <- 1 + 2*n_groups + 4*n_comps
      } else {
        # identifier, count for each group, mean for each group, test p-val for each comp, 
        # fold change for each comp, flag for each comp
        full_col <- 1 + 2*n_groups + 3*n_comps
      }
      if (!((ncol(omicsStats$Flags) == n_comps + 1) && 
            (ncol(omicsStats$P_values) == n_comps + 1) &&
            (ncol(omicsStats$Full_results) == full_col))) stop(
              "Number of columns in omicsStats dataframes is different than expected based on group_DF and comparisons attributes. Some grouping and comparison statistics might be missing from input."
            )
      
      # Check required data #
      if(is.null(omicsData$e_data) || is.null(omicsData$f_data)) stop(
        "Omicsdata requires both e_data and f_data."
      )
      
      # Check cname e_data is in e_data, cname f_data is in f_data #
      if(!(uniqedatId %in% colnames(omicsData$e_data) &&
           sampID %in% colnames(omicsData$f_data))) stop(
             paste(
               paste(uniqedatId, "column must be present in omicsData e_data."),
               paste(sampID, "column must be present in omicsData f_data.")
             )
           )
      
      # Check cname f_data is in columns of e_data #
      if(!all(unlist(omicsData$f_data[sampID]) %in% colnames(omicsData$e_data))) stop(
        paste(sampID, "column in f_data does not match column names in e_data")
      )
      
      # --omicsStats--#
    } else if (!is.null(omicsStats)){
      # Check if group_designation has been run #
      if(is.null(attr(omicsStats, "group_DF"))) stop(
        "Function group_designation() has not been run on data, 
        please run group_designation().")
      
      # Check correct object length #
      if(is.null(omicsStats$Full_results) || 
         is.null(omicsStats$Flag) || 
         is.null(omicsStats$P_value)) stop(
           "OmicsStats should contain only and all of the following dataframes: Full_results, P_value, and Flags."
         )
      
      # Check correct lengths of dataframes #
      if(nrow(omicsStats$Full_results) != nrow(omicsStats$P_value) ||
         nrow(omicsStats$Full_results) != nrow(omicsStats$Flags) ||
         nrow(omicsStats$Flags) != nrow(omicsStats$P_value)) stop(
           "Mismatched rows between omicsStats dataframes Full_results, P_value, and Flags. Check integrity of omicsStats object."
         )
      
      # Check attribute statistical test is populated #
      if(!(attr(omicsStats, "statistical_test") %in% c("combined", "anova", "gtest"))) stop(
        "OmicsStats statistical_test attribute incorrect; must be combined, anova, or gtest."
      )
      
      # Check cname e_data is in all columns #
      if(!(uniqedatId %in% colnames(omicsStats$Flags) &&
           uniqedatId %in% colnames(omicsStats$Full_results) &&
           uniqedatId %in% colnames(omicsStats$P_values))) stop(
             paste(uniqedatId, "column must be present in all omicsStats dataframes.")
           )
      
      # Check correct column number #
      n_comps <- length(attr(omicsStats, "comparisons"))
      n_groups <- length(unique(attr(omicsStats, "group_DF")$Group))
      if (attr(omicsStats, "statistical_test") == "combined"){
        # identifier, count for each group, mean for each group, p-val g for each comp, 
        # p-val t for each comp, fold change for each comp, flag for each comp
        full_col <- 1 + 2*n_groups + 4*n_comps
      } else {
        # identifier, count for each group, mean for each group, test p-val for each comp, 
        # fold change for each comp, flag for each comp
        full_col <- 1 + 2*n_groups + 3*n_comps
      }
      if (!((ncol(omicsStats$Flags) == n_comps + 1) && 
            (ncol(omicsStats$P_values) == n_comps + 1) &&
            (ncol(omicsStats$Full_results) == full_col))) stop(
              "Number of columns in omicsStats dataframes is different than expected based on group_DF and comparisons attributes. Some grouping and comparison statistics might be missing from input."
            )
      
      # --omicsData--#
    } else {
      # Check if group_designation has been run #
      if(is.null(attr(omicsData, "group_DF"))) stop(
        "Function group_designation() has not been run on data, 
        please run group_designation().")
      
      # Check required data #
      if(is.null(omicsData$e_data) || is.null(omicsData$f_data)) stop(
        "Omicsdata requires both e_data and f_data."
      )
      
      # Check cname e_data is in e_data, cname f_data is in f_data #
      if(!(uniqedatId %in% colnames(omicsData$e_data) &&
           sampID %in% colnames(omicsData$f_data))) stop(
             paste(
               paste(uniqedatId, "column must be present in omicsData e_data."),
               paste(sampID, "column must be present in omicsData f_data.")
             )
           )
      
      # Check cname f_data is in columns of e_data #
      if(!all(unlist(omicsData$f_data[sampID]) %in% colnames(omicsData$e_data))) stop(
        paste(sampID, "column in f_data does not match column names in e_data")
      )
    }
  }
}



#' @name format_plot
#' @rdname format_plot
#' @title Plot pairwise comparisons and data values in trellData object
#'
#' @description Plot pairwise comparisons and data values in trellData object. Customizable for plot types, y axis limits, paneling variable (what overall group is plotted on each graph arrangement), as well as desired variables for the y and x axis.
#'
#' @param trellData An object of class "trellData" generated from \code{\link{as.trellData}}.
#' @param p_val Specifies p-value for setting graph border colors
#' @param panel_variable Specifies what to divide trelliscope panels by, must be a column in trellData. Defaults to cnames$edata_cname of trellData.
#' @param plot_text For comparisons only: TRUE/FALSE for p-value text above data
#' @param interactive Should the plots be rendered as plotly objects?
#' @param plot_type types of plots 
#' @param y_limits y-limits
#' @param ... further arguments
#'
#' @author Rachel Richardson
#' @export
format_plot <- function(trellData, ...) {
  .format_plot(trellData,  ...)
}

.format_plot <- function(trellData, 
                         panel_variable = NULL, 
                         plot_type = NULL, 
                         p_val = 0.05, 
                         plot_text = FALSE, 
                         y_limits = NULL,
                         interactive = FALSE) {
 
  ## Input checks
  
  if(!is.null(plot_type) && !all(plot_type %in% c(
                           "abundance_global", 
                           "foldchange_global", 
                           "missing_bar",
                           "abundance_heatmap",
                           "foldchange_heatmap",
                           "presence_heatmap",
                           "foldchange_bar",
                           "abundance_boxplot"
                           ))
     ) stop(
       "plot_type specified is not supported. Must be one of the following: abundance_boxplot, foldchange_bar, abundance_global, foldchange_global, missing_bar, abundance_heatmap, foldchange_heatmap, presence_heatmap, foldchange_bar, abundance_boxplot"
       )
  
  if((length(interactive) != 1 || class(interactive) != "logical")) stop(
    "interactive must be a logical of length 1."
  )
  
  if((length(plot_text) != 1 || class(plot_text) != "logical")) stop(
    "plot_text must be a logical of length 1."
  )
  
  if(!is.null(p_val) && (length(p_val) != 1 || class(p_val) != "numeric")) stop(
    "p_val must be a numeric of length 1."
  )
  
  if(!inherits(trellData, "trellData")) stop(
    "trellData input must be of class trellData"
  )
  
  
  
  # list_y_limits handles y_limits checks
  
  ## End checks
  
  if (is.null(panel_variable)){
    panel_variable <- pmartR::get_edata_cname(trellData)
  }
   
  
  if (panel_variable == pmartR::get_edata_cname(trellData) && stringr::str_detect(plot_type, "heatmap")) {
    stop(
      paste("Heatmaps require a panel_variable that is not the same as edata cname (Default panel_variable is set to edata cname). Current edata cname:", pmartR::get_edata_cname(trellData))
    )
  } 
  
  ## Input digestion
  graphlims <- list_y_limits(plot_type = plot_type, y_limits = y_limits)
  
  
  generate_plots <- function(panel){
    
    
  
    #### Abundance Global ####
    if("abundance_global" %in% plot_type){
      
      lims <- graphlims[["abundance_global"]]
      
      if (!is.null(lims$scale) && lims$scale == "free") {
        warning( 
          "Scale option 'free' is not valid with abundance_global plots. Fixed plots will be generated based on global values.")
      }
      
      sampID <- pmartR:::get_fdata_cname(trellData)
      
      if ("Group_DF" %in% colnames(trellData$data_value)){
        group <- "Group_DF"
      } else {
        group <- "Group"
      }
      
      abundance <- grep("abundance", colnames(trellData$data_value), value = TRUE)
      
      # reorder levels #
      plot_data <- trellData$data_value[order(trellData$data_value[,group]),]
      plot_data[[sampID]] <- factor(plot_data[[sampID]], levels=unique(plot_data[[sampID]]), ordered=TRUE)
      
      plot_base <- ggplot2::ggplot() +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(
          x = plot_data[[sampID]],
          y = plot_data[[abundance]],
          color = plot_data[[group]],
          fill = plot_data[[group]]),
          alpha = 0.75,
          na.rm = TRUE) +
        ggplot2::labs(x = group, y = abundance, fill = group, color = group) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
      
      
      data <- trellData$data_value
      data <- plot_data

        rows <- as.character(data[[panel_variable]]) == panel
        df1 <- data[rows,]
        
        if(!is.null(lims$max) || 
           !is.null(lims$min) ||
           !is.null(lims$range)){
          
          increment <- set_increment(df1[[abundance]])
          setlims1 <- set_ylimits(df1[[abundance]],
                                 increment,
                                 y_range = lims$range,
                                 y_max = lims$max,
                                 y_min = lims$min,
                                 include_zero = FALSE)
        }
          
          plot_out <- plot_base +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 330,
                                                               hjust = 0,
                                                               vjust = 0.5)) +
            ggplot2::labs(x = NULL, fill = NULL, color = NULL) +
            
            ggplot2::geom_segment(
              x = as.numeric(factor(df1[[sampID]])) - 0.25,
              xend = as.numeric(factor(df1[[sampID]])) + 0.25,
              y = df1[[abundance]],
              yend = df1[[abundance]],
              ggplot2::aes(text = paste(paste(abundance, ":", sep = ""), signif(as.numeric(df1[[abundance]]), 4))),
              na.rm = TRUE) +
            
            ggplot2::geom_point(
              x = as.numeric(factor(df1[[sampID]])),
              y = df1[[abundance]],
              ggplot2::aes(text = paste(paste(abundance, ":", sep = ""), signif(as.numeric(df1[[abundance]]), 4))),
              color = NA,
              na.rm = TRUE)
          
          if(length(unique(trellData$data_value[[sampID]])) > 8){
            plot_out <- plot_out + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5))
          }
          if(length(unique(trellData$data_value[[sampID]])) > 20){
            plot_out <- plot_out + ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
              ggplot2::xlab(sampID)
          }
          
          if(exists("setlims1")){
            plot_out <- plot_out + ggplot2::coord_cartesian(ylim = setlims1)
          }
      
    }
    
    #### Foldchange Global ####
    if("foldchange_global" %in% plot_type){
      
      lims <- graphlims[["foldchange_global"]]
  
      if (!is.null(lims$scale) && lims$scale == "free") {
        warning( 
          "Scale option 'free' is not valid with foldchange_global plots. Fixed plots will be generated based on global values.")
        lims$scale <- "fixed"
      }
      
      sampID <- pmartR:::get_fdata_cname(trellData)
      comps <- "Comparison"
      
      foldchange <- grep("Fold_change", colnames(trellData$comp_stats), value = TRUE)
      
      plot_base <- ggplot2::ggplot() +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(x = trellData$comp_stats[[comps]], 
                                           y = trellData$comp_stats[[foldchange]], 
                                           color = trellData$comp_stats[[comps]],
                                           fill = trellData$comp_stats[[comps]]), 
                              alpha = 0.75, na.rm = TRUE) +
        ggplot2::labs(x = comps, y = foldchange, fill = comps, color = comps) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                       legend.position = "none")
      
      
      data <- trellData$comp_stats
      
      if(!(panel_variable %in% colnames(trellData$comp_stats))) stop(
        paste(
          "Panel variable for foldchange_global must be selected from either edata or emeta columns. The following columns are valid:", colnames(trellData$comp_stats)
      ))
      
        rows <- as.character(data[[panel_variable]]) == panel
        df2 <- data[rows,]
        
        if(!is.null(lims$max) || 
           !is.null(lims$min) ||
           !is.null(lims$range)){
          
          increment <- set_increment(df2[[foldchange]])
          setlims2 <- set_ylimits(df2[[foldchange]],
                                 increment,
                                 y_range = lims$range,
                                 y_max = lims$max,
                                 y_min = lims$min,
                                 include_zero = FALSE)
        }
    
          plot_out <- plot_base + 
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 330,
                                                               hjust = 0,
                                                               vjust = 0.5)) +
            ggplot2::labs(x = NULL)   +
            
            ggplot2::geom_segment(
              x = as.numeric(factor(df2[[comps]])) - 0.25,
              xend = as.numeric(factor(df2[[comps]])) + 0.25,
              y = df2[[foldchange]],
              yend = df2[[foldchange]],
              ggplot2::aes(text = #paste(
                             paste(paste(foldchange, ":", sep = ""), signif(as.numeric(df2[[foldchange]]), 4))),
              # paste("Rank:", df2[["rank"]]),
              # sep = "\n")),
              na.rm = TRUE) +
            
            ggplot2::geom_point(
              x = as.numeric(factor(df2[[comps]])),
              y = df2[[foldchange]],
              ggplot2::aes(text = paste(paste(foldchange, ":", sep = ""), signif(as.numeric(df2[[foldchange]]), 4))),
              color = NA,
              na.rm = TRUE)
          
          
          if(exists("setlims2")){
            plot_out <- plot_out + ggplot2::coord_cartesian(ylim = setlims2)
          }
    }
    
    #### Missing Bar ####
    if("missing_bar" %in% plot_type){
      
      lims <- graphlims[["missing_bar"]]
      
      if (!is.null(lims$scale) || !is.null(lims$range)) {
        warning( 
          "Scale options and range are not valid with missing_bar plots, scale will not be applied. Please use max and min or list(NULL) for missing_bar y_limits.  Refer to examples in ?trelliVis().")
      }
      
      if (!is.null(lims$min) && (lims$min < 0 || lims$min > 1)) {
        stop( 
          "Minimum and maximum for missing_bar is only supported for proportions; select values between 0 and 1.")
      }
      
      if (!is.null(lims$max) && (lims$max < 0 || lims$max > 1)) {
        stop( 
          "Minimum and maximum for missing_bar is only supported for proportions; select values between 0 and 1.")
      }
      
      
      if(!is.null(trellData$summary_stats)){
        
        if ("Group_DF" %in% colnames(trellData$summary_stats)){
          group <- "Group_DF"
        } else {
          group <- "Group"
        }
        
          rows <- as.character(trellData$summary_stats[[panel_variable]]) == panel
          df3 <- trellData$summary_stats[rows,]
          totals <- data.frame(pmartR::get_group_table(trellData), stringsAsFactors = FALSE)
          
          colnames(totals) <- c(group, "Total_Group_Counts")
          df3 <- suppressWarnings(dplyr::left_join(df3, totals, by = group))
          
          df3$Missing <-  df3[["Total_Group_Counts"]] - df3[["Count"]]
          df3$Non_Missing <- df3[["Count"]]
          
          df3 <- reshape2::melt(df3, 
                               id.vars = colnames(df3)[!stringr::str_detect(colnames(df3), "Missing")], 
                               value.name = "y",
                               variable.name = "Data")
          
          
          groupgraphs <- suppressWarnings(purrr::map(unique(df3[[group]]), function(df3gr){
            
            rows <- df3[[group]] == df3gr
            df3sm <- df3[rows,c("y", "Data", group, "Total_Group_Counts")]
            
            miss <- df3sm$Data == "Missing"
            df3sm$y[miss] <- sum(df3sm$y[miss])
            df3sm$y[!miss] <- sum(df3sm$y[!miss])
            
            dup <- nrow(df3sm)/2
            df3sm[["Total_Group_Counts"]] <- df3sm[["Total_Group_Counts"]]*dup
            
            df3sm <- unique.array(df3sm)
            
            df3sm <- dplyr::arrange(df3sm, Data, y)
            
            df3sm[["Data"]] <- as.character(df3sm[["Data"]])
            df3sm[["Data"]] <- factor(df3sm[["Data"]], levels = c("Missing", "Non_Missing"))
            df3sm[["y"]] <- as.integer(df3sm[["y"]])
            df3sm[["Total_Group_Counts"]] <- as.integer(df3sm[["Total_Group_Counts"]])
            
            setlims3 <- NULL
            
            if(is.null(lims$max)){
              lims$max <- 1
            }
            if(is.null(lims$min)){
              lims$min <- 0
            }
            
            increment <- set_increment(as.numeric(df3sm[["y"]]))
            setlims3 <- set_ylimits(as.numeric(df3sm[["y"]]),
                                   increment,
                                   y_max = lims$max,
                                   y_min = lims$min,
                                   include_zero = FALSE)
            
            texty1 <- paste(
              paste(
                paste(paste(group, ":", sep = ""),
                      df3sm[[group]]),
                paste("Count:", df3sm[["y"]]), 
                sep = "\n"),
              paste("Status:", df3sm[["Data"]]),
              sep = "\n")
              
              plot <- ggplot2::ggplot() +
                ggplot2::theme_bw() +
                ggplot2::geom_col(
                  ggplot2::aes(
                               x = df3sm[[group]],
                               y = df3sm[["y"]],
                               fill = df3sm[["Data"]],
                               order = df3sm[["Data"]],
                               text = texty1),
                  position = "fill") +
                ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
                ggplot2::scale_y_continuous(
                  sec.axis = ggplot2::sec_axis(
                    ~.*max(df3sm[["Total_Group_Counts"]]),
                    name = "", 
                    breaks = c(seq(from = 0,
                                   to = max(df3sm[["Total_Group_Counts"]]) - round(max(df3sm[["Total_Group_Counts"]])/10),
                                   by = round(max(df3sm[["Total_Group_Counts"]])/10)+1),
                               max(df3sm[["Total_Group_Counts"]])
                    )),
                  expand = c(0,0)) +
                ggplot2::scale_x_discrete(expand = c(0,0)) + 
                ggplot2::coord_cartesian(ylim = setlims3)
            
            return(plot)
            
          }))
          
          if (interactive){
              plotlys <- purrr::map(groupgraphs, function(plot) plotly::ggplotly(plot, tooltip = "text"))
              plot_out <- plotly::subplot(plotlys, margin = 0.1)
          } else {
              plot <- ggpubr::ggarrange(plotlist = groupgraphs, common.legend = TRUE)
              plot_out <- ggpubr::annotate_figure(plot, 
                                              right = ggpubr::text_grob("Count", size = 8, rot = 270), 
                                              left = ggpubr::text_grob("Proportion", size = 8, rot = 90))
          }
  
      } else {
        
        if ("Group_DF" %in% colnames(trellData$data_values)){
          group <- "Group_DF"
        } else {
          group <- "Group"
        }
         rows <- as.character(trellData$data_values[[panel_variable]]) == panel
          df3 <- trellData$data_values[rows,]
          
          Count <- data.frame(table(df3[complete.cases(df3),]$Group))
          totals <- data.frame(pmartR::get_group_table(trellData), stringsAsFactors = FALSE)
          n_rep <- length(unique(df3[[pmartR::get_edata_cname(trellData)]]))
          
          totals$Freq <- totals$Freq * n_rep
          
          colnames(Count) <- c(group, "Count")
          colnames(totals) <- c(group, "Total_Group_Counts")
          df3 <- suppressWarnings(dplyr::left_join(df3, Count, by = group))
          df3 <- suppressWarnings(dplyr::left_join(df3, totals, by = group))
          
          df3$Missing <-  df3[["Total_Group_Counts"]] - df3[["Count"]]
          df3$Non_Missing <- df3[["Count"]]
          
          df3 <- reshape2::melt(df3, 
                               id.vars = colnames(df3)[!stringr::str_detect(colnames(df3), "Missing")],
                               value.name = "y",
                               variable.name = "Data")
          
          
          groupgraphs <- purrr::map(unique(df3[[group]]), function(df3gr){
            rows <- df3[[group]] == df3gr
            df3sm <- df3[rows, c("y", "Data", group, "Total_Group_Counts")]
            df3sm <- unique.array(df3sm )
            
            df3sm <- dplyr::arrange(df3sm, Data, y)
            
            df3sm[["Data"]] <- as.character(df3sm[["Data"]])
            df3sm[["Data"]] <- factor(df3sm[["Data"]], levels = c("Missing", "Non_Missing"))
            
            if(any(is.na(df3sm[["y"]]))) {
              df3sm[["y"]][is.na(df3sm[["y"]])] <- 0
              df3sm[["y"]][df3sm[["Data"]] == "Missing"] <- df3sm[["Total_Group_Counts"]][df3sm[["Data"]] == "Missing"]
            }
            
            setlims3 <- NULL
            
            if(is.null(lims$max)){
              lims$max <- 1
            }
            if(is.null(lims$min)){
              lims$min <- 0
            }
  
            increment <- set_increment(as.numeric(df3sm[["y"]]))
            setlims3 <- set_ylimits(as.numeric(df3sm[["y"]]),
                                   increment,
                                   y_max = lims$max,
                                   y_min = lims$min,
                                   include_zero = FALSE)
            
            texty1 <- paste(
              paste(
                paste(paste(group, ":", sep = ""),
                      df3sm[[group]]),
                paste("Count:", df3sm[["y"]]), 
                sep = "\n"),
              paste("Status:", df3sm[["Data"]]),
              sep = "\n")
            
            plot <- ggplot2::ggplot() +
              ggplot2::theme_bw() +
              ggplot2::geom_col(ggplot2::aes(x = df3sm[[group]],
                                             y = df3sm[["y"]],
                                             fill = df3sm[["Data"]],
                                             order = df3sm[["Data"]],
                                             text = texty1),
                                position = "fill") +
              ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
              ggplot2::scale_y_continuous(
                sec.axis = ggplot2::sec_axis(
                  ~.*max(df3sm[["Total_Group_Counts"]]),
                  name = "", 
                  breaks = c(seq(from = 0,
                                 to = max(df3sm[["Total_Group_Counts"]]) - round(max(df3sm[["Total_Group_Counts"]])/10),
                                 by = round(max(df3sm[["Total_Group_Counts"]])/10)+1),
                             max(df3sm[["Total_Group_Counts"]])
                  )),
                expand = c(0,0)) +
              ggplot2::scale_x_discrete(expand = c(0,0)) +
              ggplot2::coord_cartesian(ylim = setlims3)
            
            return(plot)
            
          })
          
          
          if (interactive){
              plotlys <- purrr::map(groupgraphs, function(plot) plotly::ggplotly(plot, tooltip = "text"))
              plot_out <- plotly::subplot(plotlys, margin = 0.05)
          } else {
              plot <- ggpubr::ggarrange(plotlist = groupgraphs, common.legend = TRUE)
              plot_out <- ggpubr::annotate_figure(plot, 
                                                  right = ggpubr::text_grob("Count", size = 8, rot = 270), 
                                                  left = ggpubr::text_grob("Proportion", size = 8, rot = 90))
      }
      }
    }
    
    #### Abundance Heatmap ####
    if("abundance_heatmap" %in% plot_type){
      
      lims <- graphlims[["abundance_heatmap"]]
      
      ## Move outside function
      
      # if (!is.null(lims[[1]])) {
      #   warning( 
      #     "y_limits are not supported with heatmaps and will not be used.  Refer to examples in ?trelliVis().")
      # }
      
      
      if(attr(trellData, "meta_info") && 
         panel_variable != pmartR::get_edata_cname(trellData)){
        
        if ("Group_DF" %in% colnames(trellData$data_values)){
          group <- "Group_DF"
        } else {
          group <- "Group"
        }
        
        abundance <- grep("abundance", colnames(trellData$data_values), value = T)
        
        pal <- grDevices::colorRampPalette(c("red", "yellow"))
        
          df4 <- trellData$data_values[as.character(trellData$data_values[[panel_variable]]) == panel,]
          
          df4[[pmartR::get_fdata_cname(trellData)]] <- ordered(
            df4[[pmartR::get_fdata_cname(trellData)]], 
            levels = rev(sort(unique(df4[[pmartR::get_fdata_cname(trellData)]]))))
          
          df4[[pmartR::get_edata_cname(trellData)]] <- ordered(
            df4[[pmartR::get_edata_cname(trellData)]], 
            levels = rev(sort(unique(df4[[pmartR::get_edata_cname(trellData)]]))))
          
          texty2 <- paste(
            paste(
              paste(paste(pmartR::get_edata_cname(trellData), ":", sep = ""),
                    df4[[pmartR::get_edata_cname(trellData)]]),
              paste(paste(pmartR::get_fdata_cname(trellData), ":", sep = ""),
                    df4[[pmartR::get_fdata_cname(trellData)]]), sep = "\n"),
            paste(paste(abundance, ":", sep = ""), signif(as.numeric(df4[[abundance]]), 4)),
            sep = "\n")
            
            plot_out <- ggplot2::ggplot() +
              ggplot2::geom_tile(ggplot2::aes(text = texty2,
                                              fill = df4[[abundance]],
                                              x = df4[[pmartR::get_fdata_cname(trellData)]],
                                              y = df4[[pmartR::get_edata_cname(trellData)]])) +
              ggplot2::scale_fill_gradientn(abundance, colours = pal(50)) +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                             axis.ticks = ggplot2::element_blank(),
                             plot.title = ggplot2::element_text(size=2)) +
              ggplot2::theme_bw() +
              ggplot2::scale_x_discrete(expand = c(0, 0)) +
              ggplot2::scale_y_discrete(expand = c(0, 0)) +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 330, hjust = 0, size = 6),
                             axis.ticks = ggplot2::element_blank(),
                             legend.title = ggplot2::element_text(size = 8)) +
              ggplot2::xlab("")
            
            if (length(unique(df4[[pmartR::get_edata_cname(trellData)]])) > 35){
              plot_out <- plot_out + ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
                ggplot2::ylab(pmartR::get_edata_cname(trellData))
            } else {
              plot_out <- plot_out + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6)) +
                ggplot2::ylab("")
            }
            
            if (length(unique(df4[[pmartR::get_fdata_cname(trellData)]])) > 35){
              plot_out <- plot_out + ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
                ggplot2::xlab(pmartR::get_fdata_cname(trellData))
            } else {
              plot_out <- plot_out + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6)) +
                ggplot2::xlab("")
            }
      }
    }
    
    #### Foldchange Heatmap ####
    if("foldchange_heatmap" %in% plot_type){
      
      
      lims <- graphlims[["foldchange_heatmap"]]
      
      if (!is.null(lims[[1]])) {
        warning( 
          "y_limits are not supported with heatmaps and will not be used.  Refer to examples in ?trelliVis().")
      }
      
      if(attr(trellData, "meta_info") && panel_variable != pmartR::get_edata_cname(trellData)){
        
        foldchange <- grep("change", colnames(trellData$comp_stats), value = T)
        
        pal <- grDevices::colorRampPalette(c("red", "yellow"))
        
        df5 <- trellData$comp_stats[as.character(trellData$comp_stats[[panel_variable]]) == panel,]
          
          texty3 <- paste(
            paste(
              paste(paste(pmartR::get_edata_cname(trellData), ":", sep = ""),
                  df5[[pmartR::get_edata_cname(trellData)]]),
              paste("Comparison:", df5[["Comparison"]]), sep = "\n"),
            paste(paste(foldchange, ":", sep = ""), signif(as.numeric(df5[[foldchange]]), 4)),
            sep = "\n")
            
            plot_out <- ggplot2::ggplot(
              df5, 
              ggplot2::aes(text = texty3,
                           x = ordered(!!rlang::sym("Comparison"),
                                       levels = rev(sort(unique(!!rlang::sym("Comparison"))))),
                           y = ordered(!!rlang::sym(pmartR::get_edata_cname(trellData)),
                                       levels = rev(sort(unique(!!rlang::sym(pmartR::get_edata_cname(trellData)))))))) +
              ggplot2::geom_tile(ggplot2::aes(fill = !!rlang::sym(foldchange))) +
              ggplot2::scale_fill_gradientn(foldchange, colours = pal(50)) +
              ggplot2::xlab("") +  ggplot2::ylab("") +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                             axis.ticks = ggplot2::element_blank(),
                             plot.title = ggplot2::element_text(size=2)) +
              ggplot2::theme_bw() + 
              ggplot2::scale_x_discrete(expand = c(0, 0)) +
              ggplot2::scale_y_discrete(expand = c(0, 0)) +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 330, hjust = 0, size = 6), 
                             axis.ticks = ggplot2::element_blank(),
                             legend.title = ggplot2::element_text(size = 8))
            
            if (length(unique(df5[[pmartR::get_edata_cname(trellData)]])) > 35){
              plot_out <- plot_out + ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
                ggplot2::ylab(pmartR::get_edata_cname(trellData))
            } else {
              plot_out <- plot_out + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6)) +
                ggplot2::ylab("")
            }
            
            if (length(unique(df5[[pmartR::get_fdata_cname(trellData)]])) > 35){
              plot_out <- plot_out + ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
                ggplot2::xlab(pmartR::get_fdata_cname(trellData))
            } else {
              plot_out <- plot_out + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6)) +
                ggplot2::xlab("")
            }
      }
    }
    
    #### Presence Heatmap ####
    if("presence_heatmap" %in% plot_type){
      
      lims <- graphlims[["presence_heatmap"]]
      
      if (!is.null(lims[[1]])) {
        warning( 
          "y_limits are not supported with heatmaps and will not be used.  Refer to examples in ?trelliVis().")
      }
      
      
      if(attr(trellData, "meta_info") && panel_variable != pmartR::get_edata_cname(trellData)){
        # if(!is.null(trellData$data_values)){
          
          abundance <- grep("abundance", colnames(trellData$data_values), value = T)
          pal <- grDevices::colorRampPalette(c("red", "yellow"))
            
            df6 <- trellData$data_values[as.character(trellData$data_values[[panel_variable]]) == panel,]
            
            df6 <- tidyr::nest(df6, -c(!!rlang::sym(pmartR::get_edata_cname(trellData)),
                                     !!rlang::sym(pmartR::get_fdata_cname(trellData))))
            df6 <- dplyr::mutate(df6, Biomolecule_Presence = purrr::map_lgl(df6$data, function(df6dat){
              !all(is.na(df6dat[[abundance]]))
            }))
            
            df6[["Biomolecule_Presence"]] <- gsub(TRUE, "Present", df6[["Biomolecule_Presence"]])
            df6[["Biomolecule_Presence"]] <- gsub(FALSE, "Absent", df6[["Biomolecule_Presence"]])
            
            
            texty4 <- paste(
              paste(
                paste(paste(pmartR::get_edata_cname(trellData), ":", sep = ""),
                      df6[[pmartR::get_edata_cname(trellData)]]),
                paste(paste(pmartR::get_fdata_cname(trellData), ":", sep = ""),
                      df6[[pmartR::get_fdata_cname(trellData)]]), sep = "\n"),
              paste("Biomolecule_Presence:", df6[["Biomolecule_Presence"]]),
              sep = "\n")
            
            
      
              
              plot_out <- ggplot2::ggplot(
                df6,
                ggplot2::aes(text = texty4,
                             x = ordered(!!rlang::sym(pmartR::get_fdata_cname(trellData)),
                                         levels = rev(sort(unique(!!rlang::sym(pmartR::get_fdata_cname(trellData)))))),
                             y = ordered(!!rlang::sym(pmartR::get_edata_cname(trellData)),
                                         levels = rev(sort(unique(!!rlang::sym(pmartR::get_edata_cname(trellData)))))))) +
                ggplot2::geom_tile(ggplot2::aes(fill = Biomolecule_Presence)) +
                ggplot2::scale_fill_manual(NULL, values = c("Present" = "blue", "Absent" = "grey")) +
                ggplot2::xlab("") +  ggplot2::ylab("") + ggplot2::labs(color = NULL) +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                               axis.ticks = ggplot2::element_blank(),
                               plot.title = ggplot2::element_text(size=2)) +
                ggplot2::theme_bw() +
                ggplot2::scale_x_discrete(expand = c(0, 0)) +
                ggplot2::scale_y_discrete(expand = c(0, 0)) +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 330, hjust = 0, size = 6),
                               axis.ticks = ggplot2::element_blank())
              
              if (length(unique(df6[[pmartR::get_edata_cname(trellData)]])) > 35){
                plot_out <- plot_out + ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
                  ggplot2::ylab(pmartR::get_edata_cname(trellData))
              } else {
                plot_out <- plot_out + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6)) +
                  ggplot2::ylab("")
              }
              
              if (length(unique(df6[[pmartR::get_fdata_cname(trellData)]])) > 35){
                plot_out <- plot_out + ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
                  ggplot2::xlab(pmartR::get_fdata_cname(trellData))
              } else {
                plot_out <- plot_out + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6)) +
                  ggplot2::xlab("")
              }
              
      }
    }
    
    #### Abundance Boxplot ####
    if("abundance_boxplot" %in% plot_type){
      
      lims <- graphlims[["abundance_boxplot"]]
      
      if ("Group_DF" %in% colnames(trellData$data_values)) {
        group_df_name <- "Group_DF"
      } else {
        group_df_name <- "Group"
      }
      
      abundance <- grep("abundance",
                        colnames(trellData$data_values),
                        value = TRUE)
      
      plotter <- trellData$data_values
      
      if(!is.null(lims$scale) && lims$scale == "fixed"){
        
        increment <- set_increment(plotter[[abundance]])
        setlims4 <- set_ylimits(plotter[[abundance]],
                               increment,
                               y_range = lims$range,
                               y_max = lims$max,
                               y_min = lims$min,
                               include_zero = FALSE)
      }
      
      nestedData <- plotter[as.character(plotter[[panel_variable]]) == panel,]
      
      ## Set hover, excluding the panel_variable ##
      nestedData_value <- add_plot_features(trellData,
                                            nestedData,
                                            p_val = p_val,
                                            panel_variable = panel_variable,
                                            value_panel_y_axis = abundance)  ### Take out of function later 
      
      nestedData_value <- unique.data.frame(nestedData_value[c(
        group_df_name,
        abundance,
        "text")])
      
      # Make ggplots #
      if(!exists("setlims4")){
        
        increment <- set_increment(nestedData_value[[abundance]])
        setlims4 <- set_ylimits(nestedData_value[[abundance]],
                               increment,
                               y_range = lims$range,
                               y_max = lims$max,
                               y_min = lims$min,
                               include_zero = FALSE)
      }
        
        plot_out <- ggplot2::ggplot() +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position='none',
                         axis.text.x = ggplot2::element_text(angle = 330, hjust = 0)) +
          ggplot2::xlab(group_df_name) + 
          ggplot2::ylab(abundance) +
          ggplot2::labs(color = "", fill = "") +
          
          ggplot2::geom_point(
            ggplot2::aes(
              text = gsub(", ", "\n", nestedData_value[["text"]]),
              x = as.character(nestedData_value[[group_df_name]]),
              y = nestedData_value[[abundance]],
              fill = nestedData_value[[group_df_name]]),
            position = "jitter",
            size = 2,
            color = "black",
            na.rm = TRUE, 
            alpha = 0.7) +
          
          ggplot2::geom_boxplot(
            ggplot2::aes(
              x = as.character(nestedData_value[[group_df_name]]),
              y = nestedData_value[[abundance]],
              fill = nestedData_value[[group_df_name]]),
            alpha = 0.2,
            position = "dodge2",
            na.rm = TRUE) +
          ggplot2::labs(x = NULL) + 
          
          ggplot2::coord_cartesian(ylim = setlims4)
    }
    
    
    #### Foldchange bar ####
    if("foldchange_bar" %in% plot_type){
      
      lims <- graphlims[["foldchange_bar"]]
      
      if ("Group_DF" %in% colnames(trellData$summary_stats)) {
        group_df_name <- "Group_DF"
      } else {
        group_df_name <- "Group"
      }
      
      plotter <- trellData$comp_stats
      
      if(!is.null(panel_variable) && !(panel_variable %in% colnames(trellData$comp_stats))) stop(
        paste(
          c("Panel variable for foldchange_bar must be selected from either edata or emeta columns. The following columns are valid:", colnames(trellData$comp_stats)), 
          collapse = " "
        ))
      
      if(!is.null(lims$scale) && lims$scale == "fixed"){
        
        setter <- plotter[["Fold_change"]] +
          (0.55 * plotter[["Fold_change"]])
        
        increment <- set_increment(setter)
        setlims5 <- set_ylimits(setter,
                               increment,
                               y_range = lims$range,
                               y_max = lims$max,
                               y_min = lims$min,
                               include_zero = TRUE)
      }
      
      nestedData <- plotter[as.character(plotter[[panel_variable]]) == panel,]
        
        ## Add border color, hover text, and label text to dataframe for plotting ##
        nestedData_comps <- add_plot_features(trellData,
                                              nestedData,
                                              p_val = p_val,
                                              panel_variable = panel_variable,
                                              comps_panel_y_axis = "Fold_change"   ### Remove from function later
                                              )
        
        nestedData_comps <- unique.data.frame(nestedData_comps[c(
          pmartR::get_edata_cname(trellData),
          "Fold_change",
          "text",
          "Comparison",
          "bord",
          "labels")])
        
        # Make ggplots #
        
        if(!exists("setlims5")){
          
          setter <- nestedData_comps[["Fold_change"]] +
            (0.55 * nestedData_comps[["Fold_change"]])
          
          increment <- set_increment(setter)
          setlims5 <- set_ylimits(setter,
                                 increment,
                                 y_range = lims$range,
                                 y_max = lims$max,
                                 y_min = lims$min,
                                 include_zero = TRUE)
        }
        
          plot_out <- ggplot2::ggplot() + 
            ggplot2::theme_bw() +
            ggplot2::geom_hline(yintercept = 0) +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 330, hjust = 0)) +
            ggplot2::xlab("Comparison") + 
            ggplot2::ylab("Fold_change") 

            plot_out <- plot_out + ggplot2::geom_col(
              ggplot2::aes(
                text = gsub(", ", "\n", nestedData_comps[["text"]]),
                x = as.character(nestedData_comps[[pmartR::get_edata_cname(trellData)]]),
                y = as.numeric(nestedData_comps[["Fold_change"]]),
                fill = as.character(nestedData_comps[["Comparison"]]),
                color = nestedData_comps[["bord"]]),
              position = "dodge2",
              size = 1,
              na.rm = TRUE) +
              ggplot2::labs(fill = "", color = "", x = NULL) +
              ggplot2::scale_color_manual(values = c("grey40" = "grey40", "black" = "black"), guide = FALSE)

            if (length(unique(nestedData_comps[[pmartR::get_edata_cname(trellData)]])) > 8){
              plot_out <- plot_out + ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 330, hjust = 0, size = 5))
            }

            if (length(unique(nestedData_comps[[pmartR::get_edata_cname(trellData)]])) > 21){
              plot_out <- plot_out + ggplot2::theme(
                axis.text.x = ggplot2::element_blank()) +
                ggplot2::xlab("Biomolecules")
            }

            if (plot_text){
              
              textcomps <- nestedData_comps
              textcomps[["Fold_change"]][is.na(textcomps[["Fold_change"]])] <- 0
              
              plot_out <-  plot_out + ggplot2::geom_text(
                ggplot2::aes(
                  x = as.character(nestedData_comps[[pmartR::get_edata_cname(trellData)]]),
                  y = (max(abs(setlims5)) -
                         (0.20 * max(abs(setlims5)))) * 
                    sign(abs(max(setlims5)) - abs(min(setlims5))),
                  label = gsub(", ", "\n", textcomps[["labels"]])),
                color = "black"
              )
            }
            
          plot_out <- plot_out + ggplot2::coord_cartesian(ylim = setlims5)
          
          if (interactive){
            plot_out <- plot_out + ggplot2::theme(legend.position = "none")
          }
        }
    
    #######
    #######
    
    if(interactive){
      
      if(any(stringr::str_detect(plot_type, "heatmap"))) {
        
        plotnames <- c("abundance_boxplot", 
                        "foldchange_bar", 
                        "abundance_global", 
                        "foldchange_global", 
                        "missing_bar",
                        "abundance_heatmap",
                        "foldchange_heatmap",
                        "presence_heatmap",
                        "foldchange_bar",
                        "abundance_boxplot")
        typenames <- plotnames[plotnames %in% plot_type]
        #Accounts for out of order input
        
        if(stringr::str_detect(typenames[1], "heatmap")){
          
          # plot_out <- plot_out + ggplot2::theme(legend.position = "none")
          plot_out <- plotly::ggplotly(plot_out, tooltip = c("text"))
          plot_out <- plot_out %>% plotly::layout(plot_bgcolor='grey',
                                                  xaxis = list(showgrid = F),
                                                  yaxis = list(showgrid = F)
                                                  )
          
        }

      }
      
        return(plotly::ggplotly(plot_out, tooltip = "text"))
      
    } else {
        return(plot_out)
    }
  }
  
  
  ##### Process plots
  
  pvs <- data.frame(unique(as.character(trellData$summary_stats[[panel_variable]])))
  
  if(nrow(pvs) == 0){
    pvs <- data.frame(unique(as.character(trellData$data_values[[panel_variable]])))
  }
  colnames(pvs) <- panel_variable
  
  pvs <- dplyr::mutate(pvs, 
                       panel = trelliscopejs::map_plot(
                         pvs[[panel_variable]], 
                         function(panel) suppressWarnings(generate_plots(panel))))
  attr(pvs, "data_class") <- attr(trellData, "data_class")
  
  return(pvs)
  
}


#' @name validate_format_plot_input
#' @rdname validate_format_plot_input
#' @title Validate inputs for omicsData and omicsStats in trelliVis processing
#' 
#' @description Checks for validity of trellData input and assigns variables where needed
#' 
#' @param trellData An object of class "trellData" generated from \code{\link{as.trellData}}.
#' @param p_val Specifies p-value for setting graph border colors
#' @param panel_variable Specifies what to divide trelliscope panels by, must be a column in trellData. Defaults to cnames$edata_cname of trellData.
#' @param plot_text For comparisons only: TRUE/FALSE for p-value text above data
#' @param plotly Should the plots be rendered as plotly objects?
#' @param ... further arguments
#' 
#' @author Rachel Richardson


# validate_format_plot_input <- function(...){
#   .validate_format_plot_input(...)
# }
# 
# .validate_format_plot_input <- function(trellData, 
#                                   comps_y_limits = NULL, comps_y_range = NULL, 
#                                   comps_y_max = NULL, comps_y_min = NULL, 
#                                   comps_include_zero = NULL,
#                                   value_y_limits = NULL, value_y_range = NULL, 
#                                   value_y_max = NULL, value_y_min = NULL,
#                                   value_include_zero = NULL,
#                                   p_val = NULL,
#                                   panel_variable = NULL, 
#                                   comps_color_variable = NULL,
#                                   comps_panel_x_axis = NULL, comps_panel_y_axis = NULL,
#                                   value_color_variable = NULL,
#                                   value_panel_x_axis = NULL, 
#                                   value_panel_y_axis = NULL,
#                                   value_plot_type = NULL, comps_plot_type = NULL,
#                                   value_plot = NULL, comps_plot = NULL, 
#                                   comps_text = NULL, plotly = NULL) {
#   
#   ##### Initial checks #####
#   
#   # Check if class is correct #
#   if(!inherits(trellData, "trellData")) stop(
#     "trellData must be of the class 'trellData'")
#   
#   if(!(is.null(p_val) && 
#        is.null(plotly) &&
#        is.null(comps_y_limits) &&
#        is.null(comps_y_range) &&
#        is.null(comps_y_max) &&
#        is.null(comps_y_min) &&
#        is.null(comps_include_zero) &&
#        is.null(comps_plot) &&
#        is.null(comps_text) &&
#        is.null(comps_plot_type) &&
#        is.null(value_y_limits) &&
#        is.null(value_y_range) &&
#        is.null(value_y_max) &&
#        is.null(value_y_min) &&
#        is.null(value_plot_type) &&
#        is.null(value_include_zero) &&
#        is.null(value_plot)
#        )) {
# 
#     
#     # Check if comp_stats is in trellData #
#     if(is.null(trellData$comp_stats) && is.null(trellData$data_values)) stop(
#       "No data values or comparison statistics in trellData to plot")
#     
#     # Check if comp_stats is in trellData (as above) #
#     if(is.na(trellData$comp_stats) && is.na(trellData$data_values)) stop(
#       "No data values or comparison statistics in trellData to plot")
#     
#     # Check check if p_val is numeric of length 1 #
#     if(!is.numeric(p_val) || (length(p_val) != 1)) stop(
#       "p_val must be a numeric of length 1")  
#     
#     # Check check if plotly is logical of length 1 #
#     if(!is.logical(plotly) || (length(plotly) != 1)) stop(
#       "plotly must be a logical (TRUE or FALSE) of length 1") 
#     
#     # Check check if comps_include_zero is logical of length 1 #
#     if(!is.logical(comps_include_zero) || (length(comps_include_zero) != 1)) stop(
#       "comps_include_zero must be a logical (TRUE or FALSE) of length 1") 
#     
#     # Check check if comps_plot is logical of length 1 #
#     if(!is.logical(comps_plot) || (length(comps_plot) != 1)) stop(
#       "comps_plot must be a logical (TRUE or FALSE) of length 1") 
#     
#     # Check check if comps_text is logical of length 1 #
#     if(!is.logical(comps_text) || (length(comps_text) != 1)) stop(
#       "comps_text must be a logical (TRUE or FALSE) of length 1") 
#     
#     # Check check if value_include_zero is logical of length 1 #
#     if(!is.logical(value_include_zero) || (length(value_include_zero) != 1)) stop(
#       "value_include_zero must be a logical (TRUE or FALSE) of length 1") 
#     
#     # Check check if value_plot is logical of length 1 #
#     if(!is.logical(value_plot) || (length(value_plot) != 1)) stop(
#       "value_plot must be a logical (TRUE or FALSE) of length 1") 
#     
#     
#     # Check if y_limits or y_range have been selected correctly #
#     if((!is.null(comps_y_limits) & !is.null(comps_y_range)) | 
#        (!is.null(value_y_limits) & !is.null(value_y_range))) stop(
#          "Input either y_limits or y_range parameters, but not both.")
#     
#     # --Comps-- #
#     # Check if only one of comps_y_max and comps_y_min has been selected with comps_y_limits or comps_y_range #
#     if(!is.null(comps_y_max) & 
#        !is.null(comps_y_min) & 
#        (!is.null(comps_y_range) | !is.null(comps_y_limits))) stop(
#          "Cannot use both comps_y_min and comps_y_max with comps_y_range 
#        or comps_y_limits parameters. Only one of comps_y_min or 
#        comps_y_max can be used.")
#     
#     # Check if comps_y_limits is in acceptable strings and length == 1 #
#     if ( !is.null(comps_y_limits)){
#       if((length(comps_y_limits) != 1)) stop(
#         "Parameter y_limits must have length = 1.")
#       if(!(comps_y_limits %in% c("fixed", "free"))) stop(
#         "Parameter y_limits must be input as either 'fixed' or 'free'.")
#     }
#     
#     # Check if comps_y_range is positive, numeric and length == 1 #
#     if (!is.null(comps_y_range)){
#       if(!is.numeric(comps_y_range)) stop(
#         "Parameter y_range must be numeric.")
#       if(length(comps_y_range) != 1) stop(
#         "Parameter y_range must have length = 1.")
#       if(!(comps_y_range > 0)) stop(
#         "Parameter y_range must be greater than zero.")
#     }
#     
#     # Check if comps_y_max is numeric and length == 1 #
#     if (!is.null(comps_y_max)){
#       if(!is.numeric(comps_y_max)) stop(
#         "Parameter y_max must be numeric.")
#       if(length(comps_y_max) != 1) stop(
#         "Parameter y_max must have length = 1.")
#     }
#     
#     # Check if comps_y_min is numeric and length == 1 #
#     if (!is.null(comps_y_min)){
#       if(!is.numeric(comps_y_min)) stop(
#         "Parameter y_min must be numeric.")
#       if(length(comps_y_min) != 1) stop(
#         "Parameter y_min must have length = 1.")
#     }
#     
#     # Check if comps_plot_type has one of the available options #
#     checkplot <- purrr::map(comps_plot_type, 
#                             function(plot) stringr::str_detect(plot, "box|col|point|scatter|bar"))
#     if (any(!unlist(checkplot))) stop(
#       "Invalid entry in comps_plot_type. Plot_type strings must contain 
#     at least one of the following: box, col, point, scatter, bar")
#     
#     # --Value-- #
#     # Check if only one of value_y_max and value_y_min has been selected with value_y_limits or value_y_range #
#     if(!is.null(value_y_max) & 
#        !is.null(value_y_min) & 
#        (!is.null(value_y_range) | !is.null(value_y_limits))) stop(
#          "Cannot use both value_y_min and value_y_max with value_y_range 
#        or value_y_limits parameters. Only one of y_min or y_max can be 
#        used.")
#     
#     # Check if value_y_limits is in acceptable strings and length == 1 #
#     if (!is.null(value_y_limits) ){
#       if((length(value_y_limits) != 1)) stop(
#         "Parameter y_limits must have length = 1.")
#       if(!(value_y_limits %in% c("fixed", "free")) ) stop(
#         "Parameter y_limits must be input as either 'fixed' or 'free'.")
#     }
#     
#     # Check if value_y_range is positive, numeric and length == 1 #
#     if (!is.null(value_y_range)){
#       if(!is.numeric(value_y_range)) stop(
#         "Parameter y_range must be numeric.")
#       if(length(value_y_range) != 1) stop(
#         "Parameter y_range must have length = 1.")
#       if(!(value_y_range > 0)) stop(
#         "Parameter y_range must be greater than zero.")
#     }
#     
#     # Check if value_y_max is numeric and length == 1 #
#     if (!is.null(value_y_max)){
#       if(!is.numeric(value_y_max)) stop(
#         "Parameter y_max must be numeric.")
#       if(length(value_y_max) != 1) stop(
#         "Parameter y_max must have length = 1.")
#     }
#     
#     # Check if value_y_min is numeric and length == 1 #
#     if (!is.null(value_y_min)){
#       if(!is.numeric(value_y_min)) stop(
#         "Parameter y_min must be numeric.")
#       if(length(value_y_min) != 1) stop(
#         "Parameter y_min must have length = 1.")
#     }
#     
#     # Check if value_plot_type has one of the available options #
#     checkplot <- purrr::map(value_plot_type, 
#                             function(plot) stringr::str_detect(plot, "box|col|point|scatter|bar"))
#     if (any(!unlist(checkplot))) stop(
#       "Invalid entry in value_plot_type. Plot_type strings must contain 
#     at least one of the following: box, col, point, scatter, bar")
#     
#   } else {
#     
#     ##### Post-Moniker and null assignment check ####
#     # --Value-- #
#     if (!is.null(trellData$data_values) && !is.null(value_panel_x_axis)) {
#       
#       # print(value_panel_x_axis)
#       # print(value_panel_y_axis)
#       # print(panel_variable)
#       # print(value_color_variable)
#       
#       # Ensure panel, color/x/y parameters are not matching #
#       if (!((value_panel_x_axis != panel_variable) && 
#             (value_panel_y_axis != panel_variable) && 
#             (panel_variable != value_color_variable))){
#         stop("Parameters for value panel_y_axis, panel_x_axis, and color_variable 
#            cannot match panel_variable. Refer to ?plot_comps for default settings 
#            or try setting all of these parameters individually.")
#       }
#       if(value_panel_x_axis == value_panel_y_axis) warning(
#         "Parameter for value panel_y_axis and panel_x_axis are identical. 
#             Refer to ?plot_comps for default settings or try setting 
#             parameters individually for different axis labels.")
#       
#       # Variable for variable-in checks #
#       allcol <- colnames(trellData$data_values)
#       allcolstr <- toString(allcol)
#       
#       # Ensure panel, color, x, and y parameters are in data_values #
#       if (!(value_panel_x_axis %in% allcol)) stop(
#         paste("Parameter value_panel_x_axis  must be in:", allcolstr))
#       if (!(value_panel_y_axis %in% allcol)) stop(
#         paste("Parameter value_panel_y_axis must  must be in:", allcolstr))
#       if (!(panel_variable %in% allcol)) stop(
#         paste("Parameter panel_variable must be in:", allcolstr))
#       if (!(value_color_variable %in% allcol)) stop(
#         paste("Parameter value_color_variable must be in:", allcolstr))
#       
#     }
#     
#     # --Comps-- #
#     
#     # Ensure summary_stats and comp_stats are both present #
#     if(((is.null(trellData$summary_stats)) & 
#         (!is.null(trellData$comp_stats))) | 
#        ((!is.null(trellData$summary_stats)) &
#         (is.null(trellData$comp_stats))) ) {
#       stop("Both summary_stats and comp_stats must be present 
#          if one or the other is in trellData.")
#     }
#     
#     if ((!is.null(trellData$summary_stats)) && 
#         (!is.null(trellData$comp_stats)) &&
#         (!is.null(comps_panel_x_axis)) ) {
#       
#       # Check if stats statistical test attribute is valid #
#       if(!(
#         attributes(trellData)$statistical_test %in% 
#         c("combined", "gtest", "anova"))) stop(
#         paste("Non-applicable statistical_test attribute in trellData object."))
#       
#       # Ensure panel, color/x/y parameters are not matching #
#       if(comps_panel_x_axis == comps_panel_y_axis) warning(
#         "Parameter for comps panel_y_axis and panel_x_axis are identical.
#     Refer to ?plot_comps for default settings or try setting parameters
#     individually for different axis labels.")
#       
#       if (!((comps_panel_x_axis != panel_variable) && 
#             (comps_panel_y_axis != panel_variable) && 
#             (panel_variable != comps_color_variable))) stop(
#               "Parameters for comps panel_y_axis, panel_x_axis, 
#             and color_variable cannot match panel_variable. Refer to ?plot_comps 
#             for default settings or try setting all of these parameters 
#             individually.")
#       
#       # Variable for variable-in checks #
#       allcol <- c(colnames(trellData$comp_stats), 
#                   colnames(trellData$summary_stats))
#       allcolstr <- toString(unique(allcol))
#       
#       # Ensure panel, color, x, and y parameters are in comp_stats or summary stats #
#       if (!(comps_panel_x_axis %in% allcol)) stop(
#         paste("Parameter comps_panel_x_axis  must be in:", allcolstr))
#       if (!(comps_panel_y_axis %in% allcol)) stop(
#         paste("Parameter comps_panel_y_axis must  must be in:", allcolstr))
#       if (!(panel_variable %in% allcol)) stop(
#         paste("Parameter panel_variable must be in:", allcolstr))
#       if (!(comps_color_variable %in% allcol)) stop(
#         paste("Parameter comps_color_variable must be in:", allcolstr))
#       
#     }
#   }
# }
#   

#' @name generate_plot_message
#' @rdname generate_plot_message
#' @title Generates a message stating the specified plotting y-limits
#' 
#' @description Generates a message stating the specified plotting y-limits
#'
#' @param trellData An object of class "trellData" generated from \code{\link{as.trellData}}.
#' @param comps_y_limits For comparisons: Set to "fixed" or "free" for automated y-axis calculating. "fixed" - axis generated based on the maximum/minimum across all plots. "free" - axis axis generated based on the maximum/minimum of individual plot.
#' @param comps_y_range For comparisons: Specify a range for the plot y-axis. Will calculate the range based on one of y_max or y_min parameters or from the median of y-values where y_max and y_min are not defined.
#' @param comps_y_max For comparisons: Sets the maximum y-value for the y-axis.
#' @param comps_y_min For comparisons: Sets the minimum y-value for the y-axis.
#' @param value_y_limits For values: Set to "fixed" or "free" for automated y-axis calculating. "fixed" - axis generated based on the maximum/minimum across all plots. "free" - axis axis generated based on the maximum/minimum of individual plot.
#' @param value_y_range For values: Specify a range for the plot y-axis. Will calculate the range based on one of y_max or y_min parameters or from the median of y-values where y_max and y_min are not defined.
#' @param value_y_max For values: Sets the maximum y-value for the y-axis.
#' @param value_y_min For values: Sets the minimum y-value for the y-axis.
#'
#' @author Rachel Richardson
#'

generate_plot_message <- function(...){
  .generate_plot_message(...)
}

.generate_plot_message <- function(trellData, 
                                   comps_y_limits = NULL, comps_y_range = NULL,
                                   comps_y_max = NULL, comps_y_min = NULL, 
                                   value_y_limits = NULL, value_y_range = NULL, 
                                   value_y_max = NULL, value_y_min = NULL){
  
  ## Input comps y limits messages, tells user the plot limit parameters ##
  #--Comps--#
  if(!is.null(trellData$summary_stats) & !is.null(trellData$comp_stats)){
    ## Input comps y limits messages, tells user the plot limit parameters ##
    if (is.null(comps_y_limits) & is.null(comps_y_range) & is.null(comps_y_max) & is.null(comps_y_min)){
      message("No specified comparison y-axis parameters. Axis y-limits will be scaled per plot, as per y_limits = 'free'.")
      comps_y_limits <- "free"
    } else if (!is.null(comps_y_limits)){
      if ((comps_y_limits == "fixed") & is.null(comps_y_max) & is.null(comps_y_min)){
        message("Specified comparison y-limit: 'fixed'. Axis y-limits will fixed for all plots based on maximum and minimum y-values.")
      } else if ((comps_y_limits == "fixed") & !is.null(comps_y_max)){
        message(paste("Specified comparison y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a maximum of y_max. Specified y_max: ", comps_y_max, sep = ""))
      } else if ((comps_y_limits == "fixed") & !is.null(comps_y_min)){
        message(paste("Specified comparison y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a minimum of y_min. Specified y_min: ", comps_y_min, sep = ""))
      } else if (comps_y_limits == "free" & is.null(comps_y_max) & is.null(comps_y_min)){
        message("Specified comparison y-limit: 'free'. Axis y-limits will be scaled per plot.")
      } else if ((comps_y_limits == "free") & !is.null(comps_y_max)){
        message(paste("Specified comparison y-limit: 'free'. Axis y-limits will be scaled per plot with a maximum of y_max. Specified y_max: ", comps_y_max, sep = ""))
      } else if ((comps_y_limits == "free") & !is.null(comps_y_min)){
        message(paste("Specified comparison y-limit: 'free'. Axis y-limits will be scaled per plot with a minimum of y_min. Specified y_min: ", comps_y_min, sep = ""))
      } 
    } else if (!is.null(comps_y_range) & is.null(comps_y_max) & is.null(comps_y_min)){
      message(paste("Specified comparison y-range: ", comps_y_range, ". Axis y-limits will range ",
                    comps_y_range," units, split over the median.", 
                    sep = ""))
    } else if (!is.null(comps_y_range) & !is.null(comps_y_min)) {
      message(paste("Specified comparison y-range: ", comps_y_range, ". Axis y-limits will range ",
                    comps_y_range," units from the y_min. Specified y_min: ", comps_y_min, sep = ""))
    } else if (!is.null(comps_y_range) & !is.null(comps_y_max)) {
      message(paste("Specified comparison y-range: ", comps_y_range, ". Axis y-limits will range ",
                    comps_y_range," units from the y_max. Specified y_max: ", comps_y_max, sep = ""))
    } else if (is.null(comps_y_min) && is.null(comps_y_range) && 
               is.null(comps_y_limits) && !is.null(comps_y_max)){
      message(paste("Specified comparison y-max: ", 
                    comps_y_max,
                    ". No range or limits specified; Axis y-limits will be scaled per plot with a maximum of y_max."))
      comps_y_limits <- "free"
    } else if (!is.null(comps_y_min) && is.null(comps_y_range) && 
               is.null(comps_y_limits) && is.null(comps_y_max)){
      message(paste("Specified comparison y-min: ", 
                    comps_y_max,
                    ". No range or limits specified; Axis y-limits will be scaled per plot with a minimum of y_min."))
      comps_y_limits <- "free"
    }
  }
  
  #--Values--#
  if(!is.null(trellData$data_values)){
    ## Input values y limits messages, tells user the plot limit parameters ##
    if (is.null(value_y_limits) & is.null(value_y_range) & is.null(value_y_max) & is.null(value_y_min)){
      message("No specified value y-axis parameters. Axis y-limits will be scaled per plot, as per y_limits = 'free'.")
      value_y_limits <- "free"
    } else if (!is.null(value_y_limits)){
      if ((value_y_limits == "fixed") & is.null(value_y_max) & is.null(value_y_min)){
        message("Specified value y-limit: 'fixed'. Axis y-limits will fixed for all plots based on maximum and minimum y-values.")
      } else if ((value_y_limits == "fixed") & !is.null(value_y_max)){
        message(paste("Specified value y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a maximum of y_max. Specified y_max: ", value_y_max, sep = ""))
      } else if ((value_y_limits == "fixed") & !is.null(value_y_min)){
        message(paste("Specified value y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a minimum of y_min. Specified y_min: ", value_y_min, sep = ""))
      } else if (value_y_limits == "free" & is.null(value_y_max) & is.null(value_y_min)){
        message("Specified value y-limit: 'free'. Axis y-limits will be scaled per plot.")
      } else if ((value_y_limits == "free") & !is.null(value_y_max)){
        message(paste("Specified value y-limit: 'free'. Axis y-limits will be scaled per plot with a maximum of y_max. Specified y_max: ", value_y_max, sep = ""))
      } else if ((value_y_limits == "free") & !is.null(value_y_min)){
        message(paste("Specified value y-limit: 'free'. Axis y-limits will be scaled per plot with a minimum of y_min. Specified y_min: ", value_y_min, sep = ""))
      } 
    } else if (!is.null(value_y_range) & is.null(value_y_max) & is.null(value_y_min)){
      message(paste("Specified value y-range: ", value_y_range, ". Axis y-limits will range ",
                    value_y_range," units, split over the median.", 
                    sep = ""))
    } else if (!is.null(value_y_range) & !is.null(value_y_min)) {
      message(paste("Specified value y-range: ", value_y_range, ". Axis y-limits will range ",
                    value_y_range," units from the y_min. Specified y_min: ", value_y_min, sep = ""))
    } else if (!is.null(value_y_range) & !is.null(value_y_max)) {
      message(paste("Specified value y-range: ", value_y_range, ". Axis y-limits will range ",
                    value_y_range," units from the y_max. Specified y_max: ", value_y_max, sep = ""))
    } else if (is.null(value_y_min) && is.null(value_y_range) && 
               is.null(value_y_limits) && !is.null(value_y_max)){
      message(paste("Specified comparison y-max: ", 
                    value_y_max,
                    ". No range or limits specified; Axis y-limits will be scaled per plot with a maximum of y_max."))
      value_y_limits <- "free"
    } else if (!is.null(value_y_min) && is.null(value_y_range) && 
               is.null(value_y_limits) && is.null(value_y_max)){
      message(paste("Specified comparison y-min: ", 
                    value_y_max,
                    ". No range or limits specified; Axis y-limits will be scaled per plot with a minimum of y_min."))
      value_y_limits <- "free"
    }
  }
  
  return(list(comps_y_limits , value_y_limits))
  
}


#' @name add_plot_features
#' @rdname add_plot_features
#' @title Adds features to dataframe for plotting
#' 
#' @description Adds border color, hover, and label text as applicable to nested trellData dataframe for plotting
#'
#' @param trellData An object of class "trellData" generated from \code{\link{as.trellData}}.
#' @param nestedData Nested trellData dataframe
#' @param p_val Specifies p-value for setting graph border colors
#' @param panel_variable Specifies what to divide trelliscope panels by, must be a column in trellData. 
#' @param comps_panel_y_axis Specifies what column should be on the y-axis, must be a column in trellData and numeric. 
#' @param value_panel_y_axis Specifies what column should be on the y-axis, must be a column in trellData and numeric. 
#'
#' @author Rachel Richardson
#'
#' @export
#'

add_plot_features <- function(...){
  .add_plot_features(...)
}

.add_plot_features <- function(trellData,
                               nestedData,
                               p_val = NULL, 
                               panel_variable = NULL,
                               comps_panel_y_axis = NULL,
                               value_panel_y_axis = NULL){
  
  colors <- c(NA, "grey40", "black")
  
  if(is.null(value_panel_y_axis)){
    
    # Set border colors based on significance #
    bord <- rep(colors[1], nrow(nestedData))
    if("P_value_G" %in% colnames(nestedData)){
      bord[nestedData$P_value_G < p_val & !is.na(nestedData$P_value_G)] <-  colors[2]
    }
    if("P_value_T" %in% colnames(nestedData)){
      bord[nestedData$P_value_T < p_val & !is.na(nestedData$P_value_T)] <- colors[3]
    }
    
    # Set hover/labels excluding the panel_variable #
    hover_want_comps <- c(attributes(trellData)$cnames$edata_cname,
                          "Group", "Count", "Comparison",
                          "P_value_G", "P_value_T", "Fold_change") 
    if (!(comps_panel_y_axis %in% hover_want_comps)){
      hover_want_comps <- c(hover_want_comps, comps_panel_y_axis)
    }
    if (panel_variable %in% hover_want_comps){
      hover_want_comps <- hover_want_comps[!(hover_want_comps %in% panel_variable)]
    }
    hover_labs_comps <- hover_want_comps[hover_want_comps %in% colnames(nestedData)]
    
    text_labs_comps <- purrr::map(1:nrow(nestedData), function(row){
      row_text <- purrr::map(hover_labs_comps, function(label){
        if (label %in% c("P_value_G", "P_value_T")) {
          return(paste(paste(label, ":", sep = ""), 
                       signif(nestedData[row,][[label]])))
        } else if (label %in% c("Fold_change", comps_panel_y_axis)){
          return(paste(paste(label, ":", sep = ""), 
                       signif(nestedData[row,][[label]])))
        } else {
          return(paste(paste(label, ":", sep = ""), 
                       nestedData[row,][[label]]))
        }
      })
      return(toString(row_text, sep = ", "))
    })
    
    label_labs_comps <- purrr::map(1:nrow(nestedData), function(row){
      row_label <- purrr::map(hover_labs_comps, function(label){
        if (label %in% c("P_value_G", "P_value_T")){
          return(paste(paste(label, ":", sep = ""), 
                       signif(nestedData[row,][[label]], 3)))
        } 
      })
      return(toString(row_label, sep = ", "))
    })
    
    label_labs_comps <- stringr::str_remove_all(unlist(label_labs_comps), "NULL, ")
    label_labs_comps <- stringr::str_remove_all(label_labs_comps, "NULL")
    
    nestedData_comps <- data.frame(nestedData,
                                   bord = bord,
                                   text = unlist(text_labs_comps),
                                   labels = label_labs_comps)
    
    return(nestedData_comps)
    
  } else {
    
    # Set hover excluding the panel_variable #
    hover_want <- c(attributes(trellData)$cnames$edata_cname, 
                    attributes(trellData)$cnames$fdata_cname,
                    "Group", grep("abundance",
                                  colnames(trellData$data_values),
                                  value = TRUE)) 
    if (!(value_panel_y_axis %in% hover_want)){
      hover_want <- c(hover_want, value_panel_y_axis)
    }
    if (panel_variable %in% hover_want){
      hover_want <- hover_want[!(hover_want %in% panel_variable)]
    }
    hover_labs <- hover_want[hover_want %in% colnames(nestedData)]
    
    text_labs <- purrr::map(1:nrow(nestedData), function(row){
      row_text <- purrr::map(hover_labs, function(label){
        if (label %in% c(grep("abundance", 
                              colnames(trellData$data_values), 
                              value = TRUE), value_panel_y_axis)) {
          return(paste(paste(label, ":", sep = ""),
                signif(nestedData[row,][[label]])))
        } else {
          return(paste(paste(label, ":", sep = ""), nestedData[row,][[label]]))
        }
      })
      
      return(toString(row_text, sep = ", "))
    })
    nestedData_value <- data.frame(nestedData, text = unlist(text_labs))
    
    return(nestedData_value)
  }
  
}


#' @name list_y_limits
#' @rdname list_y_limits
#' @title Digests user input for y_limits
#' 
#' @description Takes lists/strings/numeric input and digests as a list of y_limits for each plot type.
#'
#' @param plot_type plots to graph in trelliVis
#' @param y_limits User input y-axis limits for plots
#'
#' @author Rachel Richardson
#'
#' @export

list_y_limits <- function(plot_type, y_limits){
  .list_y_limits(plot_type, y_limits)
}

.list_y_limits <- function(plot_type, y_limits){
    
    limlist <- list()
    plotlims <- list()
    plotlims2 <- list()
    
    ## Set if null
    if(is.null(y_limits)){
      y_limits <- "fixed"
    }
    
    ## Find and assign any list elements
    lists <- purrr::map(y_limits, function(elements){
      is.list(elements)
    })
    
    ## Catch for one item list
    if(is.list(y_limits) && 
       !any(purrr::map_lgl(1:length(y_limits), function(item) is.list(y_limits[[item]])))){
      if (is.null(names(y_limits))) stop(
        "Invalid format. List inputs must be named. Refer to ?trelliVis for examples."
      )
      if (any(duplicated(names(y_limits)))) stop(
        paste("List inputs may not have duplicate names. Use a list of lists to specify different plot y_limits. Duplicate names:",
              names(y_limits)[duplicated(names(y_limits))])
      )
      if (!all(names(y_limits) %in% c("min", "max", "range", "scale"))) stop(
        "List names for y_limits are restricted to 'min', 'max', 'range', and 'scale'"
      )
      
      plotlims[1:length(plot_type)] <- list(y_limits)
    }
    
    ###### Checks for singular input
    if (length(y_limits) == 1 && (y_limits == "fixed" || y_limits == "free")){
      limlist[["scale"]] <- y_limits
      plotlims[1:length(plot_type)] <- list(limlist)
      
    } else if (is.numeric(y_limits) && length(y_limits) == 2){
      limlist[["min"]] <- min(y_limits)
      limlist[["max"]] <- max(y_limits)
      plotlims[1:length(plot_type)] <- list(limlist)
      
      #### Handles lists of y_limit input
    } else if (any(unlist(lists))){
      
      ##Check list correctness
      if (any(purrr::map_lgl(y_limits, function(elements){
        any(purrr::map_lgl(elements, function(item2) is.list(item2)))
      }))) stop(
        "List of lists of lists are not supported in y_limits.")
      
      if (length(y_limits) != length(plot_type)) stop(
        "Lists in y_limits must be the same length as plot_type. Refer to examples in ?trelliVis().")
      
      if (any(duplicated(names(y_limits)))) stop(
        paste("List inputs may not have duplicate names. Use a list of lists to specify different plot y_limits. Duplicate names:",
              names(y_limits)[duplicated(names(y_limits))])
      )  
      
      if (any(!(names(y_limits) %in% plot_type))) stop(
        paste("Names in y_limits should match plot_type. plot_type:",
              paste(plot_type, collapse = ", "))
      )
      
      y_limits <- y_limits[plot_type]
      plotlims <- y_limits[unlist(lists)]
      
      if(!is.null(plotlims) && length(plotlims) > 1){
        purrr::map(plotlims, function(listitem){
          if (length(listitem) > 1 && is.null(names(listitem))) stop(
            "Invalid format. List inputs must be named. Refer to ?trelliVis for examples."
          )
          if (any(duplicated(names(listitem)))) stop(
            paste("List inputs may not have duplicate names. Use a list of lists to specify different plot y_limits. Duplicate names:",
                  names(listitem)[duplicated(names(listitem))])
          )
          
          if (!all(names(listitem) %in% c("min", "max", "range", "scale"))) stop(
            "List names for y_limits are restricted to 'min', 'max', 'range', and 'scale'"
          )
          
        })
        
      } else if (!is.null(plotlims)){
        if (is.null(names(plotlims[[1]]))) stop(
          "Invalid format. List inputs must be named. Refer to ?trelliVis for examples."
        )
        
        if (any(duplicated(names(plotlims[[1]])))) stop(
          paste("List inputs may not have duplicate names. Use a list of lists to specify different plot y_limits. Duplicate names:",
                names(plotlims[[1]])[duplicated(names(plotlims[[1]]))])
        )
        
        if (!all(names(plotlims[[1]]) %in% c("min", "max", "range", "scale"))) stop(
          "List names for y_limits are restricted to 'min', 'max', 'range', and 'scale'"
        )
        
      }
      
      
      
      ## Find and assign non-list elements
      plotlims2 <- purrr::map(y_limits[!unlist(lists)], function(limit){
        
        listlims2 <- list()
        if (limit == "fixed" || limit == "free"){
          listlims2$scale <- limit
        } else if (is.numeric(limit) && length(limit) == 2){
          listlims2$min <- limit[1]
          listlims2$max <- limit[2]
        }
        return(listlims2)
      })
    }
    
    if(length(plotlims) == 0 && length(plotlims2) == 0) stop(
      "Invalid input. y_limits input is restricted to named lists for each plot type, strings 'fixed' or 'free', or a numeric vector of length 2 specifying maximum and minimum y_values. Refer to ?trelliVis for examples."
    )
    
    output <- c(plotlims, plotlims2)
    nms <- c(plot_type[unlist(lists)], plot_type[!unlist(lists)])
    nms <- nms[!is.na(nms)]
    names(output) <- nms
    
    return(output)
  }



#' @name set_increment
#' @rdname set_increment
#' @title Sets y-axis increment for trellData labels plotting
#' 
#' @description Sets y-axis increment for trellData plotting. Used by plot_comp.
#'
#' @param yvalues y-values for plotting
#' @param testtype consideration for different statistical tests
#'
#' @author Rachel Richardson
#'
#' @export
set_increment <- function(yvalues, ...){
  .set_increment(yvalues,  ...)
}

.set_increment <- function(yvalues, include_zero = TRUE){
  ## Initial Checks and Replacement ##
  
  # Remove NA, duplicate values #
  yvalues <- unique(yvalues[!is.na(yvalues)])
  if (length(yvalues) == 0){
    yvalues <- 0
  }
  
  # Add zero if including zero #
  if (include_zero == TRUE){
    yvalues <- c(yvalues, 0)
  }
  
  # Check if yvalues is a numeric vector #
  if(!is.vector(yvalues) || !inherits(yvalues, "numeric")) stop(
    "yvalues must be a numeric vector")
  
  # Check if include_zero is logical #
  if(!is.logical(include_zero) || !(length(include_zero) == 1) ) stop(
    "include_zero must be a length 1 logical. (TRUE/FALSE)")
  
  if(length(yvalues)==1){
    increment <- yvalues/20
  } else {
    increment <- (max(yvalues) - min(yvalues))/20
  }
  
  return(abs(increment))
}

#' @name set_ylimits
#' @rdname set_ylimits
#' @title Sets y-axis limits for trellData plotting
#' 
#' @description Sets y-axis limits for trellData plotting. Used by plot_comp.
#'
#' @param yvalues y-values for plotting
#' @param increment An increment set based on the maximum/minimum y-values
#' @param y_range Specify a range for the plot y-axis. Will calculate the range based on one of y_max or y_min parameters or from the median of y-values where y_max and y_min are not defined.
#' @param y_max Sets the maximum y-value for the y-axis.
#' @param y_min Sets the minimum y-value for the y-axis.
#'
#' @author Rachel Richardson
#'
#' @export
#' 
set_ylimits <- function(yvalues, increment, ...){
  .set_ylimits(yvalues, increment, ...)
}

.set_ylimits <- function(yvalues, increment, y_range = NULL, 
                         y_max = NULL, y_min = NULL, include_zero = TRUE){
  
  # Catch NAs #
  yvalues <- yvalues[!is.na(yvalues)]
  if (length(yvalues) == 0){
    yvalues <- 0
  }
  
  ## Initial Checks and Replacement ##
  
  # Check if yvalues is numeric vector #
  if(!is.vector(yvalues) || !inherits(yvalues, "numeric")) stop(
    "yvalues must be a numeric vector")
  
  # Check if increment is numeric length 1 #
  if(!is.numeric(increment) || !(length(increment) == 1)) stop(
    "increment must be a length 1 numeric")
  
  # Check if y_range is numeric length 1 #
  if(!is.null(y_range) && 
     (!inherits(y_range, "numeric") || !(length(y_range) == 1))) stop(
       "y_range must be a length 1 numeric")
  # Check if y_max is numeric length 1 #
  if(!is.null(y_max) &&
     (!inherits(y_max, "numeric") || !(length(y_max) == 1))) stop(
       "y_max must be a length 1 numeric")
  # Check if y_min is numeric length 1 #
  if(!is.null(y_min) &&
     (!inherits(y_min, "numeric") || !(length(y_min) == 1))) stop(
       "y_min must be a length 1 numeric")
  
  if(!is.null(y_min) && !is.null(y_max) && !is.null(y_range)) stop(
    "y_range must be NULL when y_max and y_min are assigned.")
  
  # Check if include_zero is logical #
  if(!is.logical(include_zero) || !(length(include_zero) == 1)) stop(
    "include_zero must be a length 1 logical. (TRUE/FALSE)")
  
  
  ## Set Limits ##
  # Catch for pre-set y-limits and y_range #
  if (!is.null(y_min) && !is.null(y_max)){
    maxi <- y_max
    mini <- y_min
    
    if(maxi < mini) stop ("Invalid max and min. Max < Min")
    
    return(c(mini, maxi))
  } else if( !is.null(y_range) && !is.null(y_min)){
    mini <- y_min
    maxi <- y_min + y_range
    return(c(mini, maxi))
  } else if( !is.null(y_range) && !is.null(y_max)){
    maxi <- y_max
    mini <- y_max - y_range
    return(c(mini, maxi))
  } else if (!is.null(y_range)){
    maxi <- median(yvalues) + y_range/2
    mini <- median(yvalues) - y_range/2
    return(c(mini, maxi))
  }
  
  # Set maxi and mini based on difference between max and min values #
  maxi <- max(yvalues) + 3*increment
  mini <- min(yvalues) - 3*increment
  
  # Adjust for maximums below 0 and minimums above 0 #
  if (!(maxi > 0) && (include_zero == TRUE)){
    maxi <- 3*increment
  }
  if (!(mini < 0) && (include_zero == TRUE)){
    mini <- -3*increment
  }
  
  # Adjust for specified y_min and y_max #
  if (!is.null(y_max)){
    maxi <- y_max
  }
  if (!is.null(y_min)){
    mini <- y_min
  }
  
  # Adjust for both 0 (where increment is zero and yvalues are 0) #
  if ((mini == 0) && (maxi == 0)){
    mini <- -1
    maxi  <- 1
  }
  
  ## Return limits ##
  # Minimum y value = mini, maximum y value = maxi #
  return(c(mini, maxi))
}

#' @name data_cogs
#' @rdname data_cogs
#' @title Plot pairwise comparisons and data values in trellData object
#'
#' @description Plot pairwise comparisons and data values in trellData object. Customizable for plot types, y axis limits, paneling variable (what overall group is plotted on each graph arrangement), as well as desired variables for the y and x axis.
#'
#' @param nested_plot A nested table generated from trellData using formatplot()
#' @param trellData A nested table generated from trellData using formatplot()
#' @param p_val Numeric that specifies p-value for Boolean significance cognotic. Default is 0.05.
#' @param try_URL Will attempt to link to PubChem, LipidMaps, or Uniprot based on information in edata_cname of omicsData or specified mapping_col for peptide data. Default is FALSE.
#'
#'
#' @author Rachel Richardson
#' @export


data_cogs <- function(...) {
  .data_cogs( ...)
}

.data_cogs <- function(nested_plot,
                       trellData,
                       p_val = 0.05,
                       try_URL = FALSE){
  
  ## Variable checks ##
  
  # Check check if p_val is numeric of length 1 #
  if(!is.numeric(p_val) | (length(p_val) != 1)) stop(
    "p_val must be a numeric of length 1")  
  
  # Check check if try_URL is boolean of length 1 #
  if(!is.logical(try_URL) | (length(try_URL) != 1)) stop(
    "try_URL must be a TRUE/FALSE of length 1")  
  
  ## Assign unique ID, then panel variable ##
  uniqueID <- pmartR::get_edata_cname(trellData)
  panel_variable <- colnames(nested_plot)[1]
  parent_dat <- attr(trellData, "parent_class")
  
  
  if (!is.null(trellData$comp_stats)){
    stats <- TRUE
    trell_comp <- trellData$comp_stats
    trell_summ <- trellData$summary_stats
  } else {
    stats <- FALSE
  }
  if (!is.null(trellData$data_values)){
    values <- TRUE
    trell_values <- trellData$data_values
  } else {
    values <- FALSE
  }
  
  if("pepData" %in% parent_dat && 
     !is.null(pmartR::get_emeta_cname(trellData)) &&
     values){
    peppro <- unique(trell_values[c(uniqueID, 
                                    pmartR::get_emeta_cname(trellData))])
    if(any(duplicated(peppro[[uniqueID]]))){
      degen <- peppro[[uniqueID]][duplicated(peppro[[uniqueID]])]
    } else {
      degen <- NULL
    }
  }
  
  out <- dplyr::mutate(
    nested_plot,
    cogs = trelliscopejs::map_cog(
      as.character(nested_plot[[panel_variable]]),
      function(panel){
        
        # 
        if(!is.null(pmartR::get_emeta_cname(trellData))){
          joiner <- c(uniqueID, pmartR::get_emeta_cname(trellData))
        } else {
          joiner <- uniqueID
        }
        
        ## Generate slice trelldata into panel rows and create merged df of all dfs
        if (values && stats){
          panel_values <- trell_values[as.character(trell_values[[panel_variable]]) == panel,]
          panel_comp <- trell_comp[as.character(trell_comp[[panel_variable]]) == panel,]
          panel_summ <- trell_summ[as.character(trell_summ[[panel_variable]]) == panel,]
          
          byvc  <- colnames(panel_values)[colnames(panel_values) %in% colnames(panel_comp)]
          bycs  <- colnames(panel_comp)[colnames(panel_comp) %in% colnames(panel_summ)]
          
          addOrigCogs <- dplyr::left_join(panel_values, 
                                          panel_comp, 
                                          by = byvc)
          addOrigCogs <- dplyr::left_join(addOrigCogs,
                                          panel_summ,
                                          by = unique(c(byvc, bycs, "Group")))
          
        } else if (!is.null(trellData$data_values)){
          panel_values <- trell_values[trell_values[[panel_variable]] == panel,]
          addOrigCogs <- panel_values
          
        } else {
          panel_comp <- trell_comp[trell_comp[[panel_variable]] == panel,]
          panel_summ <- trell_summ[trell_summ[[panel_variable]] == panel,]
          bycs  <- colnames(panel_comp)[colnames(panel_comp) %in% colnames(panel_summ)]
          
          addOrigCogs <- dplyr::left_join(panel_comp,
                                          panel_summ,
                                          by = unique(c(bycs, "Group")))
        }
        cogs <- addOrigCogs
        
        if(stats){
          stat_test <- attr(trellData, "statistical_test")
          if(stat_test == "combined"){
            cogs <- dplyr::mutate(
              cogs, 
              Sig_T_p_value = purrr::map_lgl(
                cogs[["P_value_T"]],
                function(value) !(value > p_val) && !is.na(value)),
              Sig_G_p_value = purrr::map_lgl(
                cogs[["P_value_G"]],
                function(value) !(value > p_val) && !is.na(value)))
          } else if (stat_test == "anova"){
            cogs <- dplyr::mutate(
              cogs, 
              Sig_T_p_value = purrr::map_lgl(
                cogs[["P_value_T"]],
                function(value) !(value > p_val) && !is.na(value)))
          } else {
            cogs <- dplyr::mutate(
              cogs, 
              Sig_G_p_value = purrr::map_lgl(
                cogs[["P_value_G"]],
                function(value) !(value > p_val) && !is.na(value)))
          }
        }
        
        if("pepData" %in% parent_dat){
          if(!is.null(degen)){
            cogs <- dplyr::mutate(
              cogs,
              Degenerate_peptide = purrr::map_lgl(
                cogs[[uniqueID]],
                function(peptide) peptide %in% degen))
          }
          
          if(!is.null(try_URL) && 
             !is.null(pmartR::get_emeta_cname(trellData))){
            cogs <- dplyr::mutate(
              cogs,
              Protein_URL = purrr::map_chr(
                cogs[[pmartR::get_emeta_cname(trellData)]],
                function(protein){
                  searchname <- stringr::str_extract(protein, "[A-Z0-9]+_[A-Z]+")
                  if (is.na(searchname)){
                    searchname <- stringr::str_extract(protein, "[A-Z0-9]{6,}")
                  }
                  if (!is.na(searchname)){
                    searchlink <- paste('https://www.uniprot.org/uniprot/', 
                                        searchname, sep = "")
                    return(searchlink)
                  } else {
                    return("Could not extract protein name")
                  }
                }))
          }
        } else if ("proData" %in% parent_dat){
          if(!is.null(try_URL)){
            cogs <- dplyr::mutate(
              cogs,
              Protein_URL = purrr::map_chr(
                cogs[[uniqueID]],
                function(protein){
                  searchname <- stringr::str_extract(protein, "[A-Z0-9]+_[A-Z]+")
                  if (is.na(searchname)){
                    searchname <- stringr::str_extract(protein, "[A-Z0-9]{6,}")
                  }
                  if (!is.na(searchname)){
                    searchlink <- paste('https://www.uniprot.org/uniprot/', 
                                        searchname, sep = "")
                    return(searchlink)
                  } else {
                    return("Could not extract protein name")
                  }
                }))
          }
        } else if ("metabData" %in% parent_dat){
          if(!is.null(try_URL)){
            cogs <- dplyr::mutate(
              cogs,
              Metabolite_URL = purrr::map_chr(
                cogs[[uniqueID]],
                function(metabolite){
                  searchlink <- paste(
                    "https://pubchem.ncbi.nlm.nih.gov/compound/", 
                    metabolite, 
                    sep = "")
                  return(searchlink)
                }))
          }
          
        } else if ("lipidData" %in% parent_dat){
          
          if(!is.null(try_URL)){
            cogs <- dplyr::mutate(
              cogs,
              Lipid_Maps_Hits = purrr::map_chr(
                cogs[[uniqueID]],
                function(lipid){
                  searchname <- gsub(" ", "", lipid)
                  searchname <- gsub(":", " ", searchname)
                  searchname <- stringr::str_remove_all(searchname, "[[:punct:]]")
                  searchname <- gsub(" ", ":", searchname)
                  searchlink <- paste(
                    'https://www.lipidmaps.org/data/structure/LMSDFuzzySearch.php?Name=',
                    searchname, 
                    sep = "")
                  return(searchlink)
                }))
          }
        }
        
        abundance <- colnames(cogs)[stringr::str_detect(colnames(cogs), "abundance")]
        
        coltrans <- function(column){
          ### Needs help from new function defining column attributes
          if(column %in% c(joiner, 
                           pmartR::get_fdata_cname(trellData),
                           abundance,
                           "Comparison", 
                           "Flag", 
                           "Count",
                           "peps_per_pro",
                           "n_peps_used",
                           "Mean",
                           "P_value_T",
                           "P_value_G",
                           "Fold_change",
                           "Sig_T_p_value",
                           "Sig_G_p_value",
                           "Protein_URL",
                           "Metabolite_URL",
                           "Lipid_Maps_Hits")){
            
            if (column == "Comparison"){
              return(trelliscopejs::cog(
                toString(unique(cogs[[column]])),
                desc = "Pairwise comparison(s)."))
              
            } else if (column == "Flag"){
              return(trelliscopejs::cog(
                toString(unique(cogs[[column]])),
                desc = "Trend of foldchange difference, where 1 and -1 indicate positive or negative quantitative changes within a pairwise comparison and 2 and -2 indicate qualititative changes within a pairwise comparison."))
              
            } else if (column == "Count") {
              x <- unique(cogs[[column]])
              if (length(x) != length(unique(cogs[["Group"]]))){
                x <- toString(rep(x, length(unique(cogs[["Group"]]))))
              } else {
                x <- toString(x)
              }
              return(trelliscopejs::cog(
                x,
                desc = "Number of sample observations in respective experimental groups."))
              
            } else if (column == "Protein_URL") {
              if(toString(unique(cogs[[column]])) == "Could not extract protein name"){
                return(trelliscopejs::cog(
                  toString(unique(cogs[[column]])), desc = "Link to protein UniProt page"))
              } else {
                return(trelliscopejs::cog_href(
                  toString(unique(cogs[[column]])),
                  desc = "Link to protein UniProt page"))
              }
              
            } else if (column == "Metabolite_URL") {
              return(trelliscopejs::cog_href(
                toString(unique(cogs[[column]])),
                desc = "Link to metabolite PubChem page."))
              
            } else if (column == "Lipid_Maps_Hits") {
              return(trelliscopejs::cog_href(
                toString(unique(cogs[[column]])),
                desc = "Link to Lipid Maps search of lipid species."))
            } else if (column == "peps_per_pro"){
              if(is.na(mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE))){
                return(trelliscopejs::cog(
                  0,
                  desc = "Number of peptides mapped to a given protein (used in pmartR::protein_quant() rollup method). Average is computed where multiple proteins are in a panel. (Value of 0 == NA) "))
              }
              return(trelliscopejs::cog(
                suppressWarnings(as.numeric(mean(cogs[[column]]))),
                desc = "Number of peptides mapped to a given protein (used in pmartR::protein_quant() rollup method). Average is computed where multiple proteins are in a panel. (Value of 0 == NA) "))
              
            } else if (column == "n_peps_used"){
              if(is.na(mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE))){
                return(trelliscopejs::cog(
                  0,
                  desc = "Number of peptides used in pmartR::protein_quant() rollup method. Average is computed where multiple proteins are in a panel. (Value of 0 == NA) "))
              }
              return(trelliscopejs::cog(
                mean(suppressWarnings(as.numeric(cogs[[column]]))),
                desc = "Number of peptides used in pmartR::protein_quant() rollup method. Average is computed where multiple proteins are in a panel. (Value of 0 == NA) "))
              
            } else if (column == "Mean"){
              return(trelliscopejs::cog(
                toString(unique(cogs[[column]])),
                desc = "Mean of abundances/intensities for experimental groups."))
              
            } else if (column == "Fold_change"){
              if(is.na(mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE))){
                return(trelliscopejs::cog(
                  0,
                  desc = "Mean of foldchange difference (per panel). (Value of 0 == NA) "))
              }
              return(trelliscopejs::cog(
                mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE),
                desc = "Mean of foldchange difference (per panel). (Value of 0 == NA) "))
              
            } else if (column == "Sig_T_p_value") {
              return(trelliscopejs::cog(
                toString(any(cogs[[column]])),
                desc = "Significant ANOVA results based on input p-value threshold (default == 0.05)."))
              
            } else if (column == "P_value_G") {
              if(is.na(mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE))){
                return(trelliscopejs::cog(
                  1,
                  desc = "Mean G-test p-values (per panel). (Value of 1 == NA) "))
              }
              return(trelliscopejs::cog(
                mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE),
                desc = "Mean G-test p-values (per panel)."))
              
            } else if (column == "P_value_T") {
              if(is.na(mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE))){
                return(trelliscopejs::cog(
                  1,
                  desc = "Mean ANOVA p-values (per panel).  (Value of 1 == NA) "))
              }
              return(trelliscopejs::cog(
                mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE),
                desc = "Mean ANOVA p-values (per panel)."))
              
            } else if (column == "Sig_G_p_value") {
              return(trelliscopejs::cog(
                toString(any(suppressWarnings(as.logical(cogs[[column]])))),
                desc = "Significant G-test results based on input p-value threshold (default == 0.05)."))
              
            } else if (column == abundance) {
              return(trelliscopejs::cog(
                mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE),
                desc = "Average abundance/intensity of values."))
              
            } else {
              return(trelliscopejs::cog(toString(unique(cogs[[column]])), desc = "User defined variable. (Categorical)"))
            }
            
          } else if (is.numeric(cogs[[column]])){
            if(is.na(mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE))){
              return(trelliscopejs::cog(
                0,
                desc = "User defined variable. (Numeric mean, Value of 0 == NA) "))
            }
            return(trelliscopejs::cog(
              mean(cogs[[column]], na.rm = TRUE), 
              desc = "User defined variable. (Numeric mean, Value of 0 == NA) "))
            
          } else if (is.logical(cogs[[column]]) || 
                     !any(is.na(suppressWarnings(as.logical(cogs[[column]]))))){
            return(trelliscopejs::cog(toString(any(cogs[[column]])), 
                                      desc = "User defined variable. (Categorical)"))
            
          } else if (is.character(cogs[[column]]) && 
                     !any(is.na(suppressWarnings(as.numeric(cogs[[column]]))))){
            if(is.na(mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE))){
              return(trelliscopejs::cog(
                0,
                desc = "User defined variable. (Numeric mean, Value of 0 == NA) "))
            }
            return(trelliscopejs::cog(mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE), 
                                      desc = "User defined variable. (Numeric mean, Value of 0 == NA) "))
            
          } else {
            return(trelliscopejs::cog(toString(unique(cogs[[column]])), 
                                      desc = "User defined variable. (Categorical)"))
          }
        }
        
        cog_out <- cogs[1,]
        for(col in colnames(cog_out)) {
          cog_out[[col]] <- coltrans(col)
        }
         
        cog_out[[panel_variable]] <- NULL
        
        return(cog_out)
      }))
  
  return(out)
  
}

#' @name trelliVis
#' @rdname trelliVis
#' @title Plot pairwise comparisons and data values in trellData object
#'
#' @description Plot pairwise comparisons and data values of omicsData and omicsStats objects. Customizable for plot types, y axis limits, paneling variable (what overall group is plotted on each graph arrangement), as well as desired variables for the y and x axis.
#'
#' @param omicsData A pmartR object of class pepData, lipidData, metabData, or proData. Can use list(pepData, proData) for associated data.
#' @param omicsStats A statistical results object produced by running \code{imd_anova} on omicsData. Can use list(pepStats, proStats) for associated data.
#' @param omicsFormat Output of as.trellData() function
#' @param p_val Numeric that specifies p-value for significance calculations. Default is 0.05.
#' @param panel_variable String: Name of column that plot panels are sorted by (e.g. each plotting arrangement has a unique identifier from panel variable). Default is emeta_cname if present, edata_cname where emeta_cname is not present.
#' @param try_URL Will attempt to link to PubChem, LipidMaps, or Uniprot based on information in edata_cname of omicsData or specified mapping_col for peptide data. Default is FALSE.
#' @param trelli_name String: name of display, or list of names where a list is provided for omicsData and omicsStats
#' @param trelli_path_out String: path to where trelliscope is stored. Default is "./TrelliDisplay"
#' @param interactive Should the plots be rendered as plotly objects?
#' @param plot_text Disable plot text
#' @param y_limits Y limits 
#' @param plot_type plots for plotting
#' @param self_contained Should display be generated in document? Defaults to FALSE
#'
#' @author Rachel Richardson
#' @export
trelliVis <- function(...) {
  .trelliVis(...)
}

.trelliVis <- function(omicsData = NULL, omicsStats = NULL,
                       omicsFormat = NULL, p_val = 0.05,
                       panel_variable = NULL, 
                       try_URL = FALSE, trelli_name = NULL,
                       trelli_path_out = "TrelliDisplay", 
                       plot_text = FALSE, interactive = FALSE,
                       y_limits = NULL, plot_type = NULL,
                       self_contained = FALSE) {
  
  
  #store_object, custom_cog_df, plot package = ggplot, rbokeh, etc, trelliscope additional arguments
  
  if(is.null(omicsFormat) && is.null(omicsData) && is.null(omicsStats)) stop(
    "At least one of omicsData, omicsStats, or omicsFormat must be populated."
  )
  
  #####
  ## Switch Stats and Omics data as appropriate ##
  #####
  
if(!is.null(omicsData) || !is.null(omicsStats)){
  
  if (is.null(omicsStats) & inherits(omicsData, 'statRes')){
    omicsStats <- omicsData
    omicsData <- NULL
  }
  
  if (!is.null(omicsStats) && 
      !is.null(omicsData) && 
      inherits(omicsStats, c("proData", "pepData", 
                             "metabData", "lipidData")) &&
      inherits(omicsData, 'statRes')){
    temp <- omicsStats
    omicsStats <- omicsData
    omicsData <- temp
    rm(temp)
  }
  
  
  
  #####
  ## Check variable type/lenghts
  #####
  
  # If lists of omicsData or trellData are used, make sure panel variables are specified for both #
  if ((class(omicsData) == "list" | 
       class(omicsStats) == "list") && 
      !is.null(panel_variable) && 
      length(panel_variable) != max(length(omicsStats), 
                                    length(omicsData))) stop(
                                      "Panel variable must be specified for each index in omicsStats/omicsData"
                                    )
  
  #####
  ## Modification for rollups that don't inherit pep info in emeta
  #####
  
  #Adds e_meta of associated peptide data to associated protein data
  if(class(omicsData) == "list" &&
     length(omicsData) != 1){
    vectorpro <- c(inherits(omicsData[1][[1]], "proData"), 
                   inherits(omicsData[2][[1]], "proData"))
    omicsData[vectorpro][[1]]$e_meta <- suppressWarnings(dplyr::left_join(omicsData[vectorpro][[1]]$e_meta,
                                                                          omicsData[!vectorpro][[1]]$e_meta))
  }
  
  
}
  # if (is.null(omicsStats) & inherits(omicsData, 'statRes')){
  #   omicsStats <- omicsData
  #   omicsData <- NULL
  # }
  # 
  # if (!is.null(omicsStats) && 
  #     !is.null(omicsData) && 
  #     inherits(omicsStats, c("proData", "pepData", 
  #                            "metabData", "lipidData")) &&
  #     inherits(omicsData, 'statRes')){
  #   temp <- omicsStats
  #   omicsStats <- omicsData
  #   omicsData <- temp
  #   rm(temp)
  # }
  
  
  
  # #####
  # ## Check variable type/lenghts
  # #####
  # 
  # # If lists of omicsData or trellData are used, make sure panel variables are specified for both #
  # if ((class(omicsData) == "list" | 
  #     class(omicsStats) == "list") && 
  #     !is.null(panel_variable) && 
  #     length(panel_variable) != max(length(omicsStats), 
  #                                   length(omicsData))) stop(
  #                                     "Panel variable must be specified for each index in omicsStats/omicsData"
  #                                   )
  
  # Check check if try_URL is boolean of length 1 #
  if(!is.logical(try_URL) | (length(try_URL) != 1)) stop(
    "try_URL must be a TRUE/FALSE of length 1")  
  
  #####
  ## Modification for rollups that don't inherit pep info in emeta
  #####
  # 
  # #Adds e_meta of associated peptide data to associated protein data
  # if(class(omicsData) == "list" &&
  #    length(omicsData) != 1){
  #   vectorpro <- c(inherits(omicsData[1][[1]], "proData"), 
  #                  inherits(omicsData[2][[1]], "proData"))
  #   omicsData[vectorpro][[1]]$e_meta <- suppressWarnings(dplyr::left_join(omicsData[vectorpro][[1]]$e_meta,
  #                                                                         omicsData[!vectorpro][[1]]$e_meta))
  # }
  
  if (is.null(omicsFormat)){
    # Re-format objects for plotting #
    trellData <- as.trellData(omicsData, omicsStats)
  } else {
    trellData <- omicsFormat
  }
  
  # If a pep/pro pair is listed, act on each item #
  if(class(trellData) == "list"){
    
    
    # Fill plot type
    if (is.null(plot_type)){
      if(is.null(trellData[[1]]$comp_stats)){
        plot_type <- "abundance_boxplot"
      } else if (is.null(trellData[[1]]$data_values)){
        plot_type <- "foldchange_bar"
      } else {
        plot_type <- list("abundance_boxplot", "foldchange_bar")
      }
    }
    
    
    if (!pmartR::get_data_norm(trellData[[1]]) &&
        (is.null(attributes(trellData[[1]])$isobaric_info$norm_info$is_normalized) || #### update helper function
         attributes(trellData[[1]])$isobaric_info$norm_info$is_normalized != TRUE)) stop(
           "Input must be normalized prior to plotting; use normalize_global or normalize_isobaric as appropriate."
         )
    
    if (!pmartR::get_data_norm(trellData[[2]]) &&
        (is.null(attributes(trellData[[2]])$isobaric_info$norm_info$is_normalized) || #### update helper function
         attributes(trellData[[2]])$isobaric_info$norm_info$is_normalized != TRUE)) stop(
           "Input must be normalized prior to plotting; use normalize_global or normalize_isobaric as appropriate."
         )
    
    ## Check trelli_name ##
    # Default trelliscope names correspond to data types entered #
    if(is.null(trelli_name)){
      trelli_name <- vector("list", 2)
      trelli_name[[1]] <- paste(attr(trellData, "data_types")[1], plot_type, sep = "_")
      trelli_name[[2]] <- paste(attr(trellData, "data_types")[2], plot_type, sep = "_")
      names(trelli_name) <- attr(trellData, "data_types")
      

      # if trelli_name is not of length list, add identifiers #
    } else if (length(trelli_name) == 1){
      label <- trelli_name
      trelli_name <- vector("list", 2)
      trelli_name[[1]] <- paste(paste(label, attr(trellData, "data_types")[1], sep = "_"), plot_type, sep = "_")
      trelli_name[[2]] <- paste(paste(label, attr(trellData, "data_types")[2], sep = "_"), plot_type, sep = "_")
      names(trelli_name) <- attr(trellData, "data_types")
      
      # if trelli_name too long, cut to length(trellData)  #
    } else if (length(trelli_name) != length(plot_types)*2) {
      message("Length of trelli_name is not equal to the number of displays generated. Only the first name will be used.")
      label <- trelli_name[1]
      trelli_name <- vector("list", 2)
      trelli_name[[1]] <- paste(label, attr(trellData, "data_types")[1], plot_type, sep = "_")
      trelli_name[[2]] <- paste(label, attr(trellData, "data_types")[2], plot_type, sep = "_")
      names(trelli_name) <- attr(trellData, "data_types")
    } else {
      label <- trelli_name
      trelli_name <- vector("list", 2)
      trelli_name[[1]] <- trelli_name[1:length(plot_types)]
      trelli_name[[2]] <- trelli_name[(length(plot_types)+1):(length(plot_types)*2)]
      names(trelli_name) <- attr(trellData, "data_types")
    }
    
    # Ensure NULLs ar the correct length #
    if (is.null(omicsData)){
      omicsData <- rep(list(NULL), length(trellData))
    }
    if (is.null(omicsStats)){
      omicsStats <- rep(list(NULL), length(trellData))
    }
    if (is.null(panel_variable)){
      panel_variable <- rep(list(NULL), length(trellData))
    }
    
    # Nest data and generate trelliscope plots #
    
    tictoc::tic("Generate plots")
    
    nested_plot <- purrr::map2(trellData, panel_variable, function(pairedplotter, pan){
      
      nest_out <- purrr::map(plot_type, function(types){
        format_plot(trellData = pairedplotter, 
                    p_val = p_val,
                    panel_variable = pan, 
                    plot_type = types,
                    plot_text = plot_text,
                    y_limits = y_limits,
                    interactive = interactive)
        })
      return(nest_out)
      })
    
    tictoc::toc()
    
    tictoc::tic("Generate auto cogs") 
    
    # Generate default cognostics #
    nest_plot_cog_list <- purrr::pmap(
      list(trellData,
          nested_plot,
          rev(trelli_name)),
      function(trell, nests, links){
        
        cogs <- suppressWarnings(data_cogs(
          nested_plot = nests[[1]],
          trellData = trell,
          p_val = p_val,
          try_URL = try_URL))
        
        pan <- colnames(cogs)[1]
        
        cog_out <- purrr::map2(nests, links, function(nest, link){
          cogs2 <- dplyr::mutate(nest, cogs = cogs$cogs)
          
          if (pmartR:::get_data_class(trell) == "pepData"){
            
            x <- trelliscopejs::cog_disp_filter(
              link,
              var = pan,
              val = cogs2[[pan]],
              desc = "Link to a display of the mapping protein."
            )
            x <- gsub(";type:select;val", ";type:regex;val", x)
            cogs2 <- dplyr::mutate(cogs2, Protein_display = x)
            
          } else {
            cogs2 <- dplyr::mutate(cogs2, Peptide_display = trelliscopejs::cog_disp_filter(
                                       link,
                                       var = pan,
                                       val = cogs2[[pan]],
                                       desc = "Link to a display with all associatated peptides."
                                     ))
          }
          return(cogs2)
        })
          
          return(cog_out)
          })
        
    tictoc::toc()
    
    # Generate trelliscope display #
    tictoc::tic("Pipe into trelliscope")
    
    out <- purrr::map2(nest_plot_cog_list, 
                trelli_name, function(display1, name1){
                  purrr::map2(display1, name1, function(display, name){
                    trelliscopejs::trelliscope(display, as.character(name), nrow = 1, ncol = 2,
                                               path = as.character(trelli_path_out), thumb = TRUE, state = list(
                                                 sort = list(trelliscopejs::sort_spec(names(display[1]), dir = "desc")), 
                                                 labels = list(names(display[1]))))
                  })
      })
    
    tictoc::toc()
    
    return(out[[1]][[1]])
    
  # Where not a pep/pro pair: #
  } else {
    
    # Fill plot type
    if (is.null(plot_type)){
      if(is.null(trellData$comp_stats)){
        plot_type <- "abundance_boxplot"
      } else if (is.null(trellData$data_values)){
        plot_type <- "foldchange_bar"
      } else {
        plot_type <- list("abundance_boxplot", "foldchange_bar")
      }
    }
    
    if (!pmartR::get_data_norm(trellData) &&
        (is.null(attributes(trellData)$isobaric_info$norm_info$is_normalized) || #### update helper function
         attributes(trellData)$isobaric_info$norm_info$is_normalized != TRUE)) stop(
           "Input must be normalized prior to plotting; use normalize_global or normalize_isobaric as appropriate."
         )
    
    # Default: Name the display after parent classes (omicsData and omicStats clases) #
    if(is.null(trelli_name)){
      trelli_name <- paste(paste(attr(trellData, "parent_class"), sep = "_"), plot_type, sep = "_")
      # Make sure trelliname is length 1 #
    } else if (length(trelli_name) != (length(plot_type))) {
      warning("trelli_name length is not equal to the number of displays generated from plot_type; using the first entry in trelli_name.")
      trelli_name <- paste(trelli_name[[1]], plot_type, sep = "_")
    }
    
    # Nest data and generate trelliscope plots #
    
    displays <- purrr::map2(plot_type, trelli_name, function(types, names){
      
      tictoc::tic("Generate plots")
      
      nested_plot <- format_plot(trellData = trellData, 
                                 p_val = p_val,
                                 plot_type = types,
                                 plot_text = plot_text,
                                 y_limits = y_limits,
                                 panel_variable = panel_variable,
                                 interactive = interactive)
      tictoc::toc()
      
      tictoc::tic("Generate auto cogs")
      
      # Generate default cognostics #
      nest_plot_cog <- suppressWarnings(data_cogs(nested_plot = nested_plot, 
                                                  trellData = trellData,
                                                  p_val = p_val, 
                                                  try_URL = try_URL))
      tictoc::toc()
      
      tictoc::tic("Pipe into trelliscope")
      
      # Generate trelliscope display #
      if(self_contained){
        out <- nest_plot_cog %>%
          trelliscopejs::trelliscope(name = as.character(names), nrow = 1, ncol = 2,
                                     self_contained = TRUE, 
                                     thumb = TRUE,
                                     state = list(
                                       sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")), 
                                       labels = list(names(nest_plot_cog[1])))
          )
      } else {
        
        out <- nest_plot_cog %>%
          trelliscopejs::trelliscope(name = as.character(names), nrow = 1, ncol = 2,
                                     path = as.character(trelli_path_out), 
                                     thumb = TRUE,
                                     state = list(
                                       sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")), 
                                       labels = list(names(nest_plot_cog[1])))
          )
        
      }
      tictoc::toc()
      return(out)
      })
    return(displays)
  }
}
