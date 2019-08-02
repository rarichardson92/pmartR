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
    data_values <- suppressWarnings(dplyr::left_join(data_values, joingroupDF))
    
    # Join with e_meta if present #
    if(!is.null(omicsData$e_meta)) {
      data_values <- suppressWarnings(dplyr::left_join(data_values, 
                               omicsData$e_meta))
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
                              omicsData$e_meta))
      summary_stats <- suppressWarnings(dplyr::left_join(summary_stats,
                              omicsData$e_meta))
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
    if(is.null(omicsStats) & is.null(omicsData)) stop(
      "as.trellData() requires at least one of the following: 
      omicsStats, omicsData")
    
    # Check that omicsData and omicsStats are the correct classes #
    if (!is.null(omicsData) & 
        !inherits(omicsData, c("proData", 
                               "pepData", 
                               "metabData", 
                               "lipidData"))) stop(
                                 "omicsData must be of class 'proData', 
                                 'pepData', 'metabData', or 'lipidData'")
    if(!is.null(omicsStats) & !inherits(omicsStats, "statRes")) stop(
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
 
  if (is.null(panel_variable)){
    panel_variable <- pmartR::get_edata_cname(trellData)
  }
  
  if (is.null(plot_type)){
    plot_type <- list("abundance_boxplot", "foldchange_bar")
  }
   
  ## Input digestion
  graphlims <- list_y_limits(plot_type = plot_type, y_limits = y_limits)
  
  All_plots <- list()
  
  if (panel_variable == pmartR::get_edata_cname(trellData) && stringr::str_detect(plot_type, "heatmap")) {
    stop(
      paste("Heatmaps require a panel_variable that is not the same as edata cname (Default panel_variable is set to edata cname). Current edata cname:", pmartR::get_edata_cname(trellData))
    )
  } 
  
  cores<- parallel::detectCores()
  cl<- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  
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
    
    plot_base <- ggplot2::ggplot() +
      ggplot2::theme_bw() +
      ggplot2::geom_boxplot(ggplot2::aes(x = trellData$data_value[[sampID]], 
                                         y = trellData$data_value[[abundance]], 
                                        color = trellData$data_value[[group]],
                                        fill = trellData$data_value[[group]]),
                            alpha = 0.75,
                            na.rm = TRUE) +
      ggplot2::labs(x = group, y = abundance, fill = group, color = group) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    
    
    data <- trellData$data_value
    # print(rank(trellData$data_value[[abundance]]))
    # data <- dplyr::mutate(data, rank = paste(rank(trellData$data_value[[abundance]]),"/", nrow(data), sep = ""))
    # print(data)
    
    
    output_df <- data.frame(unique(as.character(trellData$data_value[[panel_variable]])), stringsAsFactors = FALSE)
    colnames(output_df) <- panel_variable
    
    All_plots[["abundance_global"]] <- suppressWarnings(
      foreach::foreach(i=1:length(output_df[[panel_variable]]))%dopar%{
        
        cat(paste("abundance_global:", i, "/", length(output_df[[panel_variable]]), "plots\n"), 
            file="log.trelliVis.txt", append=TRUE)
        panel <- output_df[[panel_variable]][i]
    
      
      
    # All_plots[["abundance_global"]] <- suppressWarnings(purrr::map(output_df[[panel_variable]], function(panel){
      rows <- data[[panel_variable]] == panel
      df <- data[rows,]
      
      if(!is.null(lims$max) || 
         !is.null(lims$min) ||
         !is.null(lims$range)){
        
        increment <- set_increment(df[[abundance]])
        setlims <- set_ylimits(df[[abundance]],
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
          x = as.numeric(factor(df[[sampID]])) - 0.25,
          xend = as.numeric(factor(df[[sampID]])) + 0.25,
          y = df[[abundance]],
          yend = df[[abundance]],
          ggplot2::aes(text = paste(paste(abundance, ":", sep = ""), signif(df[[abundance]], 4))),
          na.rm = TRUE) +
        
        ggplot2::geom_point(
          x = as.numeric(factor(df[[sampID]])),
          y = df[[abundance]],
          ggplot2::aes(text = paste(paste(abundance, ":", sep = ""), signif(df[[abundance]], 4))),
          color = NA,
          na.rm = TRUE)
      
      if(length(unique(trellData$data_value[[sampID]])) > 8){
        plot_out <- plot_out + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5))
      }
      if(length(unique(trellData$data_value[[sampID]])) > 20){
        plot_out <- plot_out + ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
          ggplot2::xlab(sampID)
      }
      
      if(exists("setlims")){
        plot_out <- plot_out + ggplot2::coord_cartesian(ylim = setlims)
      }
      
      return(plot_out)
      
    })#)
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
    # print(rank(trellData$comp_stats[[foldchange]]))
    # data <- dplyr::mutate(data, rank = paste(rank(trellData$data_value[[abundance]]),"/", nrow(data), sep = ""))
    # print(data)
    
    output_df <- data.frame(unique(as.character(trellData$comp_stats[[panel_variable]])), stringsAsFactors = FALSE)
    colnames(output_df) <- panel_variable
    
    All_plots[["foldchange_global"]] <- suppressWarnings(purrr::map(output_df[[panel_variable]], function(panel){
      rows <- trellData$comp_stats[[panel_variable]] == panel
      df <- trellData$comp_stats[rows,]
      
      if(!is.null(lims$max) || 
         !is.null(lims$min) ||
         !is.null(lims$range)){
        
        increment <- set_increment(df[[foldchange]])
        setlims <- set_ylimits(df[[foldchange]],
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
          x = as.numeric(factor(df[[comps]])) - 0.25,
          xend = as.numeric(factor(df[[comps]])) + 0.25,
          y = df[[foldchange]],
          yend = df[[foldchange]],
          ggplot2::aes(text = #paste(
            paste(paste(foldchange, ":", sep = ""), signif(df[[foldchange]], 4))),
            # paste("Rank:", df[["rank"]]),
            # sep = "\n")),
          na.rm = TRUE) +
        
        ggplot2::geom_point(
          x = as.numeric(factor(df[[comps]])),
          y = df[[foldchange]],
          ggplot2::aes(text = paste(paste(foldchange, ":", sep = ""), signif(df[[foldchange]], 4))),
          color = NA,
          na.rm = TRUE)
      
      
      if(exists("setlims")){
        plot_out <- plot_out + ggplot2::coord_cartesian(ylim = setlims)
      }
      
      return(plot_out)
    }))
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
      
      output_df <- data.frame(unique(
        as.character(trellData$summary_stats[[panel_variable]]), stringsAsFactors = FALSE))
      colnames(output_df) <- panel_variable
      
      All_plots[["missing_bar"]] <- suppressWarnings(purrr::map(output_df[[panel_variable]], function(panel){
        rows <- as.character(trellData$summary_stats[[panel_variable]]) == panel
        df <- trellData$summary_stats[rows,]
        totals <- data.frame(pmartR::get_group_table(trellData), stringsAsFactors = FALSE)
        
        colnames(totals) <- c(group, "Total_Group_Counts")
        df <- dplyr::left_join(df, totals, by = group)
        
        df$Missing <-  df[["Total_Group_Counts"]] - df[["Count"]]
        df$Non_Missing <- df[["Count"]]
        
        df <- reshape2::melt(df, 
                             id.vars = colnames(df)[!stringr::str_detect(colnames(df), "Missing")], 
                             value.name = "y",
                             variable.name = "Data")
        
        
        groupgraphs <- suppressWarnings(purrr::map(unique(df[[group]]), function(dfgr){
          
          rows <- df[[group]] == dfgr
          dfsm <- df[rows,]
          
          miss <- dfsm$Data == "Missing"
          dfsm$y[miss] <- sum(dfsm$y[miss])
          dfsm$y[!miss] <- sum(dfsm$y[!miss])
          
          dup <- nrow(dfsm)/2
          dfsm[["Total_Group_Counts"]] <- dfsm[["Total_Group_Counts"]]*dup
          
          setlims <- NULL
          
          if(is.null(lims$max)){
            lims$max <- 1
          }
          if(is.null(lims$min)){
            lims$min <- 0
          }
          
          increment <- set_increment(dfsm[["y"]])
          setlims <- set_ylimits(dfsm[["y"]],
                                 increment,
                                 y_max = lims$max,
                                 y_min = lims$min,
                                 include_zero = FALSE)
          
          texty <- paste(
            paste(
              paste(paste(group, ":", sep = ""),
              unique(dfsm[[group]])),
              paste("Count:", unique(dfsm[["y"]])), 
              sep = "\n"),
            paste("Status:", unique(dfsm[["Data"]])),
            sep = "\n")
          
          plot_out <- ggplot2::ggplot() +
            ggplot2::geom_col(
              ggplot2::aes(text = texty,
              x = unique(dfsm[[group]]),
              y = unique(dfsm[["y"]]),
              fill = unique(dfsm[["Data"]])),
              position = "fill") +
            ggplot2::theme_bw() +
            ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
            ggplot2::scale_y_continuous(
              sec.axis = ggplot2::sec_axis(
                ~.*max(dfsm[["Total_Group_Counts"]]),
                name = "", 
                breaks = c(seq(from = 0,
                               to = max(dfsm[["Total_Group_Counts"]]) - round(max(dfsm[["Total_Group_Counts"]])/10),
                               by = round(max(dfsm[["Total_Group_Counts"]])/10)+1),
                           max(dfsm[["Total_Group_Counts"]])
                )),
              expand = c(0,0)) +
            ggplot2::scale_x_discrete(expand = c(0,0)) + 
            ggplot2::coord_cartesian(ylim = setlims)
          
          
          return(plot_out)
          
        }))
        
        if (interactive){

          plotlys <- purrr::map(groupgraphs, function(plot) plotly::ggplotly(plot, tooltip = "text"))
          plot <- plotly::subplot(plotlys, margin = 0.05)

        } else {
            
          plot <- ggpubr::ggarrange(plotlist = groupgraphs, common.legend = TRUE)
          plot <- ggpubr::annotate_figure(plot, 
                                          right = ggpubr::text_grob("Count", size = 8, rot = 270), 
                                          left = ggpubr::text_grob("Proportion", size = 8, rot = 90))

        }

        return(plot)
        
      }))
    } else {
      
      if ("Group_DF" %in% colnames(trellData$data_values)){
        group <- "Group_DF"
      } else {
        group <- "Group"
      }
      
      output_df <- data.frame(unique(
        as.character(trellData$data_values[[panel_variable]]), stringsAsFactors = FALSE))
      colnames(output_df) <- panel_variable
      
      All_plots[["missing_bar"]] <- suppressWarnings(purrr::map(output_df[[panel_variable]], function(panel){
        rows <- trellData$data_values[[panel_variable]] == panel
        df <- trellData$data_values[rows,]
        
        Count <- data.frame(table(df[complete.cases(df),]$Group))
        totals <- data.frame(pmartR::get_group_table(trellData), stringsAsFactors = FALSE)
        n_rep <- length(unique(df[[pmartR::get_edata_cname(trellData)]]))
        
        totals$Freq <- totals$Freq * n_rep
        
        colnames(Count) <- c(group, "Count")
        colnames(totals) <- c(group, "Total_Group_Counts")
        df <- dplyr::left_join(df, Count, by = group)
        df <- dplyr::left_join(df, totals, by = group)
        
        df$Missing <-  df[["Total_Group_Counts"]] - df[["Count"]]
        df$Non_Missing <- df[["Count"]]
        
        df <- reshape2::melt(df, 
                             id.vars = colnames(df)[!stringr::str_detect(colnames(df), "Missing")],
                             value.name = "y",
                             variable.name = "Data")
        
        
        groupgraphs <- purrr::map(unique(df[[group]]), function(dfgr){
          rows <- df[[group]] == dfgr
          dfsm <- df[rows,]
          
          setlims <- NULL
          
          if(is.null(lims$max)){
            lims$max <- 1
          }
          if(is.null(lims$min)){
            lims$min <- 0
          }

          increment <- set_increment(dfsm[["y"]])
          setlims <- set_ylimits(dfsm[["y"]],
                                 increment,
                                 y_max = lims$max,
                                 y_min = lims$min,
                                 include_zero = FALSE)
          
          plot_out <- ggplot2::ggplot() +
            ggplot2::theme_bw() +
            ggplot2::geom_col(ggplot2::aes(x = unique(dfsm[[group]]),
                                           y = unique(dfsm[["y"]]),
                                           fill = unique(dfsm[["Data"]])),
                              position = "fill") +
            ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
            ggplot2::scale_y_continuous(
              sec.axis = ggplot2::sec_axis(
                ~.*max(dfsm[["Total_Group_Counts"]]),
                name = "", 
                breaks = c(seq(from = 0,
                               to = max(dfsm[["Total_Group_Counts"]]) - round(max(dfsm[["Total_Group_Counts"]])/10),
                               by = round(max(dfsm[["Total_Group_Counts"]])/10)+1),
                           max(dfsm[["Total_Group_Counts"]])
                )),
              expand = c(0,0)) +
            ggplot2::scale_x_discrete(expand = c(0,0)) +
            ggplot2::coord_cartesian(ylim = setlims)
          
          return(plot_out)
          
        })
        
        plot <- ggpubr::ggarrange(plotlist = groupgraphs, common.legend = TRUE)
        plot <- ggpubr::annotate_figure(plot, right = "Count of samples", left = "Proportion of samples")
        
        return(plot)
      }))
    }
  }
  
  #### Abundance Heatmap ####
  
  if("abundance_heatmap" %in% plot_type){
    
    lims <- graphlims[["abundance_heatmap"]]
    
    if (!is.null(lims[[1]])) {
      warning( 
        "y_limits are not supported with heatmaps and will not be used.  Refer to examples in ?trelliVis().")
    }
    
    
    if(attr(trellData, "meta_info") && 
       panel_variable != pmartR::get_edata_cname(trellData)){
      
      if ("Group_DF" %in% colnames(trellData$data_values)){
        group <- "Group_DF"
      } else {
        group <- "Group"
      }
      
      abundance <- grep("abundance", colnames(trellData$data_values), value = T)
      
      pal <- grDevices::colorRampPalette(c("red", "yellow"))
      
      output_df <- trellData$data_values[c(abundance, 
                                           group,
                                           pmartR::get_fdata_cname(trellData),
                                           pmartR::get_edata_cname(trellData),  
                                           panel_variable)]
      
      output_df <- tidyr::nest(output_df, -!!rlang::sym(panel_variable))
      
      All_plots[["abundance_heatmap"]] <- suppressWarnings(purrr::map(output_df$data, function(df){
        
        df[[pmartR::get_fdata_cname(trellData)]] <- ordered(
          df[[pmartR::get_fdata_cname(trellData)]], 
          levels = rev(sort(unique(df[[pmartR::get_fdata_cname(trellData)]]))))
        
        df[[pmartR::get_edata_cname(trellData)]] <- ordered(
          df[[pmartR::get_edata_cname(trellData)]], 
          levels = rev(sort(unique(df[[pmartR::get_edata_cname(trellData)]]))))
        
        
        texty <- paste(
          paste(
            paste(paste(pmartR::get_edata_cname(trellData), ":", sep = ""),
                  df[[pmartR::get_edata_cname(trellData)]]),
            paste(paste(pmartR::get_fdata_cname(trellData), ":", sep = ""),
                  df[[pmartR::get_fdata_cname(trellData)]]), sep = "\n"),
          paste(paste(abundance, ":", sep = ""), signif(df[[abundance]], 4)),
          sep = "\n")

        
        heatmap <- ggplot2::ggplot() +
          ggplot2::geom_tile(ggplot2::aes(text = texty,
            fill = df[[abundance]],
            x = df[[pmartR::get_fdata_cname(trellData)]],
            y = df[[pmartR::get_edata_cname(trellData)]])) +
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
        
        if (length(unique(df[[pmartR::get_edata_cname(trellData)]])) > 35){
          heatmap <- heatmap + ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
            ggplot2::ylab(pmartR::get_edata_cname(trellData))
        } else {
          heatmap <- heatmap + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6)) +
            ggplot2::ylab("")
        }
        
        if (length(unique(df[[pmartR::get_fdata_cname(trellData)]])) > 35){
          heatmap <- heatmap + ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
            ggplot2::xlab(pmartR::get_fdata_cname(trellData))
        } else {
          heatmap <- heatmap + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6)) +
            ggplot2::xlab("")
        }
        
        return(heatmap)
      }))
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
      
      output_df <- trellData$comp_stats[c(foldchange, 
                                          "Comparison", 
                                          pmartR::get_edata_cname(trellData),  
                                          panel_variable)]
      
      output_df <- tidyr::nest(output_df, -!!rlang::sym(panel_variable))
      
      All_plots[["foldchange_heatmap"]] <- suppressWarnings(purrr::map(output_df$data, function(df){
        
        
        texty <- paste(
          paste(
            paste(paste(pmartR::get_edata_cname(trellData), ":", sep = ""),
                df[[pmartR::get_edata_cname(trellData)]]),
            paste("Comparison:", df[["Comparison"]]), sep = "\n"),
          paste(paste(foldchange, ":", sep = ""), signif(df[[foldchange]], 4)),
          sep = "\n")
        
        heatmap <- ggplot2::ggplot(
          df, 
          ggplot2::aes(text = texty,
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
        
        if (length(unique(df[[pmartR::get_edata_cname(trellData)]])) > 35){
          heatmap <- heatmap + ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
            ggplot2::ylab(pmartR::get_edata_cname(trellData))
        } else {
          heatmap <- heatmap + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6)) +
            ggplot2::ylab("")
        }
        
        if (length(unique(df[[pmartR::get_fdata_cname(trellData)]])) > 35){
          heatmap <- heatmap + ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
            ggplot2::xlab(pmartR::get_fdata_cname(trellData))
        } else {
          heatmap <- heatmap + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6)) +
            ggplot2::xlab("")
        }
        
        
        return(heatmap)
      }))
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
      if(!is.null(trellData$data_values)){
        
        abundance <- grep("abundance", colnames(trellData$data_values), value = T)
        pal <- grDevices::colorRampPalette(c("red", "yellow"))
        
        output_df <- trellData$data_values[c(abundance,
                                             pmartR::get_fdata_cname(trellData),
                                             pmartR::get_edata_cname(trellData),
                                             panel_variable)]
        output_df <- tidyr::nest(output_df, -!!rlang::sym(panel_variable))
        All_plots[["presence_heatmap"]] <- suppressWarnings(purrr::map(output_df$data, function(df){
          df <- tidyr::nest(df, -c(!!rlang::sym(pmartR::get_edata_cname(trellData)),
                                   !!rlang::sym(pmartR::get_fdata_cname(trellData))))
          df <- dplyr::mutate(df, Biomolecule_Presence = purrr::map_lgl(df$data, function(dfdat){
            !all(is.na(dfdat[[abundance]]))
          }))
          
          df[["Biomolecule_Presence"]] <- gsub(TRUE, "Present", df[["Biomolecule_Presence"]])
          df[["Biomolecule_Presence"]] <- gsub(FALSE, "Absent", df[["Biomolecule_Presence"]])
          
          
          texty <- paste(
            paste(
              paste(paste(pmartR::get_edata_cname(trellData), ":", sep = ""),
                    df[[pmartR::get_edata_cname(trellData)]]),
              paste(paste(pmartR::get_fdata_cname(trellData), ":", sep = ""),
                    df[[pmartR::get_fdata_cname(trellData)]]), sep = "\n"),
            paste("Biomolecule_Presence:", df[["Biomolecule_Presence"]]),
            sep = "\n")
          
          heatmap <- ggplot2::ggplot(
            df,
            ggplot2::aes(text = texty,
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
          
          if (length(unique(df[[pmartR::get_edata_cname(trellData)]])) > 35){
            heatmap <- heatmap + ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
              ggplot2::ylab(pmartR::get_edata_cname(trellData))
          } else {
            heatmap <- heatmap + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6)) +
              ggplot2::ylab("")
          }
          
          if (length(unique(df[[pmartR::get_fdata_cname(trellData)]])) > 35){
            heatmap <- heatmap + ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
              ggplot2::xlab(pmartR::get_fdata_cname(trellData))
          } else {
            heatmap <- heatmap + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6)) +
              ggplot2::xlab("")
          }
          
          return(heatmap)
          
        }))
      }
    }
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
    
    if(!is.null(lims$scale) && lims$scale == "fixed"){
      
      setter <- plotter[["Fold_change"]] +
        (0.55 * plotter[["Fold_change"]])
      
      increment <- set_increment(setter)
      setlims <- set_ylimits(setter,
                             increment,
                             y_range = lims$range,
                             y_max = lims$max,
                             y_min = lims$min,
                             include_zero = TRUE)
    }
    
    output_df <- plotter %>% tidyr::nest(-panel_variable)
    
    #  # #Subset large groups ########### Take out later ######################################
    # if (nrow(nestplotter) > 10){
    #   nestplotter <- nestplotter[1:10,]
    # }
    
    All_plots[["foldchange_bar"]] <- suppressWarnings(purrr::map(output_df$data, function(nestedData) {
      
      ## Add border color, hover text, and label text to dataframe for plotting ##
      nestedData_comps <- add_plot_features(trellData,
                                            nestedData,
                                            p_val = p_val,
                                            panel_variable = panel_variable,
                                            comps_panel_y_axis = "Fold_change"   ### Remove from function later
                                            )
      
      
      # Make ggplots #
      
      if(!is.null(lims$scale) && lims$scale != "fixed" ||
         !is.null(lims$max) ||
         !is.null(lims$min) ||
         !is.null(lims$range)){
        
        setter <- nestedData_comps[["Fold_change"]] +
          (0.55 * nestedData_comps[["Fold_change"]])
        
        increment <- set_increment(setter)
        setlims <- set_ylimits(setter,
                               increment,
                               y_range = lims$range,
                               y_max = lims$max,
                               y_min = lims$min,
                               include_zero = TRUE)
      }
      
      plot_comps <- ggplot2::ggplot() + 
        ggplot2::theme_bw() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::theme(
                       axis.text.x = ggplot2::element_text(angle = 330, hjust = 0)) +
        
        ggplot2::xlab("Comparison") + 
        ggplot2::ylab("Fold_change") 
      
      if (pmartR::get_edata_cname(trellData) %in% colnames(nestedData_comps)){
        
        plot_comps <- plot_comps + ggplot2::geom_col(
          ggplot2::aes(
            text = gsub(", ", "\n", nestedData_comps[["text"]]),
            x = nestedData_comps[[pmartR::get_edata_cname(trellData)]],
            y = as.numeric(nestedData_comps[["Fold_change"]]),
            fill = as.character(nestedData_comps[["Comparison"]]),
            color = nestedData_comps[["bord"]]),
          position = "dodge2",
          size = 1,
          na.rm = TRUE) +
          ggplot2::labs(fill = "", color = "", x = NULL) +
          ggplot2::scale_color_manual(values = c("grey40" = "grey40", "black" = "black"), guide = FALSE)
        
        if (length(unique(nestedData_comps[[pmartR::get_edata_cname(trellData)]])) > 8){
          plot_comps <- plot_comps + ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 330, hjust = 0, size = 5))
        }
        
        if (length(unique(nestedData_comps[[pmartR::get_edata_cname(trellData)]])) > 21){
          plot_comps <- plot_comps + ggplot2::theme(
            axis.text.x = ggplot2::element_blank()) +
            ggplot2::xlab("Biomolecules")
        }
        
      } else {
        
        plot_comps <- plot_comps + ggplot2::geom_col(
          ggplot2::aes(
            text = gsub(", ", "\n", nestedData_comps[["text"]]),
            x = as.character(nestedData_comps[["Comparison"]]),
            y = as.numeric(nestedData_comps[["Fold_change"]]),
            fill = as.character(nestedData_comps[["Comparison"]])
            ),
          color = nestedData_comps[["bord"]],
          position = "dodge2",
          size = 1,
          na.rm = TRUE) +
          ggplot2::labs(fill = "", color = "", x = NULL) +
          ggplot2::scale_color_manual(values = c("grey40" = "grey40", "black" = "black")) + 
          ggplot2::theme(legend.position = "none")
        
        if (plot_text){
          
          textcomps <- nestedData_comps
          textcomps[["Fold_change"]][is.na(textcomps[["Fold_change"]])] <- 0
          
          plot_comps <-  plot_comps + ggplot2::geom_text(
            ggplot2::aes(
              x = as.character(textcomps[["Comparison"]]),
              y = textcomps[["Fold_change"]] + 
                (0.5 * textcomps[["Fold_change"]]), 
              label = gsub(", ", "\n", textcomps[["labels"]])),
            color = "black"
          )
        }
        
      }

      if (!is.null(lims[[1]])){
        plot_comps <- plot_comps + ggplot2::coord_cartesian(ylim = setlims)
      }
      
      if (interactive){
        plot_comps <- plot_comps + ggplot2::theme(legend.position = "none")
      }
      
      return(plot_comps)
    }))
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
      setlims <- set_ylimits(plotter[[abundance]],
                             increment,
                             y_range = lims$range,
                             y_max = lims$max,
                             y_min = lims$min,
                             include_zero = FALSE)
    }
    
    output_df <- plotter %>% tidyr::nest(-panel_variable)
    
    All_plots[["abundance_boxplot"]] <- suppressWarnings(purrr::map(output_df$data, function(nestedData) {
      
      ## Set hover, excluding the panel_variable ##
      nestedData_value <- add_plot_features(trellData,
                                            nestedData,
                                            p_val = p_val,
                                            panel_variable = panel_variable,
                                            value_panel_y_axis = abundance)  ### Take out of function later 
      
      # Make ggplots #
      if((!is.null(lims$scale) && lims$scale != "fixed") ||
         !is.null(lims$max) ||
         !is.null(lims$min) ||
         !is.null(lims$range)){
        increment <- set_increment(nestedData_value[[abundance]])
        setlims <- set_ylimits(nestedData_value[[abundance]],
                               increment,
                               y_range = lims$range,
                               y_max = lims$max,
                               y_min = lims$min,
                               include_zero = FALSE)
      }
      
      plot_value <- ggplot2::ggplot() +
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
          na.rm = TRUE) +
        
        ggplot2::geom_boxplot(
          ggplot2::aes(
            x = as.character(nestedData_value[[group_df_name]]),
            y = nestedData_value[[abundance]],
            fill = nestedData_value[[group_df_name]]),
          alpha = 0.2,
          position = "dodge2",
          na.rm = TRUE) +
        ggplot2::labs(x = NULL)
      
      if (!is.null(lims[[1]])){
        plot_value <- plot_value + ggplot2::coord_cartesian(ylim = setlims)
      }
      return(plot_value)
    }))
  }
  
  
  #######
  #######
  
  parallel::stopCluster()

  
  ##### Process plots
  
  pvs <- data.frame(unique(as.character(trellData$summary_stats[[panel_variable]])))
  
  if(nrow(pvs) == 0){
    pvs <- data.frame(unique(as.character(trellData$data_values[[panel_variable]])))
  }
  
  colnames(pvs) <- panel_variable
  
  pvs <- dplyr::mutate(pvs, panel = trelliscopejs::pmap_plot(All_plots, function(abundance_global = NULL, 
                                                                          foldchange_global = NULL, 
                                                                          missing_bar = NULL, 
                                                                          abundance_heatmap = NULL, 
                                                                          foldchange_heatmap = NULL, 
                                                                          presence_heatmap = NULL, 
                                                                          abundance_boxplot = NULL, 
                                                                          foldchange_bar = NULL){
    
    plots <- list(
      abundance_global = abundance_global,
      foldchange_global = foldchange_global,
      abundance_heatmap = abundance_heatmap,
      foldchange_heatmap = foldchange_heatmap,
      presence_heatmap = presence_heatmap,
      missing_bar = missing_bar,
      abundance_boxplot = abundance_boxplot,
      foldchange_bar = foldchange_bar)
    
    plots <- plots[-which(sapply(plots, is.null))]
    
    if(interactive == FALSE){
      if (length(plot_type) > 3){
  
        out_plot <- patchwork::wrap_plots(plots)
        
      } else if (length(plot_type) == 3) {
        
        out_plot <- ggpubr::ggarrange(plotlist = plots,
                                      ncol = 2,
                                      nrow = 2
        )
        
      } else if (length(plot_type) == 2) {
        
        out_plot <- ggpubr::ggarrange(plotlist = plots,
                                      ncol = 1,
                                      nrow = 2
        )
        
      } else {
        out_plot <- ggpubr::ggarrange(plotlist = plots,
                                      ncol = 1,
                                      nrow = 1
        )
      }
    } else {
      
      if (any(stringr::str_detect(names(plots), "heatmap"))){
        heatplots <- names(plots)[stringr::str_detect(names(plots), "heatmap")]
        plots[heatplots] <- purrr::map(heatplots, function(heatitem){
          plots[[heatitem]] <- plotly::ggplotly(plots[[heatitem]] + ggplot2::theme(legend.position = "none"), 
                                                tooltip = c("text"))
          plots[[heatitem]] <- plots[[heatitem]] %>% plotly::layout(plot_bgcolor='grey',
                                                                    xaxis = list(showgrid = F),
                                                                    yaxis = list(showgrid = F))
          return(plots[[heatitem]])
        })
      }
      
      if (length(plot_type) > 6){
        out_plot <- plotly::subplot(plots, nrows = 3, margin = 0.05, titleX = TRUE, titleY = TRUE)
      } else if (length(plot_type) < 6 && length(plot_type) != 1) {
        out_plot <- plotly::subplot(plots, nrows = 2, margin = 0.05, titleX = TRUE, titleY = TRUE)
      } else {
        out_plot <- plotly::ggplotly(plots[[1]])
      }
      
      for(num in 1:length(out_plot$x$data)){
        out_plot$x$data[[num]]$text <- gsub(
          "nested(\\w|[[:punct:]])+[[:space:]]+(\\w|[[:punct:]])+<br \\/>", 
          "", 
          out_plot$x$data[[num]]$text)
        out_plot$x$data[[num]]$text <- gsub(
          "as.\\w+\\(",
          "", 
          out_plot$x$data[[num]]$text)
        out_plot$x$data[[num]]$text <- gsub(
          "nested.+", 
          "", 
          out_plot$x$data[[num]]$text)

      }
      
      
    }
    
    return(out_plot)
  }))
  
  
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
      limlist[["min"]] <- y_limits[1]
      limlist[["max"]] <- y_limits[2]
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
#' @param omicsData A pmartR object of class pepData, lipidData, metabData, or proData
#' @param omicsStats A statistical results object produced by running \code{imd_anova} on omicsData.
#' @param p_val Numeric that specifies p-value for Boolean significance cognotic. Default is 0.05.
#' @param mapping_col String: For proData - name of column with peptide information. For pepData - name of column with protein information. Default is NULL.
#' @param panel_variable String: Name of column that plot panels are sorted by (e.g. each plotting arrangement has a unique identifier from panel variable). Default is emeta_cname if present, edata_cname where emeta_cname is not present.
#' @param try_URL Will attempt to link to PubChem, LipidMaps, or Uniprot based on information in edata_cname of omicsData or specified mapping_col for peptide data. Default is FALSE.
#'
#'
#' @author Rachel Richardson
#' @export
data_cogs <- function(...) {
  .data_cogs( ...)
}

.data_cogs <- function(nested_plot = NULL, omicsData = NULL, 
                       omicsStats = NULL,
                      p_val = 0.05,
                      mapping_col = NULL, panel_variable = NULL, 
                      try_URL = FALSE){
  
  
  tictoc::tic("Initial checks")
  
  if(class(omicsData) == "list" & length(omicsData) == 1){
    omicsData <- omicsData[[1]]
  } else if (class(omicsData) == "list") stop (
    "Lists are not supported for data_cogs omicsData input"
  )
  if(class(omicsStats) == "list" & length(omicsStats) == 1){
    omicsStats <- omicsStats[[1]]
  } else if (class(omicsStats) == "list") stop (
    "Lists are not supported for data_cogs omicsStats input"
  )
  
  ## Switch Stats and Omics data as appropriate ##
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
  
  validate_omics_input(omicsData = omicsData, omicsStats = omicsStats)
  
  ## Moniker Variables ##
  e_data <- omicsData$e_data
  f_data <- omicsData$f_data
  e_meta <- omicsData$e_meta
  stats <- omicsStats$Full_results
  
  ## Variable checks ##
  
  # Check check if p_val is numeric of length 1 #
  if(!is.numeric(p_val) | (length(p_val) != 1)) stop(
    "p_val must be a numeric of length 1")  
  
  # Check check if try_URL is boolean of length 1 #
  if(!is.logical(try_URL) | (length(try_URL) != 1)) stop(
    "try_URL must be a TRUE/FALSE of length 1")  
  
  # Ensure mapping_col is length one #
  if(!is.null(mapping_col) && length(mapping_col) != 1) stop(
    "mapping_col must be a string of length 1")  
  
  # Ensure mapping_col is only called where e_meta exists #
  if(!is.null(mapping_col) && is.null(e_meta)) stop(
    "No e_meta data for mapping_col to be called on. 
    Using mapping_col requires e_meta in omicsData.")  
  
  # Ensure mapping_col is in emeta #
  if(!is.null(mapping_col) && 
     inherits(omicsData, "pepData") && 
     !(mapping_col %in% colnames(e_meta))) stop(
    "Invalid entry: mapping_col must be in the column names of pepData e_meta.")  
  
  # Warn in mapping_col is used where data is not pro/pep #
  if(!is.null(mapping_col) && 
     !inherits(omicsData, c("pepData", "proData"))){
    message(paste("Notice: mapping_col is not used for omicsData of class", class(omicsData)))
  }

  
  tictoc::toc()
  
  tictoc::tic("Fill")
  ## Fill null variables as needed ##

  if (is.null(nested_plot)){
    trellData = as.trellData(omicsData, omicsStats)
    # recursive for lists #
    if(class(trellData) == "list"){
      nested_plot <- purrr::map(trellData, format_plot)
      nestlist <- purrr::pmap(list(nested_plot, omicsData, omicsStats), data_cogs)
      return(nestlist)
    } else {
      nested_plot <- format_plot(trellData)
    }
  }

  ## Assign unique ID, then panel variable ##
  uniqueID <- attributes(omicsData)$cnames$edata_cname
  if(is.null(uniqueID)){
    uniqueID <- attributes(omicsStats)$cnames$edata_cname
  }
  panel_variable <- names(nested_plot[1])
  
  tictoc::toc()

  tictoc::tic("orig dat")
  if (!is.null(omicsStats) & !is.null(omicsData)){
    addOrigCogs <- suppressWarnings(dplyr::left_join(stats, e_data))
    if(!is.null(e_meta)){
      addOrigCogs <- suppressWarnings(dplyr::left_join(addOrigCogs, e_meta))
    }
  } else if (!is.null(omicsStats)){
    addOrigCogs <- stats
  } else {
    addOrigCogs <- e_data
    if(!is.null(e_meta)){
      addOrigCogs <- suppressWarnings(dplyr::left_join(addOrigCogs, e_meta))
    }
  }
  tictoc::toc()
  
  tictoc::tic("Additional stats/emeta")
  ## Stats Data Cogs ##
  if (!is.null(stats)){
    uniqlist <- stats[[uniqueID]]
    if (!is.null(e_meta)){
      joiner <- merge(stats, e_meta)
    } else {
      joiner <- stats
    }
  
    uniqlist <- joiner[[uniqueID]]
    panel_list <- joiner[[panel_variable]]
    
    addcogs <- data.frame(uniqlist)
    colnames(addcogs) <- uniqueID
    
    if(any(stringr::str_detect(names(joiner), "P_value_G"))){
      
      pvalg <- data.frame(t(joiner[stringr::str_detect(names(joiner), "P_value_G")]))
      colnames(pvalg) <- substring(colnames(pvalg), 2)
      Sig_G_p_value <- as.factor(uniqlist %in% uniqlist[as.numeric(names(
        dplyr::select_if(pvalg, function(x) any(x < p_val & !is.na(x)))
      ))])
      
      addcogs[["Sig_G_p_value"]] <- trelliscopejs::cog(Sig_G_p_value, desc = "Boolean for significant G test p-values,
      where TRUE indicates that the p-value is < 0.05 and is grounds for 
      rejecting the null hypothesis (independence of missing data).")
      
    }
    
    if(any(stringr::str_detect(names(joiner), "P_value_T"))){
      
      pvalt <- data.frame(t(joiner[stringr::str_detect(names(joiner), "P_value_T")]))
      colnames(pvalt) <- substring(colnames(pvalt), 2)
      Sig_T_p_value <- as.factor(uniqlist %in% uniqlist[as.numeric(names(
        dplyr::select_if(pvalt, function(x) any(x < p_val & !is.na(x)))
      ))])
      
      addcogs[["Sig_T_p_value"]] <- trelliscopejs::cog(Sig_T_p_value, desc = "Boolean for significant T test p-values,
      where TRUE indicates that the p-value is < 0.05 and is grounds for 
      rejecting the null hypothesis (no significant 
      difference between groups).")
      
    }

    # addcogs <- data.frame(
    #   uniqlist,
    #   trelliscopejs::cog(Sig_G_p_value, desc = "Boolean for significant G test p-values,
    #   where TRUE indicates that the p-value is < 0.05 and is grounds for 
    #   rejecting the null hypothesis (independence of missing data)."),
    #   trelliscopejs::cog(Sig_T_p_value, desc = "Boolean for significant T test p-values,
    #   where TRUE indicates that the p-value is < 0.05 and is grounds for 
    #   rejecting the null hypothesis (no significant 
    #   difference between groups)."))
    # colnames(addcogs) <- c(uniqueID, 'Sig_G_p_value', 'Sig_T_p_value')
    
    if (!identical(uniqlist, panel_list)){
      addcogs[[panel_variable]] <- panel_list
    }
  } else {
    # omicsData only #
    # if e_meta #
    if (!is.null(e_meta)){
      uniqlist <- e_meta[[uniqueID]]
      panel_list <- e_meta[[panel_variable]]
      if (!identical(uniqlist, panel_list) & !is.null(panel_list)){
        addcogs <-  data.frame(uniqlist, panel_list)
        colnames(addcogs) <- c(uniqueID, panel_variable)
      } else {
        addcogs <-  data.frame(uniqlist)
        colnames(addcogs) <- uniqueID
      }
    # no emeta #
    } else {
      uniqlist <- e_data[[uniqueID]]
      panel_list <- e_data[[panel_variable]]
      if (!identical(uniqlist, panel_list) & !is.null(panel_list)){
        addcogs <-  data.frame(uniqlist, panel_list)
        colnames(addcogs) <- c(uniqueID, panel_variable)
      } else {
        addcogs <-  data.frame(uniqlist)
        colnames(addcogs) <- uniqueID
      }
    }
  }
  
  tictoc::toc()
  
  
  tictoc::tic("additional cognostics")
  ## Type Specifc Cognostics ##
  
  ## pepData ##
  if ("pepData" %in% attr(nested_plot, "parent_class")){

    if (!is.null(e_meta) & !is.null(mapping_col)){
      
      # Degenerate Peptides #
      pepcol <- attr(omicsData, "cnames")$edata_cname
      mapcols <- e_meta[c(pepcol, mapping_col)]
      mapcols <- unique(mapcols)
      
      peps <- mapcols[pepcol]
      degenpep <- peps[which(duplicated(peps)==TRUE),]
      Is_degenerate <- as.factor(
        purrr::map(peps, function(pep) pep %in% degenpep)[[pepcol]])
      
        # Add to dataframe #
      addcogs <- data.frame(
        addcogs,
        as.data.frame(Is_degenerate))
      }
      
      # Try to find Protein URL #
      if ((try_URL == TRUE)){
        URLlist <- e_meta[[mapping_col]]
        searchname <- stringr::str_extract(URLlist, "[A-Z0-9]+_[A-Z]+")
        if (all(is.na(searchname))){
          searchname <- stringr::str_extract(URLlist, "[A-Z0-9]{6,}")
        }
        if (!all(is.na(searchname))){
          searchlink <- paste('https://www.uniprot.org/uniprot/', 
                              searchname, sep = "")
        }
        addcogs <- data.frame(
          addcogs,
          Protein_URL = trelliscopejs::cog_href(searchlink, 
                                 desc = "UniProt lookup using PRIDE ID 
                               or Entry name (XXXX_XXXX)"))
      }
    
    ## proData ##
  } else if ("proData" %in% attr(nested_plot, "parent_class")){
    
    # Degenerate peptides #
    if (!is.null(e_meta) & !is.null(mapping_col)){
      
      # Degenerate Peptides #
      procol <- attr(omicsData, "cnames")$edata_cname
      mapcols <- e_meta[c(procol, mapping_col)]
      mapcols <- unique(mapcols)
      
      peps <- mapcols[mapping_col]
      degenpep <- peps[which(duplicated(peps)==TRUE),]
      mapcols[mapping_col] <- purrr::map(peps, function(pep) pep %in% degenpep)[[mapping_col]]
      mapcols <- mapcols %>% nest(-procol)
      mapcols <- dplyr::mutate(mapcols, 
                        n_degenerate = purrr::map_int(mapcols$data, function(propeps){
                          sum(propeps[[mapping_col]])
                          }))
      mapcols$data <- NULL
      
      # Add to dataframe #
      addcogs <- suppressWarnings(dplyr::left_join(
        addcogs,
        mapcols))
    }
    
    # Try to find Protein URL #
    if (try_URL == TRUE){
      searchname <- stringr::str_extract(uniqlist, "[A-Z0-9]+_[A-Z]+")
      if (all(is.na(searchname))){
        searchname <- stringr::str_extract(uniqlist, "[A-Z0-9]{6,}")
      }
      if (!all(is.na(searchname))){
        searchlink <- paste('https://www.uniprot.org/uniprot/', 
                            searchname, sep = "")
      }
      addcogs <- data.frame(
        addcogs,
        Protein_URL = trelliscopejs::cog_href(searchlink, 
                               desc = "UniProt lookup using PRIDE ID 
                                 or Entry name (XXXX_XXXX)"))
    }
    
    ## lipidData ##
  } else if ("lipidData" %in% attr(nested_plot, "parent_class")){
    
    # Search LipidMaps with the lipid name #
    if (try_URL == TRUE){
      searchname <- gsub(" ", "", uniqlist)
      searchname <- gsub(":", " ", searchname) %>% 
        stringr::str_remove_all("[[:punct:]]")
      searchname <- gsub(" ", ":", searchname)
      searchlink <- paste(
        'https://www.lipidmaps.org/data/structure/LMSDFuzzySearch.php?Name=',
        unlist(searchname), 
        sep = "")
      
      addcogs <- data.frame(
        addcogs,
        Lipid_Map_Hits = trelliscopejs::cog_href(searchlink, 
                                  desc = "Lipid Maps hits for uniqueID name"))
    }
    
    ## metabData ##  
  } else if ("metabData" %in% attr(nested_plot, "parent_class")){
    
    # PubChem URL #
    if (try_URL == TRUE){
      searchlink <- paste(
        "https://pubchem.ncbi.nlm.nih.gov/compound/",
        uniqlist,
        sep = "")
      addcogs <- data.frame(
        addcogs,
        PubChem_URL = trelliscopejs::cog_href(searchlink,
                               desc = "PubChem page for uniqueID name"))
    }
  }

  tictoc::toc()
  
  tictoc::tic("joining with nested data - map")
  
  allcogs <- suppressWarnings(dplyr::left_join(addcogs, addOrigCogs) %>% tidyr::nest(-panel_variable))
  
  transcogs <- list(purrr::map(allcogs$data, function(panel_data){
    if (nrow(panel_data) >1 ){
      datalist <-  list()
      changecol <- c(rep(FALSE, ncol(panel_data)))
      for (colnum in 1:ncol(panel_data)){
        panel_dat <- unlist(panel_data[colnum])
        
        # Cat of factor columns that are not T/F #
        if (is.factor(panel_dat) && 
            !stringr::str_detect(names(panel_data[colnum]), "Sig_[TG]_p_value") &&
            !stringr::str_detect(names(panel_data[colnum]), "Is_degenerate")){
          # datalist[colnum] <- capture.output(
          #   cat(unique(levels(panel_dat)[panel_dat]), sep = ", "))
          datalist[colnum] <- toString(unique(levels(panel_dat)[panel_dat]))
          # Takes mean of numeric columns, includes numeric strings not in cnames #
        } else if ((is.numeric(panel_dat) | 
                    !any(is.na(as.numeric(stats::na.omit(panel_dat))))) &&
                   !(names(panel_dat) %in% attr(omicsStats, "cnames")) && 
                   !stringr::str_detect(names(panel_data[colnum]), "Sig_[TG]_p_value")&&
                   !stringr::str_detect(names(panel_data[colnum]), "Is_degenerate")){
          changecol[colnum] <- TRUE
          datalist[colnum] <- mean(as.numeric(stats::na.omit(panel_dat)))
          
          # Takes logicals, returns TRUE if any TRUE #
        } else if (!any(is.na(as.logical(panel_dat)))){
          datalist[colnum] <- as.character(
            any(as.logical(panel_dat)))
        } else {
          # Cats characters, other #
          # datalist[colnum] <- capture.output(
          #   cat(unique(panel_dat), sep = ", "))
          datalist[colnum] <- toString(unique(panel_dat))
        }
      }
    } else {
      return(panel_data)
    }
      panelcogs <- data.frame(lapply(rbind(datalist), unlist))
      colnames(panelcogs) <- colnames(panel_data)
      colnames(panelcogs[changecol]) <- paste(colnames(panelcogs[changecol]), "_mean")
      return(panelcogs)

    
  }
  
  ))
  tictoc::toc()
  

  tictoc::tic("rbind")
  transcogs <- do.call(dplyr::bind_rows, transcogs)
  # transcogs <- data.table::rbindlist(transcogs)
  allcogs$data <- NULL
  # allcogs <- data.table::data.table(allcogs)
  
  tictoc::toc()
  
  tictoc::tic("cbind")
  allcogs <- cbind(allcogs, transcogs)
  tictoc::toc()
  tictoc::tic("leftjoin")
  # nested_plot <- allcogs[nested_plot]
  
  nested_plot <- suppressWarnings(dplyr::left_join(nested_plot, allcogs))
  tictoc::toc()
  
  
  return(nested_plot)
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
#' @param mapping_col String: For associated proData/pepData - name of column in peptide data with protein information. Default is NULL.
#' @param panel_variable String: Name of column that plot panels are sorted by (e.g. each plotting arrangement has a unique identifier from panel variable). Default is emeta_cname if present, edata_cname where emeta_cname is not present.
#' @param try_URL Will attempt to link to PubChem, LipidMaps, or Uniprot based on information in edata_cname of omicsData or specified mapping_col for peptide data. Default is FALSE.
#' @param trelli_name String: name of display, or list of names where a list is provided for omicsData and omicsStats
#' @param trelli_path_out String: path to where trelliscope is stored. Default is "./TrelliDisplay"
#' @param interactive Should the plots be rendered as plotly objects?
#' @param plot_text Disable plot text
#' @param y_limits Y limits 
#' @param plot_type plots for plotting
#'
#' @author Rachel Richardson
#' @export
trelliVis <- function(...) {
  .trelliVis(...)
}

.trelliVis <- function(omicsData = NULL, omicsStats = NULL,
                       omicsFormat = NULL, p_val = 0.05, 
                       mapping_col = NULL, panel_variable = NULL, 
                       try_URL = FALSE, trelli_name = NULL,
                       trelli_path_out = "TrelliDisplay", 
                       plot_text = FALSE, interactive = FALSE,
                       y_limits = NULL, plot_type = NULL) {
  
  
  #store_object, custom_cog_df, plot package = ggplot, rbokeh, etc, trelliscope additional arguments
  
  #####
  ## Switch Stats and Omics data as appropriate ##
  #####
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
  # If lists of omicsData or trellData are used with mapping col, #
  # make sure panel variables are specified for both #
  if ((class(omicsData) == "list" | 
       class(omicsStats) == "list") && 
      !is.null(mapping_col) && 
      length(mapping_col) != max(length(omicsStats), 
                                    length(omicsData))) stop(
                                      "mapping_col must be specified for each index in omicsStats/omicsData"
                                    )
  
  # Check check if try_URL is boolean of length 1 #
  if(!is.logical(try_URL) | (length(try_URL) != 1)) stop(
    "try_URL must be a TRUE/FALSE of length 1")  
  
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
  
  if (is.null(omicsFormat)){
    # Re-format objects for plotting #
    trellData <- as.trellData(omicsData, omicsStats)
  } else {
    trellData <- omicsFormat
  }
  

  # If a pep/pro pair is listed, act on each item #
  if(class(trellData) == "list"){
    
    
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
      trelli_name <- attr(trellData, "data_types")
      
      # if trelli_name is not of length list, add identifiers #
    } else if (length(trelli_name) == 1){
      trelli_name <- paste(trelli_name, 
                           attr(trellData, "data_types"),
                           sep = "_")
      
      # if trelli_name too long, cut to length(trellData)  #
    } else if (length(trelli_name) > length(trellData)) {
      warning( "Length trelli_name exceeds the number of generated displays. 
            Only the first %d names will be used.", length(trellData))
      trelli_name <- trelli_name[1:length(trellData)]
      
      # if trelli_name too short and not 1, throw error  #
    } else if (length(trelli_name) < length(trellData)) stop(
      paste("Length of trelli_name != 1 and less than length of trellData; 
            please use a trelli_name of length 1 or length", length(trellData)))
    
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
    if (is.null(mapping_col)){
      mapping_col <- rep(list(NULL), length(trellData))
    }
    
    # Nest data and generate trelliscope plots #
    nested_plot <- purrr::map2(trellData, panel_variable, function(pairedplotter, pan){
      format_plot(trellData = pairedplotter, 
                  p_val = p_val,
                  panel_variable = pan, 
                  plot_type = plot_type,
                  plot_text = plot_text,
                  y_limits = y_limits,
                  interactive = interactive)
      })
    
    # Generate default cognostics #
    nest_plot_cog_list <- purrr::pmap(
      list(nested_plot,
           omicsData,
           omicsStats, 
           panel_variable,
           mapping_col),
      function(nest, dat, stat, pan, map){
        suppressWarnings(data_cogs(nested_plot = nest, 
                  omicsData = dat,  
                  omicsStats = stat,
                  p_val = p_val,
                  mapping_col = map, 
                  panel_variable = pan, 
                  try_URL = try_URL))
        })

    # Generate trelliscope display #
    purrr::map2(nest_plot_cog_list, 
              trelli_name, function(display, name){
                trelliscopejs::trelliscope(display, as.character(name), nrow = 1, ncol = 2,
                  path = as.character(trelli_path_out), thumb = TRUE, state = list(
                    sort = list(trelliscopejs::sort_spec(names(display[1]), dir = "desc")), 
                    labels = list(names(display[1]))))
    })
    
  # Where not a pep/pro pair: #
  } else {
    
    if (!pmartR::get_data_norm(trellData) &&
        (is.null(attributes(trellData)$isobaric_info$norm_info$is_normalized) || #### update helper function
         attributes(trellData)$isobaric_info$norm_info$is_normalized != TRUE)) stop(
           "Input must be normalized prior to plotting; use normalize_global or normalize_isobaric as appropriate."
         )
    
    # Default: Name the display after parent classes (omicsData and omicStats clases) #
    if(is.null(trelli_name)){
      # trelli_name <- capture.output(
      #   cat(attr(trellData, "parent_class"), sep = "_"))
      trelli_name <- paste(attr(trellData, "parent_class"), collapse = "_")
      # Make sure trelliname is length 1 #
    } else if (length(trelli_name) != 1) {
      warning("trelli_name length is greater than one where only one display is generated; using only the first entry in trelli_name.")
      trelli_name <- trelli_name[[1]]
    }
    
    # Nest data and generate trelliscope plots #
    nested_plot <- format_plot(trellData = trellData, 
                               p_val = p_val,
                               panel_variable = panel_variable, 
                               plot_type = plot_type,
                               plot_text = plot_text,
                               y_limits = y_limits,
                               interactive = interactive)

    
    # Generate default cognostics #
    nest_plot_cog <- suppressWarnings(data_cogs(nested_plot = nested_plot, 
                               omicsData = omicsData, 
                               omicsStats = omicsStats, 
                               p_val = p_val, 
                               mapping_col = mapping_col, 
                               panel_variable = panel_variable, 
                               try_URL = try_URL))
    
    # Generate trelliscope display #
    nest_plot_cog %>%
      trelliscopejs::trelliscope(name = as.character(trelli_name), nrow = 1, ncol = 2,
                  path = as.character(trelli_path_out), thumb = TRUE,
                  state = list(
                    sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")), 
                    labels = list(names(nest_plot_cog[1]))))
  }
}
