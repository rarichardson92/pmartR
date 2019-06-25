#' @name format_data
#' @rdname format_data
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
format_data <- function(...){
  .format_data(...)
}

.format_data <- suppressWarnings(function(omicsData = NULL, omicsStats = NULL){
  
  ## Checks and recursive for lists as input in format_data ##
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
#' @title Recursive call of format_data and associated checks for list inputs
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
      plotterlist <- purrr::map2(omicsData, omicsStats, format_data)
      
    } else {
      plotterlist <- format_data(omicsData[[1]], omicsStats[[1]])
      return(plotterlist)
    }
    
    #--omicsData--#
  } else if(!is.null(omicsData)){
    if (length(omicsData) != 1){
      
      validate_omics_input(omicsData = omicsData, omicsStats = omicsStats)
      classlist <- omicsData %>% purrr::map(function(omics) attr(omics, which = "class"))
      plotterlist <- purrr::map(omicsData, format_data)
      
    } else {
      plotterlist <- format_data(omicsData[[1]])
      return(plotterlist)
    }
    
    #--omisStats--#
  } else {
    if (length(omicsStats) > 1){
      
      validate_omics_input(omicsData = omicsData, omicsStats = omicsStats)
      classlist <- omicsStats %>% purrr::map(function(omics) attr(omics, which = "data_class"))
      plotterlist <- purrr::map(omicsStats, format_data)
      
    } else {
      plotterlist <- format_data(omicsStats[[1]])
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
  
  
  ## Checks and recursive for lists as input in format_data ##
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
      "format_data() requires at least one of the following: 
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
#' @param trellData An object of class "trellData" generated from \code{\link{format_data}}.
#' @param comps_y_limits For comparisons: Set to "fixed" or "free" for automated y-axis calculating. "fixed" - axis generated based on the maximum/minimum across all plots. "free" - axis axis generated based on the maximum/minimum of individual plot.
#' @param comps_y_range For comparisons: Specify a range for the plot y-axis. Will calculate the range based on one of y_max or y_min parameters or from the median of y-values where y_max and y_min are not defined.
#' @param comps_y_max For comparisons: Sets the maximum y-value for the y-axis.
#' @param comps_y_min For comparisons: Sets the minimum y-value for the y-axis.
#' @param value_y_limits For values: Set to "fixed" or "free" for automated y-axis calculating. "fixed" - axis generated based on the maximum/minimum across all plots. "free" - axis axis generated based on the maximum/minimum of individual plot.
#' @param value_y_range For values: Specify a range for the plot y-axis. Will calculate the range based on one of y_max or y_min parameters or from the median of y-values where y_max and y_min are not defined.
#' @param value_y_max For values: Sets the maximum y-value for the y-axis.
#' @param value_y_min For values: Sets the minimum y-value for the y-axis.
#' @param p_val Specifies p-value for setting graph border colors
#' @param panel_variable Specifies what to divide trelliscope panels by, must be a column in trellData. Defaults to cnames$edata_cname of trellData.
#' @param comps_panel_x_axis Specifies what column should be on the x-axis, must be a column in trellData. Default setting plots pairwise comparisons along x-axis.
#' @param comps_panel_y_axis Specifies what column should be on the y-axis, must be a column in trellData and numeric. Default setting plots fold change for combined and anova testing and counts for g-test.
#' @param comps_color_variable Specifies what column should distingush color, must be a column in trellData. Default settings is set to "Group."
#' @param value_panel_x_axis Specifies what column should be on the x-axis, must be a column in trellData. Default setting plots pairwise comparisons along x-axis.
#' @param value_panel_y_axis Specifies what column should be on the y-axis, must be a column in trellData and numeric. Default setting plots fold change for combined and anova testing and counts for g-test.
#' @param value_color_variable Specifies what column should distingush color, must be a column in trellData. Default settings is set to "Group."
#' @param value_plot_type For values: Specifies plot types for graphing; must be a list of strings where "box", "bar", or "point" is specified. Combined strings like "boxpoint" will plot both in the same graph.
#' @param comps_plot_type For comparisons: Specifies plot types for graphing; must be a list of strings where "box", "bar", or "point" is specified. Combined strings like "boxpoint" will plot both in the same graph.
#' @param comps_include_zero For comparisons: Should plots show y = 0?
#' @param value_include_zero For values: Should plots show y = 0?
#' @param value_plot For values: In the presence of data_values in trellData, should the plot be rendered? Default is TRUE.
#' @param comps_plot For comparisons: In the presence of summary_stats and comp_stats in trellData, should the plot be rendered? Default is TRUE.
#' @param comps_text For comparisons only: TRUE/FALSE for p-value text above data
#' @param plotly Should the plots be rendered as plotly objects?
#' @param ... further arguments
#'
#' @author Rachel Richardson
#' @export
format_plot <- function(trellData, ...) {
   .format_plot(trellData,  ...)
}

.format_plot <- function(trellData, 
                        comps_y_limits = NULL, comps_y_range = NULL, 
                        comps_y_max = NULL, comps_y_min = NULL, 
                        comps_include_zero = TRUE,
                        value_y_limits = NULL, value_y_range = NULL, 
                        value_y_max = NULL, value_y_min = NULL,
                        value_include_zero = FALSE,
                        p_val = 0.05,
                        panel_variable = attributes(trellData)$cnames$edata_cname, 
                        comps_color_variable = NULL,
                        comps_panel_x_axis = "Comparison", comps_panel_y_axis = NULL,
                        value_color_variable = NULL,
                        value_panel_x_axis = NULL, 
                        value_panel_y_axis = NULL,
                        value_plot_type = "boxpoint", comps_plot_type = "col",
                        value_plot = TRUE, comps_plot = TRUE, 
                        comps_text = TRUE, plotly = TRUE) {

  tictoc::tic("initial check")
  ## Initial Checks ##
  validate_format_plot_input(
    trellData,
    comps_y_limits = comps_y_limits, 
    comps_y_range = comps_y_range, 
    comps_y_max = comps_y_max, 
    comps_y_min = comps_y_min, 
    comps_include_zero = comps_include_zero,
    value_y_limits = value_y_limits, 
    value_y_range = value_y_range, 
    value_y_max = value_y_max, 
    value_y_min = value_y_min,
    value_include_zero = value_include_zero,
    p_val = p_val,
    value_plot_type = value_plot_type, 
    comps_plot_type = comps_plot_type,
    value_plot = value_plot,
    comps_plot = comps_plot,
    comps_text = comps_text,
    plotly = plotly)
  
  tictoc::toc()
  
  tictoc::tic("assign vars")
  ## Re-assign Variables (if not specified) ##
  if (is.null(panel_variable)){
    panel_variable  <- attributes(trellData)$cnames$edata_cname
  }
  
  # --Comps-- #
  
  if (!is.null(trellData$summary_stats) & 
      !is.null(trellData$comp_stats)){
    
    # Sets stats test and colors for borders #
    option <- attr(trellData, "statistical_test")
    colors <- c("NA", "darkgrey", "black")
    
    # Sets default y-values based on stats test #
    if(is.null(comps_panel_y_axis) & (option == "gtest")){
      comps_panel_y_axis <- "Count"
    } else if (is.null(comps_panel_y_axis)){
      comps_panel_y_axis <- "Fold_change"
    } else {
      comps_panel_y_axis <- comps_panel_y_axis
    }
    
    # Sets default value_color_variable to Group defined in group_DF attribute #
    if (is.null(comps_color_variable)){
      if ("Group_DF" %in% colnames(trellData$summary_stats)){
        comps_color_variable <- "Group_DF"
      } else {
        comps_color_variable <- "Group"
      }
    }
    
    # Correct for 'Group_DF' misinput # 
    if (comps_color_variable == "Group_DF" & 
        !("Group_DF" %in% colnames(trellData$summary_stats)) &
        "Group" %in% colnames(trellData$summary_stats)) {
      comps_color_variable <- "Group"
    }
    
    tictoc::toc()
    
    tictoc::tic("second val")
    
    ## Validate potential re-assignments ##
    validate_format_plot_input(trellData,
                               panel_variable = panel_variable, 
                               comps_color_variable = comps_color_variable,
                               comps_panel_x_axis = comps_panel_x_axis,
                               comps_panel_y_axis = comps_panel_y_axis)
    
    tictoc::toc()
  }
  
  # --Value-- #
  if (!is.null(trellData$data_values)) {
    
    # Set default value_panel_y_axis #
    if (is.null(value_panel_y_axis)){
      value_panel_y_axis  <-  grep("abundance", 
                                   colnames(trellData$data_values), 
                                   value = TRUE)
    }
    
    # Sets default value_color_variable to Group defined in group_DF attribute #
    if (is.null(value_color_variable)){
      if ("Group_DF" %in% colnames(trellData$data_values) &
          "Group" %in% colnames(trellData$data_values)){
        value_color_variable <- "Group_DF"
      } else {
        value_color_variable <- "Group"
      }
    }
    
    # Sets default value_panel_x_axis to Group defined in group_DF attribute #
    if (is.null(value_panel_x_axis)){
      if ("Group_DF" %in% colnames(trellData$data_values) &
          "Group" %in% colnames(trellData$data_values)){
        value_panel_x_axis <- "Group_DF"
      } else {
        value_panel_x_axis <- "Group"
      }
    }
    
    # Correct for 'Group_DF' misinput for value_color_variable or value_panel_x_axis # 
    if (value_color_variable == "Group_DF" & 
        !("Group_DF" %in% colnames(trellData$data_values)) &
        "Group" %in% colnames(trellData$data_values)) {
      value_color_variable <- "Group"
    }
    if (value_panel_x_axis == "Group_DF" & 
        !("Group_DF" %in% colnames(trellData$data_values)) &
        "Group" %in% colnames(trellData$data_values)) {
      value_panel_x_axis <- "Group"
    }
    
    tictoc::tic("second val")
    ## Validate potential re-assignments ##
    validate_format_plot_input(trellData,
                               panel_variable = panel_variable,
                               value_color_variable = value_color_variable,
                               value_panel_x_axis = value_panel_x_axis,
                               value_panel_y_axis = value_panel_y_axis)
    
    tictoc::toc()
  }
  
  
  tictoc::tic("Message")
  ## Inform user of selected plotting parameters ##
  generate_plot_message(trellData,
                        comps_y_limits = comps_y_limits, 
                        comps_y_range = comps_y_range,
                        comps_y_max = comps_y_max, 
                        comps_y_min = comps_y_min,
                        value_y_limits = value_y_limits, 
                        value_y_range = value_y_range,
                        value_y_max = value_y_max, 
                        value_y_min = value_y_min)
  
  tictoc::toc()
  
  ## Generate specified y-limits for plotting (if applicable) ##
  #--Values--#
  if(!is.null(trellData$data_values)){
    
    # Set ylims if appropriate #
    if (!is.null(value_y_limits)){
      if (value_y_limits == 'fixed'){
        increment <- set_increment(trellData$data_values[[value_panel_y_axis]], 
                                   include_zero = value_include_zero)
        ylims <- set_ylimits(trellData$data_values[[value_panel_y_axis]], 
                             increment = increment,
                             y_min = value_y_min, y_max = value_y_max, 
                             include_zero = value_include_zero)
      }
    } else if (is.null(value_y_limits) & is.null(value_y_range)){
      increment <- set_increment(trellData$data_values[[value_panel_y_axis]], 
                                 include_zero = value_include_zero)
      ylims <- set_ylimits(trellData$data_values[[value_panel_y_axis]], 
                           increment = increment,
                           y_min = value_y_min, y_max = value_y_max, 
                           include_zero = value_include_zero)
    }
  }
  
  #--Comps--#
  if(!is.null(trellData$comp_stats)){
    # Generate increment and y limits based on parameters and y-values #
    if (!is.null(comps_y_limits)){
      if (comps_y_limits == 'fixed'){
        comps_increment <- set_increment(trellData$comp_stats[[comps_panel_y_axis]], 
                                         include_zero = comps_include_zero)
        comps_ylims <- set_ylimits(trellData$comp_stats[[comps_panel_y_axis]], 
                                   increment = comps_increment,
                                   y_min = comps_y_min, y_max = comps_y_max,
                                   include_zero = comps_include_zero)
      }
    } else if (is.null(comps_y_limits) & is.null(comps_y_range)){
      comps_increment <- set_increment(trellData$comp_stats[[comps_panel_y_axis]], 
                                       include_zero = comps_include_zero)
      comps_ylims <- set_ylimits(trellData$comp_stats[[comps_panel_y_axis]], 
                                 increment = comps_increment,
                                 y_min = comps_y_min, y_max = comps_y_max,
                                 include_zero = comps_include_zero)
    }
    
    
  ## Generate nested data ##
    
    if ("Group_DF" %in% colnames(trellData$summary_stats)) {
      group_df_name <- "Group_DF"
    } else {
      group_df_name <- "Group"
    }
    
    tictoc::tic("make comp df")
    
    plotter <- tidyr::separate(trellData$comp_stats, Comparison,
                               c("comp1", "comp2"), sep = "_vs_", 
                               remove = FALSE) %>%
      reshape2::melt(id.vars = names(trellData$comp_stats),
                     value.name = group_df_name)
    plotter <- suppressWarnings(dplyr::left_join(plotter, trellData$summary_stats))
    
    if(!is.null(trellData$data_values)){
      plotter <- suppressWarnings(dplyr::left_join(trellData$data_values, plotter))
    }
    
    tictoc::toc()
    
  } else {
    plotter <- trellData$data_values
  }
  
  tictoc::tic("nesting")
  
  nestplotter <- plotter %>% tidyr::nest(-panel_variable)
  
  tictoc::toc()
  
  # #Subset large groups ########### Take out later ######################################
 if (nrow(nestplotter) > 10){
   nestplotter <- nestplotter[1:10,]
 }
  
  tictoc::tic("plotting")
  nestplotter <- nestplotter %>% 
    dplyr::mutate(panel = trelliscopejs::map_plot(data, function(nestedData) {
      
      ## Generate specified y-limits for plotting (if applicable) ##
      # For comp stats #
      if(!is.null(trellData$comp_stats) && 
         comps_plot == TRUE){
        # Generate increment and y limits based on parameters and y-values #
        if (!is.null(comps_y_limits)){
          if (comps_y_limits == 'free'){
            comps_increment <-set_increment(nestedData[[comps_panel_y_axis]],
                                            include_zero = comps_include_zero)
            comps_ylims <- set_ylimits(nestedData[[comps_panel_y_axis]], 
                                 increment = comps_increment,
                                 y_min = comps_y_min, y_max = comps_y_max,
                                 include_zero = comps_include_zero)
          }
        } else if (!is.null(comps_y_range)){
          comps_increment <- set_increment(nestedData[[comps_panel_y_axis]], 
                                           include_zero = comps_include_zero)
          comps_ylims <- set_ylimits(nestedData[[comps_panel_y_axis]],
                               increment = comps_increment, y_min = comps_y_min, 
                               y_max = comps_y_max, y_range = comps_y_range,
                               include_zero = comps_include_zero)
        }
      
        ## Add border color, hover text, and label text to dataframe for plotting ##
        nestedData_comps <- add_plot_features(trellData,
                                              nestedData,
                                              p_val = p_val,
                                              panel_variable = panel_variable,
                                              comps_panel_y_axis = comps_panel_y_axis,
                                              colors = colors)
        
      # Make ggplots #
      plot_comps_all <- purrr::map(1:length(comps_plot_type), function (typenum){
        type <- comps_plot_type[typenum]
        
        plot_comps <- ggplot2::ggplot(
          data = nestedData_comps,
          ggplot2::aes(x = as.character(nestedData_comps[[comps_panel_x_axis]]),
                       y = nestedData_comps[[comps_panel_y_axis]],
                       color = bord,
                       fill = nestedData_comps[[comps_color_variable]],
                       text = gsub(", ", "\n", text),
                       label = gsub(", ", "\n", labels)
                       )
          ) +
          ggplot2::geom_hline(yintercept = 0) +
          ggplot2::xlab(comps_panel_x_axis) + 
          ggplot2::ylab(comps_panel_y_axis) +
          ggplot2::labs(fill = "", color = "") +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, 
                                                             hjust = 1, 
                                                             vjust = 0.5), 
                         legend.position='none') +
          ggplot2::guides(color = FALSE) +
          ggplot2::coord_cartesian(ylim=comps_ylims)

        if (comps_text == TRUE) {
          # Set text spacing above/below barplot #
          if (any(is.na(nestedData[[comps_panel_y_axis]]))){
            tempfold <- nestedData[[comps_panel_y_axis]]
            tempfold[is.na(tempfold)] <- 0
            if (option == "gtest"){
              textadj <- max(tempfold) + sign(max(tempfold))*comps_increment*4
            } else {
              textadj <- tempfold + sign(tempfold)*comps_increment*4
            }
          } else {
            if (option == "gtest"){
              textadj <- max(nestedData[[comps_panel_y_axis]]) + comps_increment*4
            } else {
              textadj <- nestedData[[comps_panel_y_axis]] + 
                sign(nestedData[[comps_panel_y_axis]])*comps_increment*4
            }
          }
          
          plot_comps <- plot_comps +
            ggplot2::geom_text(ggplot2::aes(y = textadj, 
                                   group = textadj), 
                               color = "black",
                               position = ggplot2::position_dodge(width = 0.9))
          
        }
        
        typelist <- c()
        if (stringr::str_detect(type, pattern = "col|bar")){
          typelist <- c("bar")
          plot_comps <- plot_comps + ggplot2::geom_col(position = "dodge", 
                                                      size = 1
                                                      ) +
               ggplot2::scale_color_manual(values = levels(nestedData_comps$bord))
        }
        if (stringr::str_detect(type, pattern = "point|scatter")){
          typelist <- c(typelist, "scatter")
          plot_comps <- plot_comps + ggplot2::geom_point(position = "identity", 
                                                         size = 2) +
              ggplot2::scale_color_manual(values = levels(nestedData_comps$bord))
        }
        if (stringr::str_detect(type, pattern = "box")){
          if(!is.na(na.omit(nestedData_comps[[comps_panel_y_axis]])) && 
             !is.na(var(na.omit(nestedData_comps[[comps_panel_y_axis]]))) &&
             var(na.omit(nestedData_comps[[comps_panel_y_axis]])) != 0 ){
            typelist <- c(typelist, "box")
            plot_comps <- plot_comps +  
              suppressWarnings(ggplot2::geom_boxplot(alpha = 0.2, 
                                    position = "dodge2"))
          } else {
            message("Varience of y_value is zero, boxplot is not applicable. Plotting y-value as a line (geom_crossbar).")
            typelist <- c(typelist, "line")
            plot_comps <- plot_comps +  
              ggplot2::geom_crossbar(position = "dodge", ggplot2::aes(
                ymin = nestedData_comps[[comps_panel_y_axis]],
                ymax = nestedData_comps[[comps_panel_y_axis]]), color = "black")
          }
        }
        
        if (plotly == TRUE){
          # Make and return plotly for map_plot #
          comps_plotly <- plotly::ggplotly(plot_comps, tooltip = c("text"))
          for (plotter in 1:length(comps_plotly$x$data)){
            if (!(comps_plotly$x$data[[plotter]]$type %in% typelist)){
              # comps_plotly$x$data[[plotter]]$showlegend <- FALSE
              comps_plotly$x$data[[plotter]]$hovertext <- stringr::str_remove(
                comps_plotly$x$data[[plotter]]$hovertext,
                "Group: .+\nCount: .+\n")
            }

          }
          
          return(comps_plotly)
        } else {
          
          return(plot_comps)
        }
      })
      
      if (plotly == TRUE){
        arr_comps_plot <- plotly::subplot(plot_comps_all, shareY = TRUE,
                                    margin = 0.02)
      } else {
        arr_comps_plot <- ggpubr::ggarrange(plotlist = plot_comps_all)
        }
      }
      
      if(!is.null(trellData$data_values) & 
         value_plot == TRUE){
        # Generate increment and y limits based on parameters and y-values #
        if (!is.null(value_y_limits)){
          if (value_y_limits == 'free'){
            increment <-set_increment(nestedData[[value_panel_y_axis]], 
                                      include_zero = value_include_zero)
            ylims <- set_ylimits(nestedData[[value_panel_y_axis]], 
                                 increment = increment,
                                 y_min = value_y_min, y_max = value_y_max, 
                                 include_zero = value_include_zero)
          }
        } else if (!is.null(value_y_range)){
          increment <- set_increment(nestedData[[value_panel_y_axis]], 
                                     include_zero = value_include_zero)
          ylims <- set_ylimits(nestedData[[value_panel_y_axis]],
                               increment = increment, y_min = value_y_min, 
                               y_max = value_y_max, y_range = value_y_range, 
                               include_zero = value_include_zero)
        }
        
        ## Set hover, excluding the panel_variable ##
        nestedData_value <- add_plot_features(trellData,
                                              nestedData,
                                              p_val = p_val,
                                              panel_variable = panel_variable,
                                              value_panel_y_axis = value_panel_y_axis)
        
        # Make ggplots #
        plot_value_all <- purrr::map(1:length(value_plot_type), function (typenum){
          type <- value_plot_type[typenum]
          
          plot_value <- ggplot2::ggplot(data = nestedData_value,
                                    ggplot2::aes(x = as.character(nestedData_value[[value_panel_x_axis]]), 
                                                 y = nestedData_value[[value_panel_y_axis]],
                                                 fill = nestedData_value[[value_color_variable]],
                                                 text = gsub(", ", "\n", text),
                                    )
          ) +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, 
                                                               hjust = 1, 
                                                               vjust = 0.5),
                           legend.position='none') +
            ggplot2::xlab(value_panel_x_axis) + 
            ggplot2::ylab(value_panel_y_axis) +
            ggplot2::labs(color = "", fill = "") +
            ggplot2::coord_cartesian(ylim=ylims)
          
          typelist <- c()
          
          if (stringr::str_detect(type, pattern = "col|bar")){
            typelist <- c("bar")
            plot_value <- plot_value + ggplot2::geom_col(position = "dodge", size = 1)
          }
          if (stringr::str_detect(type, pattern = "point|scatter")){
            typelist <- c(typelist, "scatter")
            if(plotly == TRUE){
              plot_value <- plot_value + ggplot2::geom_point(position = "identity", 
                                                             size = 2, 
                                                             color = "black")
            } else {
              plot_value <- plot_value + 
                ggplot2::geom_point(position = "identity",
                                    size = 2, 
                                    ggplot2::aes(colour = nestedData_value[[value_color_variable]]))
            }
          }
          if (stringr::str_detect(type, pattern = "box")){
            if(!is.na(na.omit(nestedData_value[[value_panel_y_axis]])) &&
               !is.na(var(na.omit(nestedData_value[[value_panel_y_axis]]))) &&
               var(na.omit(nestedData_value[[value_panel_y_axis]])) != 0 ){
              typelist <- c(typelist, "box")
              omitna <- nestedData_value[!is.na(nestedData_value[[value_panel_y_axis]]),]
              plot_value <- plot_value + suppressWarnings(
                ggplot2::geom_boxplot(data = omitna, alpha = 0.2,
                                      position = "dodge2", 
                                      ggplot2::aes(
                                        y = omitna[[value_panel_y_axis]],
                                        x = omitna[[value_panel_x_axis]],
                                        fill = omitna[[value_color_variable]],
                                        text = NA)
                                      ))
            } else {
              typelist <- c(typelist, "line")
              plot_value <- plot_value +  
                ggplot2::geom_crossbar(position = "dodge2", ggplot2::aes(
                  ymin = nestedData_value[[value_panel_y_axis]],
                  ymax = nestedData_value[[value_panel_y_axis]],
                  color = nestedData_value[[value_color_variable]]))
            }
          }
          
          if (plotly == TRUE){
            value_plotly <- plotly::ggplotly(plot_value, tooltip = c("text"))
            return(value_plotly)
          } else {
            return(plot_value)
            }
        })
        
        if (plotly == TRUE){
          arr_value_plot <- plotly::subplot(plot_value_all, shareY = TRUE,
                                    margin = 0.02)
        } else {
          arr_value_plot <- ggpubr::ggarrange(plotlist = plot_value_all)
        }
      }
      
      if(plotly == TRUE){
        if(!is.null(trellData$data_values) & !is.null(trellData$comp_stats)){
          all_plotly <- plotly::subplot(list(arr_value_plot, arr_comps_plot), nrows = 2, 
                            titleX = TRUE, titleY = TRUE,
                            margin = 0.1)
        } else if (!is.null(trellData$data_values)) {
          all_plotly <- plotly::subplot(list(arr_value_plot), titleX = TRUE, titleY = TRUE, 
                            margin = 0.1)
        } else {
          all_plotly <- plotly::subplot(list(arr_comps_plot), titleX = TRUE, titleY = TRUE, 
                            margin = 0.1)
        }
        return(all_plotly)
      } else {
        if(!is.null(trellData$data_values) & !is.null(trellData$comp_stats)){

          all_ggplot <- ggpubr::ggarrange(plotlist = list(arr_value_plot, arr_comps_plot), nrow = 2)
        } else if (!is.null(trellData$data_values)) {
          all_ggplot <- arr_value_plot
        } else {
          all_ggplot <- arr_comps_plot
        }
        return(all_ggplot)
      }

    }
  )
)
  tictoc::toc()
  # Return nested table #
  
  #### Don't actually need data
  nestplotter <- nestplotter[,c(1,3)]
  
  attr(nestplotter, "parent_class") <- attr(trellData, "parent_class")
  return(nestplotter)
}



#' @name validate_format_plot_input
#' @rdname validate_format_plot_input
#' @title Validate inputs for omicsData and omicsStats in trelliVis processing
#' 
#' @description Checks for validity of trellData input and assigns variables where needed
#'
#' @param trellData An object of class "trellData" generated from \code{\link{format_data}}.
#' @param comps_y_limits For comparisons: Set to "fixed" or "free" for automated y-axis calculating. "fixed" - axis generated based on the maximum/minimum across all plots. "free" - axis axis generated based on the maximum/minimum of individual plot.
#' @param comps_y_range For comparisons: Specify a range for the plot y-axis. Will calculate the range based on one of y_max or y_min parameters or from the median of y-values where y_max and y_min are not defined.
#' @param comps_y_max For comparisons: Sets the maximum y-value for the y-axis.
#' @param comps_y_min For comparisons: Sets the minimum y-value for the y-axis.
#' @param value_y_limits For values: Set to "fixed" or "free" for automated y-axis calculating. "fixed" - axis generated based on the maximum/minimum across all plots. "free" - axis axis generated based on the maximum/minimum of individual plot.
#' @param value_y_range For values: Specify a range for the plot y-axis. Will calculate the range based on one of y_max or y_min parameters or from the median of y-values where y_max and y_min are not defined.
#' @param value_y_max For values: Sets the maximum y-value for the y-axis.
#' @param value_y_min For values: Sets the minimum y-value for the y-axis.
#' @param p_val Specifies p-value for setting graph border colors
#' @param panel_variable Specifies what to divide trelliscope panels by, must be a column in trellData. Defaults to cnames$edata_cname of trellData.
#' @param comps_panel_x_axis Specifies what column should be on the x-axis, must be a column in trellData. Default setting plots pairwise comparisons along x-axis.
#' @param comps_panel_y_axis Specifies what column should be on the y-axis, must be a column in trellData and numeric. Default setting plots fold change for combined and anova testing and counts for g-test.
#' @param comps_color_variable Specifies what column should distingush color, must be a column in trellData. Default settings is set to "Group."
#' @param value_panel_x_axis Specifies what column should be on the x-axis, must be a column in trellData. Default setting plots pairwise comparisons along x-axis.
#' @param value_panel_y_axis Specifies what column should be on the y-axis, must be a column in trellData and numeric. Default setting plots fold change for combined and anova testing and counts for g-test.
#' @param value_color_variable Specifies what column should distingush color, must be a column in trellData. Default settings is set to "Group."
#' @param value_plot_type For values: Specifies plot types for graphing; must be a list of strings where "box", "bar", or "point" is specified. Combined strings like "boxpoint" will plot both in the same graph.
#' @param comps_plot_type For comparisons: Specifies plot types for graphing; must be a list of strings where "box", "bar", or "point" is specified. Combined strings like "boxpoint" will plot both in the same graph.
#' @param comps_include_zero For comparisons: Should plots show y = 0?
#' @param value_include_zero For values: Should plots show y = 0?
#' @param value_plot For values: In the presence of data_values in trellData, should the plot be rendered? Default is TRUE.
#' @param comps_plot For comparisons: In the presence of summary_stats and comp_stats in trellData, should the plot be rendered? Default is TRUE.
#' @param comps_text For comparisons only: TRUE/FALSE for p-value text above data
#' @param plotly Should the plots be rendered as plotly objects?
#' @param ... further arguments
#'
#' @author Rachel Richardson
#'

validate_format_plot_input <- function(...){
  .validate_format_plot_input(...)
}

.validate_format_plot_input <- function(trellData, 
                                  comps_y_limits = NULL, comps_y_range = NULL, 
                                  comps_y_max = NULL, comps_y_min = NULL, 
                                  comps_include_zero = NULL,
                                  value_y_limits = NULL, value_y_range = NULL, 
                                  value_y_max = NULL, value_y_min = NULL,
                                  value_include_zero = NULL,
                                  p_val = NULL,
                                  panel_variable = NULL, 
                                  comps_color_variable = NULL,
                                  comps_panel_x_axis = NULL, comps_panel_y_axis = NULL,
                                  value_color_variable = NULL,
                                  value_panel_x_axis = NULL, 
                                  value_panel_y_axis = NULL,
                                  value_plot_type = NULL, comps_plot_type = NULL,
                                  value_plot = NULL, comps_plot = NULL, 
                                  comps_text = NULL, plotly = NULL) {
  
  ##### Initial checks #####
  
  # Check if class is correct #
  if(!inherits(trellData, "trellData")) stop(
    "trellData must be of the class 'trellData'")
  
  if(!(is.null(p_val) && 
       is.null(plotly) &&
       is.null(comps_y_limits) &&
       is.null(comps_y_range) &&
       is.null(comps_y_max) &&
       is.null(comps_y_min) &&
       is.null(comps_include_zero) &&
       is.null(comps_plot) &&
       is.null(comps_text) &&
       is.null(comps_plot_type) &&
       is.null(value_y_limits) &&
       is.null(value_y_range) &&
       is.null(value_y_max) &&
       is.null(value_y_min) &&
       is.null(value_plot_type) &&
       is.null(value_include_zero) &&
       is.null(value_plot)
       )) {

    
    # Check if comp_stats is in trellData #
    if(is.null(trellData$comp_stats) && is.null(trellData$data_values)) stop(
      "No data values or comparison statistics in trellData to plot")
    
    # Check if comp_stats is in trellData (as above) #
    if(is.na(trellData$comp_stats) && is.na(trellData$data_values)) stop(
      "No data values or comparison statistics in trellData to plot")
    
    # Check check if p_val is numeric of length 1 #
    if(!is.numeric(p_val) || (length(p_val) != 1)) stop(
      "p_val must be a numeric of length 1")  
    
    # Check check if plotly is logical of length 1 #
    if(!is.logical(plotly) || (length(plotly) != 1)) stop(
      "plotly must be a logical (TRUE or FALSE) of length 1") 
    
    # Check check if comps_include_zero is logical of length 1 #
    if(!is.logical(comps_include_zero) || (length(comps_include_zero) != 1)) stop(
      "comps_include_zero must be a logical (TRUE or FALSE) of length 1") 
    
    # Check check if comps_plot is logical of length 1 #
    if(!is.logical(comps_plot) || (length(comps_plot) != 1)) stop(
      "comps_plot must be a logical (TRUE or FALSE) of length 1") 
    
    # Check check if comps_text is logical of length 1 #
    if(!is.logical(comps_text) || (length(comps_text) != 1)) stop(
      "comps_text must be a logical (TRUE or FALSE) of length 1") 
    
    # Check check if value_include_zero is logical of length 1 #
    if(!is.logical(value_include_zero) || (length(value_include_zero) != 1)) stop(
      "value_include_zero must be a logical (TRUE or FALSE) of length 1") 
    
    # Check check if value_plot is logical of length 1 #
    if(!is.logical(value_plot) || (length(value_plot) != 1)) stop(
      "value_plot must be a logical (TRUE or FALSE) of length 1") 
    
    
    # Check if y_limits or y_range have been selected correctly #
    if((!is.null(comps_y_limits) & !is.null(comps_y_range)) | 
       (!is.null(value_y_limits) & !is.null(value_y_range))) stop(
         "Input either y_limits or y_range parameters, but not both.")
    
    # --Comps-- #
    # Check if only one of comps_y_max and comps_y_min has been selected with comps_y_limits or comps_y_range #
    if(!is.null(comps_y_max) & 
       !is.null(comps_y_min) & 
       (!is.null(comps_y_range) | !is.null(comps_y_limits))) stop(
         "Cannot use both comps_y_min and comps_y_max with comps_y_range 
       or comps_y_limits parameters. Only one of comps_y_min or 
       comps_y_max can be used.")
    
    # Check if comps_y_limits is in acceptable strings and length == 1 #
    if ( !is.null(comps_y_limits)){
      if((length(comps_y_limits) != 1)) stop(
        "Parameter y_limits must have length = 1.")
      if(!(comps_y_limits %in% c("fixed", "free"))) stop(
        "Parameter y_limits must be input as either 'fixed' or 'free'.")
    }
    
    # Check if comps_y_range is positive, numeric and length == 1 #
    if (!is.null(comps_y_range)){
      if(!is.numeric(comps_y_range)) stop(
        "Parameter y_range must be numeric.")
      if(length(comps_y_range) != 1) stop(
        "Parameter y_range must have length = 1.")
      if(!(comps_y_range > 0)) stop(
        "Parameter y_range must be greater than zero.")
    }
    
    # Check if comps_y_max is numeric and length == 1 #
    if (!is.null(comps_y_max)){
      if(!is.numeric(comps_y_max)) stop(
        "Parameter y_max must be numeric.")
      if(length(comps_y_max) != 1) stop(
        "Parameter y_max must have length = 1.")
    }
    
    # Check if comps_y_min is numeric and length == 1 #
    if (!is.null(comps_y_min)){
      if(!is.numeric(comps_y_min)) stop(
        "Parameter y_min must be numeric.")
      if(length(comps_y_min) != 1) stop(
        "Parameter y_min must have length = 1.")
    }
    
    # Check if comps_plot_type has one of the available options #
    checkplot <- purrr::map(comps_plot_type, 
                            function(plot) stringr::str_detect(plot, "box|col|point|scatter|bar"))
    if (any(!unlist(checkplot))) stop(
      "Invalid entry in comps_plot_type. Plot_type strings must contain 
    at least one of the following: box, col, point, scatter, bar")
    
    # --Value-- #
    # Check if only one of value_y_max and value_y_min has been selected with value_y_limits or value_y_range #
    if(!is.null(value_y_max) & 
       !is.null(value_y_min) & 
       (!is.null(value_y_range) | !is.null(value_y_limits))) stop(
         "Cannot use both value_y_min and value_y_max with value_y_range 
       or value_y_limits parameters. Only one of y_min or y_max can be 
       used.")
    
    # Check if value_y_limits is in acceptable strings and length == 1 #
    if (!is.null(value_y_limits) ){
      if((length(value_y_limits) != 1)) stop(
        "Parameter y_limits must have length = 1.")
      if(!(value_y_limits %in% c("fixed", "free")) ) stop(
        "Parameter y_limits must be input as either 'fixed' or 'free'.")
    }
    
    # Check if value_y_range is positive, numeric and length == 1 #
    if (!is.null(value_y_range)){
      if(!is.numeric(value_y_range)) stop(
        "Parameter y_range must be numeric.")
      if(length(value_y_range) != 1) stop(
        "Parameter y_range must have length = 1.")
      if(!(value_y_range > 0)) stop(
        "Parameter y_range must be greater than zero.")
    }
    
    # Check if value_y_max is numeric and length == 1 #
    if (!is.null(value_y_max)){
      if(!is.numeric(value_y_max)) stop(
        "Parameter y_max must be numeric.")
      if(length(value_y_max) != 1) stop(
        "Parameter y_max must have length = 1.")
    }
    
    # Check if value_y_min is numeric and length == 1 #
    if (!is.null(value_y_min)){
      if(!is.numeric(value_y_min)) stop(
        "Parameter y_min must be numeric.")
      if(length(value_y_min) != 1) stop(
        "Parameter y_min must have length = 1.")
    }
    
    # Check if value_plot_type has one of the available options #
    checkplot <- purrr::map(value_plot_type, 
                            function(plot) stringr::str_detect(plot, "box|col|point|scatter|bar"))
    if (any(!unlist(checkplot))) stop(
      "Invalid entry in value_plot_type. Plot_type strings must contain 
    at least one of the following: box, col, point, scatter, bar")
    
  } else {
    
    ##### Post-Moniker and null assignment check ####
    # --Value-- #
    if (!is.null(trellData$data_values)) {
      
      # Ensure panel, color/x/y parameters are not matching #
      if (!((value_panel_x_axis != panel_variable) & 
            (value_panel_y_axis != panel_variable) & 
            (panel_variable != value_color_variable))){
        stop("Parameters for value panel_y_axis, panel_x_axis, and color_variable 
           cannot match panel_variable. Refer to ?plot_comps for default settings 
           or try setting all of these parameters individually.")
      }
      if(value_panel_x_axis == value_panel_y_axis) warning(
        "Parameter for value panel_y_axis and panel_x_axis are identical. 
            Refer to ?plot_comps for default settings or try setting 
            parameters individually for different axis labels.")
      
      # Variable for variable-in checks #
      allcol <- colnames(trellData$data_values)
      allcolstr <- toString(allcol)
      
      # Ensure panel, color, x, and y parameters are in data_values #
      if (!(value_panel_x_axis %in% allcol)) stop(
        paste("Parameter value_panel_x_axis  must be in:", allcolstr))
      if (!(value_panel_y_axis %in% allcol)) stop(
        paste("Parameter value_panel_y_axis must  must be in:", allcolstr))
      if (!(panel_variable %in% allcol)) stop(
        paste("Parameter panel_variable must be in:", allcolstr))
      if (!(value_color_variable %in% allcol)) stop(
        paste("Parameter value_color_variable must be in:", allcolstr))
      
    }
    
    # --Comps-- #
    
    # Ensure summary_stats and comp_stats are both present #
    if(((is.null(trellData$summary_stats)) & 
        (!is.null(trellData$comp_stats))) | 
       ((!is.null(trellData$summary_stats)) &
        (is.null(trellData$comp_stats))) ) {
      stop("Both summary_stats and comp_stats must be present 
         if one or the other is in trellData.")
    }
    
    if ((!is.null(trellData$summary_stats)) & 
        (!is.null(trellData$comp_stats))) {
      
      # Check if stats statistical test attribute is valid #
      if(!(
        attributes(trellData)$statistical_test %in% 
        c("combined", "gtest", "anova"))) stop(
        paste("Non-applicable statistical_test attribute in trellData object."))
      
      # Ensure panel, color/x/y parameters are not matching #
      if(comps_panel_x_axis == comps_panel_y_axis) warning(
        "Parameter for comps panel_y_axis and panel_x_axis are identical.
    Refer to ?plot_comps for default settings or try setting parameters
    individually for different axis labels.")
      
      if (!((comps_panel_x_axis != panel_variable) && 
            (comps_panel_y_axis != panel_variable) && 
            (panel_variable != comps_color_variable))) stop(
              "Parameters for comps panel_y_axis, panel_x_axis, 
            and color_variable cannot match panel_variable. Refer to ?plot_comps 
            for default settings or try setting all of these parameters 
            individually.")
      
      # Variable for variable-in checks #
      allcol <- c(colnames(trellData$comp_stats), 
                  colnames(trellData$summary_stats))
      allcolstr <- toString(unique(allcol))
      
      # Ensure panel, color, x, and y parameters are in comp_stats or summary stats #
      if (!(comps_panel_x_axis %in% allcol)) stop(
        paste("Parameter comps_panel_x_axis  must be in:", allcolstr))
      if (!(comps_panel_y_axis %in% allcol)) stop(
        paste("Parameter comps_panel_y_axis must  must be in:", allcolstr))
      if (!(panel_variable %in% allcol)) stop(
        paste("Parameter panel_variable must be in:", allcolstr))
      if (!(comps_color_variable %in% allcol)) stop(
        paste("Parameter comps_color_variable must be in:", allcolstr))
      
    }
  }
}
  

#' @name generate_plot_message
#' @rdname generate_plot_message
#' @title Generates a message stating the specified plotting y-limits
#' 
#' @description Generates a message stating the specified plotting y-limits
#'
#' @param trellData An object of class "trellData" generated from \code{\link{format_data}}.
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
}


#' @name add_plot_features
#' @rdname add_plot_features
#' @title Adds features to dataframe for plotting
#' 
#' @description Adds border color, hover, and label text as applicable to nested trellData dataframe for plotting
#'
#' @param trellData An object of class "trellData" generated from \code{\link{format_data}}.
#' @param nestedData Nested trellData dataframe
#' @param p_val Specifies p-value for setting graph border colors
#' @param panel_variable Specifies what to divide trelliscope panels by, must be a column in trellData. 
#' @param comps_panel_y_axis Specifies what column should be on the y-axis, must be a column in trellData and numeric. 
#' @param value_panel_y_axis Specifies what column should be on the y-axis, must be a column in trellData and numeric. 
#' @param colors Specifies colors to use for border.
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
                               value_panel_y_axis = NULL,
                               colors = NULL){
  
  if(is.null(value_panel_y_axis)){
    
    # Set border colors based on significance #
    bord <- rep(colors[1], length(nestedData$P_value_T))
    bord[nestedData$P_value_G < p_val & !is.na(nestedData$P_value_G)] <-  colors[2]
    bord[nestedData$P_value_T < p_val & !is.na(nestedData$P_value_T)] <- colors[3]
    
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
  maxi <- max(yvalues) + 5*increment
  mini <- min(yvalues) - 5*increment
  
  # Adjust for maximums below 0 and minimums above 0 #
  if (!(maxi > 0) && (include_zero == TRUE)){
    maxi <- 5*increment
  }
  if (!(mini < 0) && (include_zero == TRUE)){
    mini <- -5*increment
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
#' @param trellData An object of class "trellData" generated from \code{\link{format_data()}}.
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
    trellData = format_data(omicsData, omicsStats)
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
    pvalg <- data.frame(t(joiner[stringr::str_detect(names(joiner), "P_value_G")]))
    pvalt <- data.frame(t(joiner[stringr::str_detect(names(joiner), "P_value_T")]))
    colnames(pvalg) <- substring(colnames(pvalg), 2)
    colnames(pvalt) <- substring(colnames(pvalt), 2)
    Sig_G_p_value <- as.factor(uniqlist %in% uniqlist[as.numeric(names(
      select_if(pvalg, function(x) any(x < p_val & !is.na(x)))
    ))])
    Sig_T_p_value <- as.factor(uniqlist %in% uniqlist[as.numeric(names(
      select_if(pvalt, function(x) any(x < p_val & !is.na(x)))
    ))])
    
    addcogs <- data.frame(
      uniqlist,
      trelliscopejs::cog(Sig_G_p_value, desc = "Boolean for significant G test p-values,
      where TRUE indicates that the p-value is < 0.05 and is grounds for 
      rejecting the null hypothesis (independence of missing data)."),
      trelliscopejs::cog(Sig_T_p_value, desc = "Boolean for significant T test p-values,
      where TRUE indicates that the p-value is < 0.05 and is grounds for 
      rejecting the null hypothesis (no significant 
      difference between groups)."))
    colnames(addcogs) <- c(uniqueID, 'Sig_G_p_value', 'Sig_T_p_value')
    
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
#' @param omicsFormat Output of format_data() function
#' @param p_val Numeric that specifies p-value for significance calculations. Default is 0.05.
#' @param mapping_col String: For associated proData/pepData - name of column in peptide data with protein information. Default is NULL.
#' @param panel_variable String: Name of column that plot panels are sorted by (e.g. each plotting arrangement has a unique identifier from panel variable). Default is emeta_cname if present, edata_cname where emeta_cname is not present.
#' @param try_URL Will attempt to link to PubChem, LipidMaps, or Uniprot based on information in edata_cname of omicsData or specified mapping_col for peptide data. Default is FALSE.
#' @param trelli_name String: name of display, or list of names where a list is provided for omicsData and omicsStats
#' @param trelli_path_out String: path to where trelliscope is stored. Default is "./TrelliDisplay"
#' @param comps_y_limits For omicsStats: Set to "fixed" or "free" for automated y-axis calculating. "fixed" - axis generated based on the maximum/minimum across all plots. "free" - axis axis generated based on the maximum/minimum of individual plot.
#' @param comps_y_range For omicsStats: Specify a range for the plot y-axis. Will calculate the range based on one of y_max or y_min parameters or from the median of y-values where y_max and y_min are not defined.
#' @param comps_y_max For omicsStats: Sets the maximum y-value for the y-axis.
#' @param comps_y_min For omicsStats: Sets the minimum y-value for the y-axis.
#' @param comps_panel_x_axis For omicsStats: Specifies what column should be on the x-axis. Default setting plots pairwise comparisons along x-axis.
#' @param comps_panel_y_axis For omicsStats: Specifies what column should be on the y-axis and numeric. Default setting plots fold change for combined and anova testing and counts for g-test.
#' @param comps_color_variable For omicsStats: Specifies what column should distingush color. Default settings is set to "Group."
#' @param comps_plot_type For omicsStats: Specifies plot types for graphing; must be a list of strings where "box", "bar", or "point" is specified. Combined strings like "boxpoint" will plot both in the same graph.
#' @param comps_include_zero For omicsStats: Should plots show y = 0?
#' @param comps_text For omicsStats only: Should p-value text be displayed with plot points/bars/boxes?
#' @param comps_plot For omicsStats: In the presence of summary_stats and comp_stats in trellData, should the plot be rendered? Default is TRUE.
#' @param value_y_limits For omicsData: Set to "fixed" or "free" for automated y-axis calculating. "fixed" - axis generated based on the maximum/minimum across all plots. "free" - axis axis generated based on the maximum/minimum of individual plot.
#' @param value_y_range For omicsData: Specify a range for the plot y-axis. Will calculate the range based on one of y_max or y_min parameters or from the median of y-values where y_max and y_min are not defined.
#' @param value_y_max For omicsData: Sets the maximum y-value for the y-axis.
#' @param value_y_min For omicsData: Sets the minimum y-value for the y-axis.
#' @param value_panel_x_axis For omicsData: Specifies what column should be on the x-axis. Default setting plots pairwise comparisons along x-axis.
#' @param value_panel_y_axis For omicsData: Specifies what column should be on the y-axis and numeric. Default setting plots fold change for combined and anova testing and counts for g-test.
#' @param value_color_variable For omicsData: Specifies what column should distingush color. Default settings is set to "Group."
#' @param value_plot_type For omicsData: Specifies plot types for graphing; must be a list of strings where "box", "bar", or "point" is specified. Combined strings like "boxpoint" will plot both in the same graph.
#' @param value_include_zero For omicsData: Should plots show y = 0?
#' @param value_plot For omicsData: In the presence of data_values in trellData, should the plot be rendered? Default is TRUE.
#' @param plotly Should the plots be rendered as plotly objects?
#'
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
                       comps_y_limits = NULL, comps_y_range = NULL, 
                       comps_y_max = NULL, comps_y_min = NULL, 
                       comps_include_zero = TRUE, value_y_limits = NULL, 
                       value_y_range = NULL, value_y_max = NULL, 
                       value_y_min = NULL, value_include_zero = FALSE, 
                       comps_color_variable = NULL,
                       comps_panel_x_axis = "Comparison", 
                       comps_panel_y_axis = NULL, value_color_variable = NULL, 
                       value_panel_x_axis = NULL, value_panel_y_axis = NULL, 
                       value_plot_type = "boxpoint", comps_plot_type = "col", 
                       value_plot = TRUE, comps_plot = TRUE, comps_text = TRUE,
                       plotly = TRUE) {
  
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
    trellData <- format_data(omicsData, omicsStats)
  } else {
    trellData <- omicsFormat
  }
  
  # If a pep/pro pair is listed, act on each item #
  if(class(trellData) == "list"){
    
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
    
    # #Adds e_meta of associated peptide data to associated protein data
    # if(class(omicsData) == "list" &&
    #    length(omicsData) != 1){
    #   vectorpro <- c(inherits(omicsData[1][[1]], "proData"), 
    #                  inherits(omicsData[2][[1]], "proData"))
    #   omicsData[vectorpro][[1]]$e_meta <- suppressWarnings(dplyr::left_join(omicsData[vectorpro][[1]]$e_meta,
    #                                                 omicsData[!vectorpro][[1]]$e_meta))
    # }
    
    # Nest data and generate trelliscope plots #
    nested_plot <- purrr::map2(trellData, panel_variable, function(pairedplotter, pan){
      format_plot(pairedplotter, 
                  comps_y_limits = comps_y_limits, 
                  comps_y_range = comps_y_range, 
                  comps_y_max = comps_y_max, 
                  comps_y_min = comps_y_min, 
                  comps_include_zero = comps_include_zero,
                  value_y_limits = value_y_limits, 
                  value_y_range = value_y_range, 
                  value_y_max = value_y_max, 
                  value_y_min = value_y_min,
                  value_include_zero = value_include_zero,
                  p_val = p_val,
                  panel_variable = pan, 
                  comps_color_variable = comps_color_variable,
                  comps_panel_x_axis = comps_panel_x_axis, 
                  comps_panel_y_axis = comps_panel_y_axis,
                  value_color_variable = value_color_variable,
                  value_panel_x_axis = value_panel_x_axis, 
                  value_panel_y_axis = value_panel_y_axis,
                  value_plot_type = value_plot_type, 
                  comps_plot_type = comps_plot_type,
                  value_plot = value_plot,
                  comps_plot = comps_plot,
                  comps_text = comps_text,
                  plotly = plotly)
      })
    
    # Generate default cognostics #
    nest_plot_cog_list <- purrr::pmap(
      list(nested_plot,
           omicsData,
           omicsStats, 
           panel_variable,
           mapping_col),
      function(nest, dat, stat, pan, map){
        data_cogs(nested_plot = nest, 
                  omicsData = dat,  
                  omicsStats = stat,
                  p_val = p_val,
                  mapping_col = map, 
                  panel_variable = pan, 
                  try_URL = try_URL)
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
    nested_plot <- format_plot(trellData,
                               comps_y_limits = comps_y_limits, 
                               comps_y_range = comps_y_range, 
                               comps_y_max = comps_y_max, 
                               comps_y_min = comps_y_min, 
                               comps_include_zero = comps_include_zero,
                               value_y_limits = value_y_limits, 
                               value_y_range = value_y_range, 
                               value_y_max = value_y_max, 
                               value_y_min = value_y_min,
                               value_include_zero = value_include_zero,
                               p_val = p_val,
                               panel_variable = panel_variable, 
                               comps_color_variable = comps_color_variable,
                               comps_panel_x_axis = comps_panel_x_axis, 
                               comps_panel_y_axis = comps_panel_y_axis,
                               value_color_variable = value_color_variable,
                               value_panel_x_axis = value_panel_x_axis, 
                               value_panel_y_axis = value_panel_y_axis,
                               value_plot_type = value_plot_type, 
                               comps_plot_type = comps_plot_type,
                               value_plot = value_plot,
                               comps_plot = comps_plot,
                               comps_text = comps_text,
                               plotly = plotly)
    
    # Generate default cognostics #
    nest_plot_cog <- data_cogs(nested_plot = nested_plot, 
                               omicsData = omicsData, 
                               omicsStats = omicsStats, 
                               p_val = p_val, 
                               mapping_col = mapping_col, 
                               panel_variable = panel_variable, 
                               try_URL = try_URL)
    
    # Generate trelliscope display #
    nest_plot_cog %>%
      trelliscopejs::trelliscope(name = as.character(trelli_name), nrow = 1, ncol = 2,
                  path = as.character(trelli_path_out), thumb = TRUE,
                  state = list(
                    sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")), 
                    labels = list(names(nest_plot_cog[1]))))
  }
}
