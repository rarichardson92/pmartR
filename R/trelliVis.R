#' @name as.trellData
#' @rdname as.trellData
#' @title Convert Omics data and pairwise statistics to a plotting object
#'
#' @description Converts a ResObject and its respective OmicsData into an easily plottable object. Called automatically in trelliVis() function.
#'
#' @param omicsData A pmartR object of class pepData, lipidData, metabData, or proData
#' @param omicsStats A statistical results object produced by running \code{imd_anova} on omicsData.
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
#' @seealso \link[pmartR]{trelliVis}
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' library(pmartR)
#' data("metab_object")
#' mymetabData <- pmartR::edata_transform(metab_object, "log10")
#' mymetabData <- pmartR::group_designation(metab_object, "Condition")
#' mymetabData <- pmartR::normalize_global(metab_object, "all", "median", apply_norm = TRUE, backtransform = TRUE)
#' as.trellData(omicsData = mymetabData)
#' 
#' mymetabStats <- pmartR::imd_anova(mymetabData, test_method = "combined")
#' as.trellData(omicsData = mymetabData, omicsStats = mymetabStats)
#' 
#' data("isobaric_object")
#' mypepData <- pmartR::edata_transform(isobaric_object, "log10")
#' mypepData <- pmartR::normalize_isobaric(mypepData, apply_norm = T) # For isobaric normalization
#' mypepData <- pmartR::group_designation(mypepData, "Group")
#' mypepData <- pmartR::normalize_global(mypepData, "all", "median", apply_norm = TRUE, backtransform = TRUE)
#' myproRoll <- pmartR::protein_quant(mypepData, "rrollup")
#' mypepStats <- pmartR::imd_anova(mypepData, test_method = "combined")
#' myproStats <- pmartR::imd_anova(myproRoll, test_method = "combined")
#' pmartR::as.trellData(omicsData = list(mypepData, myproRoll), omicsStats = list(mypepStats, myproStats))
#' 
#'}
#'
#' @author Rachel Richardson
#'
#' @export
#' 
as.trellData <- function(omicsData = NULL, omicsStats = NULL){
  .as.trellData(omicsData, omicsStats)
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
        dplyr::left_join(omicsData$f_data, 
                         by = sampID))
    
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
    data_values <- suppressWarnings(
      dplyr::left_join(data_values, 
                       joingroupDF,
                       by = dplyr::intersect(colnames(data_values), colnames(joingroupDF))))
    
    
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
#' For printing an S3 object of type 'trellData'; Displays truncated dataframes for each element in trellData.
#' 
#'@rdname print-trellData
#'@param trellData The trellData object to print.
#'
#'@examples
#'
#' dontrun{
#' library(pmartRdata)
#' library(pmartR)
#' data("metab_object")
#' mymetabData <- pmartR::edata_transform(metab_object, "log10")
#' mymetabData <- pmartR::group_designation(metab_object, "Condition")
#' mymetabData <- pmartR::normalize_global(metab_object, "all", "median", apply_norm = TRUE, backtransform = TRUE)
#' print(pmartR::as.trellData(mymetabData))
#' }
#' 
#' @seealso \link[pmartR]{as.trellData}
#' 
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
#' @param omicsData A pmartR object of class pepData, lipidData, metabData, or proData
#' @param omicsStats A statistical results object produced by running \code{imd_anova} on omicsData.
#'
#' @seealso \link[pmartR]{as.trellData}
#'
#' @author Rachel Richardson
#'

recursive_format <- function(omicsData = NULL, omicsStats = NULL){
  .recursive_format(omicsData, omicsStats)
}

.recursive_format <- function(omicsData = NULL, omicsStats = NULL){

  #--Both--#
  if (!is.null(omicsData) && !is.null(omicsStats)) {
    if(length(omicsData) == 1 && length(omicsStats) == 1) {
      plotterlist <- as.trellData(omicsData[[1]], omicsStats[[1]])
      return(plotterlist)
    } else {
      validate_omics_input(omicsData = omicsData, omicsStats = omicsStats)
      classlist <- omicsData %>% purrr::map(function(omics) attr(omics, which = "class"))
      plotterlist <- purrr::map2(omicsData, omicsStats, as.trellData)
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
#' @title Validate inputs for omicsData and omicsStats in as.trellData processing.
#' 
#' @description Checks for validity of omicsData and omicsStats inputs.
#'
#' @param omicsData A pmartR object of class pepData, lipidData, metabData, or proData
#' @param omicsStats A statistical results object produced by running \code{imd_anova} on omicsData.
#'
#' @seealso \link[pmartR]{as.trellData}
#'
#' @author Rachel Richardson
#'

validate_omics_input <- function(omicsData = NULL, omicsStats = NULL){
  .validate_omics_input(omicsData, omicsStats)
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
    if (!is.null(omicsData) && !is.null(omicsStats)) {
      
      # Check that the list inputs are equal in length (1 and 1, or 2 and 2) #
      if(length(omicsData) != length(omicsStats)) stop(
          "List length does not match; lists for omicsData amd omicsStats should contain
          the same protien set(s) and peptide set(s)."
      )
      if(length(omicsData) != 1 && length(omicsStats) != 1) {
        # # Check that the list inputs are equal in length (1 and 1, or 2 and 2) #
        # if (length(omicsData) != length(omicsStats)) stop(
        #   "List length does not match; lists for omicsData amd omicsStats should contain
        #   the same protien set(s) and peptide set(s).")
        
        # Check that omicsData input is a list() not a c() of omics data #
        if (any(unlist(purrr::map(c(omicsData), is.data.frame)))) stop(
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
    
    # if((class(omicsData) == "list" && class(omicsStats) != "NULL") || 
    #    (class(omicsData) != "NULL" && class(omicsStats) == "list")) stop(
    #      "Dual input using lists in omicsData and omicsStats must be lists of equal length."
    #    )

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
      if(is.null(attr(omicsStats, "statistical_test")) || 
         !(attr(omicsStats, "statistical_test") %in% c("combined", "anova", "gtest"))) stop(
        "OmicsStats statistical_test attribute incorrect or NULL; must be combined, anova, or gtest."
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
      
      # # Check required data # ##### redundant with Biomolecules in omicsStats do not match biomolecules in omicsData, column does not match between omicsData and omicsStats errors
      # if(is.null(omicsData$e_data) || is.null(omicsData$f_data)) stop( 
      #   "Omicsdata requires both e_data and f_data."
      # )
      
      # # Check cname e_data is in e_data, cname f_data is in f_data #    ##### redundant with Biomolecules in omicsStats do not match biomolecules in omicsData, column does not match between omicsData and omicsStats errors
      # if(!(uniqedatId %in% colnames(omicsData$e_data) && 
      #      sampID %in% colnames(omicsData$f_data))) stop(
      #        paste(
      #          paste(uniqedatId, "column must be present in omicsData e_data."),
      #          paste(sampID, "column must be present in omicsData f_data.")
      #        )
      #      )
      
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
      if(is.null(attr(omicsStats, "statistical_test"))|| 
         !(attr(omicsStats, "statistical_test") %in% c("combined", "anova", "gtest"))) stop(
        "OmicsStats statistical_test attribute incorrect or NULL; must be combined, anova, or gtest."
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


#' @name custom_fn
#' @rdname custom_fn
#' @title Test function for format_plot(), custom_plot argument
#'
#' @description Test function for format_plot(), custom_plot argument. Plots a ggplot where x = 1:10, y = 1:10.
#'
#' @param data data subset from a trellData object.
#' 
#' @seealso \link[pmartR]{format_plot}
#' 
#' @examples
#'
#' dontrun{
#' library(pmartRdata)
#' library(pmartR)
#' library(ggplot2)
#' data("metab_object")
#' mymetabData <- pmartR::edata_transform(metab_object, "log10")
#' mymetabData <- pmartR::group_designation(metab_object, "Condition")
#' mymetabData <- pmartR::normalize_global(metab_object, "all", "median", apply_norm = TRUE, backtransform = TRUE)
#' 
#' custom_fn <- function(data){
#' ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = 1:10, y = 1:10))
#' }
#' 
#' pmartR:::format_plot(pmartR::as.trellData(mymetabData), custom_plot = "custom_fn")
#' }
#' 
#' @author Rachel Richardson
#' 
#' 
custom_fn <- function(data){
  ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = 1:10, y = 1:10))
}

#' @name custom_fn_cog
#' @rdname custom_fn_cog
#' @title Test function for data_cog(), custom_cog argument
#'
#' @description Test function for data_cog(), custom_cog argument. Returns a tibble with one row and 3 columns.
#'
#' @param data data subset from a trellData object.
#' 
#' @examples
#'
#' dontrun{
#' library(pmartRdata)
#' library(pmartR)
#' library(ggplot2)
#' data("metab_object")
#' mymetabData <- pmartR::edata_transform(metab_object, "log10")
#' mymetabData <- pmartR::group_designation(metab_object, "Condition")
#' mymetabData <- pmartR::normalize_global(metab_object, "all", "median", apply_norm = TRUE, backtransform = TRUE)
#' 
#' custom_fn <- function(data){
#' ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = 1:10, y = 1:10))
#' }
#' 
#' custom_fn_cog <- function(data){
#' tibble::tibble(test1 = 1, test2 = 2, test3 = 3)
#' }
#' 
#' format <- pmartR::as.trellData(mymetabData)
#' plot <- pmartR:::format_plot(format, custom_plot = "custom_fn")
#' pmartR:::data_cogs(plot, format, custom_cog = "custom_fn_cog")
#' 
#' }
#' 
#' @seealso \link[pmartR]{data_cogs}
#' 
#' @author Rachel Richardson
custom_fn_cog <- function(data){
  tibble::tibble(test1 = 1, test2 = 2, test3 = 3)
}

#' @name format_plot
#' @rdname format_plot
#' @title Plot pairwise comparisons and data values in trellData object
#'
#' @description Plot pairwise comparisons and data values in trellData object. Customizable for plot types, y axis limits, paneling variable (what overall group is plotted on each graph arrangement), as well as desired variables for the y and x axis. Run in main trelliVis function.
#'
#' @param trellData An object of class "trellData" generated from \code{\link{as.trellData}}.
#' @param p_val Specifies p-value for setting graph border colors
#' @param panel_variable Specifies what to divide trelliscope panels by, must be a column in trellData. Defaults to cnames$edata_cname of trellData.
#' @param plot_text For plot_type == foldchange_bar only: TRUE/FALSE for p-value text above data
#' @param interactive Should the plots be rendered as plotly objects?
#' @param plot_type Types of plots. Restricted to "abundance_boxplot", "abundance_global", "abundance_heatmap", "foldchange_bar", "foldchange_global", "foldchange_heatmap", "missing_bar", "presence_heatmap"
#' @param y_limits y-limits for plots. Accepts scale, max, min, and range (see details and examples)
#' @param custom_plot User defined plotting function to be executed on specified data subsets. Other format_plot specifications do not apply to this plot. Should return a single plot per function call. Veiwing the data using as.trellData is highly encouraged to facillitate function development.
#'
#' @seealso \link[pmartR]{as.trellData}
#' @seealso \link[pmartR]{trelliVis}
#'
#' @details Descriptions of plot_type values and y-limits are as follows:
#' \tabular{ll}{
#' abundance_boxplot \tab Boxplots generated from trellData abundance values. Only available if omicsData was passed in as.trellData to generate trellData object. \cr
#' \tab \cr
#' abundance_global \tab  Biomolecule-specific abundance values compared to global abundances across all biomolecules. Only available if omicsData was passed in as.trellData to generate trellData object. \cr
#' \tab \cr
#' abundance_heatmap \tab  Heatmap of biomolecule abundances in with a mapping variable (e_meta) across samples. Only available if omicsData with e_meta was passed in as.trellData to generate trellData object. \cr
#' \tab \cr
#' foldchange_bar \tab Bar graphs generated from trellData foldchange values. Only available if omicsStats was passed in as.trellData to generate trellData object. \cr
#' \tab\cr
#' foldchange_global \tab Biomolecule-specific foldchange values compared to global foldchanges across all biomolecules. Only available if omicsStats was passed in as.trellData to generate trellData object. \cr
#' \tab \cr
#' foldchange_heatmap \tab  Heatmap of biomolecule foldchange values in with a mapping variable (e_meta) across samples. Only available if omicsStats AND omicsData was passed in as.trellData to generate trellData object. \cr
#' \tab \cr
#' missing_bar \tab Bar graph of the proportion of samples missing/present for each panel_variable. \cr
#' \tab\cr
#' presence_heatmap \tab Heatmap of biomolecule presence/absence in with a mapping variable (e_meta) across samples. Only available if omicsData with e_meta was passed in as.trellData to generate trellData object. \cr
#' }
#' Valid y_limits entries are as follows:
#' \tabular{ll}{
#' scale \tab Options include "free" or "fixed", where "free" allows each panel to auto-scale based on values and "fixed" uses the same scaling across all panels.\cr
#' \tab \cr
#' min \tab Minimum value on the y-axis. Where max argument is provided, must be less than max. \cr
#' \tab \cr
#' max \tab Maximum value on the y-axis. Where min argument is provided, must be greater than min. \cr
#' \tab \cr
#' range \tab A numeric defining the range of y-axis limits, centered on the median of values plotted, OR from a min/max value (if provided). \cr
#' }
#'
#' @examples
#'
#' dontrun{
#' library(pmartRdata)
#' library(pmartR)
#' library(ggplot2)
#' data("metab_object")
#' mymetabData <- pmartR::edata_transform(metab_object, "log10")
#' mymetabData <- pmartR::group_designation(metab_object, "Condition")
#' mymetabData <- pmartR::normalize_global(metab_object, "all", "median", apply_norm = TRUE, backtransform = TRUE)
#' format      <- pmartR::as.trellData(mymetabData)
#' 
#' pmartR:::format_plot(format, plot_type = "abundance_boxplot", y_limits = "free")
#' 
#' pmartR:::format_plot(format, plot_type = c("abundance_boxplot", "missing_bar"), y_limits = "fixed")
#' 
#' pmartR:::format_plot(format, plot_type = c("abundance_boxplot", "missing_bar"), y_limits = list(abundance_boxplot = list(min = 3, max = 7), missing_bar = list(max = 0.7)))
#' 
#' pmartR:::format_plot(format, plot_type = "abundance_boxplot", y_limits = list(abundance_boxplot = list(min = 3, max = 7)))
#' 
#' pmartR:::format_plot(format, plot_type = "abundance_boxplot", y_limits = list(abundance_boxplot = list(min = 3, max = 7)))
#' 
#' pmartR:::format_plot(format, plot_type = "abundance_boxplot", y_limits = list(abundance_boxplot = list(range = 7)))
#' 
#' data("isobaric_object")
#' mypepData  <- pmartR::edata_transform(isobaric_object, "log10")
#' mypepData  <- pmartR::normalize_isobaric(mypepData, apply_norm = T) # For isobaric normalization
#' mypepData  <- pmartR::group_designation(mypepData, "Group")
#' mypepData  <- pmartR::normalize_global(mypepData, "all", "median", apply_norm = TRUE, backtransform = TRUE)
#' myproRoll  <- pmartR::protein_quant(mypepData, "rrollup")
#' mypepStats <- pmartR::imd_anova(mypepData, test_method = "combined")
#' myproStats <- pmartR::imd_anova(myproRoll, test_method = "combined")
#' 
#' format2    <- pmartR::as.trellData(omicsData = mypepData, omicsStats = mypepStats)
#' format3    <- pmartR::as.trellData(omicsData = list(mypepData, myproRoll), omicsStats = list(mypepStats, myproStats))
#' 
#' custom_fn <- function(data){
#' ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = 1:10, y = 1:10))
#' }
#' 
#' pmartR:::format_plot(format2, plot_type = "foldchange_bar", panel_variable = "Protein", plot_text = TRUE, p_val = 0.001, interactive = TRUE)
#' pmartR:::format_plot(format3, plot_type = "abundance_boxplot", custom_plot = "custom_fn", panel_variable = c("Protein", "Protein"))
#' 
#' }
#'
#' @author Rachel Richardson
#'
format_plot <- function(trellData, 
                        panel_variable = NULL, 
                        plot_type = NULL, 
                        p_val = 0.05, 
                        plot_text = FALSE, 
                        y_limits = NULL,
                        interactive = FALSE,
                        custom_plot = NULL) {
  .format_plot(trellData, 
               panel_variable, 
               plot_type, 
               p_val, 
               plot_text, 
               y_limits,
               interactive,
               custom_plot)
}

.format_plot <- function(trellData, 
                         panel_variable = NULL, 
                         plot_type = NULL, 
                         p_val = 0.05, 
                         plot_text = FALSE, 
                         y_limits = NULL,
                         interactive = FALSE,
                         custom_plot = NULL) {
 
  ## Input checks
 
  if (is.null(panel_variable)){
    panel_variable <- attr(trellData, "cname")$edata_cname
  }
  
  if (is.null(p_val)){
    p_val <- 0.05
  }
   
  validate_format_plot_input(
    trellData = trellData, 
    panel_variable = panel_variable, 
    plot_type = plot_type, 
    p_val = p_val, 
    plot_text = plot_text, 
    interactive = interactive,
    custom_plot = custom_plot
  )
  
  # list_y_limits handles y_limits checks

  ## Input digestion
  if(!is.null(plot_type)){
    graphlims <- list_y_limits(plot_type = plot_type, y_limits = y_limits)
  } else if(!is.null(y_limits)){
    warning("y-limits are not applied to custom plots.")
  }
  
  
  generate_plots <- function(panel){
    
    #### Custom Plot ####
    if(!is.null(custom_plot)){
      trellData2 <- trellData
      if(!is.null(trellData2$data_values)){
        trellData2$data_values <- trellData2$data_values[as.character(trellData2$data_values[[panel_variable]]) == panel,]
      }
      if(!is.null(trellData2$comp_stats)){
        trellData2$comp_stats <- trellData2$comp_stats[as.character(trellData2$comp_stats[[panel_variable]]) == panel,]
        trellData2$summary_stats <- trellData2$summary_stats[as.character(trellData2$summary_stats[[panel_variable]]) == panel,]
      }
      plot_out <- eval(parse(text = paste0(custom_plot, "(trellData2)")))
      }
    #### Abundance Global ####
    if("abundance_global" %in% plot_type){
      
      lims <- graphlims[["abundance_global"]]
      
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
      
      # if(!(panel_variable %in% colnames(trellData$comp_stats))) stop(
      #   paste(
      #     "Panel variable for foldchange_global must be selected from either edata or emeta columns. The following columns are valid:", colnames(trellData$comp_stats)
      # ))
      
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
      
      # if (!is.null(lims$scale) || !is.null(lims$range)) {
      #   warning( 
      #     "Scale options and range are not valid with missing_bar plots, scale will not be applied. Please use max and min or list(NULL) for missing_bar y_limits.  Refer to examples in ?trelliVis().")
      # }
      # 
      # if (!is.null(lims$min) && (lims$min < 0 || lims$min > 1)) {
      #   stop( 
      #     "Minimum and maximum for missing_bar is only supported for proportions; select values between 0 and 1.")
      # }
      # 
      # if (!is.null(lims$max) && (lims$max < 0 || lims$max > 1)) {
      #   stop( 
      #     "Minimum and maximum for missing_bar is only supported for proportions; select values between 0 and 1.")
      # }
      
      
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
      # }
    }
    
    #### Foldchange Heatmap ####
    if("foldchange_heatmap" %in% plot_type){
      
      
      lims <- graphlims[["foldchange_heatmap"]]
      
      # if (!is.null(lims[[1]])) {
      #   warning( 
      #     "y_limits are not supported with heatmaps and will not be used.  Refer to examples in ?trelliVis().")
      # }
        
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
    
    #### Presence Heatmap ####
    if("presence_heatmap" %in% plot_type){
      
      lims <- graphlims[["presence_heatmap"]]
      
      # if (!is.null(lims[[1]])) {
      #   warning( 
      #     "y_limits are not supported with heatmaps and will not be used.  Refer to examples in ?trelliVis().")
      # }
          
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
                                            y_axis = abundance)
      
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
                                              y_axis = "Fold_change"
                                              )
        
        keep <- colnames(nestedData_comps) %in% c(
          pmartR::get_edata_cname(trellData),
          "Fold_change",
          "text",
          "Comparison",
          "bord",
          "labels")
        
        nestedData_comps <- unique.data.frame(nestedData_comps[keep])
        
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
    
    if(interactive && !is.null(plot_type)){
      
      if(any(stringr::str_detect(plot_type, "heatmap"))) {
        
        # plotnames <- c("abundance_boxplot", 
        #                 "foldchange_bar", 
        #                 "abundance_global", 
        #                 "foldchange_global", 
        #                 "missing_bar",
        #                 "abundance_heatmap",
        #                 "foldchange_heatmap",
        #                 "presence_heatmap",
        #                 "foldchange_bar",
        #                 "abundance_boxplot")
        # typenames <- plotnames[plotnames %in% plot_type]
        # #Accounts for out of order input
        # 
        # if(stringr::str_detect(typenames[1], "heatmap")){

          # plot_out <- plot_out + ggplot2::theme(legend.position = "none")
          plot_out <- plotly::ggplotly(plot_out, tooltip = c("text"))
          plot_out <- plot_out %>% plotly::layout(plot_bgcolor='grey',
                                                  xaxis = list(showgrid = F),
                                                  yaxis = list(showgrid = F)
                                                  )
          
        # }

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
#' @param interactive Should the plots be rendered as plotly objects?
#' @param plot_type types of plots, specified in format_plot
#' @param y_limits y-limits, specified in format_plot
#' @param custom_plot User defined plotting function to be executed on specified data subsets. Other format_plot specifications do not apply to this plot. Should return a single plot per function call. Veiwing the data using as.trellData is highly encouraged to facillitate function development.
#'
#' @seealso \link[pmartR]{format_plot}
#'
#' @author Rachel Richardson


validate_format_plot_input <- function(trellData,
                                       plot_type,
                                       p_val,
                                       panel_variable,
                                       plot_text,
                                       custom_plot,
                                       interactive){
  .validate_format_plot_input(trellData,
                              plot_type,
                              p_val,
                              panel_variable,
                              plot_text,
                              custom_plot,
                              interactive)
}

.validate_format_plot_input <- function(
  trellData,
  plot_type,
  p_val,
  panel_variable,
  plot_text,
  custom_plot,
  interactive){

  ##### Initial checks #####

  # Check if class is correct #
  if(!inherits(trellData, "trellData")) stop(
    "trellData must be of the class 'trellData'")

  # Check if data is in trellData #
  if((is.null(trellData$comp_stats) || is.na(trellData$comp_stats)) &&
     (is.null(trellData$data_values) || is.na(trellData$data_values))) stop(
    "No data values or comparison statistics in trellData to plot")

  # # Check if data is in trellData (as above) #
  # if(is.na(trellData$comp_stats) && is.na(trellData$data_values)) stop(
  #   "No data values or comparison statistics in trellData to plot")
  
  # Check if p_val is numeric of length 1 #
  if(!is.numeric(p_val) || length(p_val) != 1) stop(
    "p_val must be a numeric of length 1")

  # Check if plotly is logical of length 1 #
  if(!is.logical(interactive) || (length(interactive) != 1)) stop(
    "interactive must be a logical (TRUE or FALSE) of length 1")

  # Check if plot_text is logical of length 1 #
  if(!is.logical(plot_text) || (length(plot_text) != 1)) stop(
    "plot_text must be a logical (TRUE or FALSE) of length 1")

  # Check if custom_plot is a string of length 1 #
  if(!is.null(custom_plot) && (!is.character(custom_plot) || (length(custom_plot) != 1))) stop(
    "custom_plot must be a character string (TRUE or FALSE) of length 1. Example: custom_plot = 'custom_fn' ")
  
  # Check if custom_plot function exists #
  if(!is.null(custom_plot) && !exists(custom_plot)) stop(
    paste("custom_plot function not found! Input custom_plot function name:", custom_plot))
  
  # Check plot_type input
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
    "plot_type specified is not supported. Must be a character string in the following: abundance_boxplot, foldchange_bar, abundance_global, foldchange_global, missing_bar, abundance_heatmap, foldchange_heatmap, presence_heatmap, foldchange_bar, abundance_boxplot"
  )
  
  # Check plot_type length #
  if(!is.null(plot_type) && (length(plot_type) != 1 )) stop(
    "Invalid plot_type input. Refer to examples in ?pmartR:::format_plot."
  )
  
  # Check if plotting info is available #
  if(is.null(plot_type) && is.null(custom_plot)) stop(
    "No plotting information found in plot_type or custom_plot."
  )
  
  # Check heatmap requirements #
  if (!is.null(plot_type) && 
      panel_variable == pmartR::get_edata_cname(trellData) && 
      stringr::str_detect(plot_type, "heatmap")) {
    stop(
      paste("Heatmaps require a panel_variable that is not the same as edata cname (Default panel_variable is set to edata cname). Current edata cname:", pmartR::get_edata_cname(trellData))
    )
  } 
  
  temp <- Reduce("|", list(unlist(lapply(trellData, is.data.frame)), 
                       unlist(lapply(trellData, is.null))))
  
  # Check structure of trellData
  if(!all(temp)) stop(
    "All trellData fields must be NULL or data.frames."
  )
  if((!is.null(trellData$comp_stats) && is.null(trellData$summary_stats) ) || 
     (is.null(trellData$comp_stats)  && !is.null(trellData$summary_stats))) stop(
    "Where comp_stats is generated, summary_stats are required and vice versa."
  )
  if(!is.null(trellData$comp_stats) && 
     (is.null(attr(trellData, "statistical_test")) || 
      !(attr(trellData, "statistical_test") %in% c("combined", "anova", "gtest")))) stop(
    "Statistical_test attribute must be set to combined, anova, or gtest."
  )
  
  # Check if panel variable is correct #
  allcol <- lapply(trellData, colnames)
  if (!all(panel_variable %in% unlist(allcol))) stop(
    paste("Panel_variable is not valid. String must be in the following: ", toString(unique(unlist(allcol))))
  )
  
  # Check plot_type input
  if(is.null(trellData$comp_stats) && 
     !is.null(plot_type) &&
     plot_type %in% c(
    "foldchange_global", 
    "foldchange_heatmap",
    "foldchange_bar"
  )) stop(
    "Statistical comparisons (omicsStats) are required for plot_type specified (foldchange_global, foldchange_heatmap, or foldchange_bar). Refer to pmatR::imd_anova() for computing statistics."  )
  
  # Check plot_type input
  if(is.null(trellData$data_values) && 
     !is.null(plot_type) &&
     plot_type %in% c(
    "abundance_global", 
    "abundance_heatmap",
    "presence_heatmap",
    "abundance_boxplot"
  )) stop(
    "Data values (omicsData) are required for plot_type specified (abundance_global, abundance_heatmap, presence_heatmap, abundance_boxplot).")
  
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
#' @param y_axis Specifies what column should be on the y-axis, must be a column in trellData and numeric. 
#'
#' @seealso \link[pmartR]{format_plot}
#'
#' @author Rachel Richardson
#'
#'
#'

add_plot_features <- function(trellData,
                              nestedData,
                              p_val = NULL, 
                              panel_variable = NULL,
                              y_axis){
  .add_plot_features(trellData,
                     nestedData,
                     p_val, 
                     panel_variable,
                     y_axis)
}

.add_plot_features <- function(trellData,
                               nestedData,
                               p_val = NULL, 
                               panel_variable = NULL,
                               y_axis){
  
  colors <- c(NA, "grey40", "black")
  
  if(y_axis == "Fold_change"){
    
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
    if (!(y_axis %in% hover_want_comps)){
      hover_want_comps <- c(hover_want_comps, y_axis)
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
        } else if (label %in% c("Fold_change", y_axis)){
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
    if (!(y_axis %in% hover_want)){
      hover_want <- c(hover_want, y_axis)
    }
    if (panel_variable %in% hover_want){
      hover_want <- hover_want[!(hover_want %in% panel_variable)]
    }
    hover_labs <- hover_want[hover_want %in% colnames(nestedData)]
    
    text_labs <- purrr::map(1:nrow(nestedData), function(row){
      row_text <- purrr::map(hover_labs, function(label){
        if (label %in% c(grep("abundance", 
                              colnames(trellData$data_values), 
                              value = TRUE), y_axis)) {
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
#' @seealso \link[pmartR]{format_plot}
#'
#' @author Rachel Richardson
#'
#'

list_y_limits <- function(plot_type, y_limits){
  .list_y_limits(plot_type, y_limits)
}

.list_y_limits <- function(plot_type, y_limits){
  
  if(!(plot_type %in% c(
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
    "plot_type specified is not supported. Must be a character string in the following: abundance_boxplot, foldchange_bar, abundance_global, foldchange_global, missing_bar, abundance_heatmap, foldchange_heatmap, presence_heatmap, foldchange_bar, abundance_boxplot"
  )
    
    limlist <- list()
    plotlims <- list()
    plotlims2 <- list()
    
    ## Set if null
    if(is.null(y_limits)){
      y_limits <- "fixed"
      nullflag <- TRUE
    } else {
      nullflag <- FALSE
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
      
      plotlims <- list(y_limits)
      names(plotlims) <- plot_type
      
      ### Singular inputs
    } else if (y_limits == "fixed" || y_limits == "free"){
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
      
      if (any(duplicated(names(y_limits)))) stop(
        paste("List inputs may not have duplicate names. Use a list of lists to specify different plot y_limits. Duplicate names:",
              names(y_limits)[duplicated(names(y_limits))])
      )  
      
      if (any(!(plot_type %in% names(y_limits)))) stop(
        paste("List of list inputs must be named. Names in y_limits should match plot_type. plot_type:",
              paste(plot_type, collapse = ", "))
      )
      
      y_limits2 <- y_limits[plot_type]
      plotlims <- y_limits2[unlist(lists)]
      
      if(!is.null(plotlims) && length(plotlims) > 1){
        purrr::map(plotlims, function(listitem){
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
      nonlist <- names(unlist(lists))[!as.logical(unlist(lists))]
      plotlims2 <- purrr::map(nonlist, function(nonlistname){
        limit <- y_limits[[nonlistname]]
        listlims2 <- list()
        if (limit == "fixed" || limit == "free"){
          listlims2$scale <- limit
        } else if (is.numeric(limit) && length(limit) == 2){
          listlims2$min <- min(limit)
          listlims2$max <- max(limit)
        } else stop(
          "For mixed input of lists and non-lists, non-list input is restricted to c(min, max), 'free', or 'fixed'"
        )
        return(listlims2)
      })
    }
    
    if(length(plotlims) == 0 && length(plotlims2) == 0) stop(
      "Invalid input. y_limits input is restricted to named lists for each plot type, strings 'fixed' or 'free', or a numeric vector of length 2 specifying maximum and minimum y_values. Refer to ?trelliVis for examples."
    )
    
    if(!exists("nonlist")){
      nonlist <- plot_type[!unlist(lists)]
    }
    
    output <- c(plotlims, plotlims2)
    nms <- c(plot_type[unlist(lists)], nonlist)
    nms <- nms[!is.na(nms)]
    names(output) <- nms
    
    check <- output[[plot_type]]
    
    
    
    if("min" %in% names(check) && 
       (length(check[["min"]]) != 1 ||
       !is.numeric(check[["min"]]))) stop(
      "Min must be a numeric of length 1."
    )
    if("max" %in% names(check) && 
       (length(check[["max"]]) != 1 ||
        !is.numeric(check[["max"]]))) stop(
          "Max must be a numeric of length 1."
        )
    
    if("range" %in% names(check) && 
       (length(check[["range"]]) != 1 ||
        !is.numeric(check[["range"]]) ||
        !(check[["range"]] > 0))) stop(
          "Range must be a numeric of length 1 that is greater than 0."
        )
    if("scale" %in% names(check) && 
       (length(check[["scale"]]) != 1 ||
        !is.character(check[["scale"]]) ||
        !(check[["scale"]] %in% c("free", "fixed")))) stop(
          "Scale must be a character of length 1 in 'free' or 'fixed'."
        )
    
    if("min" %in% names(check) && 
       "max" %in% names(check) && 
       !(check[["min"]] < check[["max"]])) stop(
         "Set min must be less than set max."
       )
    if("range" %in% names(check) && 
       "scale" %in% names(check) &&
       ("max" %in% names(check) || "min" %in% names(check))) stop(
         "Range y-limits are not supported with scale y-limits and set min/max."
       )
    if("min" %in% names(check) && 
       "max" %in% names(check) && 
       "range" %in% names(check)) stop(
         "Range y-limits are not supported with set min and max y-limits."
       )
    if("min" %in% names(check) && 
       "max" %in% names(check) && 
       "scale" %in% names(check)) stop(
         "Scale y-limits are not supported with set min and max y-limits."
       )
    
    if (plot_type %in% c("abundance_heatmap", "foldchange_heatmap", "presence_heatmap") && !nullflag) {
      warning( 
        "y_limits are not supported with heatmaps and will not be used.  Refer to examples in ?trelliVis().")
    }
    
    if (plot_type %in% c("abundance_global", "foldchange_global") && !is.null(check$scale) && check$scale == "free") {
      warning( 
        "Scale option 'free' is not valid with global plots. Fixed plots will be generated based on global values.")
      output[[plot_type]]$scale <- "fixed"
    }
    
    if (plot_type == "missing_bar" && (!is.null(check$scale) || !is.null(check$range))) {
      warning( 
        "Scale options and range are not valid with missing_bar plots, scale will not be applied. Please use max and min or list(NULL) for missing_bar y_limits.  Refer to examples in ?trelliVis().")
    }
    
    if (plot_type == "missing_bar" && (!is.null(check$min) && (check$min < 0 || check$min > 1))) {
      stop( 
        "Minimum and maximum for missing_bar is only supported for proportions; select values between 0 and 1.")
    }
    
    if (plot_type == "missing_bar" && (!is.null(check$max) && (check$max < 0 || check$max > 1))) {
      stop( 
        "Minimum and maximum for missing_bar is only supported for proportions; select values between 0 and 1.")
    }
    
    return(output)
  }

#' @name set_increment
#' @rdname set_increment
#' @title Sets y-axis increment for trellData labels plotting
#' 
#' @description Sets y-axis increment for trellData plotting. Used by plot_comp.
#'
#' @param yvalues y-values for plotting
#' @param include_zero Should zero be included in plot scaling?
#'
#' @seealso \link[pmartR]{format_plot}
#' @seealso \link[pmartR]{list_y_limits}
#' @seealso \link[pmartR]{set_ylimits}
#'
#' @author Rachel Richardson
#'
#'
set_increment <- function(yvalues, include_zero = TRUE){
  .set_increment(yvalues, include_zero)
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
#' @param include_zero Should zero be included regardless of scaling?
#'
#' @seealso \link[pmartR]{format_plot}
#' @seealso \link[pmartR]{list_y_limits}
#' @seealso \link[pmartR]{set_increment}
#'
#' @author Rachel Richardson
#'
#'
#' 
set_ylimits <- function(yvalues, 
                        increment, 
                        y_range = NULL, 
                        y_max = NULL, 
                        y_min = NULL, 
                        include_zero = TRUE){
  .set_ylimits(yvalues, 
               increment, 
               y_range, 
               y_max, 
               y_min, 
               include_zero)
}

.set_ylimits <- function(yvalues, 
                         increment, 
                         y_range = NULL, 
                         y_max = NULL, 
                         y_min = NULL, 
                         include_zero = TRUE
                         ){
  
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
    
    # if(maxi < mini) stop ("Invalid max and min. Max < Min")
    # 
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
#' @description Plot pairwise comparisons and data values in trellData object. Customizable for plot types, y axis limits, paneling variable (what overall group is plotted on each graph arrangement), as well as desired variables for the y and x axis. Called in trelliVis function.
#'
#' @param nested_plot A nested table generated from trellData using formatplot()
#' @param trellData A nested table generated from trellData using formatplot()
#' @param p_val Numeric that specifies p-value for Boolean significance cognotic. Default is 0.05.
#' @param try_URL Will attempt to link to PubChem, LipidMaps, or Uniprot based on information in specified col. Default is NULL.
#' @param custom_cog The name of a user generated function with cognostics. Function should act on a subset of data and output a dataframe (or tibble or equivelent) with one row (summary of rows).
#'
#' @seealso \link[pmartR]{as.trellData}
#' @seealso \link[pmartR]{format_plot}
#' @seealso \link[pmartR]{trelliVis}
#'
#' @examples
#'
#' dontrun{
#' library(pmartRdata)
#' library(pmartR)
#' library(ggplot2)
#' data("metab_object")
#' mymetabData <- pmartR::edata_transform(metab_object, "log10")
#' mymetabData <- pmartR::group_designation(metab_object, "Condition")
#' mymetabData <- pmartR::normalize_global(metab_object, "all", "median", apply_norm = TRUE, backtransform = TRUE)
#' format      <- pmartR::as.trellData(mymetabData)
#' plot        <- pmartR:::format_plot(format, plot_type = "abundance_boxplot", y_limits = "free")
#' 
#' custom_fn_cog <- function(data){
#' tibble::tibble(test1 = 1, test2 = 2, test3 = 3)
#' }
#' 
#' pmartR::data_cogs(plot, format, try_URL = TRUE, custom_cog = "custom_fn_cog")
#' }
#'
#' @author Rachel Richardson
#'


data_cogs <- function(nested_plot,
                      trellData,
                      p_val = 0.05,
                      try_URL = NULL,
                      custom_cog = NULL) {
  .data_cogs(nested_plot,
             trellData,
             p_val,
             try_URL,
             custom_cog)
}

.data_cogs <- function(nested_plot,
                       trellData,
                       p_val = 0.05,
                       try_URL = NULL,
                       custom_cog = NULL){
  
  ## Variable checks, also checked at start of trelliVis
  ## Extensive testing for trellData executed in format plot validation function
  ## nested_plot should be piped from format_plot directly with trellData and nested_plot
  ## custom cog tested breifly at the beginning of trelliVis
  
  # Check check if p_val is numeric of length 1 #
  if(!is.numeric(p_val) || length(p_val) != 1) stop(
    "p_val must be a numeric of length 1")  
  
  # Check check if try_URL is string of length 1 #
  if(!is.null(try_URL) && (!is.character(try_URL) || length(try_URL) != 1)) stop(
    "try_URL must be a string of length 1")  
  
  # Check check if nested_plot is dataframe from format_plot  #
  if(!is.data.frame(nested_plot) || ncol(nested_plot) != 2 || colnames(nested_plot)[2] != "panel") stop(
    "nested_plot must be a data.frame passed from format_plot; requires ncol == 2, and colnames(nested_plot)[2] == 'panel'")  
  
  # Check check if trellData correct class #
  if(!inherits(trellData, "trellData")) stop(
    "trellData must be class trellData.")
  
  # Check check if custom cog is a string of length 1 #
  if(!is.null(custom_cog) && (!is.character(custom_cog) || length(custom_cog) != 1)) stop(
    "custom_cog argument must be a string of length 1, defining the name of the custom function.")  
  
  
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
          
          addOrigCogs <- dplyr::left_join(panel_values, 
                                          panel_comp, 
                                          by = intersect(colnames(panel_values),
                                                         colnames(panel_comp)))
          addOrigCogs <- dplyr::left_join(addOrigCogs,
                                          panel_summ,
                                          by = intersect(colnames(addOrigCogs),
                                                         colnames(panel_summ)))
          
        } else if (!is.null(trellData$data_values)){
          panel_values <- trell_values[trell_values[[panel_variable]] == panel,]
          addOrigCogs <- panel_values
          
        } else {
          panel_comp <- trell_comp[trell_comp[[panel_variable]] == panel,]
          panel_summ <- trell_summ[trell_summ[[panel_variable]] == panel,]
          # bycs  <- colnames(panel_comp)[colnames(panel_comp) %in% colnames(panel_summ)]
          
          addOrigCogs <- dplyr::left_join(panel_comp,
                                          panel_summ,
                                          by = intersect(colnames(panel_summ),
                                                         colnames(panel_comp)))
        }
        cogs <- addOrigCogs
        
        if(stats){
          stat_test <- attr(trellData, "statistical_test")
          if(stat_test == "combined"){
            cogs <- dplyr::mutate(
              cogs, 
              Sig_T_p_value = purrr::map_lgl(
                as.numeric(cogs[["P_value_T"]]),
                function(value) !(value > p_val) && !is.na(value)),
              Sig_G_p_value = purrr::map_lgl(
                as.numeric(cogs[["P_value_G"]]),
                function(value) !(value > p_val) && !is.na(value)))
          } else if (stat_test == "anova"){
            cogs <- dplyr::mutate(
              cogs, 
              Sig_T_p_value = purrr::map_lgl(
                as.numeric(cogs[["P_value_T"]]),
                function(value) !(value > p_val) && !is.na(value)))
          } else {
            cogs <- dplyr::mutate(
              cogs, 
              Sig_G_p_value = purrr::map_lgl(
                as.numeric(cogs[["P_value_G"]]),
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
          
          if(!is.null(try_URL)){
            cogs <- dplyr::mutate(
              cogs,
              Protein_URL = purrr::map_chr(
                cogs[[try_URL]],
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
                cogs[[try_URL]],
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
                cogs[[try_URL]],
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
                cogs[[try_URL]],
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
                  desc = "Link to protein UniProt page",
                  default_label = TRUE))
              }
              
            } else if (column == "Metabolite_URL") {
              return(trelliscopejs::cog_href(
                toString(unique(cogs[[column]])),
                desc = "Link to metabolite PubChem page.",
                default_label = TRUE))
              
            } else if (column == "Lipid_Maps_Hits") {
              return(trelliscopejs::cog_href(
                toString(unique(cogs[[column]])),
                desc = "Link to Lipid Maps search of lipid species.", 
                default_label = TRUE))
            } else if (column == "peps_per_pro"){
              return(trelliscopejs::cog(
                mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE),
                desc = "Number of peptides mapped to a given protein (used in pmartR::protein_quant() rollup method). Average is computed where multiple proteins are in a panel. (Value of 0 == NA) ",
                default_label = TRUE))
              
            } else if (column == "n_peps_used"){
              return(trelliscopejs::cog(
                mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE),
                desc = "Number of peptides used in pmartR::protein_quant() rollup method. Average is computed where multiple proteins are in a panel. (Value of 0 == NA) "))
              
            } else if (column == "Mean"){
              return(trelliscopejs::cog(
                toString(unique(cogs[[column]])),
                desc = "Mean of abundances/intensities for experimental groups where statistical tests could be performed."))
              
            } else if (column == "Fold_change"){
              return(trelliscopejs::cog(
                mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE),
                desc = "Mean of foldchange difference (per panel). (Value of 0 == NA) "))
              
            } else if (column == "Sig_T_p_value") {
              return(trelliscopejs::cog(
                toString(any(cogs[[column]])),
                desc = "Significant ANOVA results based on input p-value threshold (default == 0.05)."))
              
            } else if (column == "P_value_G") {
              return(trelliscopejs::cog(
                mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE),
                desc = "Mean G-test p-values (per panel)."))
              
            } else if (column == "P_value_T") {
              return(trelliscopejs::cog(
                mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE),
                desc = "Mean ANOVA p-values (per panel)."))
              
            } else if (column == "Sig_G_p_value") {
              return(trelliscopejs::cog(
                toString(any(suppressWarnings(as.logical(cogs[[column]])))),
                desc = "Significant G-test results based on input p-value threshold (default == 0.05)."))
              
            } else if (length(abundance) > 0 && column == abundance) {
              return(trelliscopejs::cog(
                mean(suppressWarnings(as.numeric(cogs[[column]])), na.rm = TRUE),
                desc = "Average abundance/intensity of values."))
              
            } else {
              return(trelliscopejs::cog(toString(unique(cogs[[column]])), desc = "User defined variable. (Categorical)"))
            }
            
          } else if (is.numeric(cogs[[column]])){
            return(trelliscopejs::cog(
              mean(cogs[[column]], na.rm = TRUE), 
              desc = "User defined variable. (Numeric mean, Value of 0 == NA) "))
            
          } else if (is.logical(cogs[[column]]) || 
                     !any(is.na(suppressWarnings(as.logical(cogs[[column]]))))){
            return(trelliscopejs::cog(toString(any(cogs[[column]])), 
                                      desc = "User defined variable. (Categorical)"))
            
          } else if (is.character(cogs[[column]]) && 
                     !any(is.na(suppressWarnings(as.numeric(cogs[[column]]))))){
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
        
        if(!is.null(custom_cog)){ ## Run custom function on cogs using function name
         x <- eval(parse(text = paste0(custom_cog, "(cogs)")))
         cog_out <- cbind(cog_out,x)
        }
        
        return(cog_out)
      }))
  
  return(out)
  
}

#' @name trelliVis
#' @rdname trelliVis
#' @title Automated trelliscopejs plotting from omicsData
#'
#' @description Generates trelliscopejs displays using pairwise comparisons and data values of omicsData and omicsStats objects. Customizable for plot types, y axis limits, paneling variable (what overall group is plotted on each graph arrangement), as well as desired variables for the y and x axis.
#'
#' @param omicsData A pmartR object of class pepData, lipidData, metabData, or proData. Can use list(pepData, proData) for associated data.
#' @param omicsStats A statistical results object produced by running \code{imd_anova} on omicsData. Can use list(pepStats, proStats) for associated data.
#' @param omicsFormat Output of as.trellData() function
#' @param p_val Numeric that specifies p-value for significance calculations. Default is 0.05.
#' @param panel_variable String: Name of column that plot panels are sorted by (e.g. each plotting arrangement has a unique identifier from panel variable). Default is emeta_cname if present, edata_cname where emeta_cname is not present.
#' @param try_URL Will attempt to link to PubChem, LipidMaps, or Uniprot based on information in specified column. Default is NULL.
#' @param trelli_name String: name of display, or list of names where a list is provided for omicsData and omicsStats
#' @param trelli_path_out String: path to where trelliscope is stored. Default is "./TrelliDisplay"
#' @param interactive Should the plots be rendered as plotly objects?
#' @param plot_text Disable plot text
#' @param y_limits Y limits 
#' @param plot_type plots for plotting
#' @param self_contained Should display be generated in document? Defaults to FALSE
#' @param custom_cog The name of a user generated function with cognostics. Function should act on a subset of data and output a dataframe (or tibble or equivelent) with one row (summary of rows).
#' @param custom_plot User defined plotting function to be executed on specified data subsets. Other format_plot specifications do not apply to this plot. Should return a single plot per function call. Veiwing the data using as.trellData is highly encouraged to facillitate function development.
#' @param display When FALSE, will return arguments to be passed into trelliscopejs::trelliscope() as a list without generating a display.
#' @param ... Additional arguments for trelliscope() function; trelliVis supports arguments jsonp, split_sig, auto_cog, height, width, desc, md_desc. Argument panel_col is currently not supported and may produce errors. Arguments state, group, ncol, nrow, and thumb are preset (display = FALSE is useful for modifying these). Refer to ?trelliscopejs::trelliscope()
#'
#' @details Descriptions of plot_type values and y-limits are as follows:
#' \tabular{ll}{
#' abundance_boxplot \tab Boxplots generated from trellData abundance values. Only available if omicsData was passed in as.trellData to generate trellData object. \cr
#' \tab \cr
#' abundance_global \tab  Biomolecule-specific abundance values compared to global abundances across all biomolecules. Only available if omicsData was passed in as.trellData to generate trellData object. \cr
#' \tab \cr
#' abundance_heatmap \tab  Heatmap of biomolecule abundances in with a mapping variable (e_meta) across samples. Only available if omicsData with e_meta was passed in as.trellData to generate trellData object. \cr
#' \tab \cr
#' foldchange_bar \tab Bar graphs generated from trellData foldchange values. Only available if omicsStats was passed in as.trellData to generate trellData object. \cr
#' \tab\cr
#' foldchange_global \tab Biomolecule-specific foldchange values compared to global foldchanges across all biomolecules. Only available if omicsStats was passed in as.trellData to generate trellData object. \cr
#' \tab \cr
#' foldchange_heatmap \tab  Heatmap of biomolecule foldchange values in with a mapping variable (e_meta) across samples. Only available if omicsStats AND omicsData was passed in as.trellData to generate trellData object. \cr
#' \tab \cr
#' missing_bar \tab Bar graph of the proportion of samples missing/present for each panel_variable. \cr
#' \tab\cr
#' presence_heatmap \tab Heatmap of biomolecule presence/absence in with a mapping variable (e_meta) across samples. Only available if omicsData with e_meta was passed in as.trellData to generate trellData object. \cr
#' }
#' Valid y_limits entries are as follows:
#' \tabular{ll}{
#' scale \tab Options include "free" or "fixed", where "free" allows each panel to auto-scale based on values and "fixed" uses the same scaling across all panels.\cr
#' \tab \cr
#' min \tab Minimum value on the y-axis. Where max argument is provided, must be less than max. \cr
#' \tab \cr
#' max \tab Maximum value on the y-axis. Where min argument is provided, must be greater than min. \cr
#' \tab \cr
#' range \tab A numeric defining the range of y-axis limits, centered on the median of values plotted, OR from a min/max value (if provided). \cr
#' }
#'
#' @seealso \link[pmartR]{as.trellData}
#' @seealso \link[trelliscopejs]{trelliscope}
#'
#'
#' @examples
#'
#' dontrun{
#' library(pmartRdata)
#' library(pmartR)
#' library(ggplot2)
#' data("metab_object")
#' mymetabData <- pmartR::edata_transform(metab_object, "log10")
#' mymetabData <- pmartR::group_designation(metab_object, "Condition")
#' mymetabData <- pmartR::normalize_global(metab_object, "all", "median", apply_norm = TRUE, backtransform = TRUE)
#' 
#' pmartR::trelliVis(omicsData = mymetabData, plot_type = "abundance_boxplot", y_limits = "fixed")
#' pmartR::trelliVis(omicsFormat = pmartR::as.trellData(mymetabData), plot_type = "abundance_boxplot", y_limits = "fixed")
#' pmartR::trelliVis(mymetabData, plot_type = c("abundance_boxplot", "missing_bar"), y_limits = list(abundance_boxplot = list(min = 3, max = 7), missing_bar = list(max = 0.7)))
#' pmartR::trelliVis(mymetabData, y_limits = list(abundance_boxplot = list(min = 3, range = 7))))
#' pmartR::trelliVis(mymetabData, width = 700, height = 200, ncol = 2)
#' 
#' data("isobaric_object")
#' mypepData  <- pmartR::edata_transform(isobaric_object, "log10")
#' mypepData  <- pmartR::normalize_isobaric(mypepData, apply_norm = T) # For isobaric normalization
#' mypepData  <- pmartR::group_designation(mypepData, "Group")
#' mypepData  <- pmartR::normalize_global(mypepData, "all", "median", apply_norm = TRUE, backtransform = TRUE)
#' myproRoll  <- pmartR::protein_quant(mypepData, "rrollup")
#' mypepStats <- pmartR::imd_anova(mypepData, test_method = "combined")
#' myproStats <- pmartR::imd_anova(myproRoll, test_method = "combined")
#' 
#' pmartR::trelliVis(omicsData = mypepData, omicsStats = mypepStats, try_URL = TRUE, self_contained = TRUE)
#' pmartR::trelliVis(omicsData = list(mypepData, myproRoll), omicsStats = list(mypepStats, myproStats), trelli_name = "test", trelli_path_out = "./Here")
#' 
#' custom_fn <- function(data){
#' ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = 1:10, y = 1:10))
#' }
#' 
#' custom_fn_cog <- function(data){
#' tibble::tibble(test1 = 1, test2 = 2, test3 = 3)
#' }
#' 
#' pmartR::trelliVis(omicsData = mypepData, omicsStats = mypepStats, plot_type = "foldchange_bar", panel_variable = "Protein", plot_text = TRUE, p_val = 0.001, interactive = TRUE)
#' pmartR:::format_plot(omicsData = list(mypepData, myproRoll), omicsStats = list(mypepStats, myproStats), plot_type = "abundance_boxplot", custom_plot = "custom_fn", costom_cog = "custom_fn_cog", panel_variable = c("Protein", "Protein"))
#' 
#' }
#'
#'
#' @author Rachel Richardson
#' @export
trelliVis <- function(omicsData = NULL, omicsStats = NULL,
                      omicsFormat = NULL, p_val = 0.05,
                      panel_variable = NULL, 
                      try_URL = NULL, trelli_name = NULL,
                      trelli_path_out = "TrelliDisplay", 
                      plot_text = FALSE, interactive = FALSE,
                      y_limits = NULL, plot_type = NULL,
                      self_contained = FALSE,
                      custom_cog = NULL,
                      custom_plot = NULL,
                      display = TRUE,
                      ...) {
  
  .trelliVis(omicsData, omicsStats,
             omicsFormat, p_val,
             panel_variable, 
             try_URL, trelli_name,
             trelli_path_out, 
             plot_text, interactive,
             y_limits, plot_type,
             self_contained,
             custom_cog,
             custom_plot,
             display,
             ...)
}

.trelliVis <- function(omicsData = NULL, omicsStats = NULL,
                       omicsFormat = NULL, p_val = 0.05,
                       panel_variable = NULL, 
                       try_URL = NULL, trelli_name = NULL,
                       trelli_path_out = "TrelliDisplay", 
                       plot_text = FALSE, interactive = FALSE,
                       y_limits = NULL, plot_type = NULL,
                       self_contained = FALSE,
                       custom_cog = NULL,
                       custom_plot = NULL, 
                       display = TRUE,
                       ...) {
  
  
  #####
  ## Initial checks ##
  #####
  
  #store_object, custom_cog_df, plot package = ggplot, rbokeh, etc, trelliscope additional arguments
  
  if(is.null(omicsFormat) && is.null(omicsData) && is.null(omicsStats)) stop(
    "At least one of omicsData, omicsStats, or omicsFormat must be populated."
  )
  
  # Check user inputs for bad paths
  if(!is.null(trelli_name)){
    if(any(stringr::str_detect(trelli_name, "[^A-z0-9_\\/\\.\\\\]+"))) warning(
      "Caution: Non-word characters (outside of periods and slashes) detected in trelliname that may cause errors in file creation. Ensure that trelliname is a permissable file name for your system."
    )
  }

  if(any(stringr::str_detect(trelli_path_out, "[^A-z0-9_\\/\\.\\\\]+"))) warning(
    "Caution: Non-word characters (outside of periods and slashes) detected in trelli_path_out that may cause errors in file creation. Ensure that trelli_path_out is a permissable folder name for your system."
    )
  
  # Check user input for ... entry  
  if (any(c("state", "self_contained", "path", "name", "panel_col", "group", "ncol", "nrow") %in% names(list(...)))) stop(
    paste("The following trelliscopejs arguments are currently not supported for user input: ", toString(c("ncol", "nrow", "state", "self_contained", "path", "name", "panel_col", "group")))
  )
  
  
  
  #####
  ## Switch Stats and Omics data as appropriate
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
  if ((class(omicsData) == "list" || 
       class(omicsStats) == "list") && 
      !is.null(panel_variable) && 
      length(panel_variable) != max(length(omicsStats), 
                                    length(omicsData))) stop(
                                      "Panel variable must be specified for each index in omicsStats/omicsData (e.g. c('Peptide', 'Protein'))"
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
  
  # Check check if try_URL is string of length 1 #
  if(!is.null(try_URL) && (!is.character(try_URL) || length(try_URL) != 1)) stop(
    "try_URL must be a string of length 1")  
  
  if (is.null(omicsFormat)){
    # Re-format objects for plotting #
    trellData <- as.trellData(omicsData, omicsStats)
  } else {
    trellData <- omicsFormat
  }
  
  # If a pep/pro pair is listed, act on each item #
  if(class(trellData) == "list"){
    
    # Check linking function in trelliscopejs
    if(is.null(trelliscopejs::cog_disp_filter)) stop(
      "Function trelliscopejs::cog_disp_filter() not found! This function is required for linked displays. Consider updating your trelliscopejs package or installing from the dev branch using devtools::install_github('hafen/trelliscopejs@dev'). (The list input of omicsData and/or omicsStats can be processed as unlinked displays in trelliVis for circumstances where updates are not available)"
    )
    
    # Validate Custom cognostics if provided
    if(!is.null(custom_cog)){
      purrr::map(trellData, function(trell){
        
        if(is.null(panel_variable)){
          tester <- pmartR::get_edata_cname(trell)
        } else {
          tester <- panel_variable
        }
        
        name <- as.character(trell[[1]][floor(runif(1, 1, nrow(trell[[1]]))), tester])
        subset_dv <- trell$data_values[as.character(trell$data_values[[tester]]) == name,]
        subset_ss <- trell$summary_stats[as.character(trell$summary_stats[[tester]]) == name,]
        subset_cs <- trell$comp_stats[as.character(trell$comp_stats[[tester]]) == name,]
        
        subsets <- list(subset_dv, subset_ss, subset_cs)
        subsets <- subsets[as.numeric(which(!unlist(lapply(subsets, is.null))))]
        if(length(subsets) == 3){
          subset <- merge(merge(subsets[[1]], subsets[[2]]), subsets[[3]])
        } else if (length(subsets) == 2){
          subset <- merge(subsets[[1]], subsets[[2]])
        } else {
          subset <- subsets[[1]]
        }
        
        x <- eval(parse(text = paste0(custom_cog, "(subset)")))
          if(!is.data.frame(x) || nrow(x) != 1) stop(
            paste("Validation failed. Custom cog function on data subset is required to yield a data.frame with nrow == 1. Subset tested on:", tester, "==", name))
      })
      
    }
    
    # Validate Custom plot if provided
    if(!is.null(custom_plot)){
      purrr::map(trellData, function(trell){
        
        if(is.null(panel_variable)){
          tester <- pmartR::get_edata_cname(trell)
        } else {
          tester <- panel_variable
        }
        
        name <- as.character(trell[[1]][floor(runif(1, 1, nrow(trell[[1]]))), tester])
        subset_dv <- trell$data_values[as.character(trell$data_values[[tester]]) == name,]
        subset_ss <- trell$summary_stats[as.character(trell$summary_stats[[tester]]) == name,]
        subset_cs <- trell$comp_stats[as.character(trell$comp_stats[[tester]]) == name,]
        
        subsets <- list(subset_dv, subset_ss, subset_cs)
        subsets <- subsets[as.numeric(which(!unlist(lapply(subsets, is.null))))]
        if(length(subsets) == 3){
          subset <- merge(merge(subsets[[1]], subsets[[2]]), subsets[[3]])
        } else if (length(subsets) == 2){
          subset <- merge(subsets[[1]], subsets[[2]])
        } else {
          subset <- subsets[[1]]
        }
        
        x <- eval(parse(text = paste0(custom_plot, "(subset)")))
        print(length(x))
        if(!inherits(x, c("ggplot", "plotly", "rbokeh")))  stop(
          paste("Validation failed. Custom plot function on data subset is required to yield a single plot with class of ggplot, plotly, or rbokeh. Subset tested on:", tester, "==", name)
        )
        
      })
      
    }
    
    # Fill plot type
    if (is.null(plot_type) && is.null(custom_plot)){
      if(is.null(trellData[[1]]$comp_stats)){
        plot_type <- "abundance_boxplot"
      } else if (is.null(trellData[[1]]$data_values)){
        plot_type <- "foldchange_bar"
      } else {
        plot_type <- list("abundance_boxplot", "foldchange_bar")
      }
    }
    
    
    # Check for normalization of both items in list
    if (is.null(pmartR::get_data_norm(trellData[[1]])) || !pmartR::get_data_norm(trellData[[1]])){
      
      if("isobaricpepData" %in% get_data_class(trellData[[1]]) &&
         (is.null(attributes(trellData[[1]])$isobaric_info$norm_info$is_normalized) ||
          attributes(trellData[[1]])$isobaric_info$norm_info$is_normalized != TRUE)) stop(
            "Input must be normalized prior to plotting. Options include normalize_isobaric as appropriate for isobaricpepData, and normalize_global for any data type. (List input 1)")
      
      if(!("isobaricpepData" %in% get_data_class(trellData[[1]]))) stop(
        "Input must be normalized prior to plotting; please run normalize_global. (List input 1)")
      
    }
    
    if (is.null(pmartR::get_data_norm(trellData[[2]])) || !pmartR::get_data_norm(trellData[[2]])){
      
      if("isobaricpepData" %in% get_data_class(trellData[[2]]) &&
         (is.null(attributes(trellData[[2]])$isobaric_info$norm_info$is_normalized) ||
          attributes(trellData[[2]])$isobaric_info$norm_info$is_normalized != TRUE)) stop(
            "Input must be normalized prior to plotting. Options include normalize_isobaric as appropriate for isobaricpepData, and normalize_global for any data type. (List input 2)")
      
      if(!("isobaricpepData" %in% get_data_class(trellData[[2]]))) stop(
        "Input must be normalized prior to plotting; please run normalize_global. (List input 2)")
      
    }
    
    ## Check trelli_name ##
    # Default trelliscope names correspond to data types entered #
    if(!is.null(plot_type)){
      
      if(is.null(trelli_name)){
        trelli_name <- vector("list", 2)
        trelli_name[[1]] <- paste(attr(trellData, "data_types")[1], plot_type, sep = "_")
        trelli_name[[2]] <- paste(attr(trellData, "data_types")[2], plot_type, sep = "_")
        if(!is.null(custom_plot)){
          trelli_name[[1]] <- c(trelli_name[[1]], paste("Custom", attr(trellData, "data_types")[1], sep = "_"))
          trelli_name[[2]] <- c(trelli_name[[2]], paste("Custom", attr(trellData, "data_types")[2], sep = "_"))
        }
        names(trelli_name) <- attr(trellData, "data_types")
        
        
        # if trelli_name is not of length list, add identifiers #
      } else if (length(trelli_name) == 1){
        label <- trelli_name
        trelli_name <- vector("list", 2)
        trelli_name[[1]] <- paste(paste(label, attr(trellData, "data_types")[1], sep = "_"), plot_type, sep = "_")
        trelli_name[[2]] <- paste(paste(label, attr(trellData, "data_types")[2], sep = "_"), plot_type, sep = "_")
        if(!is.null(custom_plot)){
          trelli_name[[1]] <- c(trelli_name[[1]], paste("Custom", attr(trellData, "data_types")[1], sep = "_"))
          trelli_name[[2]] <- c(trelli_name[[2]], paste("Custom", attr(trellData, "data_types")[2], sep = "_"))
        }
        names(trelli_name) <- attr(trellData, "data_types")
        
        # if trelli_name too long, cut to length(trellData)  #
      } else if (((length(trelli_name) != length(plot_type)*2 + 2) && !is.null(custom_plot)) || 
                 (length(trelli_name) != length(plot_type)*2  && is.null(custom_plot)) ) {
        message("Length of trelli_name is not equal to the number of displays generated. Only the first name will be used.")
        label <- trelli_name[1]
        trelli_name <- vector("list", 2)
        trelli_name[[1]] <- paste(label, attr(trellData, "data_types")[1], plot_type, sep = "_")
        trelli_name[[2]] <- paste(label, attr(trellData, "data_types")[2], plot_type, sep = "_")
        if(!is.null(custom_plot)){
          trelli_name[[1]] <- c(trelli_name[[1]], paste("Custom", attr(trellData, "data_types")[1], sep = "_"))
          trelli_name[[2]] <- c(trelli_name[[2]], paste("Custom", attr(trellData, "data_types")[2], sep = "_"))
        }
        names(trelli_name) <- attr(trellData, "data_types")
      } else {
        label <- trelli_name
        trelli_name <- vector("list", 2)
        trelli_name[[1]] <- trelli_name[1:length(plot_type)]
        trelli_name[[2]] <- trelli_name[(length(plot_type)+1):(length(plot_type)*2)]
        names(trelli_name) <- attr(trellData, "data_types")
      }
      
    } else {
      
      if(is.null(trelli_name)){
        trelli_name <- vector("list", 2)
        trelli_name[[1]] <- paste("Custom", paste(attr(trellData, "data_types")[1], collapse = "_"), sep = "_")
        trelli_name[[2]] <- paste("Custom", paste(attr(trellData, "data_types")[2], collapse = "_"), sep = "_")
        names(trelli_name) <- c(paste(attr(trellData, "data_types")[1], collapse = "_"), 
                                paste(attr(trellData, "data_types")[2], collapse = "_"))
        
        # if trelli_name is not of length list, add identifiers #
      } else if (length(trelli_name) == 1){
        label <- trelli_name
        trelli_name <- vector("list", 2)
        trelli_name[[1]] <- c(trelli_name[[1]], paste(label, attr(trellData, "data_types")[1], sep = "_"))
        trelli_name[[2]] <- c(trelli_name[[2]], paste(label, attr(trellData, "data_types")[2], sep = "_"))
        names(trelli_name) <- attr(trellData, "data_types")
        
        # if trelli_name too long, cut to length(trellData)  #
      } else if (length(trelli_name) != 2) {
        message("Length of trelli_name is not equal to the number of displays generated. Only the first name will be used.")
        label <- trelli_name[1]
        trelli_name <- vector("list", 2)
        trelli_name[[1]] <- c(trelli_name[[1]], paste(label, attr(trellData, "data_types")[1], sep = "_"))
        trelli_name[[2]] <- c(trelli_name[[2]], paste(label, attr(trellData, "data_types")[2], sep = "_"))
        names(trelli_name) <- attr(trellData, "data_types")
      } else {
        label <- trelli_name
        trelli_name <- vector("list", 2)
        trelli_name[[1]] <- trelli_name[1]
        trelli_name[[2]] <- trelli_name[2]
        names(trelli_name) <- attr(trellData, "data_types")
      }
      
    }
    
    # Ensure NULLs or the correct length #
    if (is.null(omicsData)){
      omicsData <- rep(list(NULL), length(trellData))
    }
    if (is.null(omicsStats)){
      omicsStats <- rep(list(NULL), length(trellData))
    }
    if (is.null(panel_variable)){
      panel_variable <- rep(list(NULL), length(trellData))
    } else {
      if(!(all(panel_variable %in% c(colnames(trellData[[1]]$summary_stats),
                               colnames(trellData[[1]]$comp_stats),
                               colnames(trellData[[1]]$data_values))))) stop(
                                 paste("Panel_variable input is not present in input data columns (input list item 1). Panel_variable =", toString(panel_variable))
                               )
      if(!(all(panel_variable %in% c(colnames(trellData[[2]]$summary_stats),
                                 colnames(trellData[[2]]$comp_stats),
                                 colnames(trellData[[2]]$data_values))))) stop(
                                   paste("Panel_variable input is not present in input data columns (input list item 2). Panel_variable =", toString(panel_variable))
                                 )
    }
    
    if(!is.null(try_URL) && !(try_URL %in% c(colnames(trellData[[1]]$summary_stats),
                                   colnames(trellData[[1]]$comp_stats),
                                   colnames(trellData[[1]]$data_values)))) stop(
                                     paste("try_URL input is not present in input data columns (input list item 1). try_URL =", toString(try_URL))
                                   )
    if(!is.null(try_URL) && !(try_URL %in% c(colnames(trellData[[2]]$summary_stats),
                                   colnames(trellData[[2]]$comp_stats),
                                   colnames(trellData[[2]]$data_values)))) stop(
                                     paste("try_URL input is not present in input data columns (input list item 2). try_URL =", toString(try_URL))
                                   )
    
    
    # Nest data and generate trelliscope plots #
    
    tictoc::tic("Generate plots")
    
    group <- purrr::map(trellData, function(item) paste(pmartR:::get_data_class(item), collapse = "_"))
    
      nested_plot <- purrr::map2(trellData, panel_variable, function(pairedplotter, pan){
        
        if(!is.null(plot_type)){
          nest_out <- purrr::map(plot_type, function(types){
            format_plot(trellData = pairedplotter, 
                        p_val = p_val,
                        panel_variable = pan, 
                        plot_type = types,
                        plot_text = plot_text,
                        y_limits = y_limits,
                        interactive = interactive)
          })
        
          if(!is.null(custom_plot)){
            nest_out[[length(nest_out) + 1]] <- format_plot(trellData = pairedplotter, 
                                                            p_val = p_val,
                                                            panel_variable = pan, 
                                                            custom_plot = custom_plot,
                                                            plot_text = plot_text,
                                                            y_limits = y_limits,
                                                            interactive = interactive)
          }
        } else {
          nest_out <- list()
          nest_out[[1]] <- format_plot(trellData = pairedplotter,
                                       p_val = p_val,
                                       panel_variable = pan,
                                       custom_plot = custom_plot,
                                       plot_text = plot_text,
                                       y_limits = y_limits,
                                       interactive = interactive)
        }
        return(nest_out)
      })
    
    tictoc::toc()
    
    tictoc::tic("Generate cogs") 
    
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
          try_URL = try_URL,
          custom_cog = custom_cog))
        
        pan <- colnames(cogs)[1]
        
        cog_out <- purrr::map2(nests, links, function(nest, link){
          cogs2 <- dplyr::mutate(nest, cogs = cogs$cogs)
          
          if ("pepData" %in% pmartR:::get_data_class(trell)){
            
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
    if(display) tictoc::tic("Pipe into trelliscope")

    out <- purrr::pmap(list(nest_plot_cog_list, 
                trelli_name, group), function(display1, name1, grp){
                  purrr::map2(display1, name1, function(display2, name){
                    
                    
                    if(!display){
                      list(display_object = display2, 
                                  name = as.character(name), 
                                  self_contained = self_contained, 
                                  group = grp, 
                                  thumb = TRUE, 
                                  state = list(
                                    sort = list(trelliscopejs::sort_spec(names(display2[1]), dir = "desc")),
                                    labels = list(names(display2[1]))),
                                  ...)
                    } else {
                    
                    if(self_contained){
                      
                        trelliscopejs::trelliscope(display2, as.character(name), nrow = 1, ncol = 2,
                                                   self_contained = TRUE, 
                                                   group = grp,
                                                   thumb = TRUE, state = list(
                                                     sort = list(trelliscopejs::sort_spec(names(display2[1]), dir = "desc")), 
                                                     labels = list(names(display2[1]))), 
                                                   ...)
                        
                      } else {
                        
                        trelliscopejs::trelliscope(display2, as.character(name), nrow = 1, ncol = 2,
                                                   path = as.character(trelli_path_out), 
                                                   group = grp,
                                                   thumb = TRUE, state = list(
                                                     sort = list(trelliscopejs::sort_spec(names(display2[1]), dir = "desc")), 
                                                     labels = list(names(display2[1]))), ...)
                        
                      }
                    }
                  })
      })
    
    if(!display){
      return(out)
    } else {
      tictoc::toc()
      return(out[[1]][[1]])
    }
    
  # Where not a pep/pro pair: #
  } else {
    
    group <- paste(pmartR:::get_data_class(trellData), collapse = "_")
    
    #Check panel variable
    if(!is.null(panel_variable) && !(panel_variable %in% c(colnames(trellData$summary_stats),
                               colnames(trellData$comp_stats),
                               colnames(trellData$data_values)))) stop(
                                 paste("Panel_variable input is not present in input data columns.  Panel_variable =", toString(panel_variable))
                               )
    
    #Check try_url
    if(!is.null(try_URL) && !(try_URL %in% c(colnames(trellData$summary_stats),
                                                           colnames(trellData$comp_stats),
                                                           colnames(trellData$data_values)))) stop(
                                                             paste("try_URL input is not present in input data columns.  try_URL =", toString(try_URL))
                                                           )

    # Check custom cog
    if(!is.null(custom_cog)){
        
      if(is.null(panel_variable)){
        tester <- pmartR::get_edata_cname(trellData)
      } else {
        tester <- panel_variable
      }
      
      name <- as.character(trellData[[1]][floor(runif(1, 1, nrow(trellData[[1]]))), tester])
      subset_dv <- trellData$data_values[as.character(trellData$data_values[[tester]]) == name,]
      subset_ss <- trellData$summary_stats[as.character(trellData$summary_stats[[tester]]) == name,]
      subset_cs <- trellData$comp_stats[as.character(trellData$comp_stats[[tester]]) == name,]
      
      subsets <- list(subset_dv, subset_ss, subset_cs)
      subsets <- subsets[as.numeric(which(!unlist(lapply(subsets, is.null))))]
      if(length(subsets) == 3){
        subset <- merge(merge(subsets[[1]], subsets[[2]]), subsets[[3]])
      } else if (length(subsets) == 2){
        subset <- merge(subsets[[1]], subsets[[2]])
      } else {
        subset <- subsets[[1]]
      }
      
      x <- eval(parse(text = paste0(custom_cog, "(subset)")))
      if(!is.data.frame(x) || nrow(x) != 1) stop(
        paste("Validation failed. Custom cog function on data subset is required to yield a data.frame with nrow == 1. Subset tested on:", tester, "==", name))
    
    }
    
    # Validate Custom plot if provided
    if(!is.null(custom_plot)){
        
      if(is.null(panel_variable)){
        tester <- pmartR::get_edata_cname(trellData)
      } else {
        tester <- panel_variable
      }
      
      name <- as.character(trellData[[1]][floor(runif(1, 1, nrow(trellData[[1]]))), tester])
      subset_dv <- trellData$data_values[as.character(trellData$data_values[[tester]]) == name,]
      subset_ss <- trellData$summary_stats[as.character(trellData$summary_stats[[tester]]) == name,]
      subset_cs <- trellData$comp_stats[as.character(trellData$comp_stats[[tester]]) == name,]
      
      subsets <- list(subset_dv, subset_ss, subset_cs)
      subsets <- subsets[as.numeric(which(!unlist(lapply(subsets, is.null))))]
      if(length(subsets) == 3){
        subset <- merge(merge(subsets[[1]], subsets[[2]]), subsets[[3]])
      } else if (length(subsets) == 2){
        subset <- merge(subsets[[1]], subsets[[2]])
      } else {
        subset <- subsets[[1]]
      }
      
      x <- eval(parse(text = paste0(custom_plot, "(subset)")))
      if(!inherits(x, c("ggplot", "plotly", "rbokeh")))  stop(
        paste("Validation failed. Custom plot function on data subset is required to yield a single plot with class of ggplot, plotly, or rbokeh. Subset tested on:", tester, "==", name)
      )
    }
    
    # Fill plot type
    if (is.null(plot_type) && is.null(custom_plot)){
      if(is.null(trellData$comp_stats)){
        plot_type <- "abundance_boxplot"
      } else if (is.null(trellData$data_values)){
        plot_type <- "foldchange_bar"
      } else {
        plot_type <- list("abundance_boxplot", "foldchange_bar")
      }
    }
    
    if (is.null(pmartR::get_data_norm(trellData)) || !pmartR::get_data_norm(trellData)){
      
      if("isobaricpepData" %in% get_data_class(trellData) &&
         (is.null(attributes(trellData)$isobaric_info$norm_info$is_normalized) ||
         attributes(trellData)$isobaric_info$norm_info$is_normalized != TRUE)) stop(
           "Input must be normalized prior to plotting. Options include normalize_isobaric as appropriate for isobaricpepData, and normalize_global for any data type.")
      
      if(!("isobaricpepData" %in% get_data_class(trellData))) stop(
        "Input must be normalized prior to plotting; please run normalize_global.")
      
    }
    
    # Default: Name the display after parent classes (omicsData and omicStats clases) #
    
    if(!is.null(plot_type)){
      if(is.null(trelli_name)){
        trelli_name <- paste(paste(attr(trellData, "parent_class"), collapse = "_"), plot_type, sep = "_")
        # Make sure trelliname is length 1 #
      } else if (length(trelli_name) == 1){
        label <- trelli_name
        trelli_name <- paste(label, plot_type, sep = "_")
        if(!is.null(custom_plot)){
          trelli_name <- c(trelli_name, paste(label, "Custom", sep = "_"))
        }
      } else if (((length(trelli_name) != (length(plot_type))) && is.null(custom_plot)) || 
                 (!is.null(custom_plot) && length(trelli_name) != (length(plot_type)) + 1)) {
        warning("trelli_name length is not equal to the number of displays generated from plot_type; using the first entry in trelli_name.")
        label <- trelli_name[[1]]
        trelli_name <- paste(label, plot_type, sep = "_")
        if(!is.null(custom_plot)){
          trelli_name <- c(trelli_name, paste(label, "Custom", sep = "_"))
        }
      }
    } else {
      if(is.null(trelli_name)){
        trelli_name <- paste(paste(attr(trellData, "parent_class"), collapse = "_"), "Custom", sep = "_")
        # Make sure trelliname is length 1 #
      } else if (length(trelli_name) != 1) {
        warning("trelli_name length is not equal to the number of displays generated from plot_type; using the first entry in trelli_name.")
        trelli_name <- paste(trelli_name[[1]], "Custom", sep = "_")
      }
    }
    
    # Nest data and generate trelliscope plots #
    
    if(!is.null(plot_type)){
      
      displays <- purrr::map2(plot_type, trelli_name[1:length(plot_type)], function(types, names){
        
        tictoc::tic("Generate plots")
        
        nested_plot <- format_plot(trellData = trellData, 
                                   p_val = p_val,
                                   plot_type = types,
                                   plot_text = plot_text,
                                   y_limits = y_limits,
                                   panel_variable = panel_variable,
                                   interactive = interactive)
        tictoc::toc()
        
        tictoc::tic("Generate cogs")
        
        # Generate default cognostics #
        nest_plot_cog <- suppressWarnings(data_cogs(nested_plot = nested_plot, 
                                                    trellData = trellData,
                                                    p_val = p_val, 
                                                    try_URL = try_URL,
                                                    custom_cog = custom_cog))
        tictoc::toc()
        
        # Generate trelliscope display #
        
        if(!display){
          out <- list(display_object = nest_plot_cog, 
                      name = as.character(names), 
                      self_contained = self_contained, 
                      group = group, 
                      thumb = TRUE, 
                      state = list(
                        sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")),
                        labels = list(names(nest_plot_cog[1]))
                      ),
                      ...)
          return(out)
        }
        
        tictoc::tic("Pipe into trelliscope")
        
        if(self_contained){
          out <- nest_plot_cog %>%
            trelliscopejs::trelliscope(name = as.character(names), nrow = 1, ncol = 2,
                                       self_contained = TRUE, 
                                       thumb = TRUE,
                                       state = list(
                                         sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")), 
                                         labels = list(names(nest_plot_cog[1]))),
                                       ...
            )
        } else {
          
          
          if(!display){
            out <- list(display_object = nest_plot_cog, 
                        name = as.character(names), 
                        self_contained = self_contained, 
                        group = group, 
                        thumb = TRUE, 
                        state = list(
                          sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")),
                          labels = list(names(nest_plot_cog[1]))
                        ),
                        ...)
            return(out)
          }
          
          
          out <- nest_plot_cog %>%
            trelliscopejs::trelliscope(name = as.character(names), nrow = 1, ncol = 2,
                                       path = as.character(trelli_path_out), 
                                       thumb = TRUE,
                                       group = group,
                                       state = list(
                                         sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")),
                                         labels = list(names(nest_plot_cog[1]))
                                         ),
                                       ...
            )
        }
        tictoc::toc()
        return(out)
      })
      
      if(!is.null(custom_plot)){
        
          tictoc::tic("Generate plots")
          
          nested_plot <- format_plot(trellData = trellData, 
                                     p_val = p_val,
                                     custom_plot = custom_plot,
                                     plot_text = plot_text,
                                     y_limits = y_limits,
                                     panel_variable = panel_variable,
                                     interactive = interactive)
          tictoc::toc()
          
          tictoc::tic("Generate cogs")
          
          # Generate default cognostics #
          nest_plot_cog <- suppressWarnings(data_cogs(nested_plot = nested_plot, 
                                                      trellData = trellData,
                                                      p_val = p_val, 
                                                      try_URL = try_URL,
                                                      custom_cog = custom_cog))
          tictoc::toc()
          
          # Generate trelliscope display #
          
          if(!display){
            out <- list(display_object = nest_plot_cog, 
                        name = as.character(trelli_name[length(trelli_name)]), 
                        self_contained = self_contained, 
                        group = group, 
                        thumb = TRUE, 
                        state = list(
                          sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")), 
                          labels = list(names(nest_plot_cog[1]))),
                        ...)
          } else {
            
            tictoc::tic("Pipe into trelliscope")
          
            if(self_contained){
              out <- nest_plot_cog %>%
                trelliscopejs::trelliscope(name = as.character(trelli_name[length(trelli_name)]), nrow = 1, ncol = 2,
                                           self_contained = TRUE, 
                                           thumb = TRUE,
                                           state = list(
                                             sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")), 
                                             labels = list(names(nest_plot_cog[1]))),
                                           ...
                )
            } else {
              
              out <- nest_plot_cog %>%
                trelliscopejs::trelliscope(name = as.character(trelli_name[length(trelli_name)]), nrow = 1, ncol = 2,
                                           path = as.character(trelli_path_out), 
                                           thumb = TRUE,
                                           group = group,
                                           state = list(
                                             sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")), 
                                             labels = list(names(nest_plot_cog[1]))),
                                           ...
                )
            }
          }
          
          displays[[length(displays) + 1]] <- out
          
          if(!display) tictoc::toc()
        }
        
    } else {
      tictoc::tic("Generate plots")

      nested_plot <- format_plot(trellData = trellData, 
                                 p_val = p_val,
                                 custom_plot = custom_plot,
                                 plot_text = plot_text,
                                 y_limits = y_limits,
                                 panel_variable = panel_variable,
                                 interactive = interactive)
      tictoc::toc()
      
      tictoc::tic("Generate cogs")
      
      # Generate default cognostics #
      nest_plot_cog <- suppressWarnings(data_cogs(nested_plot = nested_plot, 
                                                  trellData = trellData,
                                                  p_val = p_val, 
                                                  try_URL = try_URL,
                                                  custom_cog = custom_cog))
      tictoc::toc()
      
      # Generate trelliscope display #
      
      if(!display){
        out <- list(display_object = nest_plot_cog, 
                    name = as.character(trelli_name[length(trelli_name)]), 
                    self_contained = self_contained, 
                    group = group, 
                    thumb = TRUE, 
                    state = list(
                      sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")), 
                      labels = list(names(nest_plot_cog[1]))),
                    ...)
      } else {
      
        tictoc::tic("Pipe into trelliscope")
        
        if(self_contained){
          out <- nest_plot_cog %>%
            trelliscopejs::trelliscope(name = as.character(trelli_name), nrow = 1, ncol = 2,
                                       self_contained = TRUE, 
                                       thumb = TRUE,
                                       state = list(
                                         sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")), 
                                         labels = list(names(nest_plot_cog[1]))),
                                       ...
            )
        } else {
          
          out <- nest_plot_cog %>%
            trelliscopejs::trelliscope(name = as.character(trelli_name), nrow = 1, ncol = 2,
                                       path = as.character(trelli_path_out), 
                                       thumb = TRUE,
                                       group = group,
                                       state = list(
                                         sort = list(trelliscopejs::sort_spec(names(nest_plot_cog[1]), dir = "desc")), 
                                         labels = list(names(nest_plot_cog[1]))),
                                       ...
            )
        }
      }
      
      displays <- out
      if(display) tictoc::toc()
      
    }
    
    return(displays)
  }
}
