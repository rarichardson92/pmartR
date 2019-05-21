#' @name as.omicsPlotter
#' @rdname as.omicsPlotter
#' @title Convert Omics data and pairwise statistics to a plotting object
#'
#' @description Converts a ResObject and its respective OmicsData into an easily plottable object.
#'
#' @param omicsData A pmartR object of class pepData, lipidData, metabData, or proData
#' @param omicsStats A statistical results object produced by running \code{imd_anova} on omicsData.
#' @param ... further arguments
#'
#' @details Objects of class 'omicsPlotter' inherit attributes from omicsStats and omicsData. These attributes can be changed from their default value by manual specification. A list of these attributes as well as their default values are as follows:
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
as.omicsPlotter <- function(...){
  .as.omicsPlotter(...)
}

.as.omicsPlotter <- function(omicsData = NULL, omicsStats = NULL, check.names = TRUE){
  
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
  
  ## Initial Checks ##
  # Check that omicsData and omicsStats are the correct classes #
  if (!is.null(omicsData) & !inherits(omicsData, c("proData", "pepData", "metabData", "lipidData"))) stop("omicsData must be of class 'proData', 'pepData', 'metabData', or 'lipidData'")
  if(!is.null(omicsStats) & !inherits(omicsStats, "statRes")) stop("omicsStats must be of the class 'statRes'")
  # Check if omicsData and omicsStats have the same cname attributes. #
  if(!is.null(omicsStats) & !is.null(omicsData)){
    # Check if omicsData and omicsStats have the same cname attributes. #
    if (!(identical(attributes(omicsData)$cnames, attributes(omicsStats)$cnames))) stop("Non-matching cname attributes in omicsStats and omicsData. Check that omicsStats is correctly derived from omicsData.")
    # Check that the biomolecule unique ID column exists in omicsData and omicsStats #
    if(!(uniqedatId %in% names(omicsStats$Full_results))) stop(paste("Column ", uniqedatId," not found in omicsStats. Requires compatible identifiers.", sep = ""))
    # Check if stats biomolecules are (at least) a subset of omicsData biomolecules. #
    if(!all(omicsStats$Full_results[[uniqedatId]] %in% omicsData$e_data[[uniqedatId]])) stop(paste("Biomolecules in omicsStats do not match biomolecules in omicsData.", sep = ""))
    # Check if stats sampleIDs are (at least) a subset of omicsData f_data sample IDs. #
    if(!all(attributes(omicsStats)$group_DF[[sampID]] %in% omicsData$f_data[[sampID]])) stop(paste(sampID, "column does not match between omicsData and omicsStats.", sep = ""))
  } 
  # Make sure at least one of omicsData or omicsStats is present #
  if(is.null(omicsStats) & is.null(omicsData)) stop ("as.omicsPlotter() requires at least one of the following: omicsStats, omicsData")
  
  if (!is.null(omicsData)){
    ## Manipulate Dataframes ##
    # Original data values from omicsData #
    # 1) Melt to assign e_data values to each sample ID
    # 2) Combine with groups in f_data
    # 3) Save attributes to be used
    data_values <- suppressWarnings(
      reshape2::melt(omicsData$e_data, 
                     id.vars = uniqedatId,
                     value.name = paste(attr(omicsData, "data_info")$data_scale,
                                        "_abundance", sep = ""),
                     variable.name = sampID) %>%
        dplyr::left_join(omicsData$f_data, by = sampID))
    
    inheritDatAtt <- attributes(omicsData)[grep("filters|^names|group_DF", 
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
        attributes(omicsStats)$comparisons %>%
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
      paste(attributes(omicsStats)$comparisons, collapse = '|')
      )]
    
    # Re-sort for any grouping variables #
    # 1) If there is a group_DF specified,
    # 2) Extract unique groups in group_DF, and for each group
    #  3) Bind extracted with unique IDs from omicsStats and a Group column
    #  4) Rename columns by removing group designation
    # 5) Bind rows from each group into one dataframe 
    if (!is.null(attributes(omicsStats)$group_DF)){
      groups <-  unique(attributes(omicsStats)$group_DF$Group)
      summary_stats <- suppressWarnings(
        dplyr::bind_rows(purrr::map(groups, function(group){
          df <- cbind(summary_stats[[uniqedatId]],
                      rep(group, nrow(summary_stats)),
                      summary_stats[,stringr::str_detect(
                        names(summary_stats), 
                        paste(group, "$", sep = ""))])
          trimname <- stringr::str_remove_all(colnames(df)[3:ncol(df)],
                                              paste("_", group, sep = ""))
          colnames(df) <- c(uniqedatId, "Group", trimname)
          return(df)
          }
        )))
    }
    
    # Save attributes to be used #
    inheritStatAtt <- attributes(omicsStats)[grep("^names", 
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
  
  # inheritDatAtt <- attributes(omicsData)[grep("filters|^names|group_DF", 
  #                             names(attributes(omicsData)), 
  #                             invert = T)]
  # inheritStatAtt <- attributes(omicsStats)[grep("^names", 
  #                                             names(attributes(omicsStats)), 
  #                                             invert = T)]
  
  attributes(res) <- c(attributes(res), inheritDatAtt, inheritStatAtt)
  class(res) <- "omicsPlotter" 
  
  if (!is.null(omicsData) & !is.null(omicsStats)){
    attr(res, "parent_class") <- c(class(omicsStats), class(omicsData))
  } else if (!is.null(omicsData)){
    attr(res, "parent_class") <- class(omicsData)
  } else {
    attr(res, "parent_class") <- class(omicsStats)
  }

  return(res)
}


#' print.omicsPlotter
#' 
#' For printing an S3 object of type 'omicsPlotter':
#' 
#'@rdname print-omicsPlotter
#'@export
#'
print.omicsPlotter<- function(omicsPlotter){
  if(!inherits(omicsPlotter, "omicsPlotter")) stop("omicsPlotter object must be of the class 'omicsPlotter'")
  
  data_values <- as.data.frame(lapply(omicsPlotter$data_values, as.character), stringsAsFactors = FALSE, check.names = attr(omicsPlotter, "check.names"))
  summary_stats <- as.data.frame(lapply(omicsPlotter$summary_stats, as.character), stringsAsFactors = FALSE, check.names = attr(omicsPlotter, "check.names"))
  comp_stats <- as.data.frame(lapply(omicsPlotter$comp_stats, as.character), stringsAsFactors = FALSE, check.names = attr(omicsPlotter, "check.names"))
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

#' @name set_increment
#' @rdname set_increment
#' @title Sets y-axis increment for Omicsplotter labels plotting
#' 
#' @description Sets y-axis increment for Omicsplotter plotting. Used by plot_comp.
#'
#' @param yvalues y-values for plotting
#' @param testtype consideration for different statistical tests
#'
#' @author Rachel Richardson
#'
#' @export
set_increment <- function(yvalues, testtype){
  .set_increment(yvalues, testtype)
}

.set_increment <- function(yvalues, testtype){
  ## Initial Checks and Replacement ##
  # Catch NAs, replace with 0 #
  if (any(is.na(yvalues))){
    yvalues[is.na(yvalues)] <- 0
  }
  # Check for numeric data #
  if(!inherits(yvalues, "numeric")) stop("yvalues must be of the class 'numeric'")
  # Check if a vector #
  if(!is.vector(yvalues)) stop("yvalues must be a numeric vector")
  # Check if test type is correct #
  if(!(testtype %in% c("anova", "gtest", "combined"))) stop("testtype must be 'anova', 'gtest', or 'combined'")

  ## Set increment based on test type ##
  if (testtype != "gtest"){
    if(length(yvalues)==1){
      increment <- abs(yvalues)/20
    } else {
      increment <- (max(yvalues) - min(yvalues))/20
    }
  } else {
    increment <- max(yvalues)/20
  }
  return(increment)
}

#' @name set_ylimits
#' @rdname set_ylimits
#' @title Sets y-axis limits for Omicsplotter plotting
#' 
#' @description Sets y-axis limits for Omicsplotter plotting. Used by plot_comp.
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
                         y_max = NULL, y_min = NULL){
  
  ## Initial Checks and Replacement ##
  # Check for numeric data #
  if(!inherits(yvalues, "numeric")) stop("yvalues must be of the class 'numeric'")
  # Check if a vector #
  if(!is.vector(yvalues)) stop("yvalues must be a numeric vector")
  # Catch NAs, replace with 0 #
  if (any(is.na(yvalues))){
    yvalues[is.na(yvalues)] <- 0
  }
  
  ## Set Limits ##
  # Catch for pre-set y-limits and y_range #
  if (!is.null(y_min) & !is.null(y_max)){
    maxi <- y_max
    mini <- y_min
    return(c(mini, maxi))
  } else if( !is.null(y_range) & !is.null(y_min)){
    mini <- y_min
    maxi <- y_min + y_range
    return(c(mini, maxi))
  } else if( !is.null(y_range) & !is.null(y_max)){
    maxi <- y_max
    mini <- y_max - y_range
    return(c(mini, maxi))
  } else if (!is.null(y_range)){
    maxi <- median(yvalues) + y_range/2
    mini <- median(yvalues) - y_range/2
    return(c(mini, maxi))
  }
  
  # Set maxi and mini based on difference between max and min values #
  maxi <- max(yvalues) + 4*increment
  mini <- min(yvalues) - 4*increment
  
  # Adjust for maximums below 0 and minimums above 0 #
  if (max(yvalues) < 0){
    maxi <- 4*increment
  }
  if (min(yvalues) > 0){
    mini <- -4*increment
  }
  
  # Adjust for specified y_min and y_max #
  if (!is.null(y_min)){
    mani <- y_max
  }
  if (!is.null(y_min)){
    mini <- y_min
  }
  
  ## Return limits ##
  # Minimum y value = mini, maximum y value = maxi #
  return(c(mini, maxi))
}


#' @name plot_comps
#' @rdname plot_comps
#' @title Plot pairwise comparisons of omicsPlotter object
#'
#' @description Plot pairwise comparisons of omicsPlotter. Plots differently depending on statistical tests ran, located in attributes.
#'
#' @param omicsPlotter An object of class "omicsPlotter" generated from \code{\link{as.omicsPlotter}}.
#' @param y_limits Set to "fixed" or "free" for automated y-axis calculating. "fixed" - axis generated based on the maximum/minimum across all plots. "free" - axis axis generated based on the maximum/minimum of individual plot.
#' @param y_range Specify a range for the plot y-axis. Will calculate the range based on one of y_max or y_min parameters or from the median of y-values where y_max and y_min are not defined.
#' @param y_max Sets the maximum y-value for the y-axis.
#' @param y_min Sets the minimum y-value for the y-axis.
#' @param p_val Specifies p-value for setting graph border colors
#' @param panel_variable Specifies what to divide trelliscope panels by, must be a column in omicsPlotter. Defaults to cnames$edata_cname of omicsPlotter.
#' @param panel_x_axis Specifies what column should be on the x-axis, must be a column in omicsPlotter. Defaults setting plots pairwise comparisons along x-axis.
#' @param panel_y_axis Specifies what column should be on the y-axis, must be a column in omicsPlotter and numeric. Defaults setting plots fold change for combined and anova testing and counts for g-test.
#'
#' @author Rachel Richardson
#' @export
plot_comps <- function(omicsPlotter, ...) {
   .plot_comps(omicsPlotter,  ...)
}

.plot_comps <- function(omicsPlotter, 
                        y_limits = NULL, y_range = NULL, 
                        y_max = NULL, y_min = NULL, p_val = 0.05, 
                        panel_variable = NULL, 
                        panel_x_axis = NULL, panel_y_axis = NULL ) {

  ## Initial Checks ##
  # Check if class is correct #
  if(!inherits(omicsPlotter, "omicsPlotter")) stop("omicsPlotter must be of the class 'omicsPlotter'")
  # Check if comp_stats is in omicsPlotter #
  if(is.null(omicsPlotter$comp_stats)) stop("No comparisons in omicsPlotter to plot")
  # Check if stats statistical test attribute is valid #
  if(!(attributes(omicsPlotter)$statistical_test %in% c("combined", "gtest", "anova"))) stop(paste("Non-applicable statistical_test attribute in omicsPlotter object."))
  # Check check if p_val is numeric of length 1 #
  if(!is.numeric(p_val) | (length(p_val) != 1)) stop("p_val must be a numeric of length 1")  
  
  # Check if y_limits or y_range have been selected correctly #
  if(!is.null(y_limits) & !is.null(y_range)) stop("Input either y_limits or y_range parameters, but not both.")
  # Check if only one of y_max and y_min has been selected with y-limits or y_range #
  if(!is.null(y_max) & !is.null(y_min) & (!is.null(y_range) | !is.null(y_limits))) stop("Cannot use both y_min and y_max with y_range or y_limits parameters. Only one of y_min or y_max can be used.")
  
  # Check if y_limits is in acceptable strings and length == 1 #
  if (!is.null(y_limits)){
    if(!(y_limits %in% c("fixed", "free")) )stop("Parameter y_limits must be input as either 'fixed' or 'free'.")
    if(length(y_limits) != 1) stop("Parameter y_limits must have length = 1.")
  }
  # Check if y_range is numeric and length == 1 #
  if (!is.null(y_range)){
    if(!is.numeric(y_range)) stop("Parameter y_range must be numeric.")
    if(length(y_range) != 1) stop("Parameter y_range must have length = 1.")
  }
  # Check if y_max is numeric and length == 1 #
  if (!is.null(y_max)){
    if(!is.numeric(y_max)) stop("Parameter y_max must be numeric.")
    if(length(y_max) != 1) stop("Parameter y_max must have length = 1.")
  }
  # Check if y_min is numeric and length == 1 #
  if (!is.null(y_min)){
    if(!is.numeric(y_min)) stop("Parameter y_min must be numeric.")
    if(length(y_min) != 1) stop("Parameter y_min must have length = 1.")
  }
  
  ## Moniker Variables ##
  option <- attributes(omicsPlotter)$statistical_test
  colors <- c("NA", "darkgrey", "black")
  if(is.null(panel_variable)){
    panel_variable <- attributes(omicsPlotter)$cnames$edata_cname
  } else {
    panel_variable <- panel_variable
  }
  if(is.null(panel_x_axis)){
    panel_x_axis <- "Comparison"
  } else {
    panel_x_axis <- panel_x_axis
  }
  if(is.null(panel_y_axis) & (option == "gtest")){
    panel_y_axis <- "Count"
  } else if (is.null(panel_y_axis)){
    panel_y_axis <- "Fold_change"
  } else {
    panel_y_axis <- panel_y_axis
  }
  
  ## Check Monikers ##
  # Ensure panel, x/y parameters are not matching #
  if (!((panel_x_axis != panel_variable) & (panel_y_axis != panel_variable))){
    stop("Parameter panel_y_axis, panel_x_axis cannot match panel_variable. Refer to ?plot_comps for default settings or try setting all of these parameters individually.")
  }
  if(panel_x_axis == panel_y_axis){
    print("Parameter panel_y_axis and panel_x_axis are identical. Refer to ?plot_comps for default settings or try setting parameters individually for different axis labels.")
  }
  # Ensure panel,x, and y parameters are in comp_stats, include summary stats if option == gtest #
  if (!(panel_x_axis %in% colnames(omicsPlotter$comp_stats))){
    if((option == "gtest") & (!(panel_x_axis %in% colnames(omicsPlotter$summary_stats)))){
      stop("Parameter panel_x_axis must be in column names of omicsPlotter comp_stats or summary_stats.")
    } else if (option != "gtest"){
      stop("Parameter panel_x_axis must be in column names of omicsPlotter comp_stats.")
    }
  }
  if (!(panel_y_axis %in% colnames(omicsPlotter$comp_stats))){
    if((option == "gtest") & (!(panel_y_axis %in% colnames(omicsPlotter$summary_stats)))){
      stop("Parameter panel_y_axis must be in column names of omicsPlotter comp_stats or summary_stats.")
    } else if (option != "gtest"){
      stop("Parameter panel_y_axis must be in column names of omicsPlotter comp_stats.")
    }
  }
  if (!(panel_variable %in% colnames(omicsPlotter$comp_stats))){
    if((option == "gtest") & (!(panel_variable %in% colnames(omicsPlotter$summary_stats)))){
      stop("Parameter panel_variable must be in column names of omicsPlotter comp_stats or summary_stats.")
    } else if (option != "gtest"){
      stop("Parameter panel_variable must be in column names of omicsPlotter comp_stats.")
    }
  }
  
  # Ensure summary_stats are present for gtest plotting #
  if(option == "gtest" & is.null(omicsPlotter$summary_stats)) stop("G test plotting requires summary statistics in omicsPlotter. Accesed via omicsPlotter$summary_stats.")

  ## Input y limits messages, tells user the plot limit parameters ##
  if (is.null(y_limits) & is.null(y_range) & is.null(y_max) & is.null(y_min)){
    print("No specified y-axis parameters. Axis y-limits will be scaled per plot, as per y_limits = 'free'.")
    y_limits <- "free"
    } else if (!is.null(y_limits)){
      if ((y_limits == "fixed") & is.null(y_max) & is.null(y_min)){
        print("Specified y-limit: 'fixed'. Axis y-limits will fixed for all plots based on maximum and minimum y-values.")
        } else if ((y_limits == "fixed") & !is.null(y_max)){
          print(paste("Specified y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a maximum of y_max. Specified y_max: ", y_max, sep = ""))
        } else if ((y_limits == "fixed") & !is.null(y_min)){
          print(paste("Specified y-limit: 'fixed'. Axis y-limits will be fixed for all plots with a minimum of y_min. Specified y_min: ", y_min, sep = ""))
        } else if (y_limits == "free" & is.null(y_max) & is.null(y_min)){
          print("Specified y-limit: 'free'. Axis y-limits will be scaled per plot.")
        } else if ((y_limits == "free") & !is.null(y_max)){
          print(paste("Specified y-limit: 'free'. Axis y-limits will be scaled per plot with a maximum of y_max. Specified y_max: ", y_max, sep = ""))
        } else if ((y_limits == "free") & !is.null(y_min)){
          print(paste("Specified y-limit: 'free'. Axis y-limits will be scaled per plot with a minimum of y_min. Specified y_min: ", y_min, sep = ""))
        } 
    } else if (!is.null(y_range) & is.null(y_max) & is.null(y_min)){
      print(paste("Specified y-range: ", y_range, ". Axis y-limits will range ",
            y_range," units, split over the median.", 
            sep = ""))
    } else if (!is.null(y_range) & !is.null(y_min)) {
      print(paste("Specified y-range: ", y_range, ". Axis y-limits will range ",
            y_range," units from the y_min. Specified y_min: ", y_min, sep = ""))
    } else if (!is.null(y_range) & !is.null(y_max)) {
      print(paste("Specified y-range: ", y_range, ". Axis y-limits will range ",
            y_range," units from the y_max. Specified y_max: ", y_max, sep = ""))
    } 
  
  ## Combined ##           ##########################################################
  if (option == "combined"){
    # Generate increment and y limits based on parameters and y-values #
    if (!is.null(y_limits)){
      if (y_limits == 'fixed'){
        increment <- set_increment(omicsPlotter$comp_stats[[panel_y_axis]], option)
        ylims <- set_ylimits(omicsPlotter$comp_stats[[panel_y_axis]], 
                             increment = increment,
                             y_min = y_min, y_max = y_max)
      }
    } else if (is.null(y_limits) & is.null(y_range)){
      increment <- set_increment(omicsPlotter$comp_stats[[panel_y_axis]], option)
      ylims <- set_ylimits(omicsPlotter$comp_stats[[panel_y_axis]], 
                           increment = increment,
                           y_min = y_min, y_max = y_max)
    }
    
    # Nest data #
    plotter <- tidyr::separate(omicsPlotter$comp_stats, Comparison, 
                               c("comp1", "comp2"), sep = "_vs_", remove = FALSE) %>%
      reshape2::melt(id.vars = names(omicsPlotter$comp_stats), 
                     value.name = "Group") %>%
      merge(omicsPlotter$summary_stats)
    nestplotter <- plotter %>% tidyr::nest(-panel_variable)
    
    ##nestplotter <- omicsPlotter$comp_stats %>% tidyr::nest(-panel_variable)
    #Subset large groups ########### Take out later ######################################
    if (nrow(nestplotter) > 10){
      nestplotter <- nestplotter[1:10,]
    }
    # Generate plots from nested data #
    # 1) Generate an increment for adjusting y limits and text label position
    # 2) Genreate y limits and text position
    # 3) Define border colors
    # 4) Generate ggplot of data
    # 5) Pipe ggplot to plotly
    nestplotter <- nestplotter %>% 
      mutate(panel = trelliscopejs::map_plot(data, function(nestedCompStats) {
      
      # Generate increment and y limits based on parameters and y-values #
      if (!is.null(y_limits)){
        if (y_limits == 'free'){
          increment <-set_increment(nestedCompStats[[panel_y_axis]], option)
          ylims <- set_ylimits(nestedCompStats[[panel_y_axis]], 
                               increment = increment,
                               y_min = y_min, y_max = y_max)
        }
      } else if (!is.null(y_range)){
        increment <- set_increment(omicsPlotter$comp_stats[[panel_y_axis]], option)
        ylims <- set_ylimits(omicsPlotter$comp_stats[[panel_y_axis]],
                             increment = increment, y_min = y_min, 
                             y_max = y_max, y_range = y_range)
      }
      
      # Set text spacing above/below barplot #
      text_increment <-set_increment(ylims, option)
      if (any(is.na(nestedCompStats[[panel_y_axis]]))){
        tempfold <- nestedCompStats[[panel_y_axis]]
        tempfold[is.na(tempfold)] <- 0
        textadj <- tempfold + sign(tempfold)*text_increment*2
      } else {
        textadj <- nestedCompStats[[panel_y_axis]] + 
          sign(nestedCompStats[[panel_y_axis]])*text_increment*2
      }
      
      # Set border colors based on significance #
      bord <- rep(colors[1], length(nestedCompStats$P_value_T))
      bord[nestedCompStats$P_value_G < p_val & !is.na(nestedCompStats$P_value_G)] <-  colors[2]
      bord[nestedCompStats$P_value_T < p_val & !is.na(nestedCompStats$P_value_T)] <- colors[3]
      nestedCompStats <- data.frame(nestedCompStats, bord)
      
      # Set hover/labels excluding the panel_variable #
      hover_want <- c( attributes(omicsPlotter)$cnames$edata_cname,
                       "Group", "Count", "Comparison",
                       "P_value_G", "P_value_T", "Fold_change") 
      if (!(panel_y_axis %in% hover_want)){
        hover_want <- c(hover_want, panel_y_axis)
      }
      if (panel_variable %in% hover_want){
        hover_want <- hover_want[!(hover_want %in% panel_variable)]
      }
      hover_labs <- hover_want[hover_want %in% colnames(nestedCompStats)]
      text_labs <- c()
      label_labs <- c()
      
      for (row in 1:nrow(nestedCompStats)){
        row_text <- c()
        row_label <- c()
        for (label in hover_labs){
          if (label %in% c("P_value_G", "P_value_T")){
            row_label <- c(row_label, 
                           paste(paste(label, ":", sep = ""), 
                                 signif(nestedCompStats[row,][[label]], 3)))
            row_text <- c(row_text,
                          paste(paste(label, ":", sep = ""), 
                                signif(nestedCompStats[row,][[label]])))
          } else if (label %in% c("Fold_change", panel_y_axis)) {
            row_text <- c(row_text,
                          paste(paste(label, ":", sep = ""), 
                                signif(nestedCompStats[row,][[label]])))
          } else {
            row_text <- c(row_text,
                          paste(paste(label, ":", sep = ""), 
                                nestedCompStats[row,][[label]]))
          }
        }
        text_labs[row] <- capture.output(cat(row_text, sep = ", "))
        label_labs[row] <- capture.output(cat(row_label, sep = ", "))
      }
      nestedCompStats <- data.frame(nestedCompStats, 
                                    text = text_labs, 
                                    labels = label_labs)
      
      # Make ggplot #
      plot22 <- ggplot2::ggplot(data = nestedCompStats,
                  ggplot2::aes(x = as.character(nestedCompStats[[panel_x_axis]]), 
                           y = nestedCompStats[[panel_y_axis]],
                           color = bord,
                           fill = Group,
                           text = gsub(", ", "\n", text),
                           label = gsub(", ", "\n", labels)
                           )
                  ) +
        ggplot2::scale_color_manual(values = levels(nestedCompStats$bord)) +
        ggplot2::geom_col(position = "dodge", size = 1) + 
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::xlab(panel_x_axis) + 
        ggplot2::ylab(panel_y_axis) +
        ggplot2::labs(fill = "", color = "") +
        ggplot2::geom_text(aes(y = textadj, group = textadj), color = "black", 
                           position = ggplot2::position_dodge(width = 0.9)) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                           hjust = 1, 
                                                           vjust = 0.5)) +
        ggplot2::guides(color = FALSE) +
        ggplot2::coord_cartesian(ylim=ylims)
      
      # Make and return plotly for map_plot #
      plotly <- plotly::ggplotly(plot22, tooltip = c("text")) 
      for (plotter in 1:length(plotly$x$data)){
        if (plotly$x$data[[plotter]]$type != "bar"){
          plotly$x$data[[plotter]]$showlegend <- FALSE
          plotly$x$data[[plotter]]$hoveron <- "none"
          plotly$x$data[[plotter]]$hoverinfo <- "none"
        }
        name <- plotly$x$data[[plotter]]$name
        groups <- unique(nestedCompStats$Group)
        name <- groups[unlist(map(groups, 
                                  function(group) grepl(group, name)))]
        plotly$x$data[[plotter]]$name <- name
        print(plotly$x$data[[plotter]]$name)
      }
      return(plotly)
    }
  )
  )
    # Return nested table #
    return(nestplotter)
    fruits <- c("one apple", "two pears", "three bananas")

  ## ANOVA ##        
  } else if (option == "anova") {       #########################################
    # Generate increment and y limits based on parameters and y-values #
    if (!is.null(y_limits)){
      if (y_limits == 'fixed'){
        increment <-set_increment(omicsPlotter$comp_stats[[panel_y_axis]], option)
        ylims <- set_ylimits(omicsPlotter$comp_stats[[panel_y_axis]], 
                             increment = increment,
                             y_min = y_min, y_max = y_max)
      }
    } else if (is.null(y_limits) & is.null(y_range)){
      increment <- set_increment(omicsPlotter$comp_stats[[panel_y_axis]], option)
      ylims <- set_ylimits(omicsPlotter$comp_stats[[panel_y_axis]], 
                           increment = increment,
                           y_min = y_min, y_max = y_max)
    }
    
    # Nest data #
    # nestplotter <- omicsPlotter$comp_stats %>% tidyr::nest(-panel_variable) ###########
    plotter <- tidyr::separate(omicsPlotter$comp_stats, Comparison, 
                               c("comp1", "comp2"), sep = "_vs_", remove = FALSE) %>%
      reshape2::melt(id.vars = names(omicsPlotter$comp_stats), 
                     value.name = "Group") %>%
      merge(omicsPlotter$summary_stats)

    nestplotter <- plotter %>% tidyr::nest(-panel_variable)
    
    ###################### Subset ###################3 to be removed ####### #######
    if (nrow(nestplotter) > 10){
      nestplotter <- nestplotter[1:10,]
    }
    
    # Generate plots from nested data #
    # 1) Generate an increment for adjusting y limits and text label position
    # 2) Define border colors
    # 3) Generate ggplot of data
    # 4) Pipe ggplot to plotly
    nestplotter <- nestplotter %>% 
      mutate(panel = trelliscopejs::map_plot(data, function(nestedCompStats) { 
        
      # Generate increment and y limits based on parameters and y-values #
      if (!is.null(y_limits)){
        if (y_limits == 'free'){
          increment <- set_increment(nestedCompStats[[panel_y_axis]], option)
          ylims <- set_ylimits(nestedCompStats[[panel_y_axis]], 
                               increment = increment,
                               y_min = y_min, y_max = y_max)
        }
      } else if (!is.null(y_range)){
        increment <- set_increment(nestedCompStats[[panel_y_axis]], option)
        ylims <- set_ylimits(nestedCompStats[[panel_y_axis]],
                             increment = increment, y_min = y_min, 
                             y_max = y_max, y_range = y_range)
      }
        
      # Set text spacing #
      text_increment <-set_increment(ylims, option)
      if (any(is.na(nestedCompStats[[panel_y_axis]]))){
        tempfold <- nestedCompStats[[panel_y_axis]]
        tempfold[is.na(tempfold)] <- 0
        textadj <- tempfold + sign(tempfold)*text_increment*2
      } else {
        textadj <- nestedCompStats[[panel_y_axis]] + 
          sign(nestedCompStats[[panel_y_axis]])*text_increment*2
      }
      
      # Set border colors based on significance #
      bord <- rep(colors[1], length(nestedCompStats$P_value_T))
      bord[nestedCompStats$P_value_G < p_val & !is.na(nestedCompStats$P_value_G)] <-  colors[2]
      bord[nestedCompStats$P_value_T < p_val & !is.na(nestedCompStats$P_value_T)] <- colors[3]
      nestedCompStats <- data.frame(nestedCompStats, bord)
      
      # Set hover/labels excluding the panel_variable #
      hover_want <- c( attributes(omicsPlotter)$cnames$edata_cname,
                       "Group", "Count", "Comparison",
                       "P_value_G", "P_value_T", "Fold_change") 
      if (!(panel_y_axis %in% hover_want)){
        hover_want <- c(hover_want, panel_y_axis)
      }
      if (panel_variable %in% hover_want){
        hover_want <- hover_want[!(hover_want %in% panel_variable)]
      }
      hover_labs <- hover_want[hover_want %in% colnames(nestedCompStats)]
      text_labs <- c()
      label_labs <- c()
      
      for (row in 1:nrow(nestedCompStats)){
        row_text <- c()
        row_label <- c()
        for (label in hover_labs){
          if (label %in% c("P_value_G", "P_value_T")){
            row_label <- c(row_label, 
                           paste(paste(label, ":", sep = ""), 
                                 signif(nestedCompStats[row,][[label]], 3)))
            row_text <- c(row_text,
                          paste(paste(label, ":", sep = ""), 
                                signif(nestedCompStats[row,][[label]])))
          } else if (label %in% c("Fold_change", panel_y_axis)) {
            row_text <- c(row_text,
                          paste(paste(label, ":", sep = ""), 
                                signif(nestedCompStats[row,][[label]])))
          } else {
            row_text <- c(row_text,
                          paste(paste(label, ":", sep = ""), 
                                nestedCompStats[row,][[label]]))
          }
        }
        text_labs[row] <- capture.output(cat(row_text, sep = ", "))
        label_labs[row] <- capture.output(cat(row_label, sep = ", "))
      }
      nestedCompStats <- data.frame(nestedCompStats, 
                                    text = text_labs, 
                                    labels = label_labs)
      
      # Make ggplot #
      # 1) Define variables, with text passed on to plotly and labels above data
      # 2) Adjust border color for significance based on flag and color vector
      # 3) Add columns and line at 0
      # 4) Update label names
      # 5) Adjust label position to be above/below data as appropriate
      # 6) Adjust label angle
      # 7) Adjust y limits to account for label positions      
      plot22 <- ggplot2::ggplot(data = nestedCompStats, 
              ggplot2::aes(x = as.character(nestedCompStats[[panel_x_axis]]), 
                           y = nestedCompStats[[panel_y_axis]],
                           color = bord,
                           fill = Group,
                           # color = as.factor(nestedCompStats[[panel_x_axis]]),
                           # fill = as.factor(nestedCompStats[[panel_x_axis]]),
                           # label = paste("p-value T:",signif(P_value_T, 3)),
                           text = gsub(", ", "\n", text),
                           label = gsub(", ", "\n", labels)
                             # paste(
                             # paste("Group:", Group),
                             # paste("Count:", Count),
                             # paste("Comparison:", Comparison),
                             # paste(paste(panel_y_axis,":", sep = ""),
                             #       round(nestedCompStats[[panel_y_axis]], 6)),
                             # paste("T-test p-value:",signif(P_value_T, 6)),
                             # sep = "\n")
                           )
              ) +
        ggplot2::scale_color_manual(values = levels(nestedCompStats$bord)) +
        ggplot2::geom_col(position = "dodge", size = 1) + 
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::xlab(panel_x_axis) + 
        ggplot2::ylab(panel_y_axis) +
        ggplot2::labs(fill = "", color = "") +
        ggplot2::geom_text(aes(y = textadj, group = textadj), color = "black", 
                           position = ggplot2::position_dodge(width = 0.9)) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                           hjust = 1, 
                                                           vjust = 0.5)) +
        ggplot2::guides(color = FALSE) +
        ggplot2::coord_cartesian(ylim=ylims) 
      
      # Make and return plotly #
      plotly <- plotly::ggplotly(plot22, tooltip = c("text")) 
      for (plotter in 1:length(plotly$x$data)){
        if (plotly$x$data[[plotter]]$type != "bar"){
          plotly$x$data[[plotter]]$showlegend <- FALSE
          plotly$x$data[[plotter]]$hoveron <- "none"
          plotly$x$data[[plotter]]$hoverinfo <- "none"
        }
        name <- plotly$x$data[[plotter]]$name
        groups <- unique(nestedCompStats$Group)
        name <- groups[unlist(map(groups, 
                                  function(group) grepl(group, name)))]
        plotly$x$data[[plotter]]$name <- name
        print(plotly$x$data[[plotter]]$name)
      }
      return(plotly)
    }
    )
    )

    return(nestplotter)
    
  ## G test ##                                              ######################
  } else if (option == "gtest") {
    
    if (!is.null(y_limits)){
      if (y_limits == 'fixed'){
        increment <-set_increment(omicsPlotter$comp_stats[[panel_y_axis]], option)
        ylims <- set_ylimits(omicsPlotter$comp_stats[[panel_y_axis]], 
                             increment = increment,
                             y_min = y_min, y_max = y_max)
      }
    } else if (is.null(y_limits) & is.null(y_range)){
      increment <- set_increment(omicsPlotter$comp_stats[[panel_y_axis]], option)
      ylims <- set_ylimits(omicsPlotter$comp_stats[[panel_y_axis]], 
                           increment = increment,
                           y_min = y_min, y_max = y_max)
    }
    
    # Generate combined and nested data for plotting #
    # 1) Split comparisons in comp_stats into multiple columns in dataframe
    # 2) Melt new stats dataframe to isolate groups per comparison per unique ID
    # 3) Merge with summary_stats for dataframe containing compared counts only
    # 4) Nest the new dataframe by unique ID
    plotter <- tidyr::separate(omicsPlotter$comp_stats, Comparison, 
                        c("comp1", "comp2"), sep = "_vs_", remove = FALSE) %>%
      reshape2::melt(id.vars = names(omicsPlotter$comp_stats), 
                     value.name = "Group") %>%
      merge(omicsPlotter$summary_stats)
    nestplotter <- plotter %>% tidyr::nest(-panel_variable)
    
    # Generate plots from nested data #
    # 1) Generate an increment for adjusting y limits and text label position
    # 2) Define border colors
    # 3) Generate ggplot of data
    # 4) Pipe ggplot to plotly
    if (nrow(nestplotter) > 10){
      nestplotter <- nestplotter[1:10,]
    }
    nestplotter <- nestplotter %>% 
      mutate(panel = trelliscopejs::map_plot(data, function(nestedCompStats) { 
      
      # Generate increment based on maximum y (no negative counts) #
      # 1) Extract non-na y values from the data
      # 2) Assign increment as 1/20th the maximum y
      
      if (!is.null(y_limits)){
        if (y_limits == 'free'){
          increment <-set_increment(nestedCompStats[[panel_y_axis]], option)
          ylims <- set_ylimits(nestedCompStats[[panel_y_axis]], 
                               increment = increment,
                               y_min = y_min, y_max = y_max)
        }
      } else if (!is.null(y_range)){
        increment <- set_increment(nestedCompStats[[panel_y_axis]], option)
        ylims <- set_ylimits(nestedCompStats[[panel_y_axis]],
                             increment = increment, y_min = y_min, 
                             y_max = y_max, y_range = y_range)
      }
      
      ## Set text spacing ##
      text_increment <-set_increment(ylims, option)
      if (any(is.na(nestedCompStats[[panel_y_axis]]))){
        tempfold <- nestedCompStats[[panel_y_axis]]
        tempfold[is.na(tempfold)] <- 0
        textadj <- max(tempfold) + sign(max(tempfold))*text_increment*2
      } else {
        textadj <- max(nestedCompStats[[panel_y_axis]]) + text_increment*2
      }
      
      ## Set border colors based on significance ##
      bord <- rep(colors[1], length(nestedCompStats$P_value_T))
      bord[nestedCompStats$P_value_G < p_val & !is.na(nestedCompStats$P_value_G)] <-  colors[2]
      bord[nestedCompStats$P_value_T < p_val & !is.na(nestedCompStats$P_value_T)] <- colors[3]
      nestedCompStats <- data.frame(nestedCompStats, bord)
      
      # Set hover/labels excluding the panel_variable #
      hover_want <- c( attributes(omicsPlotter)$cnames$edata_cname,
                       "Group", "Count", "Comparison",
                       "P_value_G", "P_value_T", "Fold_change") 
      if (!(panel_y_axis %in% hover_want)){
        hover_want <- c(hover_want, panel_y_axis)
      }
      if (panel_variable %in% hover_want){
        hover_want <- hover_want[!(hover_want %in% panel_variable)]
      }
      hover_labs <- hover_want[hover_want %in% colnames(nestedCompStats)]
      text_labs <- c()
      label_labs <- c()
      
      for (row in 1:nrow(nestedCompStats)){
        row_text <- c()
        row_label <- c()
        for (label in hover_labs){
          if (label %in% c("P_value_G", "P_value_T")){
            row_label <- c(row_label, 
                           paste(paste(label, ":", sep = ""), 
                                 signif(nestedCompStats[row,][[label]], 3)))
            row_text <- c(row_text,
                          paste(paste(label, ":", sep = ""), 
                                signif(nestedCompStats[row,][[label]])))
          } else if (label %in% c("Fold_change", panel_y_axis)) {
            row_text <- c(row_text,
                          paste(paste(label, ":", sep = ""), 
                                signif(nestedCompStats[row,][[label]])))
          } else {
            row_text <- c(row_text,
                          paste(paste(label, ":", sep = ""), 
                                nestedCompStats[row,][[label]]))
          }
        }
        text_labs[row] <- capture.output(cat(row_text, sep = ", "))
        label_labs[row] <- capture.output(cat(row_label, sep = ", "))
      }
      nestedCompStats <- data.frame(nestedCompStats, 
                                    text = text_labs, 
                                    labels = label_labs)
      
      # Make ggplot #
      # 1) Define variables, with text passed on to plotly and labels above data
      # 2) Adjust border color for significance based on flag and color vector
      # 3) Add columns and line at 0
      # 4) Update label names
      # 5) Adjust label position to be above data
      # 6) Adjust label angle
      # 7) Adjust y limits to account for label positions  
      plot22 <- ggplot2::ggplot(data = nestedCompStats,
              ggplot2::aes(x = as.character(nestedCompStats[[panel_x_axis]]), 
                           y = nestedCompStats[[panel_y_axis]],
                           fill = Group,
                           text = gsub(", ", "\n", text),
                           label = gsub(", ", "\n", labels)
                           )
              ) +
        ggplot2::scale_color_manual(values = levels(nestedCompStats$bord))  +
        ggplot2::geom_col(position = "dodge", size = 1) +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::xlab(panel_x_axis) + 
        ggplot2::ylab(panel_y_axis) +
        ggplot2::labs(fill = "", color = "") +
        ggplot2::geom_text(aes(y = textadj, group = textadj), color = "black", 
                           position = "dodge") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                           hjust = 1, 
                                                           vjust = 0.5)) +
        ggplot2::guides(color = FALSE) +
        ggplot2::coord_cartesian(ylim=ylims)
      
      # Make and return plotly #
      plotly <- plotly::ggplotly(plot22, tooltip = c("text"))
      for (plotter in 1:length(plotly$x$data)){
        if (plotly$x$data[[plotter]]$type != "bar"){
          plotly$x$data[[plotter]]$showlegend <- FALSE
          plotly$x$data[[plotter]]$hoveron <- "none"
          plotly$x$data[[plotter]]$hoverinfo <- "none"
        }
        name <- plotly$x$data[[plotter]]$name
        groups <- unique(nestedCompStats$Group)
        name <- groups[unlist(map(groups, 
                                  function(group) grepl(group, name)))]
        plotly$x$data[[plotter]]$name <- name
        print(plotly$x$data[[plotter]]$name)
      }
      return(plotly)
      }
    )
    )
    return(nestplotter)
  } 
}

