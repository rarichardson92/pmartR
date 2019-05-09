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
as.omicsPlotter <- function(omicsData, omicsStats, ...){
  .as.omicsPlotter(omicsData, omicsStats, ...)
}

.as.omicsPlotter <- function(omicsData, omicsStats, check.names = TRUE){
  
  ## Moniker Variables ##
  
  uniqedatId <- attributes(omicsData)$cnames$edata_cname
  sampID <- attributes(omicsData)$cnames$fdata_cname
  stats <- omicsStats$Full_results
  
  ## Initial Checks ##
  
  # Check that omicsData and omicsStats are the correct classes #
  if (!inherits(omicsData, c("proData", "pepData", "metabData", "lipidData"))) stop("omicsData must be of class 'proData', 'pepData', 'metabData', or 'lipidData'")
  if(!inherits(omicsStats, "statRes")) stop("omicsStats must be of the class 'statRes'")
  
  # Check if omicsData and omicsStats have the same cname attributes. #
  if(!(identical(attributes(omicsData)$cnames, attributes(omicsStats)$cnames))) stop("Non-matching cname attributes in omicsStats and omicsData. Check that omicsStats is correctly derived from omicsData.")
  
  # Check that the biomolecule unique ID column exists in omicsData and omicsStats #
  if(!(uniqedatId %in% names(omicsStats$Full_results))) stop(paste("Column ", uniqedatId," not found in omicsStats. Requires compatible identifiers.", sep = ""))
  
  # Check if stats biomolecules are (at least) a subset of omicsData biomolecules. #
  if(!all(omicsStats$Full_results[[uniqedatId]] %in% norm_data$e_data[[uniqedatId]])) stop(paste("Biomolecules in omicsStats do not match biomolecules in omicsData.", sep = ""))
  
  # Check if stats sampleIDs are (at least) a subset of omicsData f_data sample IDs. #
  if(!all(attributes(omicsStats)$group_DF[[sampID]] %in% norm_data$f_data[[sampID]])) stop(paste(sampID, "column does not match between omicsData and omicsStats.", sep = ""))
  
  #### Group check??
  
  ## Manipulate Dataframes ##

  # Original data values from omicsData #
  # 1) Melt to assign e_data values to each sample ID
  # 2) Combine with groups in f_data
  data_values <- suppressWarnings(reshape2::melt(omicsData$e_data, 
                      id.vars = uniqedatId, 
                      variable.name = sampID) %>%
                        dplyr::left_join(omicsData$f_data, by = sampID))
  
  # Comparison statistics by pairwise comparison from omicsStats #
  # 1) For each comparison, extract comparison relevant columns from omicsStats
  #  2) Bind extracted with unique IDs from omicsStats and a comparison column
  #  3) Rename columns by removing comparison designation
  # 4) Bind rows from each comparison into one dataframe 
  comp_stats <- suppressWarnings(dplyr::bind_rows(attributes(pep_stats)$comparisons %>% 
    purrr::map(function(paircomp) {
      df <- cbind(
        stats[[uniqedatId]], 
        rep(paircomp, nrow(stats)),
        stats[stringr::str_detect(names(stats), 
                         pattern = paste(paircomp, "$", sep = ""))]
      )
      trimname <- stringr::str_remove_all(colnames(df)[3:ncol(df)], 
                                 paste("_", paircomp, sep = ""))
      colnames(df) <- c(uniqedatId, "Comparison", trimname)
      return(df)
      }
    )))

  # Comparison statistics by groups from omicsStats #
  # 1) Extract all rows w/o comparison information in the column names
  summary_stats <- stats[,!stringr::str_detect(names(stats), 
                    paste(attributes(omicsStats)$comparisons, 
                    collapse = '|'))]
  
  # Re-sort for any grouping variables #
  # 1) If there is a group_DF specified,
  # 2) Extract unique groups in group_DF, and for each group
  #  3) Bind extracted with unique IDs from omicsStats and a Group column
  #  4) Rename columns by removing group designation
  # 5) Bind rows from each group into one dataframe 
  if (!is.null(attributes(omicsStats)$group_DF)){
    index <- which(names(attributes(omicsStats)$group_DF) == sampID)
    groups <-  unique(attributes(omicsStats)$group_DF[,-index])
    summary_stats <- suppressWarnings(dplyr::bind_rows(purrr::map(groups, function(group){
  
      df <- cbind(summary_stats[[uniqedatId]], 
                  rep(group, nrow(summary_stats)),
                  summary_stats[,stringr::str_detect(names(summary_stats), 
                    paste(group, "$", sep = ""))])
      trimname <- stringr::str_remove_all(colnames(df)[3:ncol(df)], 
                                 paste("_", group, sep = ""))
      colnames(df) <- c(uniqedatId, "Group", trimname)
      
      return(df)
      }
    )))
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
  
  inheritDatAtt <- attributes(omicsData)[grep("filters|^names|group_DF", 
                              names(attributes(omicsData)), 
                              invert = T)]
  inheritStatAtt <- attributes(omicsStats)[grep("^names", 
                                              names(attributes(omicsStats)), 
                                              invert = T)]
  attributes(res) <- c(attributes(res), inheritDatAtt, inheritStatAtt)
  
  class(res) <- "omicsPlotter"
  
  return(res)
}

#' @name set_ylimits
#' @rdname set_ylimits
#' @title Sets y-axis limits for Omicsplotter plotting
#' 
#' @description Sets y-axis limits for Omicsplotter plotting. Used by plot_comp.
#'
#' @param yvalues y-values for plotting
#' @param increment An increment set based on the maximum/minimum y-values
#'
#' @author Rachel Richardson
#'
#' @export
 
set_ylimits <- function(yvalues, increment){
   .set_ylimits(yvalues, increment)
}

.set_ylimits <- function(yvalues, increment){
  
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
  
  # Set maxi and mini based on difference between max and min values #
  maxi <- max(yvalues) + 5*increment
  mini <- min(yvalues) - 5*increment
  
  # Adjust for maximums below 0 and minimums above 0 #
  if (max(yvalues) < 0){
    maxi <- 5*increment
  }
  if (min(yvalues) > 0){
    mini <- -5*increment
  }
  
  ## Return limits ##
  # Minimum y value = mini, maximum y value = maxi #
  return(c(mini, maxi))
}


#' @name plot_comps
#' @rdname plot_comps
#' @title Plot pairwise comparisons of omicsPlotter
#'
#' @description Plot pairwise comparisons of omicsPlotter. Plots differently depending on statistical tests ran, located in attributes.
#'
#' @param omicsPlotter An object of class "omicsPlotter" generated from \code{\link{as.omicsPlotter}}.
#'
#' @author Rachel Richardson
#' @export
plot_comps <- function(omicsPlotter) {
   .plot_comps(omicsPlotter)
}

.plot_comps <- function(omicsPlotter) {

  ## Initial Checks ##
  
  # Check check if class is correct #
  if(!inherits(omicsPlotter, "omicsPlotter")) stop("omicsPlotter must be of the class 'omicsPlotter'")  
    
  # Check if stats statistical test attribute is valid #
  if(!(attributes(omicsPlotter)$statistical_test %in% c("combined", "gtest", "anova"))) stop(paste("Non-applicable statistical_test attribute in omicsPlotter object."))
  
  ## Moniker Variables ##
  omicsPlotterCname <- attributes(omicsPlotter)$cnames$edata_cname
  option <- attributes(omicsPlotter)$statistical_test
  uniqedatId <- attributes(omicsPlotter)$cnames$edata_cname
  colors <- c("red", "green4", "purple")

  ## Combined ##
  # Nest data #
  if (option == "combined"){
    nestplotter <- omicsPlotter$comp_stats %>% tidyr::nest(-omicsPlotterCname)
    
  # Generate plots from nested data #
    # 1) Generate an increment for adjusting y limits and text label position
    # 2) Define border colors
    # 3) Generate ggplot of data
    # 4) Pipe ggplot to plotly
    nestplotter$data[1:10] %>% purrr::map(function(nestedCompStats) { ################
      
      # Generate increment based on difference between maximum and minimum y #
      # 1) Extract non-na y values from the data
      # 2) Where there is only one value, the increment is 1/20th of that value
      # 3) Else, the increment is 1/20th the difference betweeen max and min y
      rmna <- nestedCompStats$Fold_change[!is.na(nestedCompStats$Fold_change)]
      if(length(rmna)==1){
        increment <- abs(rmna)/20
      } else {
        increment <- (max(nestedCompStats$Fold_change) - min(nestedCompStats$Fold_change))/20
      }
      
      # Make ggplot #
      # 1) Define variables, with text passed on to plotly and labels above data
      # 2) Adjust border color for significance based on flag and color vector
      # 3) Add columns and line at 0
      # 4) Update label names
      # 5) Adjust label position to be above/below data as appropriate
      # 6) Adjust label angle
      # 7) Adjust y limits to account for label positions
      plot22 <- ggplot2::ggplot(data = nestedCompStats, 
                  ggplot2::aes(x = Comparison, 
                           y = Fold_change,
                           color = Comparison,
                           fill = Comparison,
                           label = paste(
                             paste("p-value T:",signif(P_value_T, 3)),
                             paste("p-value G:",signif(P_value_G, 3)),
                             sep = "\n"),
                           text = paste(
                             paste("Fold change:",round(Fold_change, 6)),
                             paste("T-test p-value:",signif(P_value_T, 6)),
                             paste("G-test p-value:",signif(P_value_G, 6)),
                             sep = "\n")
                           )
                  ) +
        ggplot2::scale_color_manual(values = colors[abs(nestedCompStats$Flag)+1]) +
        ggplot2::geom_col() + 
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::xlab("Pairwise Comparisons") + 
        ggplot2::ylab("Fold Change") +
        ggplot2::labs(fill = "", color = "") +
        ggplot2::geom_text(y = nestedCompStats$Fold_change +
                    sign(nestedCompStats$Fold_change)*increment*3) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                  hjust = 1, 
                                                  vjust = 0.5)) +
        ggplot2::ylim(set_ylimits(nestedCompStats$Fold_change, 
                                  increment)) #Adjust for text labels
      
      # Make and return plotly #
      plotly::ggplotly(plot22, tooltip = c("text")) 
    }
  )

  ## ANOVA ##        
  } else if (option == "anova") {
    # Nest data #
    nestplotter <- omicsPlotter$comp_stats %>% tidyr::nest(-omicsPlotterCname)
    
    # Generate plots from nested data #
    # 1) Generate an increment for adjusting y limits and text label position
    # 2) Define border colors
    # 3) Generate ggplot of data
    # 4) Pipe ggplot to plotly
    nestplotter$data[1:10] %>% purrr::map(function(nestedCompStats) {  ################
      
      # Generate increment based on difference between maximum and minimum y #
      # 1) Extract non-na y values from the data
      # 2) Where there is only one value, the increment is 1/20th of that value
      # 3) Else, the increment is 1/20th the difference betweeen max and min y
      rmna <- nestedCompStats$Fold_change[!is.na(nestedCompStats$Fold_change)]
      if(length(rmna)==1){
        increment <- abs(rmna)/20
      } else {
        increment <- (max(nestedCompStats$Fold_change) - min(nestedCompStats$Fold_change))/20
      }
      
      # Make ggplot #
      # 1) Define variables, with text passed on to plotly and labels above data
      # 2) Adjust border color for significance based on flag and color vector
      # 3) Add columns and line at 0
      # 4) Update label names
      # 5) Adjust label position to be above/below data as appropriate
      # 6) Adjust label angle
      # 7) Adjust y limits to account for label positions      
      plot22 <- ggplot2::ggplot(data = nestedCompStats, 
              ggplot2::aes(x = Comparison, 
                           y = Fold_change,
                           color = Comparison,
                           fill = Comparison,
                           label = paste("p-value T:",signif(P_value_T, 3)),
                           text = paste(
                             paste("Fold change:",round(Fold_change, 6)),
                             paste("T-test p-value:",signif(P_value_T, 6)),
                             sep = "\n")
                           )
              ) +
        ggplot2::scale_color_manual(values = colors[abs(nestedCompStats$Flag)+1]) +
        ggplot2::geom_col() + 
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::xlab("Pairwise Comparisons") + 
        ggplot2::ylab("Fold Change") +
        ggplot2::labs(fill = "", color = "") +
        ggplot2::geom_text(y = nestedCompStats$Fold_change +
                    sign(nestedCompStats$Fold_change)*increment*3) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                  hjust = 1, 
                                                  vjust = 0.5)) +
        ggplot2::ylim(set_ylimits(nestedCompStats$Fold_change,
                                  increment)) #Adjust for text labels
      
      # Make and return plotly #
      plotly::ggplotly(plot22, tooltip = c("text")) 
    }
    )
    
  ## G test ##  
  } else if (option == "gtest") {
    # Generate combined and nested data for plotting #
    # 1) Split comparisons in comp_stats into multiple columns in dataframe
    # 2) Melt new stats dataframe to isolate groups per comparison per unique ID
    # 3) Merge with summary_stats for dataframe containing compared counts only
    # 4) Nest the new dataframe by unique ID
    plotter <- tidyr::separate(omicsPlotter$comp_stats, Comparison, 
                        c("comp1", "comp2"), sep = "_vs_", remove = FALSE) %>%
      reshape2::melt(id.vars = names(omicsPlotter$comp_stats), 
                     value.name = "Group") %>%
      merge(omicsPlotter$summary_stats, by = c(uniqedatId, "Group"))
    nestplotter <- plotter %>% tidyr::nest(-uniqedatId)
    
    # Generate plots from nested data #
    # 1) Generate an increment for adjusting y limits and text label position
    # 2) Define border colors
    # 3) Generate ggplot of data
    # 4) Pipe ggplot to plotly
    nestplotter$data[1:10] %>% purrr::map(function(nestedCompStats) { ################
      
      # Generate increment based on maximum y (no negative counts) #
      # 1) Extract non-na y values from the data
      # 2) Assign increment as 1/20th the maximum y
      rmna <- nestedCompStats$Count[!is.na(nestedCompStats$Count)]
      increment <- max(rmna)/20
      
      # Make ggplot #
      # 1) Define variables, with text passed on to plotly and labels above data
      # 2) Adjust border color for significance based on flag and color vector
      # 3) Add columns and line at 0
      # 4) Update label names
      # 5) Adjust label position to be above data
      # 6) Adjust label angle
      # 7) Adjust y limits to account for label positions  
      plot22 <- ggplot2::ggplot(data = nestedCompStats, 
              ggplot2::aes(x = Comparison, 
                           y = Count,
                           color = Group,
                           fill = Group,
                           label = paste(
                             paste("p-value G:",signif(P_value_G, 3)),
                             sep = "\n"),
                           text = paste("Count for Group", 
                                        paste(Group, ":", sep = ""), 
                                        Count)
                           )
              ) +
        ggplot2::scale_color_manual(values = colors[abs(nestedCompStats$Flag)+1]) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::xlab("Pairwise Comparisons") + 
        ggplot2::ylab("Count") +
        ggplot2::labs(fill = "", color = "") +
        ggplot2::geom_text(y = max(nestedCompStats$Count) + increment*3) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                  hjust = 1, 
                                                  vjust = 0.5)) +
        ggplot2::ylim(set_ylimits(nestedCompStats$Count, 
                                  increment)) #Adjust for text labels
      
      # Make and return plotly #
      plotly::ggplotly(plot22, tooltip = c("text")) 
    }
    )
  } 
}

