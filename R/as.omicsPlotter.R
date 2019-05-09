#' Convert Omics data and pairwise statistics to a plotting object
#'@name as.omicsPlotter
#'
#' Converts a ResObject and its respective OmicsData into an easily plottable object.
#'
#' @param omicsData 
#' @param omicsStats
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
#' @author Rachel Richardson
#' @seealso \code{\link{as.pepData}}
#' @seealso \code{\link{as.lipidData}}
#' @seealso \code{\link{as.metabData}}
#'
#' @export
#' 
as.omicsPlotter <- function(omicsData, omicsStats, ...){
  .as.omicsPlotter(omicsData, omicsStats, ...)
}

## peptide data ##
.as.omicsPlotter<- function(omicsData, omicsStats, check.names = TRUE){
  
  uniqedatId <- attributes(omicsData)$cnames$edata_cname
  sampID <- attributes(omicsData)$cnames$fdata_cname
  
  # initial checks #
  
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
  
  ## Manipulate dataframes ##
  stats <- omicsStats$Full_results
  
  # Original data values from omicsData #
  data_values <- melt(omicsData$e_data, 
                      id.vars = uniqedatId, 
                      variable.name = sampID) %>%
    left_join(omicsData$f_data, by = sampID)
  
  # Comparison statistics by pairwise comparison from omicsStats #
  comp_stats <- bind_rows(attributes(pep_stats)$comparisons %>% 
    map(function(paircomp) {
      df <- cbind(
        stats[[uniqedatId]], 
        rep(paircomp, nrow(stats)),
        stats[str_detect(names(stats), 
                         pattern = paste(paircomp, "$", sep = ""))]
      )
      trimname <- str_remove_all(colnames(df)[3:ncol(df)], 
                                 paste("_", paircomp, sep = ""))
      colnames(df) <- c(uniqedatId, "Comparison", trimname)
      return(df)
      }
    ))
  
  # Comparison statistics by groups from omicsStats #
  summary_stats <- stats[,!str_detect(names(stats), 
                    paste(attributes(omicsStats)$comparisons, 
                    collapse = '|'))]
  
  # Re-sort for any grouping variables #
  if (!is.null(attributes(omicsStats)$group_DF)){
    index <- which(names(attributes(omicsStats)$group_DF) == sampID)
    groups <-  unique(attributes(omicsStats)$group_DF[,-index])
    summary_stats <- bind_rows(map(groups, function(group){
  
      df <- cbind(summary_stats[[uniqedatId]], 
                  rep(group, nrow(summary_stats)),
                  summary_stats[,str_detect(names(summary_stats), 
                    paste(group, "$", sep = ""))])
      trimname <- str_remove_all(colnames(df)[3:ncol(df)], 
                                 paste("_", group, sep = ""))
      colnames(df) <- c(uniqedatId, "Group", trimname)
      
      return(df)
      }
    ))
  }
  
  # store results #
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


##################################

load("../../subset_stats_objects.rds")      ##################loading data for function testing
omicsData <- norm_data
omicsStats <- pep_stats
### need to load prot_data for pro_stats

x <- as.omicsPlotter(omicsData, omicsStats)

class(x)
attributes(x)

x$data_values
x$summary_stats
x$comp_stats

#####################
set_ylimits <- function(foldchangevalues, increment){
  #Catch NAs
  if (any(is.na(foldchangevalues))){
    foldchangevalues[is.na(foldchangevalues)] <- 0
  }
  
  #Set maxi and mini based on difference between max and min values
  #Accounts for maximums below zero and minimums above zero
  maxi <- max(foldchangevalues) + 5*increment
  mini <- min(foldchangevalues) - 5*increment
  
  #Adjust for maximums above 0 and minimums below 0
  if (max(foldchangevalues) < 0){
    maxi <- 5*increment
  }
  if (min(foldchangevalues) > 0){
    mini <- -5*increment
  }
  return(c(mini, maxi))
}

###
plot_comps <- function(omicsPlotter) {
  
  # Check if stats statistical test attribute is valid #
  if(!(attributes(omicsPlotter)$statistical_test %in% c("combined", "gtest", "anova"))) stop(paste("Non-applicable statistical_test attribute in input."))
  
  
  omicsPlotterCname <- attributes(omicsPlotter)$cnames$edata_cname
  nestplotter <- omicsPlotter$comp_stats %>% nest(-omicsPlotterCname)
  nestplotter <- nestplotter[1:10,]     ##############
  option <- attributes(omicsPlotter)$statistical_test
  uniqedatId <- attributes(omicsPlotter)$cnames$edata_cname

## Combined ## 
  if (option == "combined"){
    nestplotter$data %>% map(function(nestedCompStats) {
      rmna <- nestedCompStats$Fold_change[!is.na(nestedCompStats$Fold_change)]
      if(length(rmna)==1){
        increment <- abs(rmna)/20
      } else {
        increment <- (max(nestedCompStats$Fold_change) - min(nestedCompStats$Fold_change))/20
      }
      
      colors <- c("red", "green4", "purple")
      plot22 <- ggplot(data = nestedCompStats, 
                       aes(x = Comparison, 
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
                       )) +
        scale_color_manual(values = colors[abs(nestedCompStats$Flag)+1]) +
        geom_col() + 
        geom_hline(yintercept = 0) +
        xlab("Pairwise Comparisons") + 
        ylab("Fold Change") +
        labs(fill = "", color = "") +
        geom_text(y = nestedCompStats$Fold_change +
                    sign(nestedCompStats$Fold_change)*increment*3) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ylim(set_ylimits(nestedCompStats$Fold_change, increment)) #Adjust for text labels
      ggplotly(plot22, tooltip = c("text")) 
    }
  )

  ## ANOVA ##        
  } else if (option == "anova") {
    nestplotter$data %>% map(function(nestedCompStats) {
      rmna <- nestedCompStats$Fold_change[!is.na(nestedCompStats$Fold_change)]
      if(length(rmna)==1){
        increment <- abs(rmna)/20
      } else {
        increment <- (max(nestedCompStats$Fold_change) - min(nestedCompStats$Fold_change))/20
      }
      
      colors <- c("red", "green4", "purple")
      plot22 <- ggplot(data = nestedCompStats, 
                       aes(x = Comparison, 
                           y = Fold_change,
                           color = Comparison,
                           fill = Comparison,
                           label = paste("p-value T:",signif(P_value_T, 3)),
                           text = paste(
                             paste("Fold change:",round(Fold_change, 6)),
                             paste("T-test p-value:",signif(P_value_T, 6)),
                             sep = "\n")
                       )) +
        scale_color_manual(values = colors[abs(nestedCompStats$Flag)+1]) +
        geom_col() + 
        geom_hline(yintercept = 0) +
        xlab("Pairwise Comparisons") + 
        ylab("Fold Change") +
        labs(fill = "", color = "") +
        geom_text(y = nestedCompStats$Fold_change +
                    sign(nestedCompStats$Fold_change)*increment*3) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ylim(set_ylimits(nestedCompStats$Fold_change, increment)) #Adjust for text labels
      ggplotly(plot22, tooltip = c("text")) 
    }
    )
    
  ## gtest ##  
  } else if (option == "gtest") {
    plotter <- separate(omicsPlotter$comp_stats, Comparison, 
                        c("comp1", "comp2"), sep = "_vs_", remove = FALSE) %>%
      melt(id.vars = names(omicsPlotter$comp_stats), value.name = "Group") %>%
      merge(omicsPlotter$summary_stats, by = c(uniqedatId, "Group"))
    nestplotter <- plotter %>% nest(-uniqedatId)
    
    nestplotter[1:10,]$data%>% map(function(nestedCompStats) {
      rmna <- nestedCompStats$Count[!is.na(nestedCompStats$Count)]
      increment <- max(rmna)/20
      
      colors <- c("red", "green4", "purple")
      plot22 <- ggplot(data = nestedCompStats, 
                       aes(x = Comparison, 
                           y = Count,
                           color = Group,
                           fill = Group,
                           label = paste(
                             paste("p-value G:",signif(P_value_G, 3)),
                             sep = "\n"),
                           text = paste("Count for Group", paste(Group, ":", sep = ""), Count)
                       )
      ) +
        scale_color_manual(values = colors[abs(nestedCompStats$Flag)+1]) +
        geom_col(position = "dodge") +
        geom_hline(yintercept = 0) +
        xlab("Pairwise Comparisons") + 
        ylab("Count") +
        labs(fill = "", color = "") +
        geom_text(y = max(nestedCompStats$Count) + increment*3) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ylim(set_ylimits(nestedCompStats$Count, increment)) #Adjust for text labels
      ggplotly(plot22, tooltip = c("text")) 
    }
    )
  } 
}

####################

y <- x
z <- x

attributes(y)$statistical_test <- "anova"
attributes(z)$statistical_test <- "gtest"

attributes(x)$statistical_test
attributes(y)$statistical_test
attributes(z)$statistical_test

plot_comps(x)
plot_comps(y)
plot_comps(z)

w <- x
attributes(w)$statistical_test <- "blah"
plot_comps(w)

