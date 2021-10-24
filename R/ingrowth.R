#' Ingrowth imputation for permanent sample plots
#'
#' Adds ingrowth records to a tree data frame
#'
#' @param treedata    Tree data frame. Contains at least plot, measurement,
#'                      and tree identifiers, and diameter values.
#' @param plotdata    Plot data frame. Contains at least a plot identifier,
#'                      and the minimum diameter threshold for the plot.
#' @param plot        Name of the plot identifier in \code{treedata} and
#'                      \code{plotdata}.
#' @param meas        Name of the measurement identifier in \code{treedata}.
#' @param tree        Name of the tree identifier in \code{treedata}.
#' @param diam        Name of the diameter values in \code{treedata}.
#' @param thresh      Name of the minimum diameter threshold in \code{treedata}.
#'
#' @return  A data frame with \code{treedata} followed by the imputed ingrowth
#'            tree records. A logical variable named \code{ingrowth} is appended,
#'            containing \code{FALSE} for the observed trees and \code{TRUE}
#'            for the ingrowth. The values in the ingrowth records for
#'            variables other than the identifiers and diameters are inherited
#'            from later measurements and are often meaningless.
#'
#' @export
#' @examples  withIngrowth <- ingrowth(treeData, plotData, plot="PlotId",
#'                                     meas = "MeasId", tree = "TreeId",
#'                                     diam = "DBH", thresh = "Thresh")
#'
#' @details  The data should not include plots with substantial ingrowth at
#'             the oldest measurement. The measurement identifier values must
#'             be capable of being sorted in time, e.g., they could be
#'             consecutive integer indices, years, or dates. The other
#'             identifiers can be of any type (numbers, strings, or factors).
#'             It is advisable to discard measurements where the imputed
#'             ingrowth basal area is a large fraction of the total.
#'
#' @references #' @references García, O. (2021)"Imputing Ingrowth in Even-aged Permanent
#'  Sample Plots". \emph{Mathematical and Computational Forestry
#'  & Natural-Resource Sciences (MCFNS)} (to appear)
#'
ingrowth <- function(treedata, plotdata, plot, meas, tree, diam, thresh){
  ingr <- treedata[NULL,]  # ingrowth, initially empty
  for(plt in unique(treedata[, plot])){  # one plot at a time
    data <- droplevels(treedata[treedata[, plot] == plt,])  # data for a plot
    dmin <- plotdata[match(plt, plotdata[, plot]), thresh]  # diameter threshold
    msrmnts <- sort(unique(data[, meas]))  # measurement IDs
    for(i in msrmnts[-length(msrmnts)]){ # current measurement, exclude last
      data_i <- data[data[, meas] == i, ]
      trees_i <- data_i[, tree]  # trees in measurement i
      trees_iplus <- trees_i  # including ingrowth (none so far)
      for(j in msrmnts[msrmnts > i]){  # future measurement
        data_j <- data[data[, meas] == j, ]
        trees_j <- data_j[, tree]  # trees in measurement j
        new <- !(trees_j %in% trees_iplus) # in j but not in i yet
        if(any(new)){  # add ingrowth
          ingr_ij <- data_j[new,]  # new ingrowth for i
          ingr_ij[, meas] <- i  # fix it
          dsumi <- sum(data_i[trees_i %in% trees_j, diam])  # over the shared
          dsumj <- sum(data_j[trees_j %in% trees_i, diam])  #   observations
          ingr_ij[, diam] <- pmin((dsumi / dsumj) * ingr_ij[, diam],
                                 dmin)  # estimated diameters
          ingr <- rbind(ingr, ingr_ij)  # cummulated ingrowth
          trees_iplus <- c(trees_iplus, trees_j[new])  # plus added ingrowth
        } # end if
      } # j
    } # i
  } # plt
  treedata$ingrowth <- FALSE  # observed
  if(nrow(ingr) > 0){
    ingr$ingrowth <- TRUE  # imputed
    treedata <- rbind(treedata, ingr)
  }
  return(treedata)
}
