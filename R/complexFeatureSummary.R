#' Summarize complex features
#' @description Summarize complex features
#' @export
summarizeComplexFeatures <- function(res){
  total_confirmed_hypotheses = length(unique(res$complex_id))
  total_features = nrow(res)
  sub=res
  sub[ , `:=`( COUNT = .N , IDX = 1:.N ) , by = complex_id ]
  sub = subset(sub,COUNT>=2)
  setkey(sub, "complex_id")
  total_hypotheses_with_multiple_features = nrow(unique(sub))
  data.table(total_confirmed_hypotheses = total_confirmed_hypotheses,
    total_features = total_features,
    total_hypotheses_with_multiple_features = total_hypotheses_with_multiple_features)
}

#' Visual complex feature summary
#' @description Visual complex feature summary
#' @export
plotComplexFeatureSummary <- function(res){
  subunit_density = ggplot(res) +
    geom_density(aes(x=n_subunits_detected))
  subunit_histogram = ggplot(res) +
    geom_histogram(aes(x=n_subunits_detected),binwidth = 1)
  completeness_density = ggplot(res) +
    geom_density(aes(x=completeness))
  completeness_histogram = ggplot(res) +
    geom_histogram(aes(x=completeness),binwidth = 0.05)
  peak_corr_density = ggplot(res) +
    geom_density(aes(x=peak_corr))
  peak_corr_histogram = ggplot(res) +
    geom_histogram(aes(x=peak_corr),binwidth = 0.05)
  multiplot(subunit_density, completeness_density, peak_corr_density, subunit_histogram, completeness_histogram, peak_corr_histogram, cols=2)
}

## http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
