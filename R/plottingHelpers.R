#' Create a reproducible color map for a set of ids by ordering
#' them alphabetically
#' @importFrom scales hue_pal
createGGplotColMap <- function(ids){
  ids <- sort(unique(ids))
  colormap <- hue_pal()(length(ids))
  names(colormap) <- ids
  colormap
}

getMWcalibration <- function(fr_ann){
  tr <- lm(log(fr_ann$molecular_weight) ~ fr_ann$id)
  intercept <- as.numeric(tr$coefficients[1])
  slope <- as.numeric(tr$coefficients[2])
  FractionToMW <- as.formula(paste0("~exp(",slope, "* . + ", intercept, ")"))
  return(FractionToMW)
}
