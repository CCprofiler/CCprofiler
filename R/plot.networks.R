#' Plot networks
#' @import data.table
#' @import igraph
#' @param complex_hypothesis data.table with complex hypotheses
#' @param dt data.table with binary interactions between a and b
#' @param col_a character string specifying the color of verteces in a
#' @param col_b character string specifying the color of verteces in b
#' @param all_b logical TRUE if all vertices in b should be colored, FALSE to only color those connected to a's with the complex_hypothesis. Example: a = bait and b= prays, if FALSE, only preys that were found by pulling on the baits within the complex hypothesis are colored.
#' @export

plot_network <- function(complex_hypothesis,dt,col_a="darkblue",col_b="darkred",all_b=FALSE){
  a <- dt$a
  b <- dt$b
  dt_filtered=subset(dt, (a %in% complex_hypothesis$protein_id) & (b %in% complex_hypothesis$protein_id))
  a_filtered=dt_filtered$a
  b_filtered=dt_filtered$b
  g <- graph.empty(directed=FALSE) + vertices(unique(complex_hypothesis$protein_id))
  g$layout <- layout_in_circle
  g[V(g), V(g)] <- TRUE
  g <- simplify(g)
  E(g)$color="grey"
  V(g)$color="grey"
  g <- add_edges(g,c(rbind(a_filtered, b_filtered)),color="darkred")
  if (all_b) {
    V(g)$color[V(g)$name %in% b] <- col_b
  } else {
    V(g)$color[V(g)$name %in% b_filtered] <- col_b
  }
  V(g)$color[V(g)$name %in% a] <- col_a
  g <- simplify(g,remove.multiple=TRUE,edge.attr.comb="last")
  lab.locs <- radian.rescale(x=1:length(unique(complex_hypothesis$protein_id)), direction=-1, start=0)
  p <- plot(g,vertex.frame.color="grey",vertex.size=15,vertex.label.dist=1,vertex.label.degree=lab.locs,vertex.label.color="black",main=complex_hypothesis$complex_name[1])
  print(p)
}

radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}
