nicePath <- function (filename) {
  tryCatch({
    return(file.path(getwd(), dirname(sys.frame(1)$ofile), filename))
  }, error = function (error) {
    return(filename)
  })
}

getBioconductorPackage <- function (packageName) {
  source("https://bioconductor.org/biocLite.R")
  biocLite(packageName)
}

# Gene set enrichment data structures and methods
# Provides OBOCollection class
if (!require(GSEABase, quietly=TRUE)) getBioconductorPackage("GSEABase")

source(nicePath("ontology-explorer.r"))

# ==============================================================================
# functions

getGraphCleanedByClusters <- function(g, communities) {
  # get edges that connect two communities:
  explorers <- as.integer(which(crossing(communities, g)))

  # take all the rest:
  proles <- seq(1, length(E(g)))[-explorers]

  # delete these edges (which do not connect clusters)
  g <- delete_edges(g, proles)
  # and then clean out unusued vertices
  g <- delete.vertices(g, which(igraph::degree(g) == 0))

  return(g)
}

# ==============================================================================

fantom <- getOBOCollection(nicePath("../data/ff-phase2-140729.obo"))

g <- getIgraph(fantom)

fg <- cluster_fast_greedy(as.undirected(g))
# plot_dendrogram(fg)
plot(fg, g, vertex.size=0.01, vertex.label=NA, edge.arrow.width=0)
plot()

g2 <- getGraphCleanedByClusters(g, fg)
