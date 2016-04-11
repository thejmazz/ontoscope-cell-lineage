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

# overloading makeVisNetwork

makeVisNetwork <- function (graph, customLayout="layout_nicely") {
  nodes <- as_data_frame(graph, what="vertices")
  # colnames(nodes) <- c("id")
  nodes <- data.frame(id=nodes$name, label=nodes$label)
  # nodes <- data.frame(id=nodes$name, label='')
  # nodes <- data.frame(id=nodes$name, label=nodes$desc)

  edges <- as_data_frame(graph, what="edges")
  visNet <<- visNetwork(nodes, edges, width = "100%") %>%
    # visHierarchicalLayout(levelSeparation=250) %>%
    # visHierarchicalLayout(direction="LR", levelSeparation=250) %>%
    visIgraphLayout(layout = customLayout) %>%
    visNodes(size=5) %>%
    visEdges(arrows="to", smooth=FALSE)
  visNet
}


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

getNameFromID <- function(fantom, id) {
  keyVals <- fantom@.kv[fantom@.kv$stanza_id == id,]
  name <- keyVals$value[keyVals$key == "name"]

  return(name)
}

# ==============================================================================
# playing with different clustering algs.

fantom <- getOBOCollection(nicePath("../data/ff-phase2-140729.obo"))

g <- getIgraph(fantom)
# TODO move into getIgraph
labels <- unlist(lapply(as.list(as_ids(V(g))), function(x) getNameFromID(fantom, x)))
g <- set_vertex_attr(g, "label", value=labels)

fg <- cluster_fast_greedy(as.undirected(g))
#plot(fg, g, vertex.size=0.01, vertex.label=NA, edge.arrow.width=0)
g2 <- getGraphCleanedByClusters(g, fg)

# takes too long
# eb <- cluster_edge_betweenness(g)

im <- cluster_infomap(g)
#plot(im, g, vertex.size=0.01, vertex.label=NA, edge.arrow.width=0)

# takes too long
#lp <- cluster_label_prop(g)

le <- cluster_leading_eigen(as.undirected(g))
#plot(le, g, vertex.size=0.01, vertex.label=NA, edge.arrow.width=0)

louv <- cluster_louvain(as.undirected(g))
# plot(louv, g, vertex.size=0.01, vertex.label=NA, edge.arrow.width=0)

# too slow
# opt <- cluster_optimal(g)

# does not work with unconnected graph
# sg <- cluster_spinglass(g)

wt <- cluster_walktrap(g)
# plot(wt, g, vertex.size=0.01, vertex.label=NA, edge.arrow.width=0)

# ==============================================================================
# Let's consider only CL terms:

CLs <- as_ids(V(g))[grep("^CL", as_ids(V(g)))]

# There are 475 of these. Is a graph with only CLs connected?

gCL <- filterByGood(g, CLs)
count_components(gCL) # 1. nice.
makeVisNetwork(gCL)

makeVisNetwork(getGraphCleanedByClusters(gCL, cluster_fast_greedy(as.undirected(gCL))))


gCL2 <- getGraphCleanedByClusters(gCL, cluster_fast_greedy(as.undirected(gCL)))
labels <- unlist(lapply(nodes$id, function(x) getNameFromID(fantom, x)))
gCL2 <- set_vertex_attr(gCL2, "label", value=labels)
