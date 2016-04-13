nicePath <- function (filename) {
  tryCatch({
    return(file.path(getwd(), dirname(sys.frame(1)$ofile), filename))
  }, error = function (error) {
    return(filename)
  })
}

# basic operations with ontologies
if (!require(ontoCAT, quietly=TRUE)) getBioconductorPackage("ontoCAT")

# R bindings to C based network analysis package igraph
if (!require(igraph, quietly=TRUE)) install.packages("igraph")

# talk to vis.js, with igraph!
if (!require(visNetwork, quietly=TRUE)) install.packages("visNetwork")

# ==============================================================================
# using OWL

cco <- getOntology(normalizePath("../data/CellCultureOntology.owl"))

length(getAllTerms(cco)) # 297

# mostly CHEBIs, and EFO, then 3 others
getRootIds(cco)

for (root in as.list(getRootIds(cco))) {
  children <- getTermChildrenById(cco, root)

  for (child in children) {
    child <- unlist(strsplit(as.character(child), ":"))[1]

    print(root)
    print(child)
  }
}

# ==============================================================================
# using CSV


CCO <- read.csv(nicePath("../data/CCONT.csv"), stringsAsFactors=FALSE)
names(CCO)[-grep("^http", names(CCO))]

edges <- CCO[c("Class.ID", "Parents")]
colnames(edges) <- c("from", "to")

# go from http://www.orpha.net/ORDO/Orphanet_1390 to Orphanet_1390
edges$from <- unlist(lapply(strsplit(as.character(edges$from), "/"), function(x) tail(x, n=1)))
edges$to <- unlist(lapply(strsplit(as.character(edges$to), "/"), function(x) tail(x, n=1)))

length(unique(edges$from))
length(unique(edges$to))

extras <- unique(edges$to[!(edges$to %in% edges$from)])
length(extras) # 2
verts <- data.frame(name=c(edges$from, extras), label=c(CCO$Preferred.Label, "EXTRA", "EXTRA"))



G <- graph_from_data_frame(edges, vertices=verts)

visNet <- makeVisNetwork(G, useLabel=TRUE, cluster=TRUE, clusterAlg=cluster_fast_greedy, clusterAsUndirected=TRUE)

# minC <- rep(-Inf, vcount(G))
# maxC <- rep(Inf, vcount(G))
# minC[1] <- maxC[1] <- 0
# co <- layout_with_fr(G, minx=minC, maxx=maxC, miny=minC, maxy=maxC)
# plot(G, layout=co, vertex.size=0.01, vertex.label=NA, edge.arrow.width=0)
#makeVisNetwork(G, customLayout=co)

# too slow..
#eb <- cluster_edge_betweenness(G)
fg <- cluster_fast_greedy(as.undirected(G))

#plot_dendrogram(fg)

# get edges that connect two communities:
explorers <- as.integer(which(crossing(fg, G)))

# take all the rest:
proles <- seq(1, length(E(G)))[-explorers]

G2 <- delete_edges(G, proles)
# clean out unusued vertices
G2 <- delete.vertices(G2, which(igraph::degree(G2) == 0))

makeVisNetwork(G2)
visNet <- makeVisNetwork(G2, useLabel=TRUE, hierarchicalLayout=TRUE, direction="LR", cluster=TRUE, clusterAlg=cluster_fast_greedy, clusterAsUndirected=TRUE)
