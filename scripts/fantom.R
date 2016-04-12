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

# interface between KEGG pathways and graph model
if (!require(KEGGgraph, quietly=TRUE)) getBioconductorPackage("KEGGgraph")

source(nicePath("ontology-explorer.r"))

# overloading makeVisNetwork

makeVisNetwork <- function (graph,
  smooth=FALSE, useLabel=TRUE,
  cluster=TRUE, clusterAlg=cluster_edge_betweenness, clusterAsUndirected=FALSE,
  hierarchicalLayout=FALSE, levelSeparation=250, direction="UD",
  igraphLayout=TRUE, layout="layout_nicely") {

  nodes <- as_data_frame(graph, what="vertices")
  if (useLabel) {
    colnames(nodes) <- c("id", "label")
  } else {
    nodes <- data.frame(id=nodes$name, label=nodes$name)
  }

  if (cluster) {
    if (clusterAsUndirected) {
      clusters <- clusterAlg(as.undirected(graph))
    } else {
      clusters <- clusterAlg(graph)
    }

    nodes$group = clusters$membership
  }

  edges <- as_data_frame(graph, what="edges")

  visNet <- visNetwork(nodes, edges, width="100%")

  if (hierarchicalLayout) {
    visNet <- visHierarchicalLayout(visNet, direction=direction, levelSeparation=levelSeparation)
  } else if (igraphLayout) {
    visNet <- visIgraphLayout(visNet, layout=layout)
  }

  visNet <- visNodes(visNet, size=5)
  visNet <- visEdges(visNet, arrows="to", smooth=smooth)

  visNet

  return(visNet)
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

edge_attr(gCL, "weight") <- edge_betweenness(gCL)

makeVisNetwork(gCL, clusterAlg=cluster_fast_greedy, clusterAsUndirected=TRUE)

makeVisNetwork(gCL, layout="layout_with_dh", clusterAlg=cluster_fast_greedy, clusterAsUndirected=TRUE)
makeVisNetwork(gCL, layout="layout_with_lgl")
makeVisNetwork(gCL, layout="layout_with_kk", clusterAlg=cluster_fast_greedy, clusterAsUndirected=TRUE)

makeVisNetwork(gCL, hierarchicalLayout=TRUE)

makeVisNetwork(gCL, clusterAlg=cluster_fast_greedy, clusterAsUndirected=TRUE)


makeVisNetwork(getGraphCleanedByClusters(gCL, cluster_fast_greedy(as.undirected(gCL))), hierarchicalLayout=TRUE)

# ==============================================================================
# hematopoietic only

verts <- as_data_frame(gCL, what="vertices")
verts[grep("hematopoietic", verts$label),]

# by inspection, count does not grow after 8 orders
hematopoietics <- make_ego_graph(gCL, "CL:0000988", order=8, mode="in")[[1]]
# and has all terms with hematopoietics
verts[grep("hematopoietic", verts$label),]$name %in% as_ids(V(hematopoietics))

makeVisNetwork(hematopoietics, hierarchicalLayout=TRUE, direction="LR")
makeVisNetwork(hematopoietics, layout="layout_with_kk")

# ==============================================================================
# working with KEGG

hsa <- parseKGML2Graph(nicePath("../data/hsa04640.xml"))
getKEGGnodeData(hsa)
# has no edges...

# so, lets make our own:
# see: http://www.genome.jp/kegg-bin/show_pathway?hsa04640+3815
hsa <- make_graph(~
                           lymphoid_stem_cell-+            lymphoid_related_dendritic_cell,
                                               pro_t_cell-+lymphoid_related_dendritic_cell,
                           lymphoid_stem_cell-+pro_t_cell,
                                                           DN3-+lymphoid_related_dendritic_cell,
                                                                                                        double_positive_cell-+gamma_phi_t_cell,
                                                                                                        double_positive_cell-+CD8_T_Cell,
                                               pro_t_cell-+DN3-+DN4-+intermediate_single_positive_cell-+double_positive_cell-+CD4_T_Cell,
                                                                                                        double_positive_cell-+Regulatory_T_Cell,
                                                                                                        double_positive_cell-+NKT_Cell,
                                               pro_t_cell-+NK_Cell_Precursor-+NK_Cell,
                           lymphoid_stem_cell-+            NK_Cell_Precursor,
                           lymphoid_stem_cell-+Pro_B_Cell-+Pre_I_B_Cell-+Pre_B_II_Cell-+Immature_B_Cell-+B_Cell,
  hematopoietic_stem_cell-+lymphoid_stem_cell,
  hematopoietic_stem_cell-+myeloid_stem_cell,
                                              CFU_GEMM-+CFU_Mast-+Mast_Cell,
                                              CFU_GEMM-+CFU_Bas-+      Myeloblast-+Basophilic_Myelocyte-+Basophil,
                                              CFU_GEMM-+CFU_E0-+       Myeloblast-+Eosinophilic_Myelocyte-+Eosionophil,
                                                                CFU_MDC-+                                 Myeloid_Related_Dendritic_Cell,
                                                                                                Monocyte-+Myeloid_Related_Dendritic_Cell,
                                                        CFU_GM-+CFU_MDC-+Monoblast-+Promonocyte-+Monocyte-+Macrophage,
                           myeloid_stem_cell-+CFU_GEMM-+CFU_GM-+CFU_G-+Myeloblast-+Neutrophilic_Myelocyte-+Neutrophil,
                                              CFU_GEMM-+BFU_E-+CFU_E-+Proerythroblast-+Erythrocyte,
                                              CFU_GEMM-+BFU_MK-+CFU_MK-+Megakaryocyte-+Platelets
)


makeVisNetwork(hsa, useLabel=FALSE, hierarchicalLayout=TRUE, direction="LR")
