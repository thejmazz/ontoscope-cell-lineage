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

# if (!require(GraphAlignment, quietly=TRUE)) getBioconductorPackage("GraphAlignment")

source(nicePath("ontology-explorer.r"))

# overloading makeVisNetwork

makeVisNetwork <- function (graph,
  smooth=FALSE, useLabel=TRUE,
  cluster=TRUE, clusterAlg=cluster_edge_betweenness, clusterAsUndirected=FALSE,
  customGroups=FALSE,
  hierarchicalLayout=FALSE, levelSeparation=250, direction="UD",
  igraphLayout=TRUE, layout="layout_nicely") {

  nodes <- as_data_frame(graph, what="vertices")
  if (useLabel) {
    nodes <- data.frame(id=nodes$name, label=nodes$label)
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

  # customGroups takes precedence over cluster
  if (customGroups) {
    nodes$group <- as_data_frame(graph, what="vertices")$group
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

# visNet <- makeVisNetwork(gCL, clusterAlg=cluster_walktrap, clusterAsUndirected=FALSE, layout="layout_with_kk")

makeVisNetwork(gCL, layout="layout_with_dh", clusterAlg=cluster_fast_greedy, clusterAsUndirected=TRUE)
makeVisNetwork(gCL, layout="layout_with_lgl")
makeVisNetwork(gCL, layout="layout_with_kk", clusterAlg=cluster_fast_greedy, clusterAsUndirected=TRUE)

makeVisNetwork(gCL, hierarchicalLayout=TRUE)

makeVisNetwork(gCL, clusterAlg=cluster_fast_greedy, clusterAsUndirected=TRUE)


makeVisNetwork(getGraphCleanedByClusters(gCL, cluster_fast_greedy(as.undirected(gCL))), hierarchicalLayout=TRUE)

# ==============================================================================
# hematopoietic only

v <- as_data_frame(gCL, what="vertices")
v[grep("hematopoietic", v$label),]

# by inspection, count does not grow after 8 orders
hematopoietics <- make_ego_graph(gCL, "CL:0000988", order=8, mode="in")[[1]]
# and has all terms with hematopoietics
v[grep("hematopoietic", v$label),]$name %in% as_ids(V(hematopoietics))

visNet <- makeVisNetwork(hematopoietics, hierarchicalLayout=TRUE, direction="LR")
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

verts <- as_ids(V(hsa))
labels <- as.character(1:length(verts))

# renaming to match fantom
verts[1]; labels[1] <- "lymphoid_stem_cell"
verts[2]; labels[2] <- "lymphoid_related_dendritic_cell"
verts[3] <- "CL:0000827"; labels[3] <- "pro-T cell" # pro_t_cell ->
verts[4] <- "CL:0000807"; labels[4] <- "DN3 thymocyte" # DN3 ->
verts[5] <- "CL:0002427"; labels[5] <- "resting double-positive thymocyte" # double_positive_cell ->
verts[6] <- "CL:0000809"; labels[6] <- "double-positive, alpha-beta thymocyte" # gamma_phi_t_cell ->
verts[7]; labels[7] <- "CD8_T_Cell"
verts[8] <- "CL:0000808"; labels[8] <- "DN4 thymocyte" # DN4 ->
verts[9] <- "CL:0000805"; labels[9] <- "immature single positive thymocyte" # intermediate_single_positive_cell ->
verts[10] <- "CL:0002431"; labels[10] <- "CD4-positive, CD8-intermediate double-positive thymocyte" # CD4_T_Cell ->
verts[11] <- "CL:0000815"; labels[11] <- "regulatory T cell" # Regulatory_T_Cell ->
verts[12] <- "CL:0000814"; labels[12] <- "mature NK T cell" # NKT_Cell ->
verts[13] <- "CL:0000914"; labels[13] <- "immature NK T cell" # NK_Cell_Precursor ->
verts[14] <- "CL:0000825"; labels[14] <- "pro-NK cell" # NK_Cell ->
verts[15] <- "CL:0000826"; labels[15] <- "pro-B cell" # Pro_B_Cell ->
verts[16] <- "CL:0000817"; labels[16] <- "precursor B cell" # Pre_I_B_Cell ->
verts[17] <- "CL:0000955"; labels[17] <- "pre-B-II cell" # Pre_B_II_Cell ->
verts[18] <- "CL:0000816"; labels[18] <- "immature B cell" # Immature_B_Cell ->
verts[19] <- "CL:0000236"; labels[19] <- "B cell" # B_Cell ->
verts[20] <- "CL:0000037"; labels[20] <- "Hematopoietic Stem Cell" # hematopoietic_stem_cell ->
verts[21]; labels[21] <- "myeloid_stem_cell"
# CFU-GEMM cells are the multipotential progenitor cells for myeloid cells
verts[22] <- "CL:0000049"; labels[22] <- "common myeloid progenitor" # CFU_GEMM ->
# CFU-Mast is a colony forming unit. It gives rise to mast cells.
verts[23] <- "CL:0000831"; labels[23] <- "mast cell progenitor" # CFU_Mast ->
verts[24] <- "CL:0000097"; labels[24] <- "mast cell" # Mast_Cell ->
verts[25] <- "CL:0002028"; labels[25] <- "basophil mast progenitor cell" # CFU_Bas ->
verts[26] <- "CL:0000835"; labels[26] <- "myeloblast" # Myeloblast ->
verts[27]; labels[27] <- "Basophilic_Myelocyte"
verts[28] <- "CL:0000767"; labels[28] <- "basophil" # Basophil ->
verts[29]; labels[29] <- "CFU_E0"
verts[30]; labels[30] <- "Eosinophilic_Myelocyte"
verts[31] <- "CL:0000771"; labels[31] <- "eosinophil" # Eosionophil ->
verts[32]; labels[32] <- "CFU_MDC"
verts[33]; labels[33] <- "Myeloid_Related_Dendritic_Cell"
verts[34] <- "CL:0000576"; labels[34] <- "monocyte" # Monocyte ->
# precursor for monoblasts and myeloblasts
verts[35]; labels[35] <- "CFU_GM"
verts[36] <- "CL:0000040"; labels[36] <- "monoblast" # Monoblast ->
verts[37] <- "CL:0000559"; labels[37] <- "promonocyte" # Promonocyte ->
verts[38] <- "CL:0000235"; labels[38] <- "macrophage" # Macrophage ->
verts[39] <- "CL:0000839"; labels[39] <- "myeloid lineage restricted progenitor cell" # CFU_G ->
verts[40]; labels[40] <- "Neutrophilic_Myelocyte"
verts[41] <- "CL:0000775"; labels[41] <- "neutrophil" # Neutrophil ->
# burst forming unit
verts[42] <- "CL:0000038"; labels[42] <- "erythroid progenitor cell" # BFU_E ->
verts[43] <- "CL:0000764"; labels[43] <- "erythroid lineage cell" # CFU_E ->
verts[44] <- "CL:0000547"; labels[44] <- "proerythroblast" # Proerythroblast ->
verts[45]; labels[45] <- "Erythrocyte"
verts[46]; labels[46] <- "BFU_MK"
verts[47] <- "CL:0000553"; labels[47] <- "megakaryocyte progenitor cell" # CFU_MK ->
verts[48] <- "CL:0000556"; labels[48] <- "megakaryocyte" # Megakaryocyte ->
verts[49]; labels[49] <- "Platelets"

hsa <- set_vertex_attr(hsa, "name", value=verts)
hsa <- set_vertex_attr(hsa, "label", value=labels)
# hsa <- set_edge_attr(hsa, "weight", value=rep(1, length(E(hsa))))
visNet <- makeVisNetwork(hsa, hierarchicalLayout=TRUE, direction="LR")

# ==============================================================================
# joining graphs: hematopoietics and hsa

hematopoietics2 <- set_vertex_attr(hematopoietics, "group", value=rep("1", length(V(hematopoietics))))
makeVisNetwork(hematopoietics, cluster=FALSE, layout="layout_with_kk")

hsa2 <- set_vertex_attr(hsa, "group", value=rep(2, length(V(hsa))))
makeVisNetwork(hsa2, cluster=TRUE, hierarchicalLayout=TRUE, direction="LR")

joined <- union(hematopoietics2, hsa2, byname=TRUE)
# cleaning: set all edge weights to 1 (there are some NAs)
joined <- set_edge_attr(joined, "weight", value=rep(1, length(E(joined))))
# Clean vertex labels
V(joined)$label_1[which(is.na(V(joined)$label_1))] <- V(joined)$label_2[which(is.na(V(joined)$label_1))]
joined <- set_vertex_attr(joined, "label", value=V(joined)$label_1)
joined <- delete_vertex_attr(joined, "label_1")
joined <- delete_vertex_attr(joined, "label_2")
# Clean groups, any NAs in group_2 should be 1s
group2nas <- which(is.na(V(joined)$group_2))
V(joined)$group_2[group2nas] <- rep(1, length(group2nas))
joined <- set_vertex_attr(joined, "group", value=V(joined)$group_2)
joined <- delete_vertex_attr(joined, "group_1")
joined <- delete_vertex_attr(joined, "group_2")


as_data_frame(joined, what="edges")
as_data_frame(joined, what="vertices")


makeVisNetwork(joined, customGroups=TRUE, useLabel=TRUE, cluster=FALSE, layout="layout_with_kk")
