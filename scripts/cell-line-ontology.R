# R bindings to C based network analysis package igraph
if (!require(igraph, quietly=TRUE)) install.packages("igraph")

nicePath <- function (filename) {
  tryCatch({
    return(file.path(getwd(), dirname(sys.frame(1)$ofile), filename))
  }, error = function (error) {
    return(filename)
  })
}

clo <- read.csv(nicePath("../data/CLO.csv"), stringsAsFactors=FALSE)

# A bunch of columns are missing from many rows
# filledCols <- c()
# for (col in names(clo)) {
#   filledCols <- c(filledCols, all(!is.na(clo[col])))
# }
#
# notNas <- names(clo)[filledCols]
# notEmpty <- c()
# for (col in as.list(notNas)) {
#   notEmpty <- c(notEmpty, !any(clo[col] == ''))
# }
#
# notNas[notEmpty]

edges <- clo[c("Preferred.Label", "Parents")]
colnames(edges) <- c("from", "to")

# Give those with no parent our own custom parent
empties <- which(edges$to == '')
for (row in as.list(empties)) {
  edges$to[row] = "/JM_00000001"
}
# check:
# edges[empties,]

# Now split the parent URL into an ID
edges$to <- unlist(lapply(strsplit(as.character(edges$to), "/"), function(x) tail(x, n=1)))

# How many unique IDs?
length(unique(edges$from)) # 38600
# and parents?
length(unique(edges$to)) # 1070

G <- graph_from_data_frame(edges)

# Not everything is ^CLO:
length(grep("^CLO", as_ids(V(G))))

plot(G, vertex.size=0.01, vertex.label=NA, edge.arrow.width=0)

as_ids(V(G))[-grep("^CLO", as_ids(V(G)))]
