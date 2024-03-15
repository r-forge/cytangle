library(igraph)
load("wp3850.Rda")

########## Graph Properties ############
## Global Metrics
any_loop(G)
any_multiple(G)
edge_connectivity(G) # aka adhesion
triangles(G)
mean_distance(G)
automorphisms(G)
articulation_points(G)
is_connected(G)
clique_num(G)
largest_cliques(G)
vertex_connectivity(G) # aka cohesion
diameter(G)
gd <- get_diameter(G)
gd$label
dyad_census(G)
triad_census(G) # see man page to interpret triads
has_eulerian_cycle(G)
has_eulerian_path(G)
gg <- girth(G) # length of the shortest circle
gg[[2]]$label
global_efficiency(G)
is_chordal(G)
is_directed(G)
is_dag(G)
radius(G)
reciprocity(G)

# Per-Vertex Metrics
sort(alpha_centrality(G))
adjacent.triangles(G)
count_triangles(G)
A <- authority_score(G)
plot(sort(A$vector), pch=16)
H  <- hub_score(G)
plot(sort(H$vector), pch=16)
P <- page_rank(G)
plot(sort(P$vector), pch=16)
E <- eigen_centrality(G)
plot(sort(E$vector), pch=16)
K <- knn(G)
plot(sort(K$knn), pch = 16)
sort(betweenness(G))
sort(closeness(G))
sort(coreness(G))
sort(degree(G))
sort(eccentricity(G))
sort(strength(G))
sort(harmonic_centrality(G))
sort(local_efficiency(G))
sort(local_scan(G))
sort(power_centrality(G))

DM <- distances(G)
dim(DM)
sym <- as.dist(DM)
plot(hc <- hclust(sym, "ward.D2"))
heatmap(DM)

LM <- laplacian_matrix(G)
dim(LM)

which(V(G)$label == "TNF")
all_shortest_paths(G, 33)
all_shortest_paths(G, 37)

# Per-Edge Metrics
edge_betweenness(G)

## Communities
club <- cluster_edge_betweenness(G)
length(club)
sizes(club)
modularity(club)
algorithm(club)
is.hierarchical(club)
plot(as.dendrogram(club))
crossing(club, G)
plot(club, G)

## Other stuff
ML <- make_line_graph(G)
plot(ML)
motifs(G, size = 3) # Don't know how to interpret this answer
motifs(G, size = 4)
count_motifs(G, size = 3)
count_motifs(G, size = 4)
count_motifs(G, size = 5)

## look at the 'contract' function. Can use this to simplfy groups/complexes
## also look at 'simplfied'
##
## max_flow
## Undirected graph only
## spectrum(G) # undirected only
##
## Really slow
## ivs(G) # really slow

############ Cycles #####################
## Online freference:
## https://stackoverflow.com/questions/55091438/r-igraph-find-all-cycles
Cycles = NULL
for(v1 in V(G)) {
    for(v2 in neighbors(G, v1, mode="out")) {
        Cycles = c(Cycles, 
            lapply(all_simple_paths(G, v2, v1, mode="out"), function(p) c(v1,p)))
    }
}
length(Cycles)
table(len <- sapply(Cycles, length))
Cycles[len == 3]
Cycles[len == 5]

UniqueCycles <- Cycles[sapply(Cycles, min) == sapply(Cycles, `[`, 1)]
length(UniqueCycles)
UniqueCycles

interpret <- function(v) {
   sapply(names(v), function(N) {
    if (N == "") return("")
    V(G)$label[V(G)$name == N]
  })
}
interpret(UniqueCycles[[1]])
lapply(UniqueCycles, interpret)

track <- function(v) {
  E(G)$MIM[v]
}
lapply(UniqueCycles, track)

bully <- function(v) {
   genes <- sapply(names(v), function(N) {
    if (N == "") return("")
    V(G)$label[V(G)$name == N]
   })
   arrows <- E(G)$MIM[v]
   data.frame(genes, arrows)
}

lapply(UniqueCycles, bully)

uu <- unique(unlist(sapply(UniqueCycles, names)))
uu <- uu[uu != ""]
SG <- subgraph(G, uu)
set.seed(54321)
L <- layout_nicely(SG)
windows(width=9, height=9)
opar <- par(mai = c(0.1, 0.1, 1, 0.1), bg = "white")
plot(SG, layout = L)
legend("bottomleft", names(atypes), lwd = 4, lty = atypes, col = acols)
legend("topright", names(punk), pch = 15, col = punk)
title("Cycles in the IGF1-AKT pathway")
par(opar)
resn <- 300
dev.copy(png, file = "IGF1-Cycles.png", width = 12*resn, height = 12*resn, res = resn, bg = "white")
dev.off()

data.frame(M=E(SG)$MIM,
           C=E(SG)$color,
           L=E(SG)$lty)
data.frame(acols, atypes)
