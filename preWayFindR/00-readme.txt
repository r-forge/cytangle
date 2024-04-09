This folder contains source material for some of the constants that
are part of the WayFindR package.

# History
We started by downloading the complete database of Homo sapiens
WikiPathways, and storing them in a subfolder. We ran two perl scripts
against the entire database:

1. `wpnode.pl` finds all `DataNode` items in the GPML files, and
   counts how many times each type is used. These results are stored
   in the file `nodeTypes.txt`. That file is in tab-spearated-values
   format with two columns. The first column is the text entry in the
   `Type` attribute of the `DataNode`, and the second column is the
   number of times we saw it.
2. `whipple.pl` does the corresponding task for edges. In a GPML file,
   an edge is actually called an `Interaction`. The kind of edge is
   recorded as a text-filled atribute called `ArrowHead` stored in the
   second `Point` in the `Graphics` item inside the `Interaction`. The
   results of this pass are stored in the file `arrowTypes.txt`, but
   we only keep the names, not the counts. (The script also produces a
   file called `allArrows.txt`, which contains more detailed counts
   stored in a matrix of pathway X arrowType.

# Dispaly Items
Next, we ran an R script against each of the output files to assign
colors (to vertices or edges) and line types (to edges).

1. `nodeColor.R` assigns colors to the nodes or vertices, and stores
   a vactor called `nodeColors` and a data frame called `nodeTypes` in
   the data directory as `nodeColor.rda`.
2. `edgeColor.R` simplifies the edge types (to remove redundancies int
   he WikiPathways type list) and assigns both colors and line types
   to the reminiang simpified list of edge types. It stores
   `edgeTypes` and `edgeColors` in the daat directory as `edgeStyle/rda`.
