library(ape)

#' Tips whose individual addition to `species` gives the largest increase in
#' unrooted PD, i.e. sum(keep.tip(tree, species)$edge.length).
#'
#' Assumes a rooted, bifurcating, ultrametric tree.
#'
#' Let h(v) be the height of node v above the tips and m the MRCA of `species`.
#'  * m is not the root: the best tips are all tips in the root child clade that
#'    contains no member of `species`. Gain = 2*h(root) - h(m).
#'  * m is the root: every covered node r with exactly one covered child offers
#'    the uncovered child clade as candidates, with gain h(r). Keep the r with
#'    the largest h(r).
#'
#' Returns a data.frame(tip, gain). Zero rows if `species` is every tip.
#'
#' @param tree The phylogenetic tree to use
#' @param species The list of existing species to use
#' @return a data.frame(tip, gain). Zero rows if `species` is every tip.


max_pd_tips <- function(tree, species, tol = 1e-8) {

  # ---- 1. Validate inputs ---------------------------------------------------
  stopifnot(inherits(tree, "phylo"))
  if (!is.rooted(tree))      stop("tree must be rooted")
  if (!is.binary(tree))      stop("tree must be bifurcating")
  if (!is.ultrametric(tree)) stop("tree must be ultrametric")

  selected_tips   <- unique(as.character(species))
  n_selected      <- length(selected_tips)
  if (n_selected == 0) stop("`species` must contain at least one tip")

  labels_not_in_tree <- setdiff(selected_tips, tree$tip.label)
  if (length(labels_not_in_tree)) {
    stop("Tips not in tree: ", paste(labels_not_in_tree, collapse = ", "))
  }

  n_tips <- Ntip(tree)
  if (n_selected == n_tips) {
    # Nothing left to add
    return(data.frame(tip = character(0), gain = numeric(0)))
  }

  # ---- 2. Node heights above the tips ---------------------------------------
  # node.depth.edgelength() gives distance from the root; on an ultrametric
  # tree, (max distance - distance) is the height above the tips.
  # Nodes are numbered 1..n_tips (tips), then n_tips + 1 (the root) upwards.
  distance_from_root <- node.depth.edgelength(tree)
  node_height        <- max(distance_from_root) - distance_from_root

  # ---- 3. Tree structure ----------------------------------------------------
  # Postorder edge ordering lists every edge below a node before the edge
  # above it, so counts can be accumulated from the tips upward in one loop.
  tree     <- reorder(tree, "postorder")
  edges    <- tree$edge                 # column 1 = parent node, column 2 = child node
  root_node <- n_tips + 1L

  # children[i, ] holds the two children of internal node (n_tips + i).
  # Sorting edges by parent puts each internal node's two edges next to each
  # other, in the order of internal node number.
  children <- matrix(edges[order(edges[, 1]), 2], ncol = 2, byrow = TRUE)
  child_1  <- children[, 1]
  child_2  <- children[, 2]

  # ---- 4. Count selected tips below every node ------------------------------
  # n_selected_below[v] = number of selected tips at or below node v.
  # A node is "covered" if this is > 0.
  n_selected_below <- numeric(n_tips + tree$Nnode)
  n_selected_below[match(selected_tips, tree$tip.label)] <- 1
  for (edge_row in seq_len(nrow(edges))) {
    parent_node <- edges[edge_row, 1]
    child_node  <- edges[edge_row, 2]
    n_selected_below[parent_node] <-
      n_selected_below[parent_node] + n_selected_below[child_node]
  }

  # ---- 5. Find the MRCA node of the selected tips -----------------------------
  # Start at the root and keep stepping into the child that contains ALL
  # selected tips.
  if (length(species)==1){
    tip_node <- which(tree$tip.label == species[1])

    # 2. Get the parent node number using getParent
    mrca_node = phytools::getParent(tree, tip_node)
  }else{
    mrca_node <- ape::getMRCA(tree, species)
  }


  # ---- 6. Work out which clades hold the best new tips ------------------------
  if (mrca_node != root_node) {
    # Case A: MRCA is below the root. The root child with no selected tips
    # (row 1 of `children` is the root) is the best attachment target.
    root_children       <- children[1, ]
    uncovered_clades    <- root_children[n_selected_below[root_children] == 0]
    best_gain           <- 2 * node_height[root_node] - node_height[mrca_node]

  } else {
    # Case B: MRCA is the root, so the whole spanning subtree is covered.
    # Find internal nodes where exactly one child contains selected tips:
    # the other child is an uncovered clade that attaches at this node.
    child_1_covered <- n_selected_below[child_1] > 0
    child_2_covered <- n_selected_below[child_2] > 0
    attachment_rows <- which(xor(child_1_covered, child_2_covered))

    # Gain from attaching at a node is that node's height (path down to the
    # new tip). Keep only the highest attachment nodes, allowing ties.
    attachment_heights <- node_height[n_tips + attachment_rows]
    best_gain          <- max(attachment_heights)
    attachment_rows    <- attachment_rows[attachment_heights >= best_gain - tol]

    uncovered_clades <- ifelse(child_1_covered[attachment_rows],
                               child_2[attachment_rows],   # child 1 covered, so child 2 is the new clade
                               child_1[attachment_rows])
  }

  # ---- 7. Convert clades (node numbers) to tip labels --------------------------
  # clade_tip_indices[[i]] lists the tip numbers below internal node (n_tips + i).
  clade_tip_indices <- prop.part(tree)
  best_tip_labels <- unlist(lapply(uncovered_clades, function(clade_node) {
    if (clade_node <= n_tips) {
      tree$tip.label[clade_node]                                   # the clade is a single tip
    } else {
      tree$tip.label[clade_tip_indices[[clade_node - n_tips]]]
    }
  }))

  data.frame(tip = best_tip_labels, gain = best_gain, stringsAsFactors = FALSE)
}

#' Brute-force reference using the keep.tip definition of PD.
#' PD of a single tip is defined as 0.
#' @param tree The phylogenetic tree to use
#' @param species The list of existing species to use
#' @return a data.frame(tip, gain).
max_pd_tips_brute <- function(tree, species, tol = 1e-8) {
  pd <- function(sp) sum(keep.tip(tree, sp)$edge.length)
  base <- pd(species)
  cand <- setdiff(tree$tip.label, species)
  if (length(cand) == 0) return(data.frame(tip = character(0), gain = numeric(0)))
  gain <- vapply(cand, function(t) pd(c(species, t)) - base, numeric(1))
  data.frame(tip = cand[gain >= max(gain)- tol], gain = max(gain),
             stringsAsFactors = FALSE)
}


#' Plot the tree with the original species and the PD-maximising additions
#' highlighted. Computes max_pd_tips() unless `result` is supplied.
#'
#' @param tree The phylogenetic tree to use
#' @param species The list of existing species to use
#' @param result   optional output of max_pd_tips(), to avoid recomputing
#' @param col_species colour for the original species
#' @param col_new     colour for the PD-maximising additions
#' @param col_other   colour for all remaining tips
#' @param legend_pos  position passed to legend(); NULL for no legend
#' @param plot_length whether to include branch lengths on graphs
#' @param ...      further arguments passed to ape's plot.phylo (e.g. type = "fan")
#' Returns the vector of tip colours invisibly.
plot_max_pd_tips <- function(tree, species, result = NULL,
                             col_species = "blue", col_new = "red",
                             col_other = "grey60", cex = 0.7,
                             legend_pos = "bottomleft",plot_length=FALSE, ...) {
  if (is.null(result)) result <- max_pd_tips(tree, species)
  species  <- unique(as.character(species))
  new_tips <- as.character(result$tip)

  tip_col <- rep(col_other, Ntip(tree))
  tip_col[tree$tip.label %in% species]  <- col_species
  tip_col[tree$tip.label %in% new_tips] <- col_new

  plot(tree, tip.color = tip_col, cex = cex, ...)



  hl <- which(tip_col != col_other)
  if (length(hl)) tiplabels(pch = 16, tip = hl, col = tip_col[hl], cex = cex * 1.5)

  if (!is.null(legend_pos)) {
    legend(legend_pos, bty = "n", pch = 16, cex = 0.8,
           col    = c(col_species, col_new, col_other),
           legend = c("Selected species", "PD-maximising additions", "Other tips"))
  }
  invisible(tip_col)
  if (plot_length){
    edgelabels(round(tree$edge.length,3),col="black",bg=NULL, font=1, offset =10)
  }
}

