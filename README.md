Given a phylogenetic tree and a collection of tips in the tree, this package provides a simple function to find which other tips in the tree would (individually) maximise the phylogenetic diversity of your collection of tips.

Here phylogenetic diversity is based on Faith's measure.

Install with `devtools::install_github("alrichardbollans/rPDtips")`

Example usage:

```R

## Gets a dataframe listing the maximising tips and the diversity gain for each (gains should be equal barring some rounding).
results = rPDtips::max_pd_tips(tree, existing_tips)

# If you want to see what these tips look like you can use a utility function:
plot_max_pd_tips(tree, existing_tips, results)

```
