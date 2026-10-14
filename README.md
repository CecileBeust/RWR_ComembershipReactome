## Random Walk With Restart on pathway co-membership cliques of Reactome (version 96)

The co-membership cliques of Reactome can be downloaded at: https://zenodo.org/records/23372223

The networks in the repository are:

- `Reactome96UnweightedComembershipClique.tsv` : The unweighted pathway co-membership clique of Reactome (all proteins linked together in a clique)
- `Reactome96WeightedComembershipClique.tsv` : The pathway co-membership clique of Reactome scored with our method (that takes into account both the hierarchy and sequentiality of pathways)
- `Reactome96WeightedComembershipClique_Cosine.tsv` : The pathway co-membership clique of Reactome scored with a Cosine similarity method.
- `Reactome96WeightedComembershipClique_ Jaccard.tsv` : The pathway co-membership clique of Reactome scored with a Jaccard similarity method.
- `Reactome96WeightedComembershipClique_NO_NSP` : The pathway co-membership clique of Reactome scored with our method, but without taking into account the sequentiality (No NextStepPathway = NO NSP)
- `Reactome96WeightedComembershipClique_SUM.tsv` : The pathway co-membership clique of Reactome scored with a variant of our method (taking the sum of hierarchical and sequential scores insted of the max)