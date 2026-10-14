import multixrank
import os
import pandas as pd
import random
import networkx as nx
import numpy as np

dico_hallmarks_ids = {
    "01": "TISSUE_INVASION_METASTASIS",
    "02": "REPROGRAMMING_ENERGY_METABOLISM",
    "03": "REPLICATIVE_IMMORTALITY",
    "04": "GENOME_INSTABILITY",
    "05": "RESISTING_CELL_DEATH",
    "06": "EVADING_GROWTH_SUPPRESSORS",
    "07": "EVADING_IMMUNE_DESTRUCTION",
    "08": "SUSTAINED_ANGIOGENESIS",
    "09": "SUSTAINING_PROLIFERATIVE_SIGNALING",
    "10": "TUMOR_PROMOTING_INFLAMMATION"
}

# dico_hallmarks_ids = {
#     "01": "TISSUE_INVASION_METASTASIS"
# }

############## UNWEIGHTED #########################
for hallmark_id, hallmark_name in dico_hallmarks_ids.items():
    multixrank_obj = multixrank.Multixrank(config=f"../Configs/0_1/SEED/REACTOME_UNWEIGHTED/{hallmark_name}/{hallmark_name}_config.yml", wdir=f"../Configs/0_1/SEED/REACTOME_UNWEIGHTED/{hallmark_name}")
    ranking_df = multixrank_obj.random_walk_rank()
    os.makedirs(f"../Outputs/0_1/SEED/REACTOME_UNWEIGHTED/{hallmark_name}", exist_ok=True)
    multixrank_obj.write_ranking(ranking_df, path=f"../Outputs/0_1/SEED/REACTOME_UNWEIGHTED/{hallmark_name}")

############### WEIGHTED #########################
for hallmark_id, hallmark_name in dico_hallmarks_ids.items():
    multixrank_obj = multixrank.Multixrank(config=f"../Configs/0_1/SEED/REACTOME_WEIGHTED/{hallmark_name}/{hallmark_name}_config.yml", wdir=f"../Configs/0_1/SEED/REACTOME_WEIGHTED/{hallmark_name}")
    ranking_df = multixrank_obj.random_walk_rank()
    os.makedirs(f"../Outputs/0_1/SEED/REACTOME_WEIGHTED/{hallmark_name}", exist_ok=True)
    multixrank_obj.write_ranking(ranking_df, path=f"../Outputs/0_1/SEED/REACTOME_WEIGHTED/{hallmark_name}")

############### WEIGHTED NO NSP #########################
for hallmark_id, hallmark_name in dico_hallmarks_ids.items():
    multixrank_obj = multixrank.Multixrank(config=f"../Configs/0_1/SEED/REACTOME_WEIGHTED_NO_NSP/{hallmark_name}/{hallmark_name}_config.yml", wdir=f"../Configs/0_1/SEED/REACTOME_WEIGHTED_NO_NSP/{hallmark_name}")
    ranking_df = multixrank_obj.random_walk_rank()
    os.makedirs(f"../Outputs/0_1/SEED/REACTOME_WEIGHTED_NO_NSP/{hallmark_name}", exist_ok=True)
    multixrank_obj.write_ranking(ranking_df, path=f"../Outputs/0_1/SEED/REACTOME_WEIGHTED_NO_NSP/{hallmark_name}")

############### WEIGHTED COSINE #########################
for hallmark_id, hallmark_name in dico_hallmarks_ids.items():
    multixrank_obj = multixrank.Multixrank(config=f"../Configs/0_1/SEED/REACTOME_WEIGHTED_COSINE/{hallmark_name}/{hallmark_name}_config.yml", wdir=f"../Configs/0_1/SEED/REACTOME_WEIGHTED_COSINE/{hallmark_name}")
    ranking_df = multixrank_obj.random_walk_rank()
    os.makedirs(f"../Outputs/0_1/SEED/REACTOME_WEIGHTED_COSINE/{hallmark_name}", exist_ok=True)
    multixrank_obj.write_ranking(ranking_df, path=f"../Outputs/0_1/SEED/REACTOME_WEIGHTED_COSINE/{hallmark_name}")

############### WEIGHTED JACCARD #########################
for hallmark_id, hallmark_name in dico_hallmarks_ids.items():
    multixrank_obj = multixrank.Multixrank(config=f"../Configs/0_1/SEED/REACTOME_WEIGHTED_JACCARD/{hallmark_name}/{hallmark_name}_config.yml", wdir=f"../Configs/0_1/SEED/REACTOME_WEIGHTED_JACCARD/{hallmark_name}")
    ranking_df = multixrank_obj.random_walk_rank()
    os.makedirs(f"../Outputs/0_1/SEED/REACTOME_WEIGHTED_JACCARD/{hallmark_name}", exist_ok=True)
    multixrank_obj.write_ranking(ranking_df, path=f"../Outputs/0_1/SEED/REACTOME_WEIGHTED_JACCARD/{hallmark_name}")

############### WEIGHTED SUM #########################
for hallmark_id, hallmark_name in dico_hallmarks_ids.items():
    multixrank_obj = multixrank.Multixrank(config=f"../Configs/0_1/SEED/REACTOME_WEIGHTED_SUM/{hallmark_name}/{hallmark_name}_config.yml", wdir=f"../Configs/0_1/SEED/REACTOME_WEIGHTED_SUM/{hallmark_name}")
    ranking_df = multixrank_obj.random_walk_rank()
    os.makedirs(f"../Outputs/0_1/SEED/REACTOME_WEIGHTED_SUM/{hallmark_name}", exist_ok=True)
    multixrank_obj.write_ranking(ranking_df, path=f"../Outputs/0_1/SEED/REACTOME_WEIGHTED_SUM/{hallmark_name}")