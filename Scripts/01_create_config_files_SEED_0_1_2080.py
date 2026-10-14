import os
import pandas as pd
import random
import networkx as nx
import glob


CURRENT_DIR = os.getcwd()
CANCER_HALLMARKS_FOLDER = os.path.join(CURRENT_DIR, "../Data/CancerHallmarks")

# Load cancer hallmarks
TISSUE_INVASION_METASTASIS = pd.read_csv(os.path.join(CANCER_HALLMARKS_FOLDER, "TissueInvasionAndMetastasis-UniProt.tsv"), sep="\t", header=0)
REPROGRAMMING_ENERGY_METABOLISM = pd.read_csv(os.path.join(CANCER_HALLMARKS_FOLDER, "ReprogrammingEnergyMetabolism-UniProt.tsv"), sep="\t", header=0)
REPLICATIVE_IMMORTALITY = pd.read_csv(os.path.join(CANCER_HALLMARKS_FOLDER, "ReplicativeImmortality-UniProt.tsv"), sep="\t", header=0)
GENOME_INSTABILITY = pd.read_csv(os.path.join(CANCER_HALLMARKS_FOLDER, "GenomeInstability-UniProt.tsv"), sep="\t", header=0)
RESISTING_CELL_DEATH = pd.read_csv(os.path.join(CANCER_HALLMARKS_FOLDER, "ResistingCellDeath-UniProt.tsv"), sep="\t", header=0)
EVADING_GROWTH_SUPPRESSORS = pd.read_csv(os.path.join(CANCER_HALLMARKS_FOLDER, "EvadingGrowthSuppressors-UniProt.tsv"), sep="\t", header=0)
EVADING_IMMUNE_DESTRUCTION = pd.read_csv(os.path.join(CANCER_HALLMARKS_FOLDER, "EvadingImmuneDestruction-UniProt.tsv"), sep="\t", header=0)
SUSTAINED_ANGIOGENESIS = pd.read_csv(os.path.join(CANCER_HALLMARKS_FOLDER, "SustainedAngiogenesis-UniProt.tsv"), sep="\t", header=0)
SUSTAINING_PROLIFERATIVE_SIGNALING = pd.read_csv(os.path.join(CANCER_HALLMARKS_FOLDER, "SustainingProliferativeSignaling-UniProt.tsv"), sep="\t", header=0)
TUMOR_PROMOTING_INFLAMMATION = pd.read_csv(os.path.join(CANCER_HALLMARKS_FOLDER, "TumorPromotingInflammation-UniProt.tsv"), sep="\t", header=0)

TISSUE_INVASION_METASTASIS_LIST_PROT = TISSUE_INVASION_METASTASIS["Entry"].to_list()
REPROGRAMMING_ENERGY_METABOLISM_LIST_PROT = REPROGRAMMING_ENERGY_METABOLISM["Entry"].to_list()
REPLICATIVE_IMMORTALITY_LIST_PROT = REPLICATIVE_IMMORTALITY["Entry"].to_list()
GENOME_INSTABILITY_LIST_PROT = GENOME_INSTABILITY["Entry"].to_list()
RESISTING_CELL_DEATH_LIST_PROT = RESISTING_CELL_DEATH["Entry"].to_list()
EVADING_GROWTH_SUPPRESSORS_LIST_PROT = EVADING_GROWTH_SUPPRESSORS["Entry"].to_list()
EVADING_IMMUNE_DESTRUCTION_LIST_PROT = EVADING_IMMUNE_DESTRUCTION["Entry"].to_list()
SUSTAINED_ANGIOGENESIS_LIST_PROT = SUSTAINED_ANGIOGENESIS["Entry"].to_list()
SUSTAINING_PROLIFERATIVE_SIGNALING_LIST_PROT = SUSTAINING_PROLIFERATIVE_SIGNALING["Entry"].to_list()
TUMOR_PROMOTING_INFLAMMATION_LIST_PROT = TUMOR_PROMOTING_INFLAMMATION["Entry"].to_list()

list_hallmarks_files = [
    TISSUE_INVASION_METASTASIS,
    REPROGRAMMING_ENERGY_METABOLISM,
    REPLICATIVE_IMMORTALITY,
    GENOME_INSTABILITY,
    RESISTING_CELL_DEATH,
    EVADING_GROWTH_SUPPRESSORS,
    EVADING_IMMUNE_DESTRUCTION,
    SUSTAINED_ANGIOGENESIS,
    SUSTAINING_PROLIFERATIVE_SIGNALING,
    TUMOR_PROMOTING_INFLAMMATION
    ]

HALLMARKS = [
    TISSUE_INVASION_METASTASIS_LIST_PROT,
    REPROGRAMMING_ENERGY_METABOLISM_LIST_PROT,
    REPLICATIVE_IMMORTALITY_LIST_PROT,
    GENOME_INSTABILITY_LIST_PROT,
    RESISTING_CELL_DEATH_LIST_PROT,
    EVADING_GROWTH_SUPPRESSORS_LIST_PROT,
    EVADING_IMMUNE_DESTRUCTION_LIST_PROT,
    SUSTAINED_ANGIOGENESIS_LIST_PROT,
    SUSTAINING_PROLIFERATIVE_SIGNALING_LIST_PROT,
    TUMOR_PROMOTING_INFLAMMATION_LIST_PROT
]

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


df = pd.read_csv("../Data/ReactomeCliques/Reactome96UnweightedComembershipClique.tsv", sep="\t", header=None)
df.columns = ["source", "target"]
Ntw = nx.from_pandas_edgelist(df, source="source", target="target")
print(len(Ntw.nodes()))
list_reactome_proteins = list(Ntw.nodes())
print(len(list_reactome_proteins))

HALLMARKS_IN_REACTOME = dict()

counter = 1
for hallmark_prot_list in HALLMARKS:
    hallmark_name = dico_hallmarks_ids[f"{counter:02d}"]
    prot_in_reactome = list()
    for prot in hallmark_prot_list:
        if prot in list_reactome_proteins:
            prot_in_reactome.append(prot)
    HALLMARKS_IN_REACTOME[hallmark_name] = prot_in_reactome
    counter += 1

print(HALLMARKS_IN_REACTOME)
print(HALLMARKS_IN_REACTOME.keys())


############ CREATE SEEDS AND TARGETS FILES ######################################
random.seed(123456)

for hallmark_name, proteins in HALLMARKS_IN_REACTOME.items():
    random.shuffle(proteins)
    split_index = int(0.2 * len(proteins))
    seeds = proteins[:split_index]
    targets = proteins[split_index:]
    os.makedirs(f"Configs/0_1/SEED_2080/REACTOME_UNWEIGHTED/{hallmark_name}", exist_ok=True)

    with open(f"Configs/0_1/SEED_2080/REACTOME_UNWEIGHTED/{hallmark_name}/{hallmark_name}_seeds.txt", "w") as f:
        for protein in seeds:
            f.write(f"{protein}\n")

    with open(f"Configs/0_1/SEED_2080/REACTOME_UNWEIGHTED/{hallmark_name}/{hallmark_name}_targets.txt", "w") as f:
        for protein in targets:
            f.write(f"{protein}\n")

    os.makedirs(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED/{hallmark_name}", exist_ok=True)

    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED/{hallmark_name}/{hallmark_name}_seeds.txt", "w") as f:
        for protein in seeds:
            f.write(f"{protein}\n")

    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED/{hallmark_name}/{hallmark_name}_targets.txt", "w") as f:
        for protein in targets:
            f.write(f"{protein}\n")

    os.makedirs(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_NO_NSP/{hallmark_name}", exist_ok=True)

    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_NO_NSP/{hallmark_name}/{hallmark_name}_seeds.txt", "w") as f:
        for protein in seeds:
            f.write(f"{protein}\n")

    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_NO_NSP/{hallmark_name}/{hallmark_name}_targets.txt", "w") as f:
        for protein in targets:
            f.write(f"{protein}\n")

    os.makedirs(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_COSINE/{hallmark_name}", exist_ok=True)

    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_COSINE/{hallmark_name}/{hallmark_name}_seeds.txt", "w") as f:
        for protein in seeds:
            f.write(f"{protein}\n")

    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_COSINE/{hallmark_name}/{hallmark_name}_targets.txt", "w") as f:
        for protein in targets:
            f.write(f"{protein}\n")

    os.makedirs(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_JACCARD/{hallmark_name}", exist_ok=True)

    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_JACCARD/{hallmark_name}/{hallmark_name}_seeds.txt", "w") as f:
        for protein in seeds:
            f.write(f"{protein}\n")

    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_JACCARD/{hallmark_name}/{hallmark_name}_targets.txt", "w") as f:
        for protein in targets:
            f.write(f"{protein}\n")

    os.makedirs(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_SUM/{hallmark_name}", exist_ok=True)
    
    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_SUM/{hallmark_name}/{hallmark_name}_seeds.txt", "w") as f:
        for protein in seeds:
            f.write(f"{protein}\n")

    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_SUM/{hallmark_name}/{hallmark_name}_targets.txt", "w") as f:
        for protein in targets:
            f.write(f"{protein}\n")


############# CREATE CONFIG FILES ######################################

for hallmark_id, hallmark_name in dico_hallmarks_ids.items():
    config_file_content = f"""seed: {hallmark_name}_seeds.txt
self_loops: 0
r: 0.1
eta: [1.0]
multiplex:
    1:
        layers:
            - ../../../../../Data/ReactomeCliques/Reactome96UnweightedComembershipClique.tsv
        delta: 0
        graph_type: [00]
        tau: [1]
"""
    with open(f"Configs/0_1/SEED_2080/REACTOME_UNWEIGHTED/{hallmark_name}/{hallmark_name}_config.yml", "w") as f:
        f.write(config_file_content)


for hallmark_id, hallmark_name in dico_hallmarks_ids.items():
    config_file_content = f"""seed: {hallmark_name}_seeds.txt
self_loops: 0
r: 0.1
eta: [1.0]
multiplex:
    1:
        layers:
            - ../../../../../Data/ReactomeCliques/Reactome96WeightedComembershipClique.tsv
        delta: 0
        graph_type: [01]
        tau: [1]
"""
    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED/{hallmark_name}/{hallmark_name}_config.yml", "w") as f:
        f.write(config_file_content)


for hallmark_id, hallmark_name in dico_hallmarks_ids.items():
    config_file_content = f"""seed: {hallmark_name}_seeds.txt
self_loops: 0
r: 0.1
eta: [1.0]
multiplex:
    1:
        layers:
            - ../../../../../Data/ReactomeCliques/Reactome96WeightedComembershipClique_NO_NSP.tsv
        delta: 0
        graph_type: [01]
        tau: [1]
"""
    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_NO_NSP/{hallmark_name}/{hallmark_name}_config.yml", "w") as f:
        f.write(config_file_content)

for hallmark_id, hallmark_name in dico_hallmarks_ids.items():
    config_file_content = f"""seed: {hallmark_name}_seeds.txt
self_loops: 0
r: 0.1
eta: [1.0]
multiplex:
    1:
        layers:
            - ../../../../../Data/ReactomeCliques/Reactome96WeightedComembershipClique_Cosine.tsv
        delta: 0
        graph_type: [01]
        tau: [1]
"""
    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_COSINE/{hallmark_name}/{hallmark_name}_config.yml", "w") as f:
        f.write(config_file_content)

for hallmark_id, hallmark_name in dico_hallmarks_ids.items():
    config_file_content = f"""seed: {hallmark_name}_seeds.txt
self_loops: 0
r: 0.1
eta: [1.0]
multiplex:
    1:
        layers:
            - ../../../../../Data/ReactomeCliques/Reactome96WeightedComembershipClique_Jaccard.tsv
        delta: 0
        graph_type: [01]
        tau: [1]
"""
    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_JACCARD/{hallmark_name}/{hallmark_name}_config.yml", "w") as f:
        f.write(config_file_content)

for hallmark_id, hallmark_name in dico_hallmarks_ids.items():
    config_file_content = f"""seed: {hallmark_name}_seeds.txt
self_loops: 0
r: 0.1
eta: [1.0]
multiplex:
    1:
        layers:
            - ../../../../../Data/ReactomeCliques/Reactome96WeightedComembershipClique_SUM.tsv
        delta: 0
        graph_type: [01]
        tau: [1]
"""
    with open(f"Configs/0_1/SEED_2080/REACTOME_WEIGHTED_SUM/{hallmark_name}/{hallmark_name}_config.yml", "w") as f:
        f.write(config_file_content)
