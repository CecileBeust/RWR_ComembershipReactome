import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import auc

hallmark_id = "01"
hallmark_name = "TISSUE_INVASION_METASTASIS"

# --- UNWEIGHTED ---
seeds_file_weighted_1 = f"../Configs/0_0_5/SEED/REACTOME_UNWEIGHTED/{hallmark_name}/{hallmark_name}_seeds.txt"
targets_file_weighted_1 = f"../Configs/0_0_5/SEED/REACTOME_UNWEIGHTED/{hallmark_name}/{hallmark_name}_targets.txt"
scores_file_weighted_1 = f"../Outputs/0_0_5/SEED/REACTOME_UNWEIGHTED/{hallmark_name}/multiplex_1.tsv"

df_scores_weighted_1 = pd.read_csv(scores_file_weighted_1, sep='\t')
with open(seeds_file_weighted_1, 'r') as f:
    seeds_1 = set(line.strip() for line in f if line.strip())
with open(targets_file_weighted_1, 'r') as f:
    targets_1 = set(line.strip() for line in f if line.strip())
df_scores_weighted_1['is_target'] = df_scores_weighted_1['node'].isin(targets_1).astype(int)
df_scores_weighted_1['is_seed'] = df_scores_weighted_1['node'].isin(seeds_1).astype(int)
df_scores_weighted_1 = df_scores_weighted_1.sort_values('score', ascending=False).reset_index(drop=True)

n_targets_1 = df_scores_weighted_1['is_target'].sum()
n_total_1 = len(df_scores_weighted_1)
n_non_targets_1 = n_total_1 - n_targets_1
cumulative_targets_1 = np.cumsum(df_scores_weighted_1['is_target'].values)
cumulative_non_targets_1 = np.cumsum((1 - df_scores_weighted_1['is_target']).values)
positions_1 = np.arange(len(df_scores_weighted_1)) + 1
positions_fraction_1 = positions_1 / len(df_scores_weighted_1)
tpr_weighted_1 = cumulative_targets_1 / n_targets_1
fpr_weighted_1 = cumulative_non_targets_1 / n_non_targets_1
auc_score_weighted_1 = auc(fpr_weighted_1, tpr_weighted_1)

# --- WEIGHTED ---
seeds_file_weighted_2 = f"../Configs/0_0_5/SEED/REACTOME_WEIGHTED/{hallmark_name}/{hallmark_name}_seeds.txt"
targets_file_weighted_2 = f"../Configs/0_0_5/SEED/REACTOME_WEIGHTED/{hallmark_name}/{hallmark_name}_targets.txt"
scores_file_weighted_2 = f"../Outputs/0_0_5/SEED/REACTOME_WEIGHTED/{hallmark_name}/multiplex_1.tsv"

df_scores_weighted_2 = pd.read_csv(scores_file_weighted_2, sep='\t')
with open(seeds_file_weighted_2, 'r') as f:
    seeds_2 = set(line.strip() for line in f if line.strip())
with open(targets_file_weighted_2, 'r') as f:
    targets_2 = set(line.strip() for line in f if line.strip())
df_scores_weighted_2['is_target'] = df_scores_weighted_2['node'].isin(targets_2).astype(int)
df_scores_weighted_2['is_seed'] = df_scores_weighted_2['node'].isin(seeds_2).astype(int)
df_scores_weighted_2 = df_scores_weighted_2.sort_values('score', ascending=False).reset_index(drop=True)

n_targets_2 = df_scores_weighted_2['is_target'].sum()
n_total_2 = len(df_scores_weighted_2)
n_non_targets_2 = n_total_2 - n_targets_2
cumulative_targets_2 = np.cumsum(df_scores_weighted_2['is_target'].values)
cumulative_non_targets_2 = np.cumsum((1 - df_scores_weighted_2['is_target']).values)
positions_2 = np.arange(len(df_scores_weighted_2)) + 1
positions_fraction_2 = positions_2 / len(df_scores_weighted_2)
tpr_weighted_2 = cumulative_targets_2 / n_targets_2
fpr_weighted_2 = cumulative_non_targets_2 / n_non_targets_2
auc_score_weighted_2 = auc(fpr_weighted_2, tpr_weighted_2)

# --- WEIGHTED NO NSP ---
seeds_file_weighted_3 = f"../Configs/0_0_5/SEED/REACTOME_WEIGHTED_NO_NSP/{hallmark_name}/{hallmark_name}_seeds.txt"
targets_file_weighted_3 = f"../Configs/0_0_5/SEED/REACTOME_WEIGHTED_NO_NSP/{hallmark_name}/{hallmark_name}_targets.txt"
scores_file_weighted_3 = f"../Outputs/0_0_5/SEED/REACTOME_WEIGHTED_NO_NSP/{hallmark_name}/multiplex_1.tsv"

df_scores_weighted_3 = pd.read_csv(scores_file_weighted_3, sep='\t')
with open(seeds_file_weighted_3, 'r') as f:
    seeds_3 = set(line.strip() for line in f if line.strip())
with open(targets_file_weighted_3, 'r') as f:
    targets_3 = set(line.strip() for line in f if line.strip())
df_scores_weighted_3['is_target'] = df_scores_weighted_3['node'].isin(targets_3).astype(int)
df_scores_weighted_3['is_seed'] = df_scores_weighted_3['node'].isin(seeds_3).astype(int)
df_scores_weighted_3 = df_scores_weighted_3.sort_values('score', ascending=False).reset_index(drop=True)

n_targets_3 = df_scores_weighted_3['is_target'].sum()
n_total_3 = len(df_scores_weighted_3)
n_non_targets_3 = n_total_3 - n_targets_3
cumulative_targets_3 = np.cumsum(df_scores_weighted_3['is_target'].values)
cumulative_non_targets_3 = np.cumsum((1 - df_scores_weighted_3['is_target']).values)
positions_3 = np.arange(len(df_scores_weighted_3)) + 1
positions_fraction_3 = positions_3 / len(df_scores_weighted_3)
tpr_weighted_3 = cumulative_targets_3 / n_targets_3
fpr_weighted_3 = cumulative_non_targets_3 / n_non_targets_3
auc_score_weighted_3 = auc(fpr_weighted_3, tpr_weighted_3)

# --- WEIGHTED COSINE ---
seeds_file_weighted_4 = f"../Configs/0_0_5/SEED/REACTOME_WEIGHTED_COSINE/{hallmark_name}/{hallmark_name}_seeds.txt"
targets_file_weighted_4 = f"../Configs/0_0_5/SEED/REACTOME_WEIGHTED_COSINE/{hallmark_name}/{hallmark_name}_targets.txt"
scores_file_weighted_4 = f"../Outputs/0_0_5/SEED/REACTOME_WEIGHTED_COSINE/{hallmark_name}/multiplex_1.tsv"

df_scores_weighted_4 = pd.read_csv(scores_file_weighted_4, sep='\t')
with open(seeds_file_weighted_4, 'r') as f:
    seeds_4 = set(line.strip() for line in f if line.strip())
with open(targets_file_weighted_4, 'r') as f:
    targets_4 = set(line.strip() for line in f if line.strip())
df_scores_weighted_4['is_target'] = df_scores_weighted_4['node'].isin(targets_4).astype(int)
df_scores_weighted_4['is_seed'] = df_scores_weighted_4['node'].isin(seeds_4).astype(int)
df_scores_weighted_4 = df_scores_weighted_4.sort_values('score', ascending=False).reset_index(drop=True)

n_targets_4 = df_scores_weighted_4['is_target'].sum()
n_total_4 = len(df_scores_weighted_4)
n_non_targets_4 = n_total_4 - n_targets_4
cumulative_targets_4 = np.cumsum(df_scores_weighted_4['is_target'].values)
cumulative_non_targets_4 = np.cumsum((1 - df_scores_weighted_4['is_target']).values)
positions_4 = np.arange(len(df_scores_weighted_4)) + 1
positions_fraction_4 = positions_4 / len(df_scores_weighted_4)
tpr_weighted_4 = cumulative_targets_4 / n_targets_4
fpr_weighted_4 = cumulative_non_targets_4 / n_non_targets_4
auc_score_weighted_4 = auc(fpr_weighted_4, tpr_weighted_4)

# --- WEIGHTED JACCARD ---
seeds_file_weighted_5 = f"../Configs/0_0_5/SEED/REACTOME_WEIGHTED_JACCARD/{hallmark_name}/{hallmark_name}_seeds.txt"
targets_file_weighted_5 = f"../Configs/0_0_5/SEED/REACTOME_WEIGHTED_JACCARD/{hallmark_name}/{hallmark_name}_targets.txt"
scores_file_weighted_5 = f"../Outputs/0_0_5/SEED/REACTOME_WEIGHTED_JACCARD/{hallmark_name}/multiplex_1.tsv"

df_scores_weighted_5 = pd.read_csv(scores_file_weighted_5, sep='\t')
with open(seeds_file_weighted_5, 'r') as f:
    seeds_5 = set(line.strip() for line in f if line.strip())
with open(targets_file_weighted_5, 'r') as f:
    targets_5 = set(line.strip() for line in f if line.strip())
df_scores_weighted_5['is_target'] = df_scores_weighted_5['node'].isin(targets_5).astype(int)
df_scores_weighted_5['is_seed'] = df_scores_weighted_5['node'].isin(seeds_5).astype(int)
df_scores_weighted_5 = df_scores_weighted_5.sort_values('score', ascending=False).reset_index(drop=True)

n_targets_5 = df_scores_weighted_5['is_target'].sum()
n_total_5 = len(df_scores_weighted_5)
n_non_targets_5 = n_total_5 - n_targets_5
cumulative_targets_5 = np.cumsum(df_scores_weighted_5['is_target'].values)
cumulative_non_targets_5 = np.cumsum((1 - df_scores_weighted_5['is_target']).values)
positions_5 = np.arange(len(df_scores_weighted_5)) + 1
positions_fraction_5 = positions_5 / len(df_scores_weighted_5)
tpr_weighted_5 = cumulative_targets_5 / n_targets_5
fpr_weighted_5 = cumulative_non_targets_5 / n_non_targets_5
auc_score_weighted_5 = auc(fpr_weighted_5, tpr_weighted_5)

# --- WEIGHTED SUM ---
seeds_file_weighted_6 = f"../Configs/0_0_5/SEED/REACTOME_WEIGHTED_SUM/{hallmark_name}/{hallmark_name}_seeds.txt"
targets_file_weighted_6 = f"../Configs/0_0_5/SEED/REACTOME_WEIGHTED_SUM/{hallmark_name}/{hallmark_name}_targets.txt"
scores_file_weighted_6 = f"../Outputs/0_0_5/SEED/REACTOME_WEIGHTED_SUM/{hallmark_name}/multiplex_1.tsv"

df_scores_weighted_6 = pd.read_csv(scores_file_weighted_6, sep='\t')
with open(seeds_file_weighted_6, 'r') as f:
    seeds_6 = set(line.strip() for line in f if line.strip())
with open(targets_file_weighted_6, 'r') as f:
    targets_6 = set(line.strip() for line in f if line.strip())
df_scores_weighted_6['is_target'] = df_scores_weighted_6['node'].isin(targets_6).astype(int)
df_scores_weighted_6['is_seed'] = df_scores_weighted_6['node'].isin(seeds_6).astype(int)
df_scores_weighted_6 = df_scores_weighted_6.sort_values('score', ascending=False).reset_index(drop=True)

n_targets_6 = df_scores_weighted_6['is_target'].sum()
n_total_6 = len(df_scores_weighted_6)
n_non_targets_6 = n_total_6 - n_targets_6
cumulative_targets_6 = np.cumsum(df_scores_weighted_6['is_target'].values)
cumulative_non_targets_6 = np.cumsum((1 - df_scores_weighted_6['is_target']).values)
positions_6 = np.arange(len(df_scores_weighted_6)) + 1
positions_fraction_6 = positions_6 / len(df_scores_weighted_6)
tpr_weighted_6 = cumulative_targets_6 / n_targets_6
fpr_weighted_6 = cumulative_non_targets_6 / n_non_targets_6
auc_score_weighted_6 = auc(fpr_weighted_6, tpr_weighted_6)

# Création d'un seul graphe
plt.figure(figsize=(10, 6))

plt.plot(positions_fraction_6, tpr_weighted_6, '--', linewidth=2.5,
        label=f'Weighted SUM AUC = {auc_score_weighted_6:.2f}', color="#ea5ce1")
plt.plot(positions_fraction_2, tpr_weighted_2, '--', linewidth=2.5,
        label=f'Weighted MAX AUC = {auc_score_weighted_2:.2f}', color="#f30d28")
plt.plot(positions_fraction_3, tpr_weighted_3, '--', linewidth=2.5,
        label=f'Weighted NO NSP AUC = {auc_score_weighted_3:.2f}', color="#ea9c30")
plt.plot(positions_fraction_4, tpr_weighted_4, '--', linewidth=2.5,
        label=f'Weighted Cosine AUC = {auc_score_weighted_4:.2f}', color="#408d09")
plt.plot(positions_fraction_5, tpr_weighted_5, '--', linewidth=2.5,
        label=f'Weighted Jaccard AUC = {auc_score_weighted_5:.2f}', color="#2f7bd3")
plt.plot(positions_fraction_1, tpr_weighted_1, '--', linewidth=2.5,
        label=f'Unweighted AUC = {auc_score_weighted_1:.2f}', color="#9e32e2")
plt.plot([0, 1], [0, 1], ':', linewidth=2, label='Random', color='#7f7f7f')

plt.xlabel('RWR scores (ranked by decreasing order)')
plt.ylabel('Targets found')
plt.title(f'{hallmark_name}, r = 0.05')
plt.legend(loc="lower right")
plt.grid(False)

plt.tight_layout()
plt.savefig("../Outputs/0_0_5/Results_SEED_0_0_5_8020.png")
plt.show()