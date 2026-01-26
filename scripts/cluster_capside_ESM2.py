#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PCA of capsid protein clusters using ESM-2
Vents and Fresh analyzed independently
Colored by cluster
Each point annotated with number of sequences
Exports sequence IDs per cluster
"""

import os
import glob
import torch
import esm
import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================
# CONFIG
# ==========================
ENVIRONMENTS = {
    "vents": "/home/hugocastelan/Documents/projects/microviriadae/sonia/mmseqs_vents/clusters_vents",
    "fresh": "/home/hugocastelan/Documents/projects/microviriadae/sonia/mmseqs_fresh/clusters_fresh"
}

MAX_SEQ_PER_CLUSTER = 200
MIN_LEN = 80
MODEL_LAYER = 6

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYX")

# ==========================
# LOAD MODEL
# ==========================
print("Loading ESM-2 model...")
model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
model.eval()
batch_converter = alphabet.get_batch_converter()

# ==========================
# EMBEDDING FUNCTION
# ==========================
def embed(seq):
    seq = seq.replace("*", "").upper()
    seq = "".join([aa for aa in seq if aa in VALID_AA])

    if len(seq) < MIN_LEN:
        return None

    data = [("seq", seq)]
    _, _, tokens = batch_converter(data)

    with torch.no_grad():
        out = model(tokens, repr_layers=[MODEL_LAYER])

    emb = out["representations"][MODEL_LAYER][0, 1:-1].mean(0)
    return emb.cpu().numpy()

# ==========================
# MAIN LOOP
# ==========================
for env, folder in ENVIRONMENTS.items():

    print(f"\nProcessing {env.upper()} clusters...")
    rows = []
    cluster_seq_rows = []

    for faa in sorted(glob.glob(os.path.join(folder, "*.faa"))):
        cluster = os.path.basename(faa).replace(".faa", "")
        embs = []
        seq_ids = []
        discarded = 0

        for i, rec in enumerate(SeqIO.parse(faa, "fasta")):
            if i >= MAX_SEQ_PER_CLUSTER:
                break

            vec = embed(str(rec.seq))
            if vec is None:
                discarded += 1
                continue

            embs.append(vec)
            seq_ids.append(rec.id)

        if len(embs) < 3:
            continue

        mean_emb = np.mean(embs, axis=0)

        row = {
            "cluster": cluster,
            "n_seq": len(embs)
        }

        for i, v in enumerate(mean_emb):
            row[f"f{i}"] = v

        rows.append(row)

        cluster_seq_rows.append({
            "cluster": cluster,
            "n_seq": len(seq_ids),
            "sequence_ids": ";".join(seq_ids)
        })

        if discarded > 0:
            print(f"{cluster}: {discarded} sequences discarded after cleaning")

    if len(rows) == 0:
        print(f"No valid clusters found for {env}")
        continue

    df = pd.DataFrame(rows)

    # ==========================
    # SAVE EMBEDDINGS
    # ==========================
    emb_out = f"{env}_clusters_embeddings.tsv"
    df.to_csv(emb_out, sep="\t", index=False)

    # ==========================
    # SAVE SEQUENCE IDs PER CLUSTER
    # ==========================
    ids_df = pd.DataFrame(cluster_seq_rows)
    ids_out = f"{env}_clusters_sequence_IDs.tsv"
    ids_df.to_csv(ids_out, sep="\t", index=False)

    # ==========================
    # PCA
    # ==========================
    X = df.filter(regex="^f")
    Xs = StandardScaler().fit_transform(X)

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(Xs)

    df["PC1"] = pcs[:, 0]
    df["PC2"] = pcs[:, 1]

    pca_out = f"{env}_clusters_PCA.tsv"
    df.to_csv(pca_out, sep="\t", index=False)

    # ==========================
    # PLOT
    # ==========================
    plt.figure(figsize=(9, 8))
    ax = sns.scatterplot(
        data=df,
        x="PC1",
        y="PC2",
        hue="cluster",
        size="n_seq",
        palette="tab20",
        alpha=0.85,
        legend=False
    )

    for _, r in df.iterrows():
        ax.text(
            r["PC1"],
            r["PC2"],
            str(r["n_seq"]),
            fontsize=8,
            ha="center",
            va="center"
        )

    plt.title(f"PCA of Capsid Protein Clusters ({env})")
    plt.tight_layout()
    plt.savefig(f"{env}_clusters_PCA.pdf")
    plt.close()

    print(f"{env} analysis completed")
    print(f"Outputs:")
    print(f"  {emb_out}")
    print(f"  {pca_out}")
    print(f"  {ids_out}")
    print(f"  {env}_clusters_PCA.pdf")

print("\nAll analyses completed")
