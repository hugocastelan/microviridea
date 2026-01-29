#!/usr/bin/env python3
# Capsid Detection using ESM 

import argparse
import torch
import pandas as pd
from Bio import SeqIO
from transformers import EsmTokenizer, EsmModel
from sklearn.ensemble import IsolationForest
import numpy as np

# --------------------------------------------------
# ESM embedding
# --------------------------------------------------
def embed(seqs, tokenizer, model, device):
    X = []
    for s in seqs:
        toks = tokenizer(s, return_tensors="pt")
        toks = {k: v.to(device) for k, v in toks.items()}
        with torch.no_grad():
            out = model(**toks)
        emb = out.last_hidden_state.mean(dim=1).squeeze().cpu().numpy()
        X.append(emb)
    return np.vstack(X)

# --------------------------------------------------
# Main
# --------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--positives", required=True, help="Known capsid FASTA") #We download sequences of capside of Microviridae from NCBI
    ap.add_argument("-i", "--input", required=True, help="ORFs FASTA")
    ap.add_argument("-o", "--output", default="capsid_predictions.tsv")
    ap.add_argument("--faa", default="capsid_only.faa")
    args = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"

    print("Loading ESM-2...")
    tokenizer = EsmTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
    model = EsmModel.from_pretrained("facebook/esm2_t6_8M_UR50D").to(device)
    model.eval()

    # -------------------------------
    # Load capsides (training)
    # -------------------------------
    caps = [(r.id, str(r.seq)) for r in SeqIO.parse(args.positives, "fasta")]
    df_pos = pd.DataFrame(caps, columns=["id", "seq"])

    print("🔬 Embedding known capsids...")
    X_pos = embed(df_pos["seq"], tokenizer, model, device)

    # -------------------------------
    # Train one-class model
    # -------------------------------
    clf = IsolationForest(
        n_estimators=300,
        contamination=0.05,
        random_state=42
    )
    clf.fit(X_pos)

    # -------------------------------
    # Load ORFs
    # -------------------------------
    orfs = [(r.id, str(r.seq)) for r in SeqIO.parse(args.input, "fasta")]
    df = pd.DataFrame(orfs, columns=["id", "seq"])
    df["length"] = df["seq"].str.len()

    print(" Embedding ORFs...")
    X = embed(df["seq"], tokenizer, model, device)

    # -------------------------------
    # Predict
    # -------------------------------
    scores = clf.decision_function(X)
    preds = clf.predict(X)  # 1 = inlier, -1 = outlier

    df["capsid_score"] = scores
    df["pred_capsid"] = (preds == 1).astype(int)

    df.to_csv(args.output, sep="\t", index=False)

    # -------------------------------
    # Export FASTA
    # -------------------------------
    with open(args.faa, "w") as f:
        for _, r in df[df["pred_capsid"] == 1].iterrows():
            f.write(f">{r.id}\n{r.seq}\n")

    print(f" Capsides predicted: {df.pred_capsid.sum()}")
    print(f" {args.output}")
    print(f" {args.faa}")

if __name__ == "__main__":
    main()
