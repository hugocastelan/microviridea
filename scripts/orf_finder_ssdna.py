#!/usr/bin/env python3
"""
orf_finder_ssdna.py — ORF caller for (circular) ssDNA viruses

Outputs:
- .faa  protein sequences
- .fna  nucleotide sequences
- .tsv  ORF coordinates table
- .gff3 genome annotation

Features:
- Accepts full IUPAC nucleotide alphabet (ACGTNRYWSKMBDHV)
- Rejects protein FASTA
- Filters ORFs with internal stop codons
- Supports circular genomes
"""

import argparse
import sys
from typing import Dict, List


# Genetic code (table 11)

def _build_table():
    aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    codons = [
        "TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA","TAG","TGT","TGC","TGA","TGG",
        "CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG",
        "ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG",
        "GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG"
    ]
    return dict(zip(codons, aa))

CODON_TABLE = _build_table()
STOP_CODONS = {c for c, a in CODON_TABLE.items() if a == "*"}
START_CODONS = {"ATG", "GTG", "TTG"}  # table 11
IUPAC_NT = set("ACGTNRYWSKMBDHV")


# Utilities

def revcomp(seq: str) -> str:
    comp = str.maketrans(
        "ACGTNRYWSKMBDHVacgtnrywskmbdhv",
        "TGCANYRWSMKVHDBtgcanyrwsmkvhdb"
    )
    return seq.translate(comp)[::-1]

def translate(seq: str) -> str:
    aa = []
    for i in range(0, len(seq) - 2, 3):
        aa.append(CODON_TABLE.get(seq[i:i+3].upper(), "X"))
    return "".join(aa)

def has_internal_stop(aa: str) -> bool:
    return "*" in aa[:-1]

def is_probably_nucleotide(seq: str, threshold: float = 0.85) -> bool:
    valid = sum(c in IUPAC_NT for c in seq)
    return (valid / len(seq)) >= threshold


# ORF finding

def find_orfs(seq: str, seq_id: str, min_aa: int, circular: bool) -> List[dict]:
    L = len(seq)
    seq2 = seq + (seq if circular else "")
    orfs = []

    for strand, s in [("+", seq), ("-", revcomp(seq))]:
        for frame in range(3):
            start = None
            for i in range(frame, len(seq2) - 2, 3):
                cod = seq2[i:i+3].upper()
                if start is None and cod in START_CODONS:
                    start = i
                elif start is not None and cod in STOP_CODONS:
                    if start < L:
                        nt_len = i + 3 - start
                        if nt_len <= L and nt_len % 3 == 0:
                            s1 = start % L
                            e1 = (i + 2) % L
                            wraps = e1 < s1
                            nuc = s[s1:] + s[:e1+1] if wraps else s[s1:e1+1]
                            aa = translate(nuc)
                            if aa:
                                aa = "M" + aa[1:]
                            if len(aa) >= min_aa and not has_internal_stop(aa):
                                if aa.endswith("*"):
                                    aa = aa[:-1]

                                start_g = s1 + 1
                                end_g = e1 + 1 + (L if wraps else 0)
                                if strand == "-":
                                    start_g, end_g = L - end_g + 1, L - start_g + 1

                                orfs.append({
                                    "seq_id": seq_id,
                                    "strand": strand,
                                    "frame": frame,
                                    "start": start_g,
                                    "end": end_g,
                                    "len_nt": len(nuc),
                                    "len_aa": len(aa),
                                    "protein": aa,
                                    "nuc_seq": nuc,
                                    "wraps": wraps
                                })
                    start = None
    return orfs


# FASTA I/O

def read_fasta(path: str) -> Dict[str, str]:
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name:
            seqs[name] = "".join(buf)
    return seqs


# OUTPUTS

def write_outputs(orfs: List[dict], prefix: str):
    faa, fna = [], []
    tsv = ["orf_id\tseq_id\tstrand\tframe\tstart\tend\tlen_nt\tlen_aa"]
    gff = ["##gff-version 3"]

    for i, o in enumerate(orfs, 1):
        oid = f"ORF_{i:05d}"
        header = f">{oid}|{o['seq_id']}|{o['strand']}|{o['start']}-{o['end']}"

        faa.append(header + "\n" + o["protein"])
        fna.append(header + "\n" + o["nuc_seq"])

        tsv.append(
            f"{oid}\t{o['seq_id']}\t{o['strand']}\t{o['frame']}\t"
            f"{o['start']}\t{o['end']}\t{o['len_nt']}\t{o['len_aa']}"
        )

        attrs = (
            f"ID={oid};product=hypothetical_protein;"
            f"length_aa={o['len_aa']}"
        )
        if o["wraps"]:
            attrs += ";Note=wraps_origin"

        gff.append(
            "\t".join([
                o["seq_id"], "orf_finder_ssdna", "CDS",
                str(o["start"]), str(o["end"]), ".",
                o["strand"], "0", attrs
            ])
        )

    with open(prefix + ".faa", "w") as f:
        f.write("\n".join(faa) + "\n")
    with open(prefix + ".fna", "w") as f:
        f.write("\n".join(fna) + "\n")
    with open(prefix + ".tsv", "w") as f:
        f.write("\n".join(tsv) + "\n")
    with open(prefix + ".gff3", "w") as f:
        f.write("\n".join(gff) + "\n")


# MAIN

def main():
    ap = argparse.ArgumentParser(description="ORF caller for ssDNA viruses (IUPAC-aware)")
    ap.add_argument("fasta")
    ap.add_argument("--out-prefix", default="orfs")
    ap.add_argument("--min-aa", type=int, default=25)
    ap.add_argument("--linear", action="store_true")
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    all_orfs = []

    for sid, seq in seqs.items():
        seq = seq.upper()
        if not is_probably_nucleotide(seq):
            print(f"[skip] {sid}: not nucleotide", file=sys.stderr)
            continue
        all_orfs.extend(find_orfs(seq, sid, args.min_aa, not args.linear))

    if not all_orfs:
        print("[warn] No ORFs found", file=sys.stderr)
    else:
        write_outputs(all_orfs, args.out_prefix)
        print(
            f"[ok] Wrote {args.out_prefix}.faa .fna .tsv .gff3",
            file=sys.stderr
        )

if __name__ == "__main__":
    main()
