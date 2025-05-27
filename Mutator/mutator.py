import pandas as pd
from Bio import SeqIO

fasta_path = "P32245.fasta"
score_path = "query_scores.csv"

# Standard amino acids
aas = list("ACDEFGHIKLMNPQRSTVWY")

# Read sequence
seq_record = next(SeqIO.parse(fasta_path, "fasta"))
seq = str(seq_record.seq)
length = len(seq)

# Read conservation scores
scores = pd.read_csv(score_path)
score_col = [col for col in scores.columns if "score" in col.lower()][0]
scores = scores[score_col].tolist()
if len(scores) != length:
    raise ValueError("Sequence length and scores do not match!")

# Build original DataFrame
origin_df = pd.DataFrame({
    "Pos": range(1, length + 1),
    "AA_origin": list(seq),
    "con": scores
})

# Generate all single mutations (19 per site, exclude self)
mutants = []
for idx, row in origin_df.iterrows():
    pos = row["Pos"]
    orig_aa = row["AA_origin"]
    con = row["con"]
    for aa in aas:
        if aa != orig_aa:
            mutants.append({
                "Pos": pos,
                "AA_origin": orig_aa,
                "con": con,
                "Mutant_AA": aa,
                "Mutation_Label": f"{orig_aa}{pos}{aa}"
            })

mutant_df = pd.DataFrame(mutants)

# Save result
mutant_df.to_csv("single_site_mutants.csv", index=False)
print(mutant_df.head(30))
