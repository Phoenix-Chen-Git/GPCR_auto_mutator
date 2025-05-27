import pandas as pd

# 1. Read scores from CSV
score_df = pd.read_csv("query_scores.csv", sep=None, engine='python')
score_dict = dict(zip(score_df['Position'], score_df['Score']))

# 2. Read & modify PDB
with open("MC4R.pdb") as fin, open("MC4R_conserv.pdb", "w") as fout:
    for line in fin:
        if line.startswith(("ATOM", "HETATM")):
            res_num = int(line[22:26])
            score = score_dict.get(res_num)
            if score is not None:
                # Insert score into B-factor (columns 61-66)
                new_line = line[:60] + f"{score:6.2f}" + line[66:]
                fout.write(new_line)
            else:
                # No score for this residue, keep original B-factor
                fout.write(line)
        else:
            fout.write(line)

print("Done! New PDB saved as MC4R_conserv.pdb")
