import sys
from collections import Counter

def calc_conservation(alignment_file, query_file, output_scores_file, output_query_scores_file):
    # Read query sequence from file
    with open(query_file) as f:
        query_seq = f.read().strip().replace("\n", "").replace(" ", "")

    # Read alignment
    with open(alignment_file) as f:
        lines = f.read().splitlines()

    # Parse FASTA
    seqs = {}
    name = None
    for line in lines:
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = ''
        else:
            seqs[name] += line

    # Transpose to columns
    aln_length = len(next(iter(seqs.values())))
    columns = [ [seqs[name][i] for name in seqs] for i in range(aln_length) ]

    # Compute conservation scores
    scores = []
    for col in columns:
        col = [x for x in col if x != '-']
        if len(col) == 0:
            scores.append(0)
        else:
            freqs = Counter(col)
            scores.append( max(freqs.values()) / len(col) )

    # Write all scores
    with open(output_scores_file, 'w') as f:
        f.write("Position,Score\n")
        for i, score in enumerate(scores, 1):
            f.write(f"{i},{score:.4f}\n")

    # Find query sequence in the alignment
    query_idx = None
    for name, seq in seqs.items():
        if query_seq in seq.replace("-", ""):
            query_idx = name
            break

    if query_idx is None:
        print("Query sequence not found in alignment.")
        return

    # Extract conservation scores for non-gap positions in the query
    query_aln = seqs[query_idx]
    query_scores = []
    for aln_char, score in zip(query_aln, scores):
        if aln_char != "-":
            query_scores.append(score)

    # Write query scores
    with open(output_query_scores_file, 'w') as f:
        f.write("Position,Score\n")
        for i, score in enumerate(query_scores, 1):
            f.write(f"{i},{score:.4f}\n")

    print(f"Query scores written to {output_query_scores_file}")

# CLI Usage
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 conservation_score_with_query_txt.py <alignment.fasta> <query.txt> <all_scores.csv> <query_scores.csv>")
        sys.exit(1)

    aln_file = sys.argv[1]
    query_file = sys.argv[2]
    all_scores_file = sys.argv[3]
    query_scores_file = sys.argv[4]

    calc_conservation(aln_file, query_file, all_scores_file, query_scores_file)
