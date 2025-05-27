from Bio import SeqIO, SeqRecord, Seq

def full_mutator(position, fasta_path, out_fasta=None):
    """
    For each sequence in a multi-FASTA, generate all 20 AA mutants at the given 1-based position.
    Each mutant's ID: {original_seq_id}|{originalAA}{pos}{mutAA}
    WT sequence included as {original_seq_id}|WT_{originalAA}{pos}{originalAA}
    """
    aas = list("ACDEFGHIKLMNPQRSTVWY")
    records = []
    for seq_record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(seq_record.seq)
        if position > len(seq):
            raise ValueError(f"Position {position} exceeds sequence length for {seq_record.id}")
        orig_aa = seq[position - 1]
        for aa in aas:
            mut_seq = list(seq)
            mut_seq[position - 1] = aa
            label = f"{seq_record.id}|{orig_aa}{position}{aa}" if aa != orig_aa else f"{seq_record.id}|WT_{orig_aa}{position}{aa}"
            records.append(SeqRecord.SeqRecord(Seq.Seq(''.join(mut_seq)), id=label, description=""))
    if out_fasta:
        SeqIO.write(records, out_fasta, "fasta")
    return records

# Example usage:
full_mutator(5, "multi.fasta", "pos5_mutants.fasta")
