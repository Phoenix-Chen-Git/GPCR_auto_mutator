from Bio import SeqIO, SeqRecord, Seq

def full_mutator(position, fasta_path, out_fasta=None):
    """
    Generate all 20 AA mutants at given 1-based position.
    position: int (1-based)
    fasta_path: fasta file with a single sequence
    out_fasta: optional output file name
    returns: list of SeqRecord (all mutants, including WT)
    """
    aas = list("ACDEFGHIKLMNPQRSTVWY")
    # Read sequence
    seq_record = next(SeqIO.parse(fasta_path, "fasta"))
    seq = str(seq_record.seq)
    orig_aa = seq[position - 1]
    records = []
    for aa in aas:
        mut_seq = list(seq)
        mut_seq[position - 1] = aa
        label = f"{orig_aa}{position}{aa}" if aa != orig_aa else f"WT_{orig_aa}{position}{aa}"
        records.append(SeqRecord.SeqRecord(Seq.Seq(''.join(mut_seq)), id=label, description=""))
    if out_fasta:
        SeqIO.write(records, out_fasta, "fasta")
    return records

# Example usage:
full_mutator(5, "P32245.fasta", "pos5_mutants.fasta")
