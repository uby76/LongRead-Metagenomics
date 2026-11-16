"""Validate MetaBinner inputs and extract MAG FASTA files from final assignments."""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


def fasta(path):
    records = {}
    name = None
    parts = []
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts)
                fields = line[1:].split()
                if not fields or fields[0] in records:
                    raise ValueError(f"Missing or duplicate FASTA ID in {path}: {line}")
                name, parts = fields[0], []
            elif name is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            else:
                parts.append(line)
    if name is not None:
        records[name] = "".join(parts)
    if not records or not all(records.values()):
        raise ValueError(f"Empty FASTA or empty contig in {path}")
    return records


def align_profile(ids, source, destination, delimiter, has_header=True):
    with Path(source).open(newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader) if has_header else None
        rows = {}
        for row in reader:
            if not row or row[0] in rows:
                raise ValueError(f"Empty or duplicate profile row in {source}")
            rows[row[0]] = row
    if rows.keys() != set(ids):
        raise ValueError(f"FASTA/profile ID mismatch in {source}")
    with Path(destination).open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter=delimiter)
        if header is not None:
            writer.writerow(header)
        writer.writerows(rows[contig] for contig in ids)


def export_bins(sequences, assignments, output, group):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", group):
        raise ValueError(f"Invalid group ID: {group}")
    bins = defaultdict(list)
    assigned = set()
    with Path(assignments).open() as handle:
        for line in handle:
            if not line.strip():
                continue
            contig, bin_id = line.rstrip("\r\n").split("\t")
            if contig not in sequences or contig in assigned:
                raise ValueError(f"Unknown or duplicate assignment: {contig}")
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", bin_id):
                raise ValueError(f"Invalid bin ID: {bin_id}")
            assigned.add(contig)
            bins[bin_id].append(contig)
    if not bins:
        raise ValueError("MetaBinner reported no final MAGs.")
    mag_dir = output / "raw_MAGs"
    mag_dir.mkdir()
    with (output / "final_contig_to_MAG.csv").open("x", newline="") as handle:
        membership = csv.writer(handle)
        membership.writerow(["MAG_ID", "contig_ID"])
        for bin_id, contigs in sorted(bins.items()):
            mag = f"{group}__bin_{bin_id}"
            with (mag_dir / f"{mag}.fna").open("x") as out:
                for contig in contigs:
                    out.write(f">{contig}\n{sequences[contig]}\n")
                    membership.writerow([mag, contig])
    with (output / "unbinned_contigs.csv").open("x", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["contig_ID", "length_bp"])
        writer.writerows((c, len(s)) for c, s in sequences.items() if c not in assigned)


if __name__ == "__main__":
    command, source, *args = sys.argv[1:]
    sequences = fasta(source)
    if command == "check-filter":
        expected = [c for c, seq in sequences.items() if len(seq) > 1000]
        if list(fasta(args[0])) != expected:
            raise ValueError("Filtered FASTA differs from the expected >1,000 bp contigs.")
    elif command == "align":
        coverage, kmer, output = args
        output = Path(output)
        align_profile(sequences, coverage, output / "coverage_profile.tsv", "\t")
        align_profile(sequences, kmer, output / "kmer_profile.csv", ",")
    elif command == "export":
        assignments, output, group = args
        export_bins(sequences, assignments, Path(output), group)
    elif command != "check-fasta":
        raise ValueError(f"Unknown command: {command}")
