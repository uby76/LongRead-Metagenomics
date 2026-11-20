"""Package raw bins and account for each input contig exactly once."""
import csv
import gzip
import sys
import tarfile
from pathlib import Path


def contig_ids(path):
    with path.open() as stream:
        for line in stream:
            if line.startswith(">"):
                yield line[1:].split()[0]


def export(sample, strategy, method, fasta, directory):
    if not directory.is_dir():
        raise ValueError(f"Missing final bin directory {directory}")
    ids = {}
    for name in contig_ids(fasta):
        if name in ids:
            raise ValueError(f"Duplicate input contig {name}")
        ids[name] = ""
    files = sorted(directory.iterdir())
    with tarfile.open("raw_MAGs.tar.gz", "w:gz", compresslevel=3) as archive:
        for path in files:
            if path.suffix not in {".fa", ".fna", ".fasta"} or not path.is_file():
                raise ValueError(f"Unexpected bin output {path}")
            bin_id = "__".join([sample, strategy, method, path.stem])
            members = list(contig_ids(path))
            if not members:
                raise ValueError(f"Empty bin {path}")
            for name in members:
                if name not in ids or ids[name]:
                    raise ValueError(f"Unknown or multiply assigned contig {name}")
                ids[name] = bin_id
            archive.add(path, arcname=bin_id + ".fna")
    with gzip.open("membership.tsv.gz", "wt") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("sample", "strategy", "method", "contig", "bin"))
        writer.writerows((sample, strategy, method, name, group) for name, group in ids.items())
    with open("summary.tsv", "w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("sample", "strategy", "method", "input_contigs", "binned_contigs", "raw_bins"))
        writer.writerow((sample, strategy, method, len(ids), sum(map(bool, ids.values())), len(files)))


if __name__ == "__main__":
    sample, strategy, method, fasta, bins = sys.argv[1:]
    export(sample, strategy, method, Path(fasta), Path(bins))
