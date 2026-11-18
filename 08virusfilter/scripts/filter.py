"""Select records from the geNomad virus FASTA; validate IDs and lengths."""
import csv
import sys
from pathlib import Path

summary, fasta, expected = sys.argv[1:]
with open(summary, newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
lengths = {r["seq_name"]: int(r["length"]) for r in rows}
keep = {r["seq_name"] for r in rows if int(r["length"]) > 10000 or
        (int(r["length"]) > 1000 and "Monodnaviria" in r["taxonomy"].split(";"))}
if len(lengths) != len(rows) or len(keep) != int(expected):
    raise ValueError("Summary disagrees with label audit or contains duplicate IDs")
seen = {}
name = None
with open(fasta) as src, open("selected.fna", "w") as out:
    for line in src:
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in seen:
                raise ValueError("Duplicate FASTA ID")
            seen[name] = 0
        else:
            seen[name] += len(line.strip())
        if name in keep:
            out.write(line)
if seen != lengths:
    raise ValueError("Virus FASTA and summary differ")
