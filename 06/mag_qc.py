#!/usr/bin/env python3
"""Run CheckM2 and domain-specific RNA annotation for one group of raw MAGs."""

import argparse
import csv
import json
import math
import shutil
import subprocess
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def write_table(path, rows):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def assess(bins_dir, bins, database, out):
    # CheckM2 evaluates all raw bins. The quality sets are exported afterward.
    with (out / "checkm2.log").open("w") as log:
        subprocess.run(
            ["checkm2", "predict", "--threads", "8", "--input", str(bins_dir),
             "--extension", "fna", "--database_path", str(database),
             "--tmpdir", str(out / "tmp"), "--output-directory", str(out / "checkm2_output")],
            stdout=log, stderr=subprocess.STDOUT, check=True,
        )
    report = out / "checkm2_output" / "quality_report.tsv"
    with report.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != len(bins) or {row["Name"] for row in rows} != bins.keys():
        raise ValueError("CheckM2 report IDs do not match the input MAG IDs.")
    shutil.copyfile(report, out / "quality_report.tsv")

    medium = out / "MQ_or_better_MAGs"
    high = out / "HQ_MAGs"
    medium.mkdir()
    high.mkdir()
    quality = []
    for row in rows:
        mag = row["Name"]
        completeness = float(row["Completeness"])
        contamination = float(row["Contamination"])
        if not all(map(math.isfinite, (completeness, contamination))):
            raise ValueError(f"Nonfinite CheckM2 result for {mag}")
        mq = completeness >= 50 and contamination <= 10
        hq = completeness >= 90 and contamination < 5
        quality.append(dict(MAG_ID=mag, completeness_percent=completeness,
                            contamination_percent=contamination,
                            passes_medium_threshold=mq, passes_high_threshold=hq))
        if mq:
            shutil.copyfile(bins[mag], medium / bins[mag].name)
        if hq:
            shutil.copyfile(bins[mag], high / bins[mag].name)
    write_table(out / "MAG_quality.csv", quality)
    for folder in (medium, high):
        with tarfile.open(out / f"{folder.name}.tar.gz", "w:gz") as archive:
            archive.add(folder, arcname=folder.name)


def annotate(mag, fasta, domain, out):
    work = out / mag
    work.mkdir()
    bacterial = domain == "bacteria"
    # Run each tool with one thread per MAG using the assigned domain model.
    gff = work / "rrna.gff"
    with gff.open("w") as output, (work / "barrnap.log").open("w") as log:
        subprocess.run(
            ["barrnap", "-k", "bac" if bacterial else "arc", "--threads", "1",
             "--lencutoff", "0.8", "--reject", "0.25", str(fasta)],
            stdout=output, stderr=log, check=True,
        )
    calls = []
    with gff.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Invalid GFF record in {gff}")
            calls.append(fields[8])
    partial_16s = sum("16S" in attrs and
                      ("partial" in attrs.lower() or "aligned only" in attrs.lower())
                      for attrs in calls)
    count_16s = sum("16S" in attrs for attrs in calls)
    rrna = dict(MAG_ID=mag, domain=domain, rRNA_16S_complete=count_16s - partial_16s,
                rRNA_16S_partial=partial_16s, rRNA_16S_count=count_16s,
                rRNA_5S_count=sum("Name=5S_rRNA" in attrs for attrs in calls),
                rRNA_23S_count=sum("Name=23S_rRNA" in attrs for attrs in calls),
                total_rRNA_count=len(calls))

    table = work / "trnas.tsv"
    with (work / "trnascan.stdout.log").open("w") as output, \
            (work / "trnascan.stderr.log").open("w") as log:
        subprocess.run(
            ["tRNAscan-SE", "-B" if bacterial else "-A", "--thread", "1",
             "-o", str(table), "-m", str(work / "trnas_summary.txt"), str(fasta)],
            stdout=output, stderr=log, check=True,
        )
    hits = []
    with table.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 9 and fields[1].isdigit() and fields[2].lstrip("-").isdigit():
                hits.append(fields)
    pseudo = sum("pseudo" in " ".join(fields).lower() for fields in hits)
    trna = dict(MAG_ID=mag, domain=domain, tRNA_count=len(hits),
                tRNA_non_pseudogene_count=len(hits) - pseudo, tRNA_pseudogene_count=pseudo)
    return rrna, trna


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bins", type=Path)
    parser.add_argument("domains", type=Path)
    parser.add_argument("database", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    bins = {path.stem: path for path in sorted(args.bins.glob("*.fna")) if path.is_file()}
    with args.domains.open() as handle:
        rows = json.load(handle)
    domains = {row["MAG_ID"]: row["domain"] for row in rows}
    if not bins or len(rows) != len(domains) or domains.keys() != bins.keys():
        raise ValueError("Provide exactly one domain assignment per input MAG.")
    if set(domains.values()) - {"bacteria", "archaea"}:
        raise ValueError("Domain assignments must be bacteria or archaea.")
    if not args.database.is_file():
        raise ValueError(f"Missing CheckM2 database: {args.database}")

    args.output.mkdir()
    assess(args.bins, bins, args.database, args.output)
    annotations = args.output / "annotations"
    annotations.mkdir()
    with ThreadPoolExecutor(max_workers=8) as pool:
        jobs = [pool.submit(annotate, mag, fasta, domains[mag], annotations)
                for mag, fasta in bins.items()]
        rrna, trna = zip(*(job.result() for job in jobs))
    write_table(args.output / "rrna_counts.csv", rrna)
    write_table(args.output / "trna_counts.csv", trna)


if __name__ == "__main__":
    main()
