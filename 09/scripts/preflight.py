"""Snapshot published CheckV and geNomad metadata without touching running jobs."""
import argparse
import csv
import io
import json
import re
import subprocess
from collections import Counter
from pathlib import Path
from urllib.parse import quote

import requests

parser = argparse.ArgumentParser()
parser.add_argument("output", type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
groups = json.loads((root / "08virusfilter/data/manifest.json").read_text())
bucket = "longshort-rawdata-wb-luminous-greens-2127"
api = f"https://storage.googleapis.com/storage/v1/b/{bucket}/o"
session = requests.Session()
session.headers["Authorization"] = "Bearer " + subprocess.check_output(
    ["gcloud", "auth", "application-default", "print-access-token"], text=True).strip()


def get(name, media=False):
    return session.get(api + "/" + quote(name, safe=""),
                       params={"alt": "media"} if media else {"fields": "size,generation"},
                       timeout=60)


def table(response):
    response.raise_for_status()
    return list(csv.DictReader(io.StringIO(response.text), delimiter="\t"))


def write(name, fields, rows):
    with (args.output / name).open("x", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(fields.split())
        w.writerows(rows)


states, mapping, counts = [], [], []
for group in groups:
    s, t, n = (group[k] for k in ("sample", "strategy", "expected"))
    base = f"jianwei_result/08virusfilter/results/{s}/{t}"
    fq = base + "/selected.fna"
    fasta = get(fq)
    qc = get(base + "/quality/quality_summary.tsv", True)
    for response in (fasta, qc):
        if response.status_code != 404:
            response.raise_for_status()
    ready = fasta.ok and qc.ok
    states.append([s,t,n,"published" if ready else "pending", "gs://"+bucket+"/"+fq])
    if not ready:
        continue
    quality = table(qc)
    source = f"jianwei_result/geNomad_contigs_20261117/results/groups/{s}/{t}/genomad_output/assembly_summary/assembly_virus_summary.tsv"
    annotations = {r["seq_name"]:r for r in table(get(source,True))}
    if len(quality) != n or len({r["contig_id"] for r in quality}) != n:
        raise ValueError(f"Unexpected membership: {s}/{t}")
    levels = Counter()
    for r in quality:
        name = r["contig_id"]
        a = annotations[name]
        if int(r["contig_length"]) != int(a["length"]):
            raise ValueError(name)
        parent = re.sub(r"\|provirus_\d+_\d+$","",name)
        mapping.append([s,t,parent,name,a["coordinates"],a["length"],a["taxonomy"],r["checkv_quality"],r["completeness"]])
        levels[r["checkv_quality"]] += 1
    counts.extend([s,t,level,value,n] for level,value in sorted(levels.items()))
args.output.mkdir(parents=True,exist_ok=False)
write("input_status.tsv","sample strategy expected status selected_fasta",states)
write("sequence_metadata.tsv","sample strategy source_contig virus_id coordinates length genomad_taxonomy checkv_quality completeness",mapping)
write("quality_by_group.tsv","sample strategy quality count denominator",counts)
print("Published groups:",sum(r[3]=="published" for r in states),"of",len(states),"records:",len(mapping))
