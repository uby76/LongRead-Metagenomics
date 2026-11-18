"""Descriptive CheckV, geNomad and candidate ARG counts in published groups."""
import argparse
import csv
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from tempfile import TemporaryDirectory
from urllib.parse import quote

import requests

parser = argparse.ArgumentParser()
parser.add_argument("snapshot", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--rgi-manifest", required=True, type=Path)
parser.add_argument("--genomad-manifest", required=True, type=Path)
args = parser.parse_args()
if args.output.exists():
    parser.error("Use a new result directory.")


def inputs(path):
    return {(r["sample"], r["strategy"]): (r["uri"],r["bytes"])
            for r in json.loads(path.read_text())}


if inputs(args.rgi_manifest) != inputs(args.genomad_manifest):
    raise ValueError("RGI and geNomad input manifests differ")
session = requests.Session()
session.headers["Authorization"] = "Bearer " + subprocess.check_output(
    ["gcloud","auth","application-default","print-access-token"],text=True).strip()
api = "https://storage.googleapis.com/storage/v1/b/longshort-rawdata-wb-luminous-greens-2127/o/"
groups = defaultdict(list)
with (args.snapshot/"sequence_metadata.tsv").open() as src:
    for row in csv.DictReader(src,delimiter="\t"):
        groups[row["sample"],row["strategy"]].append(row)
quality, taxonomy, arg = [], [], []
for (sample,strategy), viruses in sorted(groups.items()):
    intervals = defaultdict(list)
    for v in viruses:
        length = int(v["length"])
        start, end = (1,length) if v["coordinates"]=="NA" else map(int,v["coordinates"].split("-"))
        if end-start+1 != length:
            raise ValueError("Unexpected viral interval length")
        intervals[v["source_contig"]].append((start,end,v["virus_id"]))
    key = f"jianwei_result/CARD_RGI_contigs/result/groups/{sample}/{strategy}/CARD_RGI.txt"
    carriers, loci = set(), set()
    with session.get(api+quote(key,safe=""),params={"alt":"media"},stream=True,timeout=120) as response:
        response.raise_for_status()
        hits = csv.DictReader((line.decode() for line in response.iter_lines()),delimiter="\t")
        if not {"Contig","Start","Stop","Cut_Off","ORF_ID"} <= set(hits.fieldnames or []):
            raise ValueError("Unexpected RGI columns")
        for hit in hits:
            if hit["Cut_Off"] not in {"Perfect","Strict"}:
                raise ValueError("Unexpected RGI category")
            low, high = sorted(map(int,(hit["Start"],hit["Stop"])))
            if low < 1:
                raise ValueError("Unexpected RGI coordinates")
            for start,end,virus in intervals.get(hit["Contig"],()):
                if start<=low<=high<=end:
                    carriers.add(virus)
                    loci.add((hit["Contig"],low,high,hit["ORF_ID"]))
    arg.append([sample,strategy,len(viruses),len(carriers),len(loci)])
    scores = [float(v["completeness"]) for v in viruses if v["completeness"] not in ("NA","","None")]
    levels = Counter(v["checkv_quality"] for v in viruses)
    quality.append([sample,strategy,len(viruses),len(scores),len(viruses)-len(scores),
                    median(scores) if scores else "NA",levels["Complete"],levels["High-quality"]])
    for lineage,n in sorted(Counter(v["genomad_taxonomy"] for v in viruses).items()):
        taxonomy.append([sample,strategy,lineage,n])
samples = sorted({s for s,_ in groups if sum(k[0]==s for k in groups)==6})
paired = [row for row in quality if row[0] in samples]
args.output.parent.mkdir(parents=True,exist_ok=True)
with TemporaryDirectory(dir=args.output.parent) as staging:
    stage=Path(staging)/"report"
    for name,header,rows in [
        ("quality/by_group","sample strategy n completeness_available completeness_missing median_completeness complete high_quality",quality),
        ("quality/six_strategy_matched","sample strategy n completeness_available completeness_missing median_completeness complete high_quality",paired),
        ("taxonomy/genomad","sample strategy lineage n",taxonomy),
        ("arg/candidate_counts","sample strategy viral_records records_with_candidate_ARG distinct_candidate_ORFs",arg),
    ]:
        file=stage/(name+".tsv");file.parent.mkdir(exist_ok=True,parents=True)
        with file.open("x",newline="") as out:
            writer=csv.writer(out,delimiter="\t",lineterminator="\n")
            writer.writerow(header.split());writer.writerows(rows)
    stage.rename(args.output)
print(len(groups),"published groups;",len(samples),"six-strategy samples;",sum(r[3] for r in arg),"viral records with candidate ARG;",sum(r[4] for r in arg),"candidate ORFs")
