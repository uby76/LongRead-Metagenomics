"""Cluster a complete 42-group snapshot. Run on a Workbench VM with BLAST+ and gcloud."""
import argparse
import csv
import itertools
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("snapshot", type=Path)
p.add_argument("output", type=Path)
p.add_argument("--threads", type=int, default=8)
args = p.parse_args()


def read(name):
    with (args.snapshot / name).open(newline="") as src:
        return list(csv.DictReader(src, delimiter="\t"))


groups = read("input_status.tsv")
if len(groups) != 42 or len({(g["sample"],g["strategy"]) for g in groups}) != 42 or any(g["status"] != "published" for g in groups):
    raise SystemExit("All 42 distinct published groups are required.")
if args.output.exists():
    raise SystemExit("Choose a new output directory.")
for tool in ("blastn","makeblastdb","gcloud","sort"):
    if not shutil.which(tool):
        raise SystemExit("Missing executable: "+tool)
metadata = read("sequence_metadata.tsv")
records = {}
for i,row in enumerate(metadata,1):
    key = row["sample"],row["strategy"],row["virus_id"]
    if key in records:
        raise ValueError(key)
    records[key] = (f"v{i:07}",int(row["length"]))
lengths = dict(records.values())
edges = defaultdict(set)
with tempfile.TemporaryDirectory() as tmp:
    work = Path(tmp)
    pooled = work / "pooled.fna"
    seen = {}
    with pooled.open("w") as dst:
        for g in groups:
            source = work / "group.fna"
            subprocess.run(["gcloud","storage","cp",g["selected_fasta"],str(source)],check=True)
            with source.open() as src:
                for line in src:
                    if line.startswith(">"):
                        tag,_ = records[g["sample"],g["strategy"],line[1:].split()[0]]
                        if tag in seen:
                            raise ValueError(tag)
                        seen[tag] = 0
                        dst.write(">"+tag+"\n")
                    else:
                        seen[tag] += len(line.strip())
                        dst.write(line.strip()+"\n")
            source.unlink()
    if seen != lengths:
        raise ValueError("Pooled FASTA differs from the audited IDs or lengths.")
    db = str(work/"db")
    subprocess.run(["makeblastdb","-dbtype","nucl","-in",str(pooled),"-out",db],check=True)
    hits = work/"hits"
    subprocess.run(["blastn","-task","megablast","-dust","no","-query",str(pooled),"-db",db,
                    "-num_threads",str(args.threads),"-max_target_seqs",str(len(lengths)),
                    "-outfmt","6 qseqid sseqid nident length qstart qend sstart send bitscore","-out",str(hits)],check=True)
    with (work/"sorted").open("w") as out:
        subprocess.run(["sort","-k1,1","-k2,2","-k9,9nr",str(hits)],stdout=out,check=True)
    with (work/"sorted").open() as src:
        for (q,s),hsps in itertools.groupby((l.split() for l in src),lambda h:(h[0],h[1])):
            if q==s or lengths[q]>lengths[s]:
                continue
            intervals=[]
            matches=aligned=covered=0
            for h in hsps:
                a,b=sorted(map(int,h[4:6]))
                c,d=sorted(map(int,h[6:8]))
                if any((a<=y and x<=b) or (c<=v and u<=d) for x,y,u,v in intervals):
                    continue
                intervals.append((a,b,c,d))
                covered+=b-a+1
                matches+=int(h[2])
                aligned+=int(h[3])
            if aligned and 100*matches>=95*aligned and 100*covered>=85*lengths[q]:
                edges[q].add(s)
                edges[s].add(q)
    reps={}
    for tag in sorted(lengths,key=lambda t:(-lengths[t],t)):
        if tag in reps:
            continue
        reps[tag]=tag
        for member in edges[tag]:
            reps.setdefault(member,tag)
    args.output.mkdir(parents=True,exist_ok=False)
    with (args.output/"membership.tsv").open("w",newline="") as out:
        w=csv.writer(out,delimiter="\t",lineterminator="\n")
        w.writerow(["vOTU","member",*metadata[0]])
        for row in metadata:
            tag,_=records[row["sample"],row["strategy"],row["virus_id"]]
            w.writerow([reps[tag],tag,*row.values()])
    with pooled.open() as src,(args.output/"representatives.fna").open("w") as out:
        for line in src:
            if line.startswith(">"):
                tag=line[1:].strip()
                retain=reps[tag]==tag
            if retain:
                out.write(line)
