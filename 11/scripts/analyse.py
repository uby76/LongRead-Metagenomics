"""Species-level quality-ranked greedy dereplication and GTDB result validation."""
import csv
import gzip
import shutil
import sys
import tarfile
from pathlib import Path


def table(path, delimiter):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write(path, fields, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def inputs(metrics):
    rows = [r for r in table(metrics, ",") if r["passes_medium_threshold"].lower() == "true"]
    genomes = {r["MAG_ID"]: r for r in rows}
    if len(genomes) != 5390 or len(genomes) != len(rows):
        raise ValueError("Unexpected input count or duplicate MAG IDs")
    return genomes


def unpack(metrics):
    expected = set(inputs(metrics))
    Path("genomes").mkdir()
    seen = set()
    for path in sorted(Path("inputs").rglob("*.gz")):
        with tarfile.open(path, "r|gz") as archive:
            for entry in archive:
                if not entry.isfile():
                    continue
                name = Path(entry.name)
                if name.suffix != ".fna" or (name.stem not in expected or name.stem in seen):
                    raise ValueError(f"Unexpected archive member {entry.name}")
                with (Path("genomes") / name.name).open("xb") as output:
                    shutil.copyfileobj(archive.extractfile(entry), output)
                seen.add(name.stem)
    if seen != expected:
        raise ValueError(f"Missing {len(expected - seen)} MAGs")


def cluster(metrics):
    genomes = inputs(metrics)
    neighbors = {name: {} for name in genomes}
    for edge in table("pairs.tsv", "\t"):
        a, b = (Path(edge[k]).stem for k in ["Ref_file", "Query_file"])
        ani, af_a, af_b = (float(edge[k]) for k in ["ANI", "Align_fraction_ref", "Align_fraction_query"])
        if a != b and ani >= 95 and min(af_a, af_b) >= 50:
            neighbors[a][b] = (ani, af_b, af_a)
            neighbors[b][a] = (ani, af_a, af_b)

    def rank(name):
        r = genomes[name]
        comp, contam = float(r["completeness_percent"]), float(r["contamination_percent"])
        return (-(comp - 5 * contam), -comp, contam, -int(r["N50_bp"]), name)

    membership = {}
    for name in sorted(genomes, key=rank):
        if name in membership:
            continue
        membership[name] = (name, 100, 100, 100)
        for other, match in neighbors[name].items():
            membership.setdefault(other, (name, *match))
    fields = ["sample_ID", "sequencing_strategy", "MAG_ID", "domain", "completeness_percent", "contamination_percent"]
    write("membership.tsv", fields + ["representative_MAG", "ANI", "AF_member", "AF_representative"],
          [[genomes[n][k] for k in fields] + list(membership[n]) for n in sorted(genomes)])
    with tarfile.open("representatives.tar.gz", "w:gz", compresslevel=2) as archive:
        for rep in sorted({v[0] for v in membership.values()}):
            archive.add(Path("genomes") / (rep + ".fna"), arcname=rep + ".fna")


def finish():
    from Bio import Phylo, SeqIO
    expected = {p.stem for p in Path("reps").glob("*.fna")}
    taxonomy, tips = {}, {}
    for domain in ("bac120", "ar53"):
        summary = Path("classify/classify") / f"gtdbtk.{domain}.summary.tsv"
        if summary.exists():
            for row in table(summary, "\t"):
                name = row["user_genome"]
                if name in taxonomy or name not in expected:
                    raise ValueError(f"Unexpected GTDB row {name}")
                taxonomy[name] = row
        msa = Path("classify/align") / f"gtdbtk.{domain}.user_msa.fasta.gz"
        if not msa.exists():
            continue
        import subprocess
        from Bio import AlignIO
        with gzip.open(msa, "rt") as handle:
            alignment = AlignIO.read(handle, "fasta")
        ids = {r.id for r in alignment}
        if len(ids) != len(alignment) or not ids <= expected or ids.intersection(tips):
            raise ValueError("Alignment IDs disagree with representatives")
        if len(ids) < 2:
            continue
        tree = Path("results") / (domain + ".unrooted.nwk")
        SeqIO.write(alignment, "domain.faa", "fasta")
        with tree.open("w") as output:
            subprocess.run(["FastTree", "-lg", "-gamma", "domain.faa"],
                           stdout=output, check=True)
        if {t.name for t in Phylo.read(tree, "newick").get_terminals()} != ids:
            raise ValueError("Tree tips disagree with marker alignment")
        tips.update(dict.fromkeys(ids, domain))
    failed = Path("classify/identify/gtdbtk.failed_genomes.tsv")
    excluded = {line.split("\t")[0] for line in failed.read_text().splitlines()} if failed.exists() else set()
    filtered = set()
    for path in Path("classify/align").glob("*.filtered.tsv"):
        filtered.update(line.split("\t")[0] for line in path.read_text().splitlines())
    missing = expected - taxonomy.keys() - excluded - filtered
    if missing:
        raise ValueError(f"Missing GTDB result for {len(missing)} representatives")
    def record(n):
        r = taxonomy.get(n, {})
        lineage = r.get("classification", "")
        species = next((k[3:] for k in lineage.split(";") if k.startswith("s__")), "")
        return n, lineage, species, r.get("classification_method",""), r.get("warnings",""), tips.get(n,""), ("in_tree" if n in tips else "failed_identify" if n in excluded else "filtered_alignment" if n in filtered else "no_inferable_alignment")
    write("results/representative_taxonomy.tsv", ["representative_MAG","GTDB_taxonomy","species","classification_method","warnings","tree_marker_set","tree_status"], map(record, sorted(expected)))


if __name__ == "__main__":
    if sys.argv[1] == "finish":
        finish()
    else:
        {"unpack": unpack, "cluster": cluster}[sys.argv[1]](sys.argv[2])
