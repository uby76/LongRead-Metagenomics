"""Audit geNomad virus summary tables. Input: one or more TSV paths."""
import csv
import sys
from pathlib import Path

writer = csv.writer(sys.stdout, delimiter="\t")
writer.writerow(["table", "predictions", "above_10000", "monodnaviria_1001_to_10000", "eligible", "Monodnaviria_label", "monodnaviridae_label"])
for filename in sys.argv[1:]:
    total = long = short = realm = literal = 0
    with Path(filename).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            labels = {s.strip() for s in row["taxonomy"].split(";")}
            bp = int(row["length"])
            total += 1
            realm += "Monodnaviria" in labels
            literal += any(s.lower() == "monodnaviridae" for s in labels)
            long += bp > 10000
            short += 1000 < bp <= 10000 and "Monodnaviria" in labels
    writer.writerow([Path(filename).name, total, long, short, long + short, realm, literal])
