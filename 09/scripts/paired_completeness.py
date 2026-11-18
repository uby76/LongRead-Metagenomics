"""Equal-weight descriptive comparison across samples with all six strategies."""
import csv
from collections import defaultdict
from statistics import median
import sys

with open(sys.argv[1]) as handle:
    rows=list(csv.DictReader(handle,delimiter="\t"))
groups=defaultdict(list)
for row in rows:
    groups[row["strategy"]].append(row)
writer=csv.writer(sys.stdout,delimiter="\t",lineterminator="\n")
writer.writerow(["strategy","matched_samples","samples_with_completeness","median_of_sample_medians_pct"])
for strategy, entries in sorted(groups.items()):
    values=[float(r["median_completeness"]) for r in entries if r["median_completeness"]!="NA"]
    writer.writerow([strategy,len(entries),len(values),round(median(values),3) if values else "NA"])
