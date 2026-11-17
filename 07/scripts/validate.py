import csv
import json
import sys
from collections import Counter
from pathlib import Path

FIELDS = ["sample_ID", "sequencing_strategy", "perfect_hits", "strict_hits", "total_hits"]

if sys.argv[1] == "merge":
    expected = {(r["sample"], r["strategy"]) for r in json.loads(Path(sys.argv[2]).read_text())}
    rows = []
    for filename in sys.argv[3:]:
        with open(filename, newline="") as handle:
            row, = csv.DictReader(handle, delimiter="\t")
        rows.append(row)
    assert len(rows) == len(expected)
    assert {(r["sample_ID"], r["sequencing_strategy"]) for r in rows} == expected
    writer = csv.DictWriter(sys.stdout, FIELDS, delimiter="\t")
    writer.writeheader()
    writer.writerows(sorted(rows, key=lambda row: (row["sample_ID"], row["sequencing_strategy"])))
else:
    sample, strategy, prefix = sys.argv[1:]
    with open(prefix + ".json") as handle:
        assert isinstance(json.load(handle), dict)
    with open(prefix + ".txt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames and "Cut_Off" in reader.fieldnames
        counts = Counter(row["Cut_Off"] for row in reader)
    assert set(counts) <= {"Strict", "Perfect"}, counts
    writer = csv.writer(sys.stdout, delimiter="\t")
    writer.writerow(FIELDS)
    writer.writerow([sample, strategy, counts["Perfect"], counts["Strict"], sum(counts.values())])
