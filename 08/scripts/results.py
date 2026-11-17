import csv
import json
import re
import sys
from pathlib import Path

columns = ["sample_ID", "viral_predictions", "provirus_predictions", "plasmid_predictions"]
if sys.argv[1] == "prepare":
    source, local = map(Path, sys.argv[2:])
    local.mkdir()
    files = sorted(source.iterdir())
    names = {path.name for path in files}
    for item in files:
        target = item.resolve()
        # Cloud Storage holds these archived symlinks as relative path strings.
        if item.name.startswith("genomad_mini_db") and item.stat().st_size <= 64:
            try:
                link = item.read_text().strip()
            except UnicodeDecodeError:
                link = ""
            if re.fullmatch(r"genomad_db[\w.]*", link) and link in names:
                target = (source / link).resolve()
        (local / item.name).symlink_to(target)
    assert (local / "version.txt").read_text().strip() == "1.7"
elif sys.argv[1] == 'verify-input':
    import base64
    import hashlib
    digest = hashlib.md5()
    with open(sys.argv[2], 'rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    if base64.b64encode(digest.digest()).decode() != sys.argv[3]:
        raise ValueError('Input checksum differs from the audited manifest')
elif sys.argv[1] == 'merge':
    fields = ['sample_ID', 'sequencing_strategy', 'viral_predictions', 'provirus_predictions', 'plasmid_predictions']
    expected = {(r['sample'], r['strategy']) for r in json.loads(Path(sys.argv[2]).read_text())}
    rows = []
    for name in sys.argv[3:]:
        with open(name, newline='') as handle:
            entries = list(csv.DictReader(handle, delimiter='\t'))
        if len(entries) != 1:
            raise ValueError('Expected one summary row per group')
        rows.extend(entries)
    if len(rows) != len(expected) or {(r['sample_ID'], r['sequencing_strategy']) for r in rows} != expected:
        raise ValueError('Incomplete or duplicate sample/strategy groups')
    writer = csv.DictWriter(sys.stdout, fields, delimiter='\t')
    writer.writeheader()
    writer.writerows(sorted(rows, key=lambda r: (r['sample_ID'], r['sequencing_strategy'])))
else:
    sample, strategy, summary = sys.argv[1:]
    summary = Path(summary)
    def read(kind):
        with (summary / f'assembly_{kind}_summary.tsv').open(newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not reader.fieldnames or not {'seq_name', 'length'} <= set(reader.fieldnames):
                raise ValueError('Invalid geNomad summary header')
            rows = list(reader)
        expected = {row['seq_name']: int(row['length']) for row in rows}
        if len(expected) != len(rows):
            raise ValueError('Duplicate geNomad sequence IDs')
        observed = {}
        name = None
        with (summary / f'assembly_{kind}.fna').open() as handle:
            for line in handle:
                if line.startswith('>'):
                    name = line[1:].split()[0]
                    if name in observed:
                        raise ValueError('Duplicate FASTA ID')
                    observed[name] = 0
                elif line.strip():
                    if name is None:
                        raise ValueError('FASTA sequence before header')
                    observed[name] += len(line.strip())
        if observed != expected:
            raise ValueError('FASTA and summary sequence IDs/lengths differ')
        return rows
    virus, plasmid = read('virus'), read('plasmid')
    writer = csv.writer(sys.stdout, delimiter='\t')
    writer.writerow(['sample_ID', 'sequencing_strategy', 'viral_predictions', 'provirus_predictions', 'plasmid_predictions'])
    writer.writerow([sample, strategy, len(virus), sum('|provirus_' in r['seq_name'] for r in virus), len(plasmid)])
