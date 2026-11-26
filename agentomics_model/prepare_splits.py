"""Prepare the train and validation splits of the senBench dataset for train.py, as Agentomics-ML prepares them.

Each prepared split links the dataset's input/ folder and holds labels.csv with the columns id and numeric_label
(1 for the dataset's positive class, 0 for its negative class, as metadata.json names them).
"""
import argparse
import csv
import json
from pathlib import Path

SPLITS = ('train', 'validation')


def prepare_split(dataset_dir: Path, out_dir: Path, split: str, label_to_number: dict[str, int]) -> None:
    source = dataset_dir / split
    target = out_dir / split
    target.mkdir(parents=True)
    (target / 'input').symlink_to((source / 'input').resolve(), target_is_directory=True)
    with (source / 'labels.csv').open(encoding='utf-8', newline='') as f:
        rows = [(row['id'], label_to_number[row['label']]) for row in csv.DictReader(f)]
    with (target / 'labels.csv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(('id', 'numeric_label'))
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument('--dataset-dir', type=Path, required=True,
                        help='the dataset folder with metadata.json, train/ and validation/')
    parser.add_argument('--out-dir', type=Path, required=True, help='a folder that does not exist yet')
    args = parser.parse_args()
    if args.out_dir.exists():
        parser.error(f'{args.out_dir} exists')
    metadata = json.loads((args.dataset_dir / 'metadata.json').read_text(encoding='utf-8'))
    label_to_number = {metadata['negative_class']: 0, metadata['positive_class']: 1}
    for split in SPLITS:
        prepare_split(args.dataset_dir, args.out_dir, split, label_to_number)


if __name__ == '__main__':
    main()
