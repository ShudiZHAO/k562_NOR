import argparse
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def main(argv=None):
    parser = argparse.ArgumentParser(description='divide reads by barcode counts')
    parser.add_argument('--xlsx', default=os.path.join(HERE, 'barcodes_correct_miss_and_mis.xlsx'),
                        help='input xlsx')
    parser.add_argument('--sheet', default='barcodes', help='sheet name')
    parser.add_argument('--outdir', default=HERE, help='outdir')
    args = parser.parse_args(argv)

    if not os.path.isfile(args.xlsx):
        raise SystemExit('[ERROR] can not find xlsx：%s' % args.xlsx)

    frame = pd.read_excel(args.xlsx, sheet_name=args.sheet, dtype=str)
    missing = [c for c in ('read_id', 'barcode_counts') if c not in frame.columns]
    if missing:
        raise SystemExit('[ERROR] xlsx lack column %s，read column：%s' % (missing, list(frame.columns)))

    frame['read_id'] = frame['read_id'].str.strip()
    frame['barcode_counts'] = frame['barcode_counts'].astype(float).astype(int)

    out = {
        'reads_2bc.csv': frame[frame['barcode_counts'] == 2],
        'reads_gt2bc.csv': frame[frame['barcode_counts'] > 2],
    }
    os.makedirs(args.outdir, exist_ok=True)
    print('read total: %d' % frame['read_id'].nunique())
    for name, part in out.items():
        path = os.path.join(args.outdir, name)
        part.to_csv(path, index=False, encoding='utf-8-sig')
        print('%s: %d read / %d rows -> %s' % (name, part['read_id'].nunique(), len(part), path))


if __name__ == '__main__':
    main()
