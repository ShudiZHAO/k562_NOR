#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Link per-fragment consensus sequences into one sequence per assembly.

Input
-----
1. Consensus directory: ``consensus/<assembly>/<k>_consensus.fasta``
   <k> is the fragment index (0-based, sorted numerically).
2. Barcode path file: ``asm/assembly_final_with_empty.txt``
   Each record starts with ``>name``, followed by a space-separated barcode
   path. A path with N barcodes has N+1 fragments:

       fragment 0   -- left of the 1st barcode
       fragment k   -- between barcode k-1 and barcode k
       fragment N   -- right of the last barcode

   ``empty`` is a zero-length placeholder barcode (occupies a slot, no
   sequence) and is handled as the empty string ``''``.

Join rule
---------
Fragments are cut at barcode boundaries and usually already contain the fixed
sequences flanking each barcode. Therefore the linker inserts only the barcode
itself between adjacent fragments (``empty`` inserts nothing). The read
structure is::

    ... left fragment tail + barcode + right fragment head ...   # empty: direct join

Fragment ends can contain a few extra bases or miss a fixed flank because of
consensus start/end positions. The linker does not trim or pad; it only inserts
barcodes in the order given by the path.

Missing fragments
-----------------
Any missing fragment (any index in 0..N) is an error: the script prints the
missing indices and the expected fragment count, exits non-zero, and writes no
output files. The bundled example data is incomplete and will trigger this
error; real data is expected to be complete.

Usage
-----
    python link_consensus.py [consensus_dir] [barcode_file] [-o out_dir]

    # Defaults match the new data layout:
    python link_consensus.py consensus asm/assembly_final_with_empty.txt
"""

import argparse
import os
import re
import sys

from Bio import SeqIO

EMPTY_TOKENS = ('', 'empty', 'nan', 'none')
DEFAULT_CONSENSUS_DIR = 'consensus'
DEFAULT_BARCODE_FILE = os.path.join('asm', 'assembly_final_with_empty.txt')
DEFAULT_OUT_DIR = 'linked_consensus'
LINE_WIDTH = 60

# File names under consensus/<assembly>/: prefer <k>_consensus.fasta,
# then <k>.consensus.fasta, and finally <k>.fasta (does not match <k>.aln.fasta)
NESTED_PATTERNS = (
    re.compile(r'^(\d+)_consensus\.fasta$'),
    re.compile(r'^(\d+)\.consensus\.fasta$'),
    re.compile(r'^(\d+)\.fasta$'),
)
# Legacy flat layout: <assembly>_<k>.fasta / <assembly>_<k>_consensus.fasta
FLAT_PATTERN = r'^%s[._](\d+)(?:[._].*)?\.fasta$'


def normalize_barcode(token):
    """Normalize a barcode: empty/nan/none placeholders -> '', otherwise strip and uppercase."""
    token = str(token).strip()
    if token.lower() in EMPTY_TOKENS:
        return ''
    return token.upper()


def read_barcode_info(barcode_file):
    """Read barcode paths and return {assembly: [barcode, ...]}; empty is stored as ''."""
    barcode_info = {}
    current = None
    with open(barcode_file, 'r', encoding='utf-8-sig') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                fields = line[1:].split()
                name = fields[0] if fields else 'Assembly_%d' % (len(barcode_info) + 1)
                if name in barcode_info:
                    print("[WARN] duplicate assembly %s in %s, ignored" % (name, barcode_file))
                    current = None
                else:
                    current = name
                    barcode_info[name] = []
            elif current is not None:
                barcode_info[current].extend(normalize_barcode(t) for t in line.split())
    return barcode_info


def _add_candidate(found, index, priority, path):
    """Keep the highest-priority candidate (lower number wins) for each fragment."""
    if index not in found or priority < found[index][0]:
        found[index] = (priority, path)


def find_fragment_files(consensus_dir, assembly):
    """Return {fragment_index: fasta_path}, supporting nested and legacy flat layouts."""
    found = {}
    assembly_dir = os.path.join(consensus_dir, assembly)
    if os.path.isdir(assembly_dir):
        for filename in sorted(os.listdir(assembly_dir)):
            for priority, pattern in enumerate(NESTED_PATTERNS):
                match = pattern.match(filename)
                if match:
                    _add_candidate(found, int(match.group(1)), priority,
                                   os.path.join(assembly_dir, filename))
                    break
    elif os.path.isdir(consensus_dir):
        pattern = re.compile(FLAT_PATTERN % re.escape(assembly))
        for filename in sorted(os.listdir(consensus_dir)):
            lowered = filename.lower()
            if any(tag in lowered for tag in ('aln', 'scores', 'depth', 'readstats')):
                continue
            match = pattern.match(filename)
            if match:
                priority = 0 if '_consensus' in filename else 1
                _add_candidate(found, int(match.group(1)), priority,
                               os.path.join(consensus_dir, filename))
    return {index: path for index, (_priority, path) in found.items()}


def read_first_sequence(fasta_file):
    """Return the first sequence in a FASTA file, warning if it has multiple records."""
    sequence = ''
    count = 0
    for record in SeqIO.parse(fasta_file, 'fasta'):
        count += 1
        if count == 1:
            sequence = str(record.seq)
    if count > 1:
        print("[WARN] %s has %d records, only the first one is used" % (fasta_file, count))
    return sequence


def get_fragment_sequences(consensus_dir, assembly):
    """Return all fragment consensus sequences for an assembly as {index: sequence}."""
    file_map = find_fragment_files(consensus_dir, assembly)
    sequences = {}
    for index in sorted(file_map):
        path = file_map[index]
        sequence = read_first_sequence(path)
        if not sequence:
            print("[ERROR] empty consensus file: %s" % path)
            continue
        sequences[index] = sequence
    return sequences


def link_sequences(sequences, barcodes):
    """Link fragments according to a barcode path.

    Return (linked_sequence, info). If any fragment is missing, linked_sequence
    is None and info = {'missing': [...], 'indices': [...],
    'out_of_range': [...], 'junctions': n}.
    """
    expected = len(barcodes) + 1
    usable = {i: s for i, s in sequences.items() if 0 <= i < expected}
    out_of_range = sorted(set(sequences) - set(usable))
    for index in out_of_range:
        print("[WARN] fragment index %d out of range 0-%d, ignored"
              % (index, expected - 1))

    indices = sorted(usable)
    missing = [i for i in range(expected) if i not in usable]
    info = {'missing': missing, 'indices': indices, 'out_of_range': out_of_range,
            'junctions': 0}
    if missing:
        return None, info

    pieces = [usable[0]]
    for previous, current in zip(indices, indices[1:]):
        # previous -> current spans barcode slots previous..current-1
        for slot in range(previous, current):
            pieces.append(barcodes[slot])
            info['junctions'] += 1
        pieces.append(usable[current])
    return ''.join(pieces), info


def write_fasta(output_file, name, sequence, width=LINE_WIDTH):
    """Write one FASTA record with width bases per line."""
    with open(output_file, 'w') as handle:
        handle.write('>%s\n' % name)
        for start in range(0, len(sequence), width):
            handle.write(sequence[start:start + width] + '\n')


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Link fragment consensus sequences according to barcode paths')
    parser.add_argument('consensus_dir', nargs='?', default=DEFAULT_CONSENSUS_DIR,
                        help='consensus directory, layout '
                             'consensus/<assembly>/<k>_consensus.fasta '
                             '(default %s)' % DEFAULT_CONSENSUS_DIR)
    parser.add_argument('barcode_file', nargs='?', default=DEFAULT_BARCODE_FILE,
                        help='barcode path file (default %s)' % DEFAULT_BARCODE_FILE)
    parser.add_argument('-o', '--out-dir', default=DEFAULT_OUT_DIR,
                        help='output directory (default %s)' % DEFAULT_OUT_DIR)
    args = parser.parse_args(argv)

    if not os.path.isdir(args.consensus_dir):
        print('[ERROR] consensus dir not found: %s' % args.consensus_dir)
        return 1
    if not os.path.isfile(args.barcode_file):
        print('[ERROR] barcode file not found: %s' % args.barcode_file)
        return 1

    barcode_info = read_barcode_info(args.barcode_file)
    if not barcode_info:
        print('[ERROR] no assembly found in %s' % args.barcode_file)
        return 1
    print('[INFO] %d assemblies loaded from %s'
          % (len(barcode_info), args.barcode_file))

    # Validate every assembly first; if any fragment is missing, write nothing
    pending = []
    failed = 0
    for assembly, barcodes in barcode_info.items():
        expected = len(barcodes) + 1
        sequences = get_fragment_sequences(args.consensus_dir, assembly)
        linked, info = link_sequences(sequences, barcodes)
        if linked is None:
            print('[ERROR] %s: missing fragment(s) %s; expected %d fragments (0-%d)'
                  % (assembly, info['missing'], expected, expected - 1))
            failed += 1
            continue
        pending.append((assembly, linked, info, expected))

    if failed:
        print('[ERROR] %d assembly(ies) incomplete; no output written' % failed)
        return 1

    os.makedirs(args.out_dir, exist_ok=True)
    for assembly, linked, info, expected in pending:
        output_file = os.path.join(args.out_dir, '%s.fasta' % assembly)
        write_fasta(output_file, assembly, linked)
        print('[INFO] %s: %d/%d fragments, %d barcode junction(s), %d bp -> %s'
              % (assembly, len(info['indices']), expected,
                 info['junctions'], len(linked), output_file))

    print('[INFO] Done: %d linked' % len(pending))
    return 0


if __name__ == '__main__':
    sys.exit(main())

