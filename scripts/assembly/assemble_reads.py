#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import csv
import io
import os
import re
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))

MARK_SUFFIX = {
    (True, True): 'head-tail',
    (True, False): 'head-',
    (False, True): '-tail',
    (False, False): '',
}

DEFAULT_NO_JOIN = 'TGAACA'
DEFAULT_CORRECT_LIST = os.path.join(HERE, 'asm', 'correct_list.csv')
DEFAULT_REMOVE = os.path.join(HERE, 'asm', 'remove.txt')
DEFAULT_INTERVAL_TOL = 0.03   # same junction: lengths within 3% of each other
DEFAULT_INTERVAL_GAP = 500    # ... but never merge lengths more than 500 bp apart

GFA_HEADER = (
    'H\tVN:Z:1.0',
    '# segment = one read; LN = NB = number of barcodes of the read (one barcode is one',
    '#           length unit), BC = barcode path of the read',
    '# link = two reads joined in an assembly; overlap and OV are given in barcodes,',
    '#        OV:Z lists the shared barcode run (the OV:i value equals the CIGAR length)',
    '# path = one assembly; field 3 lists its reads in join order, field 4 the overlap of',
    '#        each neighbouring read pair; FT = head/tail mark, BC = assembly barcode path,',
    '#        NR = number of reads used, CR = reads contained in the assembly but not joined',
)

PATH_GFA_HEADER = (
    'H\tVN:Z:1.0',
    '# path graph: one segment per assembly of assembly.gfa, LN = NB = barcodes of the',
    '#             assembly, BC = its barcode path, GR = overlap group of the segment',
    '#             (assemblies connected by links, the group is the merged assembly)',
    '# link = two assemblies sharing a run of barcodes; the CIGAR and OV:i give the run',
    '#        length in barcodes, OV:Z the shared barcodes, PS:i / QS:i the 0-based start',
    '#        of the run in the first / second assembly (links always start with the',
    '#        lower Assembly number, a path lists its own join order), TY:Z the kind:',
    '#        dovetail (suffix of one = prefix of the other, mergeable end to end),',
    '#        contained (one assembly sits inside the other), internal (the run sits',
    '#        inside both, so the two assemblies cannot be joined linearly)',
    '# path = one merged group; field 3 lists its assemblies in join order, field 4 the',
    '#        overlap of each neighbouring pair, FT = head/tail mark, BC = merged barcode',
    '#        path, NB = barcodes, NA = assemblies in the group, CR = contained assemblies',
)


def ensure_parent(path):
    parent = os.path.dirname(os.path.abspath(path))
    if parent and not os.path.isdir(parent):
        os.makedirs(parent)


def find_column(columns, candidates, given, what):
    lowered = {str(c).strip().lower(): c for c in columns}
    if given:
        if str(given).strip().lower() not in lowered:
            raise SystemExit('[ERROR] %s column %r not found; available: %s'
                             % (what, given, list(columns)))
        return lowered[str(given).strip().lower()]
    for name in candidates:
        if name in lowered:
            return lowered[name]
    raise SystemExit('[ERROR] %s column not found; available: %s' % (what, list(columns)))


def optional_column(columns, candidates):
    lowered = {str(c).strip().lower(): c for c in columns}
    for name in candidates:
        if name in lowered:
            return lowered[name]
    return None


def split_barcodes(text):
    raw = str(text).replace('|', ',').replace(';', ',')
    return [item.strip() for item in raw.split(',') if item.strip()]


def natural_key(text):
    return [int(part) if part.isdigit() else part.lower()
            for part in re.split(r'(\d+)', str(text))]


def index_key(value):
    try:
        return (0, float(value))
    except (TypeError, ValueError):
        return (1, 0.0)


def ends_label(has_head, has_tail):
    if has_head and has_tail:
        return 'head_and_tail'
    if has_head:
        return 'head'
    if has_tail:
        return 'tail'
    return 'no'


def parse_no_join(value):
    items = {item.strip() for item in split_barcodes(value)}
    if not items or {'none', 'no', 'off', '-'} & {item.lower() for item in items}:
        return frozenset()
    return frozenset(items)


def parse_bp(value):
    """Base-pair value of an interval / barcode_start cell ('' or nan -> None)."""
    text = str(value if value is not None else '').strip()
    if not text or text.lower() in ('nan', 'none', '.', 'na'):
        return None
    try:
        number = int(round(float(text)))
    except ValueError:
        return None
    return number if number >= 0 else None


def read_interval_lengths(items):
    """Barcode-to-barcode lengths of one read: [None, d1, d2, ...] in bp.

    items: [(barcode, ends, barcode_index, interval, barcode_start), ...] in
    barcode_index order.  The interval column is used as it is; when a cell is
    missing the difference of the barcode_start cells takes over.  The first
    barcode has no length in front of it, and a read without any usable length
    returns [], so the caller can leave it out of the interval table.
    """
    lengths = []
    for pos, item in enumerate(items):
        value = item[3] if len(item) > 3 else None
        if value is None and len(item) > 4 and item[4] is not None and pos:
            before = items[pos - 1][4]
            if before is not None and item[4] > before:
                value = item[4] - before
        lengths.append(value if value else None)
    return lengths if any(lengths[1:]) else []


def read_reads(path, sep=',', id_column=None, bc_column=None, ends_column=None,
               index_column=None, min_barcodes=2, corrections=None, corr_stats=None,
               interval_column=None, intervals=None):
    """Read the barcode CSV into (read_id, [barcode, ...], ends) rows.

    When `intervals` is a dict it is filled with read_id -> [bp between barcode
    i-1 and i, ...] (None for the first barcode), taken from the interval column
    or, when that one is missing, from the differences of barcode_start.  Reads
    without a usable length are left out of the dict, everything else is
    returned exactly as before.
    """
    if not os.path.isfile(path):
        raise SystemExit('[ERROR] input file not found: %s' % path)
    with io.open(path, encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter=sep)
        columns = reader.fieldnames or []
        id_col = find_column(columns, ('read_id', 'readid', 'read', 'id'),
                             id_column, 'read_id')
        bc_col = find_column(columns, ('barcode', 'barcodes', 'barcode_path', 'path'),
                             bc_column, 'barcode')
        ends_col = (find_column(columns, (), ends_column, 'is_end') if ends_column
                    else optional_column(columns, ('is_end', 'is_ends', 'ends', 'end')))
        idx_col = (find_column(columns, (), index_column, 'barcode_index') if index_column
                   else optional_column(columns, ('barcode_index', 'index', 'bc_index')))
        int_col = (find_column(columns, (), interval_column, 'interval') if interval_column
                   else optional_column(columns, ('interval', 'barcode_interval',
                                                  'bc_interval', 'gap', 'barcode_gap',
                                                  'distance')))
        start_col = optional_column(columns, ('barcode_start', 'start', 'bc_start'))

        order, grouped = [], {}
        for record in reader:
            read_id = str(record.get(id_col) or '').strip()
            if not read_id or read_id.lower() == 'nan':
                continue
            barcode = str(record.get(bc_col) or '').strip()
            if not barcode or barcode.lower() == 'nan':
                continue
            if read_id not in grouped:
                grouped[read_id] = []
                order.append(read_id)
            ends = str(record.get(ends_col) or '').strip().lower() if ends_col else 'no'
            grouped[read_id].append((barcode, ends,
                                     record.get(idx_col) if idx_col else None,
                                     parse_bp(record.get(int_col)) if int_col else None,
                                     parse_bp(record.get(start_col)) if start_col else None))

    rows, dropped = [], 0
    for read_id in order:
        items = grouped[read_id]
        if corrections and read_id in corrections:
            items, changed = apply_barcode_corrections(
                read_id, items, corrections[read_id],
                corr_stats if corr_stats is not None else {})
            if changed and corr_stats is not None:
                corr_stats['reads_changed'] = corr_stats.get('reads_changed', 0) + 1
        if idx_col:
            items.sort(key=lambda item: index_key(item[2]))
        barcodes = [item[0] for item in items]
        if len(barcodes) < min_barcodes:
            dropped += 1
            continue
        has_head = any(item[1] in ('head', 'head_and_tail') for item in items)
        has_tail = any(item[1] in ('tail', 'head_and_tail') for item in items)
        rows.append((read_id, barcodes, ends_label(has_head, has_tail)))
        if intervals is not None:
            lengths = read_interval_lengths(items)
            if lengths:
                intervals[read_id] = lengths
    if corrections:
        stats = corr_stats if corr_stats is not None else {}
        known = set(order)
        missing = [read_id for read_id in corrections if read_id not in known]
        stats['missing_read'] = sum(len(corrections[read_id]) for read_id in missing)
        stats['missing_reads'] = missing
    return rows, id_col, bc_col, ends_col, idx_col, dropped


def index_token(value):
    """Comparable form of a barcode_index cell ('1', 1 and 1.0 give the same key)."""
    text = '' if value is None else str(value).strip()
    if not text:
        return None
    try:
        return ('n', float(text))
    except ValueError:
        return ('s', text.lower())


def read_correct_list(path, sep=','):
    """Read the barcode correction table.

    Required columns: read_id / barcode_index / origin (or barcode) / corrected
    (or new); optional: barcode_start / column.  Only rows with column == 'barcode'
    are kept.  Returns ({read_id: [entry, ...]}, stats).
    """
    stats = {'entries': 0, 'reads': 0, 'skipped': 0, 'warnings': []}
    corrections = {}
    with io.open(path, encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter=sep)
        columns = reader.fieldnames or []
        id_col = find_column(columns, ('read_id', 'readid', 'read', 'id'), None,
                             'correct_list read_id')
        idx_col = find_column(columns, ('barcode_index', 'index', 'bc_index'), None,
                              'correct_list barcode_index')
        origin_col = find_column(columns, ('origin', 'barcode', 'old', 'from'), None,
                                 'correct_list origin')
        new_col = find_column(columns, ('corrected', 'correct', 'new', 'to'), None,
                              'correct_list corrected')
        start_col = optional_column(columns, ('barcode_start', 'start'))
        target_col = optional_column(columns, ('column', 'field', 'target'))

        for record in reader:
            read_id = str(record.get(id_col) or '').strip()
            index = str(record.get(idx_col) or '').strip()
            origin = str(record.get(origin_col) or '').strip()
            corrected = str(record.get(new_col) or '').strip()
            if not any((read_id, index, origin, corrected)):
                continue
            target = (str(record.get(target_col) or '').strip().lower()
                      if target_col else 'barcode')
            if target and target != 'barcode':
                stats['skipped'] += 1
                stats['warnings'].append(
                    'read %s barcode_index %s: column %r is not supported, entry skipped'
                    % (read_id, index, target))
                continue
            if not read_id or not index or not origin or not corrected:
                stats['skipped'] += 1
                stats['warnings'].append(
                    'read %s barcode_index %s: incomplete entry (origin=%r, corrected=%r), '
                    'skipped' % (read_id, index, origin, corrected))
                continue
            entries = corrections.setdefault(read_id, [])
            token = index_token(index)
            if any(index_token(item['barcode_index']) == token for item in entries):
                stats['skipped'] += 1
                stats['warnings'].append(
                    'read %s barcode_index %s: listed twice, kept the first entry'
                    % (read_id, index))
                continue
            entries.append({
                'barcode_index': index,
                'origin': origin,
                'corrected': corrected,
                'barcode_start': (str(record.get(start_col) or '').strip()
                                  if start_col else ''),
            })
            stats['entries'] += 1
    stats['reads'] = len(corrections)
    return corrections, stats


def apply_barcode_corrections(read_id, items, entries, stats):
    """Replace the barcodes of one read according to the correct list.

    items: [(barcode, ends, barcode_index, interval, barcode_start), ...] straight
    from the input CSV (the trailing fields are kept untouched).
    An entry is applied only when the read has that barcode_index and the current
    barcode equals the recorded origin; otherwise a warning is recorded and the
    entry is skipped, so a stale correct list never rewrites a different barcode.
    """
    positions = {}
    repeated = set()
    for pos, item in enumerate(items):
        token = index_token(item[2])
        if token is None:
            continue
        if token in positions:
            repeated.add(token)
        else:
            positions[token] = pos
    changed = 0
    for entry in entries:
        label = 'read %s barcode_index %s' % (read_id, entry['barcode_index'])
        if entry.get('barcode_start'):
            label += ' (barcode_start %s)' % entry['barcode_start']
        token = index_token(entry['barcode_index'])
        if token in repeated:
            stats['missing_index'] = stats.get('missing_index', 0) + 1
            stats.setdefault('warnings', []).append(
                '%s: barcode_index is not unique in the input, entry skipped' % label)
            continue
        pos = positions.get(token)
        if pos is None:
            stats['missing_index'] = stats.get('missing_index', 0) + 1
            stats.setdefault('warnings', []).append(
                '%s: no such barcode_index in the input, entry skipped' % label)
            continue
        barcode = items[pos][0]
        if barcode != entry['origin']:
            stats['mismatch'] = stats.get('mismatch', 0) + 1
            stats.setdefault('warnings', []).append(
                '%s: input barcode is %r but the correct list expects %r, entry skipped'
                % (label, barcode, entry['origin']))
            continue
        if barcode == entry['corrected']:
            stats['same'] = stats.get('same', 0) + 1
            continue
        items[pos] = (entry['corrected'],) + tuple(items[pos][1:])
        stats['applied'] = stats.get('applied', 0) + 1
        changed += 1
    return items, changed


def read_remove_list(path):
    """Read the read ids that must be left out of the assembly.

    One id per line; blank lines and lines starting with '#' are skipped, and
    several ids per line may be separated by commas or whitespace.  The ids are
    compared with the read_id column of the input as they are written, so
    '14220_minus' and '14220' are different reads.
    Returns (ids, stats) with stats['lines'] and stats['duplicates'].
    """
    ids, stats = set(), {'lines': 0, 'duplicates': 0}
    with io.open(path, encoding='utf-8-sig') as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith('#'):
                continue
            stats['lines'] += 1
            for item in re.split(r'[,\s]+', text):
                item = item.strip()
                if not item:
                    continue
                if item in ids:
                    stats['duplicates'] += 1
                    continue
                ids.add(item)
    return ids, stats


def drop_reads(rows, remove_ids):
    """Split read rows into (kept rows, removed read ids)."""
    if not remove_ids:
        return rows, []
    kept, removed = [], []
    for row in rows:
        if str(row[0]) in remove_ids:
            removed.append(row[0])
        else:
            kept.append(row)
    return kept, removed


def load_marks(rows):
    marks = {}
    for read_id, _barcodes, ends in rows:
        flags = set()
        if ends in ('head', 'head_and_tail'):
            flags.add('head')
        if ends in ('tail', 'head_and_tail'):
            flags.add('tail')
        if flags:
            marks[read_id] = flags
    return marks


def suffix_prefix_overlap(longer, shorter):
    limit = min(len(longer), len(shorter))
    for k in range(limit, 0, -1):
        if tuple(longer[-k:]) == tuple(shorter[:k]):
            return k
    return 0


def shared_runs(left, right, min_length=1, no_join=frozenset()):
    """All maximal exact barcode runs shared by two paths: (length, start_left, start_right).

    A run is reported once, at the last cell of the run, so only runs that cannot be
    extended to the right are kept; left extension is already part of the run length.
    Runs made only of no_join barcodes are skipped, like they are for merging.
    """
    rows, cols = len(left), len(right)
    previous = [0] * (cols + 1)
    runs = set()
    for i in range(1, rows + 1):
        current = [0] * (cols + 1)
        base = left[i - 1]
        for j in range(1, cols + 1):
            if base != right[j - 1]:
                continue
            length = previous[j - 1] + 1
            current[j] = length
            if length < min_length:
                continue
            if i < rows and j < cols and left[i] == right[j]:
                continue
            start_left, start_right = i - length, j - length
            if no_join and all(bc in no_join for bc in left[start_left:i]):
                continue
            runs.add((length, start_left, start_right))
        previous = current
    kept = []
    for length, start_left, start_right in sorted(runs, reverse=True):
        if any(k2 > length and s2 <= start_left and t2 <= start_right
               and s2 + k2 >= start_left + length and t2 + k2 >= start_right + length
               for k2, s2, t2 in kept):
            continue
        kept.append((length, start_left, start_right))
    return sorted(kept, key=lambda item: (-item[0], item[1], item[2]))


def overlap_kind(left, right, length, start_left, start_right):
    """Classify the overlap of two paths: dovetail, contained or internal."""
    if start_left + length == len(left) and start_right == 0:
        return 'dovetail'
    if start_right + length == len(right) and start_left == 0:
        return 'dovetail'
    if length >= len(left) or length >= len(right):
        return 'contained'
    return 'internal'


def tag_value(tags, name, kind=None):
    """Value of a GFA tag (kind None accepts any type character)."""
    for tag in tags:
        parts = str(tag).split(':', 2)
        if len(parts) == 3 and parts[0] == name and (kind is None or parts[1] == kind):
            return parts[2]
    return None


def interval_scope_applies(overlap, min_overlap, scope):
    """Does the interval clustering decide this join?  scope: min / all / off."""
    if scope == 'off' or overlap < max(min_overlap, 1):
        return False
    if scope == 'min':
        return overlap == max(min_overlap, 1)
    return True


def interval_clusters(values, tol=DEFAULT_INTERVAL_TOL, min_gap=DEFAULT_INTERVAL_GAP):
    """Group interval lengths (bp) into clusters: [(low, high, count), ...].

    The lengths are sorted and a new cluster starts as soon as two neighbouring
    lengths are further apart than max(min_gap, tol * previous length), so the
    sequencing noise of one junction (a few hundred bp) stays in one cluster
    while a second junction using the same two barcodes (thousands of bp away)
    gets its own cluster.
    """
    numbers = sorted(int(value) for value in values if value)
    clusters = []
    for value in numbers:
        if clusters:
            low, high, count = clusters[-1]
            if value - high <= max(min_gap, tol * high):
                clusters[-1] = (low, value, count + 1)
                continue
        clusters.append((value, value, 1))
    return clusters


class IntervalClusters:
    """Clusters of the barcode-to-barcode lengths of every barcode pair.

    Two barcodes that occur together in several reads are not necessarily the same
    junction of the chromosome: the DNA between them has a length, reads of one
    junction agree on it, reads of another junction differ by thousands of bp.
    Clustering those lengths separates the junctions, so an ambiguous join
    (ABC + BCD or ABC + BCE, both sharing BC) is decided by the length the reads
    measure between B and C instead of by the order of the input.
    """

    def __init__(self, read_paths, read_intervals, tol=DEFAULT_INTERVAL_TOL,
                 min_gap=DEFAULT_INTERVAL_GAP):
        self.tol = tol
        self.min_gap = min_gap
        self.reads = {read_id: list(lengths)
                      for read_id, lengths in read_intervals.items()
                      if read_id in read_paths and lengths}
        lengths = {}
        for read_id, values in self.reads.items():
            path = read_paths[read_id]
            for pos in range(min(len(path), len(values)) - 1):
                value = values[pos + 1]
                if value:
                    lengths.setdefault((path[pos], path[pos + 1]), []).append(value)
        self.pairs = {pair: interval_clusters(values, tol, min_gap)
                      for pair, values in lengths.items()}

    def cluster(self, pair, value):
        """Cluster number of one length of a barcode pair (None when unknown)."""
        clusters = self.pairs.get(pair)
        if not clusters or not value:
            return None
        for number, (low, high, _count) in enumerate(clusters):
            limit = max(self.min_gap, self.tol * max(low, high))
            if low - limit <= value <= high + limit:
                return number
        return None

    def path_cluster(self, path, prov, pos):
        """Cluster of the pair path[pos]/path[pos+1] as voted by the reads of prov.

        Returns (cluster, length, votes); (None, None, 0) when no read of prov
        spans that pair, and then the join is not decided by interval lengths.
        """
        pair = (path[pos], path[pos + 1])
        votes = {}
        for read_id, offset in prov:
            values = self.reads.get(read_id)
            if not values:
                continue
            index = pos - offset + 1      # length in front of barcode `index`
            if index <= 0 or index >= len(values):
                continue
            value = values[index]
            number = self.cluster(pair, value) if value else None
            if number is None:
                continue
            votes.setdefault(number, []).append(value)
        if not votes:
            return None, None, 0
        number, found = sorted(votes.items(), key=lambda item: (-len(item[1]), item[0]))[0]
        return number, sorted(found)[len(found) // 2], len(found)

    def describe(self, pair, value):
        """'44771 bp (cluster 1/2: 44724-44805)' for log lines."""
        clusters = self.pairs.get(pair) or []
        number = self.cluster(pair, value)
        if number is None:
            return '%s bp (no cluster)' % value
        low, high, _count = clusters[number]
        return '%s bp (cluster %d/%d: %d-%d)' % (value, number + 1, len(clusters),
                                                 low, high)


def interval_check(clusters, left_path, left_prov, left_pos, right_path, right_prov,
                   right_pos, overlap):
    """Compare the lengths of every barcode pair both paths share.

    Returns (mismatches, matches): mismatches lists (pair, left_length,
    right_length) for pairs whose clusters differ (the two paths then measure a
    different DNA piece between the same barcodes, so they cannot be joined),
    matches lists the pairs that were compared and agreed.
    """
    mismatches, matches = [], []
    for step in range(max(overlap - 1, 0)):
        pair = (left_path[left_pos + step], left_path[left_pos + step + 1])
        left_id, left_value, _left_votes = clusters.path_cluster(
            left_path, left_prov, left_pos + step)
        right_id, right_value, _right_votes = clusters.path_cluster(
            right_path, right_prov, right_pos + step)
        if left_id is None or right_id is None:
            continue
        if left_id == right_id:
            matches.append((pair, left_value, right_value))
        else:
            mismatches.append((pair, left_value, right_value))
    return mismatches, matches


def build_overlap_assemblies(tracks, min_overlap=2, no_join=frozenset(), verbose=False,
                             stats=None, intervals=None, interval_scope='min'):
    """Join the paths end to end through their suffix/prefix barcode run.

    intervals: optional IntervalClusters.  When it is given, the shared barcode
    run must measure the same distances in both paths (see interval_check),
    otherwise the join is skipped and another candidate can be used.
    interval_scope: 'min' checks only overlaps of exactly min_overlap barcodes
    (the weakest evidence, where the same two barcodes can point to several
    continuations), 'all' checks every overlap, 'off' never checks.
    """
    if stats is None:
        stats = {}
    stats.setdefault('blocked', 0)
    stats.setdefault('blocked_examples', [])
    stats.setdefault('merges', [])
    stats.setdefault('interval_keys', set())
    stats.setdefault('interval_examples', [])
    stats.setdefault('interval_confirmed', 0)

    remaining = [(tuple(path), tuple(prov)) for path, prov in tracks]
    remaining.sort(key=lambda item: len(item[0]), reverse=True)
    merged = []
    i = 0
    while i < len(remaining):
        current, current_prov = remaining[i]
        best_j, best_ov, best_match = None, 0, ()
        blocked_ov, blocked_j = 0, None
        for j, (other, other_prov) in enumerate(remaining):
            if j == i:
                continue
            ov = suffix_prefix_overlap(current, other)
            if ov <= 0:
                continue
            region = tuple(current[-ov:])
            if no_join and all(bc in no_join for bc in region):
                if ov > blocked_ov:
                    blocked_ov, blocked_j = ov, j
                continue
            matches = ()
            if intervals is not None and interval_scope_applies(ov, min_overlap,
                                                                interval_scope):
                mismatches, matches = interval_check(
                    intervals, current, current_prov, len(current) - ov,
                    other, other_prov, 0, ov)
                if mismatches:
                    key = (current, other, ov)
                    if key not in stats['interval_keys']:
                        stats['interval_keys'].add(key)
                        if len(stats['interval_examples']) < 5:
                            stats['interval_examples'].append(
                                (current, other, ov, mismatches[0]))
                    if verbose:
                        print('[BLOCK] %s + %s: %d shared barcodes but %s, not joined'
                              % (' '.join(current), ' '.join(other), ov,
                                 ', '.join('%s %s vs %s'
                                           % ('|'.join(pair),
                                              intervals.describe(pair, left_length),
                                              intervals.describe(pair, right_length))
                                           for pair, left_length, right_length
                                           in mismatches)))
                    continue
            if ov > best_ov or (ov == best_ov and len(matches) > len(best_match)):
                best_ov, best_j, best_match = ov, j, matches
        if best_ov >= min_overlap:
            other, other_prov = remaining[best_j]
            joined = current + other[best_ov:]
            shift = len(current) - best_ov
            joined_prov = current_prov + tuple((read_id, offset + shift)
                                               for read_id, offset in other_prov)
            merged.append((joined, joined_prov))
            stats['merges'].append((best_ov, tuple(current[-best_ov:])))
            if best_match:
                stats['interval_confirmed'] += 1
            if verbose:
                extra = (', %s %s ~ %s' % ('|'.join(best_match[0][0]),
                                           best_match[0][1], best_match[0][2])
                         if best_match else '')
                print('[MERGE] %s + %s (overlap %d%s) -> %d barcodes'
                      % (' '.join(current), ' '.join(other), best_ov, extra,
                         len(joined)))
            for idx in sorted((i, best_j), reverse=True):
                remaining.pop(idx)
        else:
            if blocked_j is not None and blocked_ov >= min_overlap:
                stats['blocked'] += 1
                if len(stats['blocked_examples']) < 5:
                    stats['blocked_examples'].append(
                        (current, remaining[blocked_j][0], blocked_ov))
                if verbose:
                    print('[SKIP] %s + %s: overlap region is only %s (overlap %d), not joined'
                          % (' '.join(current), ' '.join(remaining[blocked_j][0]),
                             '|'.join(current[-blocked_ov:]), blocked_ov))
            i += 1
    return merged + remaining


def remove_contained_paths(tracks, verbose=False):
    unique, index = [], {}
    for path, prov in tracks:
        key = tuple(path)
        if key not in index:
            index[key] = len(unique)
            unique.append((key, list(prov)))
        else:
            merged_prov, known = unique[index[key]][1], {rid for rid, _off in prov}
            merged_prov.extend(item for item in prov if item[0] not in known)
    unique.sort(key=lambda item: (-len(item[0]), item[0]))

    kept = []
    for i, (path, prov) in enumerate(unique):
        text = '|' + '|'.join(path) + '|'
        dropped = None
        for j, (other, _other_prov) in enumerate(unique):
            if i == j or len(other) <= len(path):
                continue
            if text in ('|' + '|'.join(other) + '|'):
                dropped = other
                break
        if dropped is None:
            kept.append((path, tuple(prov)))
        elif verbose:
            print('[DROP] %s is contained in %s' % (' '.join(path), ' '.join(dropped)))
    return kept


def path_positions(assembly, path):
    text = '|' + '|'.join(assembly) + '|'
    key = '|' + '|'.join(path) + '|'
    hits = []
    pos = text.find(key)
    while pos >= 0:
        start = text.count('|', 0, pos)
        hits.append((start, start + len(path)))
        pos = text.find(key, pos + 1)
    return hits


def annotate_assemblies(assemblies, reads, marks):
    names = ['Assembly_%d' % (i + 1) for i in range(len(assemblies))]
    flags = {name: [False, False] for name in names}
    internal = Counter()
    for read_id, path, _ends in reads:
        tags = marks.get(read_id)
        if not tags:
            continue
        for idx, assembly in enumerate(assemblies):
            hits = path_positions(assembly, path)
            if not hits:
                continue
            if 'head' in tags:
                if any(start == 0 for start, _end in hits):
                    flags[names[idx]][0] = True
                else:
                    internal['head'] += 1
            if 'tail' in tags:
                if any(end == len(assembly) for _start, end in hits):
                    flags[names[idx]][1] = True
                else:
                    internal['tail'] += 1
    return {name: (value[0], value[1]) for name, value in flags.items()}, internal


def write_paths(assemblies, out_file, marks):
    ensure_parent(out_file)
    with io.open(out_file, 'w', encoding='utf-8', newline='') as handle:
        for i, path in enumerate(assemblies, 1):
            name = 'Assembly_%d' % i
            suffix = MARK_SUFFIX.get(marks.get(name, (False, False)), '')
            handle.write('>%s%s\n' % (name, (' ' + suffix) if suffix else ''))
            handle.write(' '.join(path) + '\n')
    return len(assemblies)


def advance_chain(total, entries, min_overlap=1, label='segments'):
    """Pick a chain of (name, start, end) entries that walks a path of `total` barcodes.

    Entries are taken from left to right, always the one reaching furthest, and only
    when it overlaps the current end by at least min_overlap (falling back to 1 and 0
    when nothing else is left). Entries that do not advance the walk are returned
    separately, they are contained in the path but not part of its chain.
    """
    entries = sorted(entries, key=lambda item: (item[1], -item[2], natural_key(item[0])))
    chain, pending, pos = [], list(entries), 0
    tiers = [max(min_overlap, 1), 1, 0]
    while pos < total:
        chosen = None
        for need in tiers:
            candidates = [item for item in pending
                          if pos - item[1] >= need and item[2] > pos]
            if candidates:
                chosen = min(candidates, key=lambda item: (-item[2], item[1]))
                break
        if chosen is None:
            raise SystemExit('[ERROR] %s cannot rebuild the path at barcode %d'
                             % (label, pos))
        chain.append(chosen)
        pos = chosen[2]
        pending.remove(chosen)
    used = {item[0] for item in chain}
    leftovers = [item for item in pending if item[0] not in used]
    return chain, leftovers


def assembly_chain(path, prov, read_paths, min_overlap=1):
    entries = []
    for read_id, offset in prov:
        barcodes = read_paths[read_id]
        if list(path[offset:offset + len(barcodes)]) != list(barcodes):
            raise SystemExit('[ERROR] read %s does not sit at offset %d of its assembly'
                             % (read_id, offset))
        entries.append((read_id, offset, offset + len(barcodes)))
    return advance_chain(len(path), entries, min_overlap, 'reads')



def write_gfa(paths_file, tracks, reads, out_file, verbose=False, min_overlap=1):
    ensure_parent(out_file)
    order, marks = parse_paths_file(paths_file)
    prov_by_path = {tuple(path): tuple(prov) for path, prov in tracks}
    read_paths = {read_id: path for read_id, path, _ends in reads}

    rows, segments, links = [], {}, set()
    for name, path in order:
        prov = prov_by_path.get(tuple(path))
        if prov is None:
            raise SystemExit('[ERROR] path.txt entry %s is not part of this run' % name)
        chain, leftovers = assembly_chain(path, prov, read_paths, min_overlap)
        overlaps = [left[2] - right[1] for left, right in zip(chain, chain[1:])]
        chained = {read_id for read_id, _start, _end in chain}
        contained = sorted({entry[0] for entry in leftovers} |
                           {read_id for read_id, read_path in read_paths.items()
                            if read_id not in chained and path_positions(path, read_path)},
                           key=natural_key)
        for read_id in chained:
            segments[read_id] = read_paths[read_id]
        for left, right in zip(chain, chain[1:]):
            ov = left[2] - right[1]
            links.add((left[0], '+', right[0], '+', '%dM' % ov, 'OV:i:%d' % ov,
                       'OV:Z:%s' % '|'.join(path[right[1]:left[2]])))
        rows.append((name, path, chain, overlaps, contained))
        if verbose:
            print('[PATH] %s: %d barcodes, %d reads, overlaps [%s]%s'
                  % (name, len(path), len(chain),
                     ','.join(str(ov) for ov in overlaps),
                     ', %d contained reads' % len(contained) if contained else ''))
            print('       reads: %s' % ' '.join(entry[0] for entry in chain))

    with io.open(out_file, 'w', encoding='utf-8', newline='') as handle:
        for line in GFA_HEADER:
            handle.write(line + '\n')
        for read_id in sorted(segments, key=natural_key):
            barcodes = segments[read_id]
            handle.write('S\t%s\t*\tLN:i:%d\tNB:i:%d\tBC:Z:%s\n'
                         % (read_id, len(barcodes), len(barcodes), '|'.join(barcodes)))
        for link in sorted(links):
            handle.write('L\t%s\n' % '\t'.join(link))
        for name, path, chain, overlaps, contained in rows:
            fields = ['P', name,
                      ','.join('%s+' % entry[0] for entry in chain),
                      ','.join('%dM' % ov for ov in overlaps) or '*']
            mark = marks.get(name, '')
            if mark:
                fields.append('FT:Z:%s' % mark)
            fields.append('BC:Z:%s' % '|'.join(path))
            fields.append('NR:i:%d' % len(chain))
            if contained:
                fields.append('CR:Z:%s' % ','.join(contained))
            handle.write('\t'.join(fields) + '\n')
    return len(segments), len(links), rows


def overlap_groups(names, edges):
    """Group assemblies that share a barcode run (connected components, largest first)."""
    parent = {name: name for name in names}

    def find(name):
        while parent[name] != name:
            parent[name] = parent[parent[name]]
            name = parent[name]
        return name

    for left, right, _length, _start_left, _start_right, _kind in edges:
        a, b = find(left), find(right)
        if a != b:
            parent[a] = b
    groups = {}
    for name in names:
        groups.setdefault(find(name), []).append(name)
    return sorted((sorted(group, key=natural_key) for group in groups.values()),
                  key=lambda group: (-len(group), natural_key(group[0])))


def write_path_gfa(assemblies, marks, out_file, min_overlap=2, no_join=frozenset(),
                   verbose=False):
    """Merge the assemblies of assembly.gfa and show how they overlap.

    Segment = one assembly (its barcode path in the BC tag), link = two assemblies
    sharing a run of barcodes, path = assemblies joined end to end through dovetail
    overlaps. Returns (segment count, edges, path rows, overlap groups).
    """
    ensure_parent(out_file)
    names = ['Assembly_%d' % (i + 1) for i in range(len(assemblies))]
    paths = {name: tuple(path) for name, (path, _prov) in zip(names, assemblies)}
    sizes = {name: len(path) for name, path in paths.items()}

    edges = []
    for i, left in enumerate(names):
        for right in names[i + 1:]:
            for length, start_left, start_right in shared_runs(
                    paths[left], paths[right], min_overlap, no_join):
                edges.append((left, right, length, start_left, start_right,
                              overlap_kind(paths[left], paths[right], length,
                                           start_left, start_right)))

    tracks = [(paths[name], ((name, 0),)) for name in names]
    while True:
        merged = build_overlap_assemblies(tracks, min_overlap, no_join)
        done = len(merged) >= len(tracks)
        tracks = merged
        if done:
            break
    grouped = {}
    for path, prov in tracks:
        grouped.setdefault(tuple(path), []).extend(prov)
    tracks = sorted(((path, tuple(prov)) for path, prov in grouped.items()),
                    key=lambda item: (-len(item[0]), item[0]))

    rows = []
    edge_keys = set(edges)
    for index, (path, prov) in enumerate(tracks, 1):
        for name, offset in prov:
            if tuple(path[offset:offset + sizes[name]]) != paths[name]:
                raise SystemExit('[ERROR] assembly %s does not sit at offset %d of its '
                                 'merged path' % (name, offset))
        entries = [(name, offset, offset + sizes[name]) for name, offset in prov]
        chain, leftovers = advance_chain(len(path), entries, min_overlap, 'assemblies')
        overlaps = []
        for left, right in zip(chain, chain[1:]):
            ov = left[2] - right[1]
            overlaps.append(ov)
            if ov <= 0:
                continue
            start_left, start_right = right[1] - left[1], 0
            if tuple(paths[left[0]][start_left:start_left + ov]) != \
                    tuple(paths[right[0]][:ov]):
                raise SystemExit('[ERROR] %s and %s are joined without a shared run'
                                 % (left[0], right[0]))
            if natural_key(left[0]) > natural_key(right[0]):
                key = (right[0], left[0], ov, start_right, start_left)
            else:
                key = (left[0], right[0], ov, start_left, start_right)
            key = key + (overlap_kind(paths[key[0]], paths[key[1]], ov,
                                      key[3], key[4]),)
            if key not in edge_keys:
                edge_keys.add(key)
                edges.append(key)
        rows.append(('Path_%d' % index, path, chain, overlaps,
                     [entry[0] for entry in leftovers]))

    edges.sort(key=lambda item: (natural_key(item[0]), natural_key(item[1]),
                                 -item[2], item[3]))
    groups = overlap_groups(names, edges)
    group_of = {name: 'Group_%d' % (i + 1)
                for i, group in enumerate(groups) for name in group}
    if verbose:
        for left, right, length, start_left, start_right, kind in sorted(
                edges, key=lambda item: (-item[2], natural_key(item[0]),
                                         natural_key(item[1]))):
            print('[LINK] %s ~ %s: %d barcodes at %d/%d (%s) %s'
                  % (left, right, length, start_left, start_right, kind,
                     '|'.join(paths[left][start_left:start_left + length])))
        for name, path, chain, overlaps, contained in rows:
            if len(chain) + len(contained) > 1:
                print('[GROUP] %s: %s (overlaps [%s])%s'
                      % (name, ' + '.join(entry[0] for entry in chain),
                         ','.join(str(ov) for ov in overlaps),
                         ', %d contained' % len(contained) if contained else ''))
        for group in groups:
            if len(group) > 1:
                print('[OVERLAP] %s: %d assemblies: %s'
                      % (group_of[group[0]], len(group), ' '.join(group)))

    with io.open(out_file, 'w', encoding='utf-8', newline='') as handle:
        for line in PATH_GFA_HEADER:
            handle.write(line + '\n')
        for name in names:
            handle.write('S\t%s\t*\tLN:i:%d\tNB:i:%d\tBC:Z:%s\tGR:Z:%s\n'
                         % (name, sizes[name], sizes[name], '|'.join(paths[name]),
                            group_of[name]))
        for left, right, length, start_left, start_right, kind in edges:
            handle.write('L\t%s\t+\t%s\t+\t%dM\tOV:i:%d\tOV:Z:%s\tPS:i:%d\tQS:i:%d'
                         '\tTY:Z:%s\n'
                         % (left, right, length, length,
                            '|'.join(paths[left][start_left:start_left + length]),
                            start_left, start_right, kind))
        for name, path, chain, overlaps, contained in rows:
            fields = ['P', name,
                      ','.join('%s+' % entry[0] for entry in chain),
                      ','.join('%dM' % ov for ov in overlaps) or '*']
            head = chain[0][1] == 0 and marks.get(chain[0][0], (False, False))[0]
            tail = chain[-1][2] == len(path) and marks.get(chain[-1][0], (False, False))[1]
            suffix = MARK_SUFFIX.get((bool(head), bool(tail)), '')
            if suffix:
                fields.append('FT:Z:%s' % suffix)
            fields.append('BC:Z:%s' % '|'.join(path))
            fields.append('NB:i:%d' % len(path))
            fields.append('NA:i:%d' % (len(chain) + len(contained)))
            if contained:
                fields.append('CR:Z:%s' % ','.join(contained))
            handle.write('\t'.join(fields) + '\n')
    return len(names), edges, rows, groups


def write_mapping(reads, assemblies, out_file):
    ensure_parent(out_file)
    mapping = {}
    for read_id, path, ends in reads:
        for i, assembly in enumerate(assemblies, 1):
            hits = path_positions(assembly, path)
            if hits:
                mapping[read_id] = ('Assembly_%d' % i, hits[0][0], hits[0][1], ends)
                break
    with io.open(out_file, 'w', encoding='utf-8', newline='') as handle:
        handle.write('read_id\tassembly\tstart_index\tend_index\tn_barcodes\tis_ends\n')
        for read_id in sorted(mapping, key=natural_key):
            name, start, end, ends = mapping[read_id]
            handle.write('%s\t%s\t%d\t%d\t%d\t%s\n'
                         % (read_id, name, start, end, end - start, ends))
    return mapping


def parse_paths_file(path):
    order, marks = [], {}
    name = None
    with io.open(path, encoding='utf-8') as handle:
        for line in handle:
            line = line.rstrip('\r\n')
            if not line:
                continue
            if line.startswith('>'):
                parts = line[1:].strip().split()
                name = parts[0] if parts else ''
                order.append((name, []))
                marks[name] = parts[1] if len(parts) > 1 else ''
            elif name is not None:
                order[-1] = (name, line.split())
    return order, marks


def parse_gfa_file(path):
    segments, links, prows = {}, set(), {}
    with io.open(path, encoding='utf-8') as handle:
        for line in handle:
            parts = line.rstrip('\r\n').split('\t')
            if not parts or parts[0] in ('H', '') or parts[0].startswith('#'):
                continue
            if parts[0] == 'S':
                segments[parts[1]] = (parts[2], parts[3:])
            elif parts[0] == 'L':
                links.add(tuple(parts[1:]))
            elif parts[0] == 'P':
                segs = [seg[:-1] for seg in parts[2].split(',') if seg]
                overlaps = parts[3].split(',') if parts[3] not in ('*', '') else []
                prows[parts[1]] = (segs, overlaps, parts[4:])
    return segments, links, prows


def check_outputs(reads, tracks, marks, paths_file, gfa_file, min_barcodes, min_overlap=1):
    problems = []
    expected_names = ['Assembly_%d' % (i + 1) for i in range(len(tracks))]
    read_paths = {read_id: path for read_id, path, _ends in reads}

    order, path_marks = parse_paths_file(paths_file)
    names = [name for name, _path in order]
    if names != expected_names:
        problems.append('path.txt names are not consecutive %s (got %d: %s)'
                        % (expected_names[:1], len(names), names[:3]))
    for i, (name, path) in enumerate(order):
        if i >= len(tracks):
            problems.append('path.txt holds %d assemblies, memory holds %d'
                            % (len(order), len(tracks)))
            break
        if tuple(path) != tuple(tracks[i][0]):
            problems.append('path.txt path of %s differs from memory' % name)
        want = MARK_SUFFIX.get(marks.get(name, (False, False)), '')
        got = path_marks.get(name, '')
        if got != want:
            problems.append('path.txt mark of %s should be %r, got %r' % (name, want, got))

    want_segments, want_links, want_paths = {}, set(), {}
    for i, (path, prov) in enumerate(tracks, 1):
        name = 'Assembly_%d' % i
        try:
            chain, leftovers = assembly_chain(path, prov, read_paths, min_overlap)
        except SystemExit as exc:
            problems.append(str(exc))
            continue
        chained = [entry[0] for entry in chain]
        chained_set = set(chained)
        overlaps = [left[2] - right[1] for left, right in zip(chain, chain[1:])]
        if any(ov < 0 for ov in overlaps):
            problems.append('assembly %s has a read pair with a negative overlap' % name)
        for read_id in chained:
            want_segments[read_id] = read_paths[read_id]
        for left, right in zip(chain, chain[1:]):
            ov = left[2] - right[1]
            want_links.add((left[0], '+', right[0], '+', '%dM' % ov, 'OV:i:%d' % ov,
                            'OV:Z:%s' % '|'.join(path[right[1]:left[2]])))
        contained = sorted({entry[0] for entry in leftovers} |
                           {read_id for read_id, read_path in read_paths.items()
                            if read_id not in chained_set and path_positions(path, read_path)},
                           key=natural_key)
        tags = []
        mark = MARK_SUFFIX.get(marks.get(name, (False, False)), '')
        if mark:
            tags.append('FT:Z:%s' % mark)
        tags.append('BC:Z:%s' % '|'.join(path))
        tags.append('NR:i:%d' % len(chain))
        if contained:
            tags.append('CR:Z:%s' % ','.join(contained))
        want_paths[name] = (chained, ['%dM' % ov for ov in overlaps], tags, contained)

    segments, links, prows = parse_gfa_file(gfa_file)
    if set(segments) != set(want_segments):
        problems.append('GFA S lines differ from the reads used for joining '
                        '(missing %d / extra %d)'
                        % (len(set(want_segments) - set(segments)),
                           len(set(segments) - set(want_segments))))
    for read_id, (seq, tags) in segments.items():
        if read_id not in read_paths:
            continue
        barcodes = read_paths[read_id]
        want = ['LN:i:%d' % len(barcodes), 'NB:i:%d' % len(barcodes),
                'BC:Z:%s' % '|'.join(barcodes)]
        if seq != '*' or tags != want:
            problems.append('GFA S line %s should carry * and %s, got %s and %s'
                            % (read_id, want, seq, tags))
    if links != want_links:
        problems.append('GFA L lines differ from the read joins (missing %d / extra %d)'
                        % (len(want_links - links), len(links - want_links)))
    if list(prows) != expected_names:
        problems.append('GFA P lines differ in count/order (%d vs %d)'
                        % (len(prows), len(tracks)))
    for name in expected_names:
        if name not in want_paths:
            continue
        chained, overlaps, tags, _contained = want_paths[name]
        segs, got_overlaps, got_tags = prows.get(name, ([], [], []))
        if segs != chained:
            problems.append('GFA P line %s lists %d reads, the join chain has %d'
                            % (name, len(segs), len(chained)))
        if got_overlaps != overlaps:
            problems.append('GFA P line %s overlaps are %s, expected %s'
                            % (name, got_overlaps, overlaps))
        if got_tags != tags:
            problems.append('GFA P line %s tags are %s, expected %s'
                            % (name, got_tags, tags))

    chained_all = {read_id for _path, prov in tracks for read_id, _offset in prov}
    listed_all = {read_id for _c, _o, _t, contained in want_paths.values()
                  for read_id in contained}
    missing = [read_id for read_id, _path, _ends in reads
               if read_id not in chained_all and read_id not in listed_all]
    if missing:
        problems.append('%d reads are neither joined nor listed as contained, e.g. %s'
                        % (len(missing), ', '.join(missing[:5])))
    uncovered = [read_id for read_id, path, _ends in reads
                 if not any(path_positions(assembly, path) for assembly, _prov in tracks)]
    if uncovered:
        problems.append('%d read paths are not covered by any assembly, e.g. %s'
                        % (len(uncovered), ', '.join(uncovered[:5])))
    if len({tuple(path) for path, _prov in tracks}) != len(tracks):
        problems.append('duplicate assembly paths found')
    short = [i for i, (path, _prov) in enumerate(tracks, 1) if len(path) < min_barcodes]
    if short:
        problems.append('%d assemblies have fewer than %d barcodes: %s'
                        % (len(short), min_barcodes, short[:5]))
    return problems


def check_path_gfa(assemblies, gfa_file, min_overlap=2, no_join=frozenset()):
    """Re-read path.gfa and check its segments, links and merged paths against the assemblies."""
    problems = []
    names = ['Assembly_%d' % (i + 1) for i in range(len(assemblies))]
    paths = {name: tuple(path) for name, (path, _prov) in zip(names, assemblies)}
    segments, links, prows = parse_gfa_file(gfa_file)

    if sorted(segments) != sorted(names):
        problems.append('path GFA holds %d assembly segments, expected %d'
                        % (len(segments), len(names)))
    groups = overlap_groups(names, [(link[0], link[2], 0, 0, 0, '')
                                    for link in links if len(link) >= 3
                                    and link[0] in paths and link[2] in paths])
    group_of = {name: 'Group_%d' % (i + 1)
                for i, group in enumerate(groups) for name in group}
    for name in names:
        seq, tags = segments.get(name, ('', []))
        want = ['LN:i:%d' % len(paths[name]), 'NB:i:%d' % len(paths[name]),
                'BC:Z:%s' % '|'.join(paths[name]), 'GR:Z:%s' % group_of[name]]
        if seq != '*' or tags != want:
            problems.append('path GFA segment %s should carry * and %s, got %s and %s'
                            % (name, want, seq, tags))

    for link in links:
        if len(link) < 5:
            problems.append('path GFA link is too short: %s' % (link,))
            continue
        left, left_orient, right, right_orient, cigar = link[:5]
        tags = link[5:]
        if left not in paths or right not in paths:
            problems.append('path GFA link uses unknown assembly %s / %s' % (left, right))
            continue
        if left_orient != '+' or right_orient != '+':
            problems.append('path GFA link %s -> %s is not +/+' % (left, right))
        try:
            length = int(str(cigar).rstrip('M'))
            start_left = int(tag_value(tags, 'PS', 'i') or 0)
            start_right = int(tag_value(tags, 'QS', 'i') or 0)
        except ValueError:
            problems.append('path GFA link %s -> %s has an unreadable overlap %r'
                            % (left, right, cigar))
            continue
        if tag_value(tags, 'OV', 'i') != str(length):
            problems.append('path GFA link %s -> %s OV:i differs from %s'
                            % (left, right, cigar))
        if length < 0 or start_left < 0 or start_right < 0 \
                or start_left + length > len(paths[left]) \
                or start_right + length > len(paths[right]):
            problems.append('path GFA link %s -> %s points outside its assemblies'
                            % (left, right))
            continue
        run = tuple(paths[left][start_left:start_left + length])
        if run != tuple(paths[right][start_right:start_right + length]):
            problems.append('path GFA link %s -> %s joins different barcodes'
                            % (left, right))
        if tag_value(tags, 'OV', 'Z') != '|'.join(run):
            problems.append('path GFA link %s -> %s OV:Z is not the shared run'
                            % (left, right))
        kind = tag_value(tags, 'TY', 'Z')
        want_kind = ('abut' if length == 0 else
                     overlap_kind(paths[left], paths[right], length,
                                  start_left, start_right))
        if kind != want_kind:
            problems.append('path GFA link %s -> %s is %r, expected %r'
                            % (left, right, kind, want_kind))
        if length and length < min_overlap:
            problems.append('path GFA link %s -> %s is shorter than %d barcodes'
                            % (left, right, min_overlap))
        if no_join and run and all(bc in no_join for bc in run):
            problems.append('path GFA link %s -> %s overlaps only in %s'
                            % (left, right, ' / '.join(sorted(no_join))))

    seen = Counter()
    for name in sorted(prows, key=natural_key):
        segs, overlaps, tags = prows[name]
        if not segs:
            problems.append('path GFA path %s lists no assembly' % name)
            continue
        unknown = [seg for seg in segs if seg not in paths]
        if unknown:
            problems.append('path GFA path %s lists unknown assembly %s'
                            % (name, ', '.join(unknown[:3])))
            continue
        rebuilt = list(paths[segs[0]])
        for i, seg in enumerate(segs[1:], 1):
            cigar = overlaps[i - 1] if i - 1 < len(overlaps) else ''
            try:
                length = int(str(cigar).rstrip('M'))
            except ValueError:
                problems.append('path GFA path %s has an unreadable overlap %r'
                                % (name, cigar))
                continue
            piece = paths[seg]
            if length:
                if length > len(piece) or list(rebuilt[-length:]) != list(piece[:length]):
                    problems.append('path GFA path %s joins %s without a shared run'
                                    % (name, seg))
                if not any(((link[0] == segs[i - 1] and link[2] == seg)
                            or (link[0] == seg and link[2] == segs[i - 1]))
                           and str(link[4]) == str(cigar) for link in links):
                    problems.append('path GFA path %s has no link %s -> %s'
                                    % (name, segs[i - 1], seg))
                rebuilt += list(piece[length:])
        barcodes = tuple((tag_value(tags, 'BC', 'Z') or '').split('|'))
        if tuple(rebuilt) != barcodes:
            problems.append('path GFA path %s BC tag does not match its assemblies' % name)
        if tag_value(tags, 'NB', 'i') != str(len(barcodes)):
            problems.append('path GFA path %s NB tag is not %d' % (name, len(barcodes)))
        contained = [item for item in (tag_value(tags, 'CR', 'Z') or '').split(',')
                     if item]
        if tag_value(tags, 'NA', 'i') != str(len(segs) + len(contained)):
            problems.append('path GFA path %s NA tag is not %d'
                            % (name, len(segs) + len(contained)))
        for item in contained:
            if item not in paths:
                problems.append('path GFA path %s lists unknown contained assembly %s'
                                % (name, item))
            elif '|' + '|'.join(paths[item]) + '|' not in '|' + '|'.join(barcodes) + '|':
                problems.append('path GFA path %s says %s is contained, but it is not'
                                % (name, item))
        seen.update(segs)
        seen.update(contained)

    repeated = [name for name, count in seen.items() if count > 1]
    if repeated:
        problems.append('path GFA lists %d assemblies in more than one merged path, e.g. %s'
                        % (len(repeated), ', '.join(sorted(repeated, key=natural_key)[:5])))
    missing = [name for name in names if name not in seen]
    if missing:
        problems.append('path GFA does not place %d assemblies, e.g. %s'
                        % (len(missing), ', '.join(missing[:5])))
    return problems


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Assemble reads of reads_gt2bc.csv by barcode path '
                    '(path.txt + assembly.gfa + path.gfa)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', default=os.path.join(HERE, 'reads_gt2bc.csv'),
                        help='input CSV (one row per barcode: read_id / barcode_index / barcode / is_end)')
    parser.add_argument('--sep', default=',', help='input delimiter')
    parser.add_argument('--id_column', help='read id column (auto: read_id)')
    parser.add_argument('--bc_column', help='barcode column (auto: barcode)')
    parser.add_argument('--ends_column', help='is_end column (auto: is_end)')
    parser.add_argument('--index_column', help='barcode index column (auto: barcode_index)')
    parser.add_argument('--correct-list', '--correct_list', dest='correct_list',
                        default=DEFAULT_CORRECT_LIST,
                        help='barcode correction table (read_id / barcode_index / origin / '
                             'corrected); a missing file disables corrections')
    parser.add_argument('--remove', default=DEFAULT_REMOVE,
                        help='read ids left out of the assembly (one per line, blank lines '
                             "and '#' comments skipped, 'none' disables); a missing file "
                             'is only a warning')
    parser.add_argument('--min_barcodes', type=int, default=2,
                        help='minimum number of barcodes per read')
    parser.add_argument('--min_overlap', type=int, default=2,
                        help='minimum shared barcodes to join two paths')
    parser.add_argument('--no_join_barcodes', default=DEFAULT_NO_JOIN,
                        help="overlaps made only of these barcodes are not joined "
                             "(comma separated, 'none' disables)")
    parser.add_argument('--interval_column', help='interval column in bp between two '
                        'neighbouring barcodes (auto: interval, then the differences '
                        'of barcode_start)')
    parser.add_argument('--interval_tol', type=float, default=DEFAULT_INTERVAL_TOL,
                        help='relative difference of two interval lengths of one barcode '
                             'pair that still counts as the same junction')
    parser.add_argument('--interval_gap', type=int, default=DEFAULT_INTERVAL_GAP,
                        help='absolute difference in bp that always starts a new '
                             'interval cluster')
    parser.add_argument('--interval_scope', choices=('min', 'all', 'off'), default='min',
                        help='overlaps whose interval lengths must agree: min = only '
                             'overlaps of exactly --min_overlap barcodes, all = every '
                             'overlap, off = never')
    parser.add_argument('--paths_out', default=os.path.join(HERE, 'path.txt'),
                        help='output path.txt')
    parser.add_argument('--gfa_out', default=os.path.join(HERE, 'assembly.gfa'),
                        help='output GFA file')
    parser.add_argument('--path_gfa_out',
                        help='merged assembly-level GFA (assemblies as segments, their '
                             'overlaps as links); default: path.gfa next to --gfa_out')
    parser.add_argument('--path_min_overlap', type=int,
                        help='minimum shared barcodes to link or merge two assemblies '
                             '(default: the value of --min_overlap)')
    parser.add_argument('--map_out', help='optional read -> assembly interval table (TSV)')
    parser.add_argument('--no_marks', action='store_true', help='do not write head/tail marks')
    parser.add_argument('--gfa_star', action='store_true',
                        help='deprecated, S lines always carry * (the barcode path '
                             'is written in the BC tag)')
    parser.add_argument('--verbose', action='store_true',
                        help='print every merge/drop and the read chain of each assembly')
    args = parser.parse_args(argv)

    corr_stats = {'entries': 0, 'reads': 0, 'skipped': 0, 'warnings': [],
                  'applied': 0, 'missing_index': 0, 'mismatch': 0, 'same': 0,
                  'reads_changed': 0, 'missing_read': 0, 'missing_reads': []}
    corrections = {}
    has_correct_list = os.path.isfile(args.correct_list)
    if has_correct_list:
        corrections, corr_stats = read_correct_list(args.correct_list)
        corr_stats.update({'applied': 0, 'missing_index': 0, 'mismatch': 0, 'same': 0,
                           'reads_changed': 0, 'missing_read': 0, 'missing_reads': []})

    read_intervals = {}
    rows, id_col, bc_col, ends_col, idx_col, dropped = read_reads(
        args.input, args.sep, args.id_column, args.bc_column, args.ends_column,
        args.index_column, args.min_barcodes, corrections, corr_stats,
        args.interval_column, read_intervals)
    if not rows:
        raise SystemExit('[ERROR] no reads with >= %d barcodes in %s'
                         % (args.min_barcodes, args.input))

    print('[INFO] input      : %s' % args.input)
    print('[INFO] columns    : read_id=%s, barcode=%s, is_end=%s, index=%s%s'
          % (id_col, bc_col, ends_col or '(none)', idx_col or '(none)',
             '' if not dropped else '; skipped %d reads with too few barcodes' % dropped))
    if not has_correct_list:
        print('[INFO] correct    : %s not found, no barcode corrections'
              % args.correct_list)
    elif not corrections:
        print('[INFO] correct    : %s holds no usable barcode entry, no corrections'
              % args.correct_list)
    else:
        unapplied = (corr_stats['skipped'] + corr_stats['missing_index']
                     + corr_stats['mismatch'] + corr_stats['missing_read'])
        print('[INFO] correct    : %s (%d rows / %d reads): %d applied to %d reads, '
              '%d skipped (index missing/ambiguous %d, barcode differs %d, '
              'read not in input %d, unusable row %d)%s'
              % (args.correct_list, corr_stats['entries'] + corr_stats['skipped'],
                 corr_stats['reads'], corr_stats['applied'], corr_stats['reads_changed'],
                 unapplied, corr_stats['missing_index'], corr_stats['mismatch'],
                 corr_stats['missing_read'], corr_stats['skipped'],
                 ', %d already corrected' % corr_stats['same'] if corr_stats['same'] else ''))
        for warning in corr_stats['warnings']:
            print('[WARN] %s' % warning)
        if idx_col is None:
            print('[WARN] input has no barcode_index column, corrections cannot be applied')
        if corr_stats['missing_reads']:
            names = corr_stats['missing_reads']
            print('[WARN] %d correct-list read(s) are not in the input: %s%s'
                  % (len(names), ', '.join(names[:10]), ' ...' if len(names) > 10 else ''))

    remove_ids = set()
    removed_count = 0
    if str(args.remove).strip().lower() in ('none', 'no', 'off', '-'):
        print('[INFO] remove     : disabled (--remove %s)' % args.remove)
    elif not os.path.isfile(args.remove):
        print('[WARN] remove     : %s not found, no reads removed' % args.remove)
    else:
        remove_ids, remove_stats = read_remove_list(args.remove)
        if not remove_ids:
            print('[INFO] remove     : %s holds no read id, no reads removed' % args.remove)
        else:
            input_ids = {str(row[0]) for row in rows}
            missing = sorted(remove_ids - input_ids, key=natural_key)
            rows, removed_reads = drop_reads(rows, remove_ids)
            removed_count = len(removed_reads)
            print('[INFO] remove     : %s (%d ids, %d listed twice): %d reads removed, '
                  '%d reads kept%s'
                  % (args.remove, len(remove_ids), remove_stats['duplicates'],
                     removed_count, len(rows),
                     ', %d ids not in the input' % len(missing) if missing else ''))
            if missing:
                print('[WARN] %d remove-list read(s) are not in the input: %s%s'
                      % (len(missing), ', '.join(missing[:10]),
                         ' ...' if len(missing) > 10 else ''))
            if not rows:
                raise SystemExit('[ERROR] every read of %s was removed by %s'
                                 % (args.input, args.remove))

    read_paths = {row[0]: row[1] for row in rows}
    interval_table = None
    if args.interval_scope == 'off':
        print('[INFO] interval   : joins are decided by the barcode order only '
              '(--interval_scope off)')
    elif not read_intervals:
        print('[WARN] interval   : no interval / barcode_start column in %s, '
              'ambiguous joins stay decided by the barcode order alone' % args.input)
    else:
        interval_table = IntervalClusters(read_paths, read_intervals,
                                          args.interval_tol, args.interval_gap)
        split = [(pair, groups) for pair, groups in interval_table.pairs.items()
                 if len(groups) > 1]
        print('[INFO] interval   : %d of %d reads carry interval lengths, %d barcode '
              'pairs, %d pairs split into several clusters'
              % (len(interval_table.reads), len(rows), len(interval_table.pairs),
                 len(split)))
        print('[INFO] cluster    : same junction within %.0f%% / %d bp, checked at %s'
              % (args.interval_tol * 100, args.interval_gap,
                 'overlaps of %d barcodes' % args.min_overlap
                 if args.interval_scope == 'min' else 'every overlap'))
        if args.verbose:
            for pair, groups in sorted(split, key=lambda item: (-len(item[1]),
                                                                item[0]))[:10]:
                print('[CLUSTER] %s|%s: %s'
                      % (pair[0], pair[1],
                         ', '.join('%d-%d bp (%d reads)' % group for group in groups)))

    marks = {} if args.no_marks else load_marks(rows)
    raw_paths = [barcodes for _read_id, barcodes, _ends in rows]
    unique_paths = len({tuple(path) for path in raw_paths})
    print('[INFO] reads      : %d (different barcode paths %d, barcode types %d)'
          % (len(raw_paths), unique_paths, len({bc for path in raw_paths for bc in path})))
    if args.gfa_star:
        print('[INFO] gfa_star   : ignored, S lines always carry * (barcode path in the BC tag)')

    assemblies = [(tuple(path), ((read_id, 0),)) for read_id, path, _ends in rows]
    no_join = parse_no_join(args.no_join_barcodes)
    if no_join:
        print('[INFO] joining    : min overlap %d barcodes; overlaps made only of %s kept apart'
              % (args.min_overlap, ' / '.join(sorted(no_join))))
    else:
        print('[INFO] joining    : min overlap %d barcodes (no_join disabled)'
              % args.min_overlap)
    stats = {'blocked': 0, 'blocked_examples': [], 'merges': [],
             'interval_keys': set(), 'interval_examples': [], 'interval_confirmed': 0}
    rounds = 0
    while True:
        blocked_before = stats['blocked']
        interval_before = len(stats['interval_keys'])
        merged = build_overlap_assemblies(assemblies, args.min_overlap, no_join,
                                          args.verbose, stats, interval_table,
                                          args.interval_scope)
        rounds += 1
        blocked_now = stats['blocked'] - blocked_before
        interval_now = len(stats['interval_keys']) - interval_before
        extra = ''
        if blocked_now:
            extra = ' (%d joins blocked: overlap only %s)' % (
                blocked_now, ' / '.join(sorted(no_join)))
        if interval_now:
            extra += ' (%d joins blocked: interval lengths differ)' % interval_now
        print('[INFO] round %-2d   : %d -> %d paths%s'
              % (rounds, len(assemblies), len(merged), extra))
        done = len(merged) >= len(assemblies)
        assemblies = merged
        if done:
            break
    if no_join and stats['blocked_examples'] and args.verbose:
        print('[INFO] blocked overlaps:')
        for left, right, ov in stats['blocked_examples']:
            print('       ... %s  +  %s ... (overlap %d)'
                  % (' '.join(left[-3:]), ' '.join(right[:3]), ov))

    bad_merges = [item for item in stats['merges']
                  if no_join and all(bc in no_join for bc in item[1])]

    before = len(assemblies)
    assemblies = remove_contained_paths(assemblies, args.verbose)
    print('[INFO] contained  : %d -> %d assemblies' % (before, len(assemblies)))
    assemblies.sort(key=lambda item: (-len(item[0]), item[0]))
    paths = [path for path, _prov in assemblies]

    flags, internal = annotate_assemblies(paths, rows, marks)
    counts = Counter(MARK_SUFFIX[flags[name]] or 'none' for name in flags)
    write_paths(paths, args.paths_out, flags)
    node_count, edge_count, gfa_rows = write_gfa(args.paths_out, assemblies, rows,
                                                 args.gfa_out, args.verbose, args.min_overlap)
    path_gfa_out = args.path_gfa_out or os.path.join(
        os.path.dirname(os.path.abspath(args.gfa_out)), 'path.gfa')
    path_min_overlap = args.path_min_overlap or args.min_overlap
    path_nodes, path_edges, path_rows, path_groups = write_path_gfa(
        assemblies, flags, path_gfa_out, path_min_overlap, no_join, args.verbose)
    mapping = write_mapping(rows, paths, args.map_out) if args.map_out else {}

    problems = check_outputs(rows, assemblies, flags, args.paths_out, args.gfa_out,
                             args.min_barcodes, args.min_overlap)
    problems.extend(check_path_gfa(assemblies, path_gfa_out, path_min_overlap, no_join))
    if bad_merges:
        problems.append('%d joins used an overlap made only of %s, e.g. overlap=%d %s'
                        % (len(bad_merges), ' / '.join(sorted(no_join)), bad_merges[0][0],
                           '|'.join(bad_merges[0][1])))
    covered = len(rows) - len([1 for _rid, path, _ends in rows
                               if not any(path_positions(assembly, path)
                                          for assembly in paths)])
    joined_reads = {entry[0] for row in gfa_rows for entry in row[2]}
    print('[INFO] ===== summary =====')
    print('[INFO] reads      : %d (>= %d barcodes)' % (len(rows), args.min_barcodes))
    if removed_count:
        print('[INFO] remove     : %d read(s) left out by %s' % (removed_count, args.remove))
    print('[INFO] assembly   : %d rounds -> %d assemblies, longest %d barcodes, covers %d/%d reads'
          % (rounds, len(paths), max(len(p) for p in paths), covered, len(rows)))
    print('[INFO] join limit : %s (%d blocked joins)'
          % ('overlaps made only of %s are kept apart' % ' / '.join(sorted(no_join))
             if no_join else 'disabled', stats['blocked']))
    if interval_table is not None:
        print('[INFO] interval   : %d candidate joins rejected by differing interval '
              'lengths, %d joins confirmed by matching ones (details with --verbose)'
              % (len(stats['interval_keys']), stats['interval_confirmed']))
        for left, right, ov, detail in stats['interval_examples'][:3]:
            pair, left_length, right_length = detail
            print('       ... %s + %s at %s (overlap %d): %s vs %s bp'
                  % (' '.join(left[-3:]), ' '.join(right[:3]), '|'.join(pair), ov,
                     left_length, right_length))
    print('[INFO] merges     : %d (overlap distribution %s)'
          % (len(stats['merges']),
             ', '.join('%d:%d' % item for item in
                       sorted(Counter(ov for ov, _region in stats['merges']).items()))))
    print('[INFO] graph      : %d read nodes / %d join edges / %d assembly paths'
          % (node_count, edge_count, len(paths)))
    print('[INFO] gfa reads  : %d reads joined into chains, %d reads contained only '
          '(listed in the CR tag)'
          % (len(joined_reads), len(rows) - len(joined_reads)))
    chain_overlaps = Counter(ov for row in gfa_rows for ov in row[3])
    print('[INFO] gfa joins  : %d read->read joins, overlap distribution %s'
          % (sum(chain_overlaps.values()),
             ', '.join('%d:%d' % item for item in sorted(chain_overlaps.items())) or 'none'))
    kinds = Counter(kind for _left, _right, _len, _sl, _sr, kind in path_edges)
    joined_paths = [row for row in path_rows if len(row[2]) + len(row[4]) > 1]
    merged_groups = [group for group in path_groups if len(group) > 1]
    print('[INFO] path gfa   : %d assembly nodes / %d overlap links / %d merged paths'
          % (path_nodes, len(path_edges), len(path_rows)))
    print('[INFO] path links : dovetail %d / contained %d / internal %d, '
          'overlap distribution %s'
          % (kinds['dovetail'], kinds['contained'], kinds['internal'],
             ', '.join('%d:%d' % item for item in
                       sorted(Counter(edge[2] for edge in path_edges).items())) or 'none'))
    print('[INFO] path merge : %d paths join several assemblies end to end (dovetail overlaps '
          '>= %d are already used by assembly.gfa); the GR tag groups all %d assemblies '
          'into %d merged groups (%d with several)'
          % (len(joined_paths), path_min_overlap, path_nodes, len(path_groups),
             len(merged_groups)))
    if merged_groups:
        print('[INFO] path group : largest %d assemblies: %s%s'
              % (len(path_groups[0]), ' '.join(path_groups[0][:6]),
                 ' ...' if len(path_groups[0]) > 6 else ''))
    for left, right, length, start_left, start_right, kind in sorted(
            path_edges, key=lambda item: (-item[2], natural_key(item[0]),
                                          natural_key(item[1])))[:3]:
        print('[INFO] best overlap: %s ~ %s, %d barcodes at %d/%d (%s) %s'
              % (left, right, length, start_left, start_right, kind,
                 '|'.join(paths[int(left.split('_')[1]) - 1]
                          [start_left:start_left + length])))
    if args.no_marks:
        print('[INFO] head/tail  : not annotated (--no_marks)')
    else:
        print('[INFO] head/tail  : head- %d / -tail %d / head-tail %d / none %d '
              '(internal head hits %d, tail hits %d, not annotated)'
              % (counts['head-'], counts['-tail'], counts['head-tail'], counts['none'],
                 internal['head'], internal['tail']))
    print('[INFO] output     : %s (%d), %s, %s%s'
          % (args.paths_out, len(paths), args.gfa_out, path_gfa_out,
             ', %s (%d reads)' % (args.map_out, len(mapping)) if args.map_out else ''))
    if problems:
        print('[ERROR] self-check found %d problem(s):' % len(problems))
        for item in problems:
            print('        - %s' % item)
        return 1
    print('[INFO] self-check : passed (read chains, overlaps, marks, assembly overlaps '
          'and path.txt agree)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

