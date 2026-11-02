#!/usr/bin/env python3
import argparse
import sys
import time

import numpy as np

CODE = np.full(256, 4, np.uint8)
CODE[np.frombuffer(b"ACGT", np.uint8)] = np.arange(4, dtype=np.uint8)
COMP = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")


def read_fasta(path):
    name, chunks = None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line[:1] == b">":
                if name is not None:
                    yield name, b"".join(chunks)
                fields = line[1:].split()
                name = fields[0].decode() if fields else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, b"".join(chunks)


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def canon_hashes(seq, k):
    c = CODE[np.frombuffer(seq, np.uint8)].astype(np.uint64)
    n = len(c) - k + 1
    if n < 1:
        return np.empty(0, np.uint64)
    c = np.minimum(c, np.uint64(3))
    h = np.zeros(n, np.uint64)
    r = np.zeros(n, np.uint64)
    two = np.uint64(2)
    three = np.uint64(3)
    for j in range(k):
        h = (h << two) | c[j:j + n]
        r = (r << two) | (three - c[k - 1 - j:k - 1 - j + n])
    return np.minimum(h, r)


def window_min(h, w):
    n = len(h)
    if n < w:
        return np.empty(0, np.uint64)
    pad = (-n) % w
    big = np.uint64(0xFFFFFFFFFFFFFFFF)
    a = np.concatenate([h, np.full(pad, big, np.uint64)])
    blk = a.reshape(len(a) // w, w)
    pre = np.minimum.accumulate(blk, axis=1).ravel()
    suf = np.minimum.accumulate(np.ascontiguousarray(blk[:, ::-1]),
                                axis=1)[:, ::-1].ravel()
    return np.minimum(suf[:n - w + 1], pre[w - 1:n])


def probes(seq, k, w):
    h = canon_hashes(seq, k)
    if not len(h):
        return h
    return np.unique(h if len(h) < w else window_min(h, w))


def build_index(rows):
    lens = np.array([len(r) for r in rows], dtype=np.int64)
    hashes = np.concatenate(rows)
    ids = np.repeat(np.arange(len(rows), dtype=np.int32), lens)
    order = np.argsort(hashes, kind="stable")
    return hashes[order], ids[order]


def candidates(pr, ih, iid, n, max_df, min_probes, max_gather):
    if not len(pr):
        return np.empty(0, dtype=np.int64)
    lo = np.searchsorted(ih, pr, "left")
    hi = np.searchsorted(ih, pr, "right")
    df = hi - lo
    ok = np.flatnonzero((df > 0) & (df <= max_df))
    if len(ok) < min_probes:
        top = np.argsort(df, kind="stable")[:min_probes]
        ok = top[df[top] > 0]
    if not len(ok):
        return np.empty(0, dtype=np.int64)
    sel, total = [], 0
    for t in ok[np.argsort(df[ok], kind="stable")]:
        total += int(df[t])
        if sel and total > max_gather:
            break
        sel.append(int(t))
    hits = np.concatenate([iid[lo[t]:hi[t]] for t in sel])
    return np.flatnonzero(np.bincount(hits, minlength=n) == len(sel))


def main():
    ap = argparse.ArgumentParser(
        description="remove reads whose whole sequence occurs inside another read")
    ap.add_argument("--fasta", default="reads_all.fasta")
    ap.add_argument("--out", default="reads_all_no_contained.fasta")
    ap.add_argument("--tsv", default="")
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--w", type=int, default=100)
    ap.add_argument("--max-df", type=int, default=200)
    ap.add_argument("--min-probes", type=int, default=20)
    ap.add_argument("--max-gather", type=int, default=400000)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    args.tsv = args.tsv or args.out + ".dropped.tsv"
    t0 = time.time()

    ids, seqs = [], []
    for name, seq in read_fasta(args.fasta):
        if args.limit and len(ids) >= args.limit:
            break
        ids.append(name)
        seqs.append(seq)
    n = len(ids)
    lens = np.array([len(s) for s in seqs], dtype=np.int64)
    print("reads %d, bp %d" % (n, int(lens.sum())), file=sys.stderr, flush=True)
    if not n:
        return 1

    prs = [probes(s, args.k, args.w) for s in seqs]
    ih, iid = build_index(prs)
    print("probes %d, %.1fs" % (len(ih), time.time() - t0),
          file=sys.stderr, flush=True)

    keep = np.ones(n, dtype=bool)
    drop = {}

    def put(key, value):
        if key not in drop or (drop[key][2] == "-" and value[2] == "+"):
            drop[key] = value

    for i in range(n):
        if i % 2000 == 0:
            print("checked %d/%d, dropped %d, %.1fs"
                  % (i, n, int((~keep).sum()), time.time() - t0),
                  file=sys.stderr, flush=True)
        if not keep[i]:
            continue
        seq, li = seqs[i], lens[i]
        if li < args.k:
            cand = np.arange(n, dtype=np.int64)
        else:
            cand = candidates(prs[i], ih, iid, n, args.max_df, args.min_probes,
                              args.max_gather)
        cand = cand[(cand != i) & (lens[cand] >= li)]
        if not len(cand):
            continue
        rseq = revcomp(seq)
        for j in cand:
            j = int(j)
            other = seqs[j]
            off = other.find(seq)
            strand = "+"
            if off < 0:
                off = other.find(rseq)
                strand = "-"
                if off < 0:
                    continue
                off = lens[j] - li - off
            if lens[j] == li:
                first, second = (i, j) if i < j else (j, i)
                keep[second] = False
                put(second, (first, "identical", strand, 0))
            else:
                keep[i] = False
                put(i, (j, "contained", strand, off))

    with open(args.fasta, "rb") as fi, open(args.out, "wb") as fo:
        idx, buf = -1, []
        for line in fi:
            if line[:1] == b">":
                if 0 <= idx < n and keep[idx]:
                    fo.write(b"".join(buf))
                idx += 1
                if idx >= n:
                    break
                buf = [line]
            elif idx >= 0:
                buf.append(line)
        else:
            if 0 <= idx < n and keep[idx]:
                fo.write(b"".join(buf))

    with open(args.tsv, "w", encoding="utf-8", newline="") as fh:
        fh.write("read_id\tread_len\tcontainer_id\tcontainer_len\t"
                 "relation\tstrand\toffset_1based\n")
        for k in sorted(drop, key=lambda x: (-lens[x], ids[x])):
            j, rel, strand, off = drop[k]
            fh.write("%s\t%d\t%s\t%d\t%s\t%s\t%d\n"
                     % (ids[k], lens[k], ids[j], lens[j], rel, strand, off + 1))

    print("input %d reads / %d bp" % (n, int(lens.sum())))
    print("dropped %d reads / %d bp"
          % (len(drop), int(lens[sorted(drop)].sum()) if drop else 0))
    print("kept %d reads / %d bp" % (int(keep.sum()), int(lens[keep].sum())))
    print("out %s" % args.out)
    print("tsv %s" % args.tsv)
    print("%.1fs" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())

