#!/usr/bin/env python3
"""
export_data.py -- derive the companion datasets in data/ from a catalog built by

    python3 polytopal_census.py census catalog.db

The catalog itself is a ~2 GB SQLite file and is not distributed here.  This
script extracts from it exactly the material quoted in the paper: the summary
counts, and the full list of isometry classes that are realized in more than one
Dynkin type -- the coincidences that motivate the main theorem -- so that a
reader can inspect every one of them without re-running the census.

Notation follows the paper: mu <= lambda, lambda is the highest weight, the
polytope is P^Phi_{mu,lambda} and the polynomial K^Phi_{lambda,mu}(q).  The
catalog columns `low`/`high` are mu/lambda respectively.

Usage:
    python3 export_data.py [catalog.db] [outdir]
"""
import sys, os, csv, json, sqlite3
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from polytopal_census import poly_str, PAPER_BOUNDS

TKEY = lambda t: (t[0], int(t[1:]))


def wt(s):
    return ' '.join(str(x) for x in json.loads(s))


def main(db_path='catalog.db', outdir='data'):
    os.makedirs(outdir, exist_ok=True)
    con = sqlite3.connect(db_path)

    # ---- one pass over the realizations ------------------------------------
    types_of = defaultdict(set)      # class -> set of Dynkin types
    count_of = Counter()             # class -> number of realizations
    first_of = {}                    # (class, type) -> first (mu, lambda)
    rows = con.execute("SELECT class_id,type,low,high FROM realizations "
                       "ORDER BY class_id,type,low,high").fetchall()
    for cid, t, lo, hi in rows:
        types_of[cid].add(t)
        count_of[cid] += 1
        first_of.setdefault((cid, t), (lo, hi))

    sizes = list(count_of.values())
    multitype = sorted(cid for cid in types_of if len(types_of[cid]) > 1)
    crossfam = [cid for cid in multitype if len({t[0] for t in types_of[cid]}) > 1]
    mt_set, cf_set = set(multitype), set(crossfam)

    # ---- summary.txt -------------------------------------------------------
    L = ["Census of  N. Libedinsky,",
         '"Polytopal invariance for Lusztig\'s q-weight multiplicities"',
         "",
         "P^Phi_{mu,lambda} = (mu + C) cap (lambda - C) cap (closure(D) - rho),",
         "K^Phi_{lambda,mu}(q) = Lusztig's q-analogue of the mu-weight",
         "multiplicity in V(lambda).",
         "",
         "Exhaustive over the pairs mu <= lambda whose fundamental-weight",
         "coordinates are bounded by:", ""]
    for t in sorted(PAPER_BOUNDS, key=TKEY):
        L.append(f"    {t:3}  <= {PAPER_BOUNDS[t]}")
    L += ["",
          f"nonempty polytopes (realizations): {sum(sizes)}",
          f"isometry classes:                  {len(sizes)}",
          f"classes with >1 member:            {sum(1 for n in sizes if n > 1)}",
          f"isometric pairs:                   {sum(n*(n-1)//2 for n in sizes)}",
          f"classes realized in >1 type:       {len(multitype)}",
          f"classes bridging >1 Dynkin family: {len(crossfam)}",
          "",
          "K-agreement: for every class with at least two members, K was computed",
          "independently from each member and the results compared; no mismatch",
          "was found (`python3 polytopal_census.py verifyk catalog.db`).",
          ""]
    with open(os.path.join(outdir, 'summary.txt'), 'w') as f:
        f.write('\n'.join(L))

    # ---- class_size_histogram.csv -----------------------------------------
    hist = Counter(sizes)
    with open(os.path.join(outdir, 'class_size_histogram.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['members_in_class', 'number_of_classes'])
        for s in sorted(hist):
            w.writerow([s, hist[s]])

    # ---- one row per class -------------------------------------------------
    info = {cid: (n, k) for cid, n, k in con.execute(
        "SELECT id,sig_n,k_json FROM classes")}

    def dump_classes(name, ids):
        path = os.path.join(outdir, name)
        with open(path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['class_id', 'vertices', 'n_types', 'types',
                        'n_realizations', 'K',
                        'one_example_mu_to_lambda_per_type'])
            for cid in ids:
                sig_n, k_json = info[cid]
                ts = sorted(types_of[cid], key=TKEY)
                K = poly_str({int(k): v for k, v in json.loads(k_json).items()}) \
                    if k_json else ''
                ex = ' ; '.join(
                    f"{t}: ({wt(first_of[(cid, t)][0])}) -> ({wt(first_of[(cid, t)][1])})"
                    for t in ts)
                w.writerow([cid, sig_n, len(ts), '+'.join(ts), count_of[cid], K, ex])
        print(f"  {path}  ({len(ids)} rows)")

    dump_classes('classes_multitype.csv', multitype)
    dump_classes('classes_crossfamily.csv', crossfam)

    # ---- one row per realization of a multi-type class ---------------------
    path = os.path.join(outdir, 'realizations_multitype.csv')
    n = 0
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['class_id', 'crosses_families', 'type', 'mu', 'lambda'])
        for cid, t, lo, hi in rows:
            if cid in mt_set:
                w.writerow([cid, int(cid in cf_set), t, wt(lo), wt(hi)])
                n += 1
    print(f"  {path}  ({n} rows)")
    print(f"  {os.path.join(outdir, 'summary.txt')}")
    print(f"  {os.path.join(outdir, 'class_size_histogram.csv')}")
    con.close()


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'catalog.db',
         sys.argv[2] if len(sys.argv) > 2 else 'data')
