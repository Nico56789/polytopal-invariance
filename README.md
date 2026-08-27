# polytopal-invariance

Companion code and data for

> N. Libedinsky, *Polytopal invariance for Lusztig's q-weight multiplicities*.

To the data of an irreducible root system Φ and dominant weights μ ≤ λ the paper
attaches the convex polytope

    P^Φ_{μ,λ} = (μ + C) ∩ (λ − C) ∩ (D̄ − ρ),

where C is the cone spanned by the positive roots and D̄ the closed dominant Weyl
chamber, and proves that the isometry class of P^Φ_{μ,λ} determines the
generalized Kostka–Foulkes polynomial K<sup>Φ</sup><sub>λμ</sub>(q) — Lusztig's
q-analogue of the μ-weight multiplicity in the irreducible highest-weight
representation V(λ). This repository contains the script used for the census
quoted in the introduction of the paper, together with the derived data that
census produced.

**Notation.** λ is the highest weight and μ the lower one, exactly as in the
paper. In the code and in the catalog the two are called `high` (= λ) and
`low` (= μ); weights are always given by their fundamental-weight coordinates.

## What the script computes

* **Polytopes.** Vertices of P^Φ_{μ,λ} by exact enumeration of the defining
  inequalities, in exact rational arithmetic (Python `Fraction`), in the
  normalization where long roots have squared length 2, so that distances are
  comparable across Dynkin types.
* **Isometry.** Two polytopes are put in the same class iff their vertex sets
  admit a distance-preserving bijection — exact congruence of the matrices of
  squared distances, found by backtracking, again over the rationals. Such a
  bijection extends to an affine isometry of the polytopes.
* **Polynomials.** K<sup>Φ</sup><sub>λμ</sub>(q) by Lusztig's definition: the
  alternating Weyl sum of the q-analogue of Kostant's partition function.

The geometry is exact rational arithmetic throughout. The partition function and
the Weyl sum are exact *integer* computations, carried in numpy `int64` for
speed; since silent wraparound would be indistinguishable from a wrong answer,
both paths carry runtime guards (`_check_int64`, and a per-step sign-and-
magnitude check inside the partition-function DP) that raise `OverflowError`
rather than return a wrapped result. The only floating-point object in the file,
a Cholesky factor, is never used by the census path.

## The census of the paper

Types A1–A5, B2–B4, C3, C4, D4, D5, E6, F4, G2, exhaustive over the pairs
μ ≤ λ with fundamental-weight coordinates bounded by

| bound | types |
|-------|-----------------------------|
| 12    | A1, A2, B2, G2              |
| 8     | A3                          |
| 6     | B3, C3                      |
| 3     | A4, B4, C4, D4, F4          |
| 1     | A5, D5, E6                  |

Result: **488 981** nonempty polytopes in **215 458** isometry classes; **64 061**
classes with more than one member (35 995 576 isometric pairs); K computed
independently for every member of every such class, with **zero mismatches**;
**2 262** classes realized in more than one Dynkin type, **1 328** bridging two
different Dynkin families.

`isocheck` is a second, independent look at the class partition, by a different
algorithm (Weisfeiler–Leman colour refinement + guided matching). Its
false-negative guard — that no two distinct classes share a distance signature,
so the signature is a complete invariant on this catalog — is exhaustive. Its
false-positive guard recomputes members of each multi-member class and confirms
them against the class representative; by default it does this for two members
per class, and the number is a command-line argument.

## Contents

```
polytopal_census.py   the census program (Python 3 + numpy, self-contained)
export_data.py        rebuilds data/ from a catalog produced by `census`
data/
  summary.txt                  the counts above, as printed from the catalog
  classes_multitype.csv        the 2262 classes realized in >1 Dynkin type
  classes_crossfamily.csv      the 1328 of those bridging >1 Dynkin family
  realizations_multitype.csv   every (type, μ, λ) realizing one of those classes
  class_size_histogram.csv     number of classes with n members, for each n
```

Each row of `classes_multitype.csv` is one isometry class: its number of
vertices, the Dynkin types realizing it, how many pairs (μ, λ) realize it, the
common polynomial K, and one example pair per type, written `μ -> λ`. Every
class in that file is a coincidence of the kind the main theorem explains: the
same polytope, hence the same K, occurring in unrelated root systems.

The full catalog (a ~2 GB SQLite database) is not distributed here; it is
rebuilt from scratch by `polytopal_census.py census`, and `export_data.py`
regenerates `data/` from it.

## Usage

Requires Python 3 and numpy.

```
# one polynomial, arguments in the order  TYPE  MU  LAMBDA
# (the D5/E6 example of the paper: μ = ϖ₃+ϖ₅, λ = ϖ₁+ϖ₂+ϖ₃+ϖ₄):
python3 polytopal_census.py pair D5 0,0,1,0,1 1,1,1,1,0
python3 polytopal_census.py pair E6 0,0,0,1,1,1 1,1,1,1,0,0

# small census (types A1,A2,A3,B2,G2, coordinates <= 4):
python3 polytopal_census.py demo

# the full census of the paper (several CPU-days; SQLite catalog on disk,
# resumable, commits incrementally):
python3 polytopal_census.py census catalog.db

# on an existing catalog:
python3 polytopal_census.py verifyk  catalog.db     # recompute and recheck all K per class
python3 polytopal_census.py isocheck catalog.db 2   # independent isometry cross-validation
python3 polytopal_census.py numbers  catalog.db     # the summary numbers quoted in the paper
python3 export_data.py               catalog.db     # regenerate data/
```

The first two commands print

```
K^D5_{lambda=(1, 1, 1, 1, 0), mu=(0, 0, 1, 0, 1)}(q) =
    1q^2 + 5q^3 + 13q^4 + 19q^5 + 20q^6 + 18q^7 + 12q^8 + 6q^9 + 3q^10 + 1q^11
```

for both — one polytope, two root systems of different families, one polynomial.
That class is `214481` in `data/classes_crossfamily.csv`.

## Licence

MIT, see `LICENSE`. If you use this code or data, please cite the paper.
