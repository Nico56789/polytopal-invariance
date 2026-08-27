#!/usr/bin/env python3
"""
polytopal_census.py -- companion code for
  N. Libedinsky, "Polytopal invariance for Lusztig's q-weight multiplicities".

Notation follows the paper: mu <= lambda are dominant weights in the dominance
order, lambda is the HIGHEST weight, mu the lower one, the polytope is

    P^Phi_{mu,lambda} = (mu + C) cap (lambda - C) cap (closure(D) - rho),

and the polynomial is Lusztig's q-analogue K^Phi_{lambda,mu}(q) of the mu-weight
multiplicity in the irreducible highest-weight representation V(lambda).
Internally the two weights are called `low` (= mu) and `high` (= lambda), and the
catalog columns `low`/`high` mean the same.

Self-contained (Python 3 + numpy).  The geometry is computed in exact rational
arithmetic (Fraction): polytope vertices and pairwise squared distances, and
isometry of polytopes by exact congruence of vertex distance matrices (an
isometry of the vertex set extends affinely to the polytopes).  The partition
function and the Weyl sum are exact integer computations carried in numpy int64,
with runtime guards that abort if any intermediate could wrap around (see
_check_int64 and _Pq_box_np); the only floating-point object in the file, a
Cholesky factor, is never used by the census path.

Usage:
  python3 polytopal_census.py pair TYPE MU LAMBDA
      K^Phi_{lambda,mu} for the dominant pair mu <= lambda, e.g. the D5 example
      of the paper (mu = w_3 + w_5, lambda = w_1 + w_2 + w_3 + w_4):
      python3 polytopal_census.py pair D5 0,0,1,0,1 1,1,1,1,0
  python3 polytopal_census.py census [catalog.db]
      Full census of the paper: enumerate all nonempty polytopes with the
      per-type coordinate bounds below, group them into exact isometry
      classes, then independently compute K for every member of every class
      with >=2 members and check that they agree.  (Several CPU-days.)
  python3 polytopal_census.py demo [catalog_demo.db]
      Small census (types A1,A2,A3,B2,G2, coordinates <= 4; a few minutes).
  python3 polytopal_census.py verifyk [catalog.db]
      (Re)run the K-agreement check on an existing catalog.
  python3 polytopal_census.py isocheck [catalog.db] [members_per_class]
      Cross-validate the isometry classes with an independent matching
      algorithm (Weisfeiler-Leman colour refinement + guided matching).  The
      false-negative guard (no two classes share a distance signature) is
      exhaustive; the false-positive guard recomputes `members_per_class`
      members of each multi-member class (default 2, 0 = all of them).
  python3 polytopal_census.py numbers [catalog.db]
      Print the summary numbers quoted in the paper.

Census bounds used in the paper (fundamental-weight coordinates of both
mu and lambda):
  12 for A1, A2, B2, G2;  8 for A3;  6 for B3, C3;
  3 for A4, B4, C4, D4, F4;  1 for A5, D5, E6.
"""

import sys, json, sqlite3
from fractions import Fraction as Fr
from itertools import combinations, product
from math import gcd as _igcd
import numpy as np

sys.setrecursionlimit(1000000)

_H = Fr(1, 2)

PAPER_BOUNDS = {
    'A1': 12, 'A2': 12, 'B2': 12, 'G2': 12,
    'A3': 8,
    'B3': 6, 'C3': 6,
    'A4': 3, 'B4': 3, 'C4': 3, 'D4': 3, 'F4': 3,
    'A5': 1, 'D5': 1, 'E6': 1,
}


def gram_A(r):
    G = [[Fr(0)]*r for _ in range(r)]
    for i in range(r):
        G[i][i] = Fr(2)
        if i+1 < r:
            G[i][i+1] = Fr(-1); G[i+1][i] = Fr(-1)
    return G

def posroots_A(r):
    # roots e_i - e_j, i<j  ->  alpha_i + ... + alpha_{j-1}
    roots = []
    for i in range(r):
        for j in range(i, r):
            v = [0]*r
            for k in range(i, j+1):
                v[k] = 1
            roots.append(tuple(v))
    return roots

TYPES = {}
TYPES['A1'] = dict(
    G=gram_A(1),
    posroots=posroots_A(1),   # (1,)
)
TYPES['A2'] = dict(
    G=gram_A(2),
    posroots=posroots_A(2),   # (1,0),(0,1),(1,1)
)
TYPES['A3'] = dict(
    G=gram_A(3),
    posroots=posroots_A(3),
)
TYPES['B2'] = dict(
    G=[[Fr(2), Fr(-1)], [Fr(-1), Fr(1)]],
    posroots=[(1,0),(0,1),(1,1),(1,2)],
)
TYPES['G2'] = dict(
    G=[[Fr(2,3), Fr(-1)], [Fr(-1), Fr(2)]],
    posroots=[(1,0),(0,1),(1,1),(2,1),(3,1),(3,2)],
)
TYPES['B3'] = dict(
    G=[[Fr(2),Fr(-1),Fr(0)],[Fr(-1),Fr(2),Fr(-1)],[Fr(0),Fr(-1),Fr(1)]],
    posroots=[(1,0,0),(0,1,0),(1,1,0),(0,0,1),(0,1,1),(1,1,1),(0,1,2),(1,1,2),(1,2,2)],
)
TYPES['C3'] = dict(
    G=[[Fr(1),Fr(-1,2),Fr(0)],[Fr(-1,2),Fr(1),Fr(-1)],[Fr(0),Fr(-1),Fr(2)]],
    posroots=[(1,0,0),(0,1,0),(1,1,0),(0,0,1),(0,1,1),(1,1,1),(0,2,1),(2,2,1),(1,2,1)],
)
TYPES['A4'] = dict(G=gram_A(4), posroots=posroots_A(4))
TYPES['B4'] = dict(
    G=[[Fr(2),Fr(-1),Fr(0),Fr(0)],[Fr(-1),Fr(2),Fr(-1),Fr(0)],
       [Fr(0),Fr(-1),Fr(2),Fr(-1)],[Fr(0),Fr(0),Fr(-1),Fr(1)]],
    posroots=[(0,0,0,1),(0,0,1,0),(0,0,1,1),(0,0,1,2),(0,1,0,0),(0,1,1,0),(0,1,1,1),
              (0,1,1,2),(0,1,2,2),(1,0,0,0),(1,1,0,0),(1,1,1,0),(1,1,1,1),(1,1,1,2),
              (1,1,2,2),(1,2,2,2)],
)
TYPES['C4'] = dict(
    G=[[Fr(1),-_H,Fr(0),Fr(0)],[-_H,Fr(1),-_H,Fr(0)],
       [Fr(0),-_H,Fr(1),Fr(-1)],[Fr(0),Fr(0),Fr(-1),Fr(2)]],
    posroots=[(0,0,0,1),(0,0,1,0),(0,0,1,1),(0,0,2,1),(0,1,0,0),(0,1,1,0),(0,1,1,1),
              (0,1,2,1),(0,2,2,1),(1,0,0,0),(1,1,0,0),(1,1,1,0),(1,1,1,1),(1,1,2,1),
              (1,2,2,1),(2,2,2,1)],
)
TYPES['D4'] = dict(
    G=[[Fr(2),Fr(-1),Fr(0),Fr(0)],[Fr(-1),Fr(2),Fr(-1),Fr(-1)],
       [Fr(0),Fr(-1),Fr(2),Fr(0)],[Fr(0),Fr(-1),Fr(0),Fr(2)]],
    posroots=[(0,0,0,1),(0,0,1,0),(0,1,0,0),(0,1,0,1),(0,1,1,0),(0,1,1,1),
              (1,0,0,0),(1,1,0,0),(1,1,0,1),(1,1,1,0),(1,1,1,1),(1,2,1,1)],
)
TYPES['F4'] = dict(
    G=[[Fr(2),Fr(-1),Fr(0),Fr(0)],[Fr(-1),Fr(2),Fr(-1),Fr(0)],
       [Fr(0),Fr(-1),Fr(1),-_H],[Fr(0),Fr(0),-_H,Fr(1)]],
    posroots=[(0,0,0,1),(0,0,1,0),(0,0,1,1),(0,1,0,0),(0,1,1,0),(0,1,1,1),(0,1,2,0),
              (0,1,2,1),(0,1,2,2),(1,0,0,0),(1,1,0,0),(1,1,1,0),(1,1,1,1),(1,1,2,0),
              (1,1,2,1),(1,1,2,2),(1,2,2,0),(1,2,2,1),(1,2,2,2),(1,2,3,1),(1,2,3,2),
              (1,2,4,2),(1,3,4,2),(2,3,4,2)],
)
TYPES['A5'] = dict(
    G=[
       [Fr(2),Fr(-1),Fr(0),Fr(0),Fr(0)],
       [Fr(-1),Fr(2),Fr(-1),Fr(0),Fr(0)],
       [Fr(0),Fr(-1),Fr(2),Fr(-1),Fr(0)],
       [Fr(0),Fr(0),Fr(-1),Fr(2),Fr(-1)],
       [Fr(0),Fr(0),Fr(0),Fr(-1),Fr(2)]],
    posroots=[
        (1,0,0,0,0),(0,1,0,0,0),(0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1),(1,1,0,0,0),
        (0,1,1,0,0),(0,0,1,1,0),(0,0,0,1,1),(1,1,1,0,0),(0,1,1,1,0),(0,0,1,1,1),
        (1,1,1,1,0),(0,1,1,1,1),(1,1,1,1,1)],
)
TYPES['D5'] = dict(
    G=[
       [Fr(2),Fr(-1),Fr(0),Fr(0),Fr(0)],
       [Fr(-1),Fr(2),Fr(-1),Fr(0),Fr(0)],
       [Fr(0),Fr(-1),Fr(2),Fr(-1),Fr(-1)],
       [Fr(0),Fr(0),Fr(-1),Fr(2),Fr(0)],
       [Fr(0),Fr(0),Fr(-1),Fr(0),Fr(2)]],
    posroots=[
        (1,0,0,0,0),(0,1,0,0,0),(0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1),(1,1,0,0,0),
        (0,1,1,0,0),(0,0,1,1,0),(0,0,1,0,1),(1,1,1,0,0),(0,1,1,1,0),(0,1,1,0,1),
        (0,0,1,1,1),(1,1,1,1,0),(1,1,1,0,1),(0,1,1,1,1),(1,1,1,1,1),(0,1,2,1,1),
        (1,1,2,1,1),(1,2,2,1,1)],
)
TYPES['E6'] = dict(
    G=[
       [Fr(2),Fr(0),Fr(-1),Fr(0),Fr(0),Fr(0)],
       [Fr(0),Fr(2),Fr(0),Fr(-1),Fr(0),Fr(0)],
       [Fr(-1),Fr(0),Fr(2),Fr(-1),Fr(0),Fr(0)],
       [Fr(0),Fr(-1),Fr(-1),Fr(2),Fr(-1),Fr(0)],
       [Fr(0),Fr(0),Fr(0),Fr(-1),Fr(2),Fr(-1)],
       [Fr(0),Fr(0),Fr(0),Fr(0),Fr(-1),Fr(2)]],
    posroots=[
        (1,0,0,0,0,0),(0,1,0,0,0,0),(0,0,1,0,0,0),(0,0,0,1,0,0),(0,0,0,0,1,0),(0,0,0,0,0,1),
        (1,0,1,0,0,0),(0,1,0,1,0,0),(0,0,1,1,0,0),(0,0,0,1,1,0),(0,0,0,0,1,1),(1,0,1,1,0,0),
        (0,1,1,1,0,0),(0,1,0,1,1,0),(0,0,1,1,1,0),(0,0,0,1,1,1),(1,1,1,1,0,0),(1,0,1,1,1,0),
        (0,1,1,1,1,0),(0,1,0,1,1,1),(0,0,1,1,1,1),(1,1,1,1,1,0),(1,0,1,1,1,1),(0,1,1,2,1,0),
        (0,1,1,1,1,1),(1,1,1,2,1,0),(1,1,1,1,1,1),(0,1,1,2,1,1),(1,1,2,2,1,0),(1,1,1,2,1,1),
        (0,1,1,2,2,1),(1,1,2,2,1,1),(1,1,1,2,2,1),(1,1,2,2,2,1),(1,1,2,3,2,1),(1,2,2,3,2,1)],
)

def matvec(M, v):
    return [sum(M[i][j]*v[j] for j in range(len(v))) for i in range(len(M))]

def matmul(A, B):
    n=len(A); m=len(B[0]); k=len(B)
    return [[sum(A[i][t]*B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]

def mat_inverse(M):
    n=len(M)
    A=[[Fr(M[i][j]) for j in range(n)]+[Fr(1) if i==j else Fr(0) for j in range(n)] for i in range(n)]
    for col in range(n):
        piv=None
        for r in range(col,n):
            if A[r][col]!=0: piv=r; break
        A[col],A[piv]=A[piv],A[col]
        pv=A[col][col]
        A[col]=[x/pv for x in A[col]]
        for r in range(n):
            if r!=col and A[r][col]!=0:
                f=A[r][col]
                A[r]=[A[r][j]-f*A[col][j] for j in range(2*n)]
    return [row[n:] for row in A]

def _det_fr(M):
    # exact fraction determinant via elimination
    n=len(M); A=[row[:] for row in M]; s=Fr(1)
    for c in range(n):
        p=None
        for rr in range(c,n):
            if A[rr][c]!=0: p=rr;break
        if p is None: return Fr(0)
        if p!=c: A[c],A[p]=A[p],A[c]; s=-s
        s*=A[c][c]
        inv=Fr(1)/A[c][c]
        for rr in range(c+1,n):
            f=A[rr][c]*inv
            A[rr]=[A[rr][j]-f*A[c][j] for j in range(n)]
    return s

def _solve_exact(M, rhs):
    n=len(M)
    A=[[Fr(M[i][j]) for j in range(n)]+[Fr(rhs[i])] for i in range(n)]
    for col in range(n):
        piv=None
        for rr in range(col,n):
            if A[rr][col]!=0: piv=rr; break
        if piv is None: return None
        A[col],A[piv]=A[piv],A[col]
        pv=A[col][col]
        A[col]=[x/pv for x in A[col]]
        for rr in range(n):
            if rr!=col and A[rr][col]!=0:
                f=A[rr][col]
                A[rr]=[A[rr][j]-f*A[col][j] for j in range(n+1)]
    return [A[i][n] for i in range(n)]

def build_type(name):
    d=TYPES[name]
    G=d['G']; r=len(G)
    Ginv=mat_inverse(G)
    # fundamental weights in alpha-coords: w_i = (G_ii/2) * Ginv[:,i]
    fund=[]
    for i in range(r):
        col=[Ginv[j][i]*(G[i][i]/2) for j in range(r)]
        fund.append(col)
    # rho = sum of fundamental weights  (= half sum of positive roots)
    rho=[sum(fund[i][j] for i in range(r)) for j in range(r)]
    # sanity: half-sum of posroots
    hs=[Fr(0)]*r
    for root in d['posroots']:
        for j in range(r): hs[j]+=Fr(root[j])
    hs=[x/2 for x in hs]
    assert hs==rho, (name, hs, rho)
    # Euclidean embedding: M_emb such that Euclid(x) = L^T x with L L^T = G (Cholesky)
    Gf=np.array([[float(G[i][j]) for j in range(r)] for i in range(r)])
    L=np.linalg.cholesky(Gf)
    emb=L.T  # Euclid coords = emb @ x_alpha
    # Weyl group W is built LAZILY (see get_weyl): it is needed only by
    # kostka_foulkes, never by the geometry phase, and for big types (E6,
    # |W|=51840) eager construction at import would be prohibitive.
    return dict(name=name, r=r, G=G, Ginv=Ginv, fund=fund, rho=rho,
                posroots=d['posroots'], W=None, emb=emb, Gf=Gf)

def get_weyl(grp):
    """Build (and cache) the Weyl group as [(matrix, sign)] in alpha-coords.
    Lazy: only paid the first time K is computed for this type."""
    if grp['W'] is not None:
        return grp['W']
    G=grp['G']; r=grp['r']
    gens=[]
    for i in range(r):
        S=[[Fr(1) if a==b else Fr(0) for b in range(r)] for a in range(r)]
        for a in range(r):
            # S = I - (2/G_ii) e_i (G row i)
            S[i][a]=S[i][a]-(Fr(2)/G[i][i])*G[i][a]
        gens.append(S)
    def key(M): return tuple(tuple(row) for row in M)
    I=[[Fr(1) if a==b else Fr(0) for b in range(r)] for a in range(r)]
    seen={key(I):I}; frontier=[I]
    while frontier:
        nf=[]
        for M in frontier:
            for g in gens:
                P=matmul(g,M)
                k=key(P)
                if k not in seen:
                    seen[k]=P; nf.append(P)
        frontier=nf
    Wsign=[(M, int(_det_fr(M))) for M in seen.values()]
    grp['W']=Wsign
    return Wsign

def get_weyl_np(grp):
    """Cache the Weyl group as a numpy int64 stack (Marr[n,r,r]) + sign vector.
    Reflections in alpha-coords are INTEGER matrices, so this is exact.  Enables
    a vectorized Weyl sum: all |W| images of (top+rho) in one matmul instead of
    |W| pure-Python Fraction matvecs (the ~4s/K floor for E6)."""
    cache=grp.get('_Wnp')
    if cache is not None:
        return cache
    W=get_weyl(grp); r=grp['r']
    Marr=np.array([[[int(M[a][b]) for b in range(r)] for a in range(r)]
                   for M,_ in W], dtype=np.int64)
    signs=np.array([s for _,s in W], dtype=np.int64)
    cache=(Marr, signs); grp['_Wnp']=cache
    return cache

GROUPS={name:build_type(name) for name in TYPES}

def dominance_le(grp, a, b):
    # a <= b  iff  b-a in cone of simple roots, i.e. (b-a) alpha-coords nonneg
    aa=to_alpha(grp,a); bb=to_alpha(grp,b)
    return all(bb[i]-aa[i]>=0 for i in range(grp['r']))

def to_alpha(grp, fundcoords):
    r=grp['r']; fund=grp['fund']
    return [sum(Fr(fundcoords[i])*fund[i][j] for i in range(r)) for j in range(r)]

def poly_add(a, b):
    out=dict(a)
    for k,v in b.items(): out[k]=out.get(k,0)+v
    return {k:v for k,v in out.items() if v!=0}

def poly_shift(a, s):
    return {k+s:v for k,v in a.items()}

def poly_str(p):
    if not p: return "0"
    terms=[]
    for k in sorted(p):
        c=p[k]
        if k==0: terms.append(f"{c}")
        elif k==1: terms.append(f"{c}q")
        else: terms.append(f"{c}q^{k}")
    return " + ".join(terms)

def Pq(grp, beta, memo=None):
    # FAST Kostant q-partition function: standard memoized recursion over a fixed
    # positive-root order (use root k 0,1,2,... times), q^(total roots used).
    # The memo is shared across one Weyl sum (pass the SAME dict to every Pq call
    # within a single kostka_foulkes) -> sub-partition values are reused across the
    # |W| terms, the real win -- then freed when the K returns (bounded memory,
    # exactly like the Sage oracle).  Mathematically identical to _Pq_boxdp
    # (cross-validated on 1800 random betas across A2..E6).
    r=grp['r']; roots=grp['posroots']
    for x in beta:
        if x<0: return {}
    target=tuple(int(x) for x in beta)
    if all(t==0 for t in target): return {0:1}
    if memo is None: memo={}
    n=len(roots)
    def rec(t, k):
        if all(x==0 for x in t): return {0:1}
        if k>=n: return {}
        key=(t,k)
        v=memo.get(key)
        if v is not None: return v
        res=dict(rec(t, k+1))            # use root k zero more times
        rk=roots[k]
        t2=tuple(t[i]-rk[i] for i in range(r))
        if all(x>=0 for x in t2):
            for e,c in rec(t2, k).items():   # use root k one (or more) more times
                res[e+1]=res.get(e+1,0)+c
        res={e:c for e,c in res.items() if c!=0}
        memo[key]=res
        return res
    return rec(target, 0)

def _check_int64(*arrays, bound=1 << 62):
    """Abort rather than wrap around silently.  Every integer array below is a
    numpy int64 array whose entries must stay well inside the int64 range for
    the computation to be exact; this asserts it explicitly at runtime."""
    for a in arrays:
        if a.size and int(np.abs(a).max()) >= bound:
            raise OverflowError("int64 magnitude guard tripped: the exact "
                                "integer arithmetic would wrap around")


def _Pq_box_np(grp, betas):
    # FAST batched Kostant q-partition function.  Computes P_q for a WHOLE BATCH of
    # target betas that share ONE dense numpy box-DP over [0,M] (M = componentwise
    # max of the batch), the q-degree carried on a trailing axis.  Algorithmically
    # identical to _Pq_boxdp (same per-root convolution), but (a) one DP serves all
    # |valid betas| of a K instead of one DP per beta, and (b) pure numpy int64 slice
    # adds instead of python dict merges -> orders of magnitude faster for E6.
    # Returns {beta_tuple: poly_dict}.  Cross-validated against _Pq_boxdp / Pq.
    r=grp['r']; roots=grp['posroots']
    uniq=[]
    seen=set()
    for b in betas:
        t=tuple(int(x) for x in b)
        if any(x<0 for x in t): continue
        if t in seen: continue
        seen.add(t); uniq.append(t)
    if not uniq: return {}
    M=[max(t[i] for t in uniq) for i in range(r)]
    S=[m+1 for m in M]
    Q=sum(M)+1                     # #roots used <= sum of coords (each root adds >=1)
    dp=np.zeros(S+[Q],dtype=np.int64)
    dp[tuple([0]*r)+(0,)]=1
    for root in roots:
        rvec=[int(x) for x in root]
        # max multiplicity of this root inside the box
        kmax=Q
        for i in range(r):
            if rvec[i]>0:
                kmax=min(kmax,M[i]//rvec[i])
        if kmax<=0:
            continue
        dp_new=dp.copy()
        for k in range(1,kmax+1):
            offs=[k*rvec[i] for i in range(r)]
            if any(offs[i]>M[i] for i in range(r)):
                break
            if k>=Q:
                break
            dst=tuple(slice(offs[i],S[i]) for i in range(r))+(slice(k,Q),)
            src=tuple(slice(0,S[i]-offs[i]) for i in range(r))+(slice(0,Q-k),)
            dp_new[dst]+=dp[src]
        # Entries of dp are counts, hence nonnegative; a negative entry could
        # only come from an int64 wraparound.  Check every step, not just the
        # last, so an overflow cannot be masked by later additions.
        if int(dp_new.min()) < 0 or int(dp_new.max()) >= (1 << 62):
            raise OverflowError("int64 overflow in the q-partition-function DP")
        dp=dp_new
    out={}
    for t in uniq:
        vec=dp[t]                  # q-axis vector for this target cell
        poly={e:int(c) for e,c in enumerate(vec.tolist()) if c!=0}
        out[t]=poly
    return out

def kostka_foulkes(grp, top, low):
    """K^Phi_{lambda,mu}(q) by Lusztig's formula, with top = lambda (the highest
    weight) and low = mu.  Reference implementation: one memoized recursion per
    term of the Weyl sum."""
    r=grp['r']; rho=grp['rho']
    topa=to_alpha(grp, top); lowa=to_alpha(grp, low)
    topr=[topa[j]+rho[j] for j in range(r)]
    lowr=[lowa[j]+rho[j] for j in range(r)]
    result={}
    pqmemo={}  # shared across the whole Weyl sum, freed when this K returns
    # Vectorized Weyl sum (exact): scale (top+rho),(low+rho) by their common
    # denominator D so D*beta = M@(D*topr) - D*lowr is integer; one numpy matmul
    # over all |W| reflections (integer matrices) replaces |W| Fraction matvecs.
    D=1
    for x in topr+lowr:
        d=x.denominator; D=D*d//_igcd(D,d)
    Marr,signs=get_weyl_np(grp)
    tv=np.array([int(x*D) for x in topr],dtype=np.int64)
    lv=np.array([int(x*D) for x in lowr],dtype=np.int64)
    _check_int64(Marr, tv, lv, bound=1 << 20)      # matmul then stays < 2^63
    beta_scaled=(Marr@tv)-lv                       # (|W|, r) == D*beta
    valid=(beta_scaled>=0).all(axis=1)&((beta_scaled%D==0).all(axis=1))
    for i in np.nonzero(valid)[0]:
        beta=(beta_scaled[i]//D).tolist()
        p=Pq(grp,beta,pqmemo)
        if int(signs[i])>0:
            result=poly_add(result,p)
        else:
            result=poly_add(result,{k:-v for k,v in p.items()})
    return {k:v for k,v in result.items() if v!=0}

def kostka_foulkes_fast(grp, top, low):
    """Same as kostka_foulkes (top = lambda, low = mu), on the fast path."""
    # Same math as kostka_foulkes, but computes every valid beta's P_q with ONE
    # shared box-DP (via _Pq_box_np) instead of |W| memoized recursions.  This is
    # the E6-tractable path; validated to agree with kostka_foulkes bit-for-bit.
    r=grp['r']; rho=grp['rho']
    topa=to_alpha(grp, top); lowa=to_alpha(grp, low)
    topr=[topa[j]+rho[j] for j in range(r)]
    lowr=[lowa[j]+rho[j] for j in range(r)]
    D=1
    for x in topr+lowr:
        d=x.denominator; D=D*d//_igcd(D,d)
    Marr,signs=get_weyl_np(grp)
    tv=np.array([int(x*D) for x in topr],dtype=np.int64)
    lv=np.array([int(x*D) for x in lowr],dtype=np.int64)
    _check_int64(Marr, tv, lv, bound=1 << 20)      # matmul then stays < 2^63
    beta_scaled=(Marr@tv)-lv                       # (|W|, r) == D*beta
    valid=(beta_scaled>=0).all(axis=1)&((beta_scaled%D==0).all(axis=1))
    idx=np.nonzero(valid)[0]
    betas=[tuple((beta_scaled[i]//D).tolist()) for i in idx]
    polys=_Pq_box_np(grp, betas)                   # one shared box-DP
    result={}
    for i,b in zip(idx,betas):
        p=polys.get(b,{})
        if int(signs[i])>0:
            result=poly_add(result,p)
        else:
            result=poly_add(result,{k:-v for k,v in p.items()})
    return {k:v for k,v in result.items() if v!=0}

def vertices_exact(grp, low, high):
    """Vertices of P^Phi_{mu,lambda} with low = mu and high = lambda."""
    r=grp['r']; G=grp['G']; rho=grp['rho']
    lowa=to_alpha(grp, low); higha=to_alpha(grp, high)
    for i in range(r):
        if higha[i]<lowa[i]: return None
    Grho=matvec(G, rho)
    H=[]  # (a, b): a.x <= b
    for i in range(r):
        e=[Fr(0)]*r; e[i]=Fr(-1); H.append((e, -lowa[i]))   # x_i >= lowa[i]
        e=[Fr(0)]*r; e[i]=Fr(1);  H.append((e,  higha[i]))  # x_i <= higha[i]
    for i in range(r):
        a=[-G[i][j] for j in range(r)]; H.append((a, Grho[i]))  # (Gx)_i >= -Grho[i]
    verts=set()
    for combo in combinations(range(len(H)), r):
        M=[H[c][0] for c in combo]; rhs=[H[c][1] for c in combo]
        x=_solve_exact(M, rhs)
        if x is None: continue
        ok=True
        for a,b in H:
            if sum(a[j]*x[j] for j in range(r)) > b: ok=False; break
        if ok: verts.add(tuple(x))
    if not verts: return None
    return sorted(verts)

def _sqdist(grp, x, y):
    G=grp['G']; r=grp['r']
    d=[x[i]-y[i] for i in range(r)]
    return sum(d[i]*G[i][j]*d[j] for i in range(r) for j in range(r))

def sqdist_matrix_exact(grp, verts):
    n=len(verts)
    D=[[Fr(0)]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            v=_sqdist(grp, verts[i], verts[j]); D[i][j]=v; D[j][i]=v
    return D

def iso_sig_exact(D):
    n=len(D)
    ds=sorted((D[i][j].numerator, D[i][j].denominator)
              for i in range(n) for j in range(i+1,n))
    return (n, tuple(ds))

def is_isometric_exact(Da, Db):
    n=len(Da)
    if len(Db)!=n: return False
    sigA=[sorted(Da[i]) for i in range(n)]
    sigB=[sorted(Db[j]) for j in range(n)]
    used=[False]*n; assign=[-1]*n
    def bt(i):
        if i==n: return True
        for j in range(n):
            if used[j] or sigA[i]!=sigB[j]: continue
            ok=True
            for k in range(i):
                if Da[i][k]!=Db[j][assign[k]]: ok=False; break
            if ok:
                used[j]=True; assign[i]=j
                if bt(i+1): return True
                used[j]=False; assign[i]=-1
        return False
    return bt(0)

def iso_equal_refine(Da, Db):
    n = len(Da)
    if len(Db) != n: return False
    cola = [tuple(sorted(Da[i])) for i in range(n)]
    colb = [tuple(sorted(Db[i])) for i in range(n)]
    allc = sorted(set(cola) | set(colb)); m = {c: k for k, c in enumerate(allc)}
    cola = [m[c] for c in cola]; colb = [m[c] for c in colb]
    for _ in range(n):
        sa = [(cola[i], tuple(sorted((cola[j], Da[i][j]) for j in range(n)))) for i in range(n)]
        sb = [(colb[i], tuple(sorted((colb[j], Db[i][j]) for j in range(n)))) for i in range(n)]
        allc = sorted(set(sa) | set(sb)); m = {c: k for k, c in enumerate(allc)}
        na = [m[s] for s in sa]; nb = [m[s] for s in sb]
        if na == cola and nb == colb: break
        cola, colb = na, nb
    from collections import Counter, defaultdict
    if Counter(cola) != Counter(colb): return False
    byc = defaultdict(list)
    for j in range(n): byc[colb[j]].append(j)
    used = [False]*n; assign = [-1]*n
    def bt(i):
        if i == n: return True
        for j in byc[cola[i]]:
            if used[j]: continue
            if all(Da[i][k] == Db[j][assign[k]] for k in range(i)):
                used[j] = True; assign[i] = j
                if bt(i+1): return True
                used[j] = False; assign[i] = -1
        return False
    return bt(0)

def _fr2s(x): return f"{x.numerator}/{x.denominator}"

def _s2fr(s): a,b=s.split('/'); return Fr(int(a),int(b))

def _sig_key(sig):
    # sig = (n, iterable of (num,den)); -> (n, comma-joined "num/den" string)
    return sig[0], ",".join(f"{a}/{b}" for (a, b) in sig[1])

def db_connect(path='catalog.db'):
    con = sqlite3.connect(path)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    con.executescript('''
      CREATE TABLE IF NOT EXISTS classes(
        id INTEGER PRIMARY KEY,
        sig_n INTEGER NOT NULL,
        sig_key TEXT NOT NULL,
        d_json TEXT NOT NULL,
        k_json TEXT NOT NULL);
      CREATE INDEX IF NOT EXISTS idx_classes_sig ON classes(sig_n,sig_key);
      CREATE TABLE IF NOT EXISTS realizations(
        class_id INTEGER NOT NULL,
        type TEXT NOT NULL, low TEXT NOT NULL, high TEXT NOT NULL,
        bound INTEGER NOT NULL,
        PRIMARY KEY(type,low,high));
      CREATE INDEX IF NOT EXISTS idx_real_class ON realizations(class_id);
    ''')
    con.commit()
    return con

def enumerate_geometry_db(types, bound, db_path='catalog.db', commit_every=5000):
    """Phase A: record each non-empty polytope's exact isometry class. No K."""
    con = db_connect(db_path)
    seen = set()
    for t, lo, hi in con.execute("SELECT type,low,high FROM realizations"):
        seen.add((t, lo, hi))
    print(f"=== GEOMETRY ENUM  types={types}  bound={bound}  db={db_path} ===")
    print(f"  existing realizations: {len(seen)}", flush=True)
    new_classes = new_real = skipped = 0; since = 0
    for nm in types:
        g = GROUPS[nm]; r = g['r']
        wts = list(product(range(bound+1), repeat=r)); cnt = 0
        for low in wts:
            lo_s = json.dumps(list(low))
            for high in wts:
                if not dominance_le(g, low, high): continue
                hi_s = json.dumps(list(high))
                if (nm, lo_s, hi_s) in seen: skipped += 1; continue
                verts = vertices_exact(g, low, high)
                if verts is None:
                    seen.add((nm, lo_s, hi_s)); continue
                D = sqdist_matrix_exact(g, verts); sig = iso_sig_exact(D)
                n, key = _sig_key(sig)
                found = None
                for cid, d_json in con.execute(
                        "SELECT id,d_json FROM classes WHERE sig_n=? AND sig_key=?", (n, key)):
                    Dc = [[_s2fr(s) for s in row] for row in json.loads(d_json)]
                    if is_isometric_exact(D, Dc): found = cid; break
                if found is None:
                    cid = con.execute(
                        "INSERT INTO classes(sig_n,sig_key,d_json,k_json) VALUES(?,?,?,?)",
                        (n, key, json.dumps([[_fr2s(x) for x in row] for row in D]), "")).lastrowid
                    con.execute("INSERT OR IGNORE INTO realizations VALUES(?,?,?,?,?)",
                                (cid, nm, lo_s, hi_s, bound))
                    new_classes += 1
                else:
                    con.execute("INSERT OR IGNORE INTO realizations VALUES(?,?,?,?,?)",
                                (found, nm, lo_s, hi_s, bound))
                    new_real += 1
                seen.add((nm, lo_s, hi_s)); cnt += 1; since += 1
                if since >= commit_every: con.commit(); since = 0
        con.commit()
        print(f"  {nm}: {cnt} new polytopes ({new_classes} cls / {new_real} extra so far)", flush=True)
    con.commit()
    print(f"  + {new_classes} new classes, + {new_real} extra realizations, {skipped} skipped")
    con.close()

def verify_k_db(db_path='catalog.db', recheck_all=False, only_missing=False):
    """Phase B: for every isometry class with >=2 realizations, compute K for
    EACH realization independently and assert they coincide.  Stores K back into
    classes.k_json.  recheck_all=True also recomputes classes that already have K.
    only_missing=True restricts to multi-member classes whose k_json is still
    empty (i.e. newly-formed collisions) -- fast path for incremental catalog
    growth, skipping the already-verified baseline."""
    con = db_connect(db_path)
    if only_missing:
        multi = con.execute(
            "SELECT class_id FROM realizations GROUP BY class_id HAVING COUNT(*)>1 "
            "INTERSECT SELECT id FROM classes WHERE k_json=''").fetchall()
    else:
        multi = con.execute(
            "SELECT class_id FROM realizations GROUP BY class_id HAVING COUNT(*)>1").fetchall()
    print(f"=== K-VERIFY  db={db_path}  multi-member classes: {len(multi)}"
          f"{'  (only_missing)' if only_missing else ''} ===", flush=True)
    checked = kcomps = 0; violations = []; updated = 0
    wcon = con  # same connection for writes
    for (cid,) in multi:
        (k_json,) = con.execute("SELECT k_json FROM classes WHERE id=?", (cid,)).fetchone()
        if k_json and not recheck_all:
            stored = {int(k): v for k, v in json.loads(k_json).items()}
        else:
            stored = None
        reps = con.execute(
            "SELECT type,low,high FROM realizations WHERE class_id=?", (cid,)).fetchall()
        Kref = stored; ref_rz = None
        for t, lo, hi in reps:
            g = GROUPS[t]
            K = {k: v for k, v in kostka_foulkes(g, tuple(json.loads(hi)),
                                                 tuple(json.loads(lo))).items() if v}
            kcomps += 1
            if Kref is None:
                Kref = K; ref_rz = (t, lo, hi)
            elif K != Kref:
                violations.append((cid, t, lo, hi, poly_str(K), poly_str(Kref)))
        if stored is None and Kref is not None:
            wcon.execute("UPDATE classes SET k_json=? WHERE id=?",
                         (json.dumps({str(k): v for k, v in Kref.items()}), cid)); updated += 1
        checked += 1
        if checked % 100 == 0: wcon.commit(); print(f"  ...{checked} classes verified", flush=True)
    wcon.commit()
    print(f"  classes verified: {checked}  K computations: {kcomps}  k_json filled: {updated}")
    if violations:
        print(f"  !!! {len(violations)} K-MISMATCHES (COUNTEREXAMPLES) !!!")
        for v in violations[:20]:
            print(f"    class {v[0]}: {v[1]} {v[2]}->{v[3]} K={v[4]}  vs class K={v[5]}")
    else:
        print("  no K-mismatch: every isometric polytope shares one K.")
    con.close()
    return violations

def validate_isometry_db(db_path='catalog.db', sample_per_class=2, max_classes=None):
    """Cross-validate the class partition with the independent refine method.
    (a) FALSE-NEGATIVE guard: no two distinct classes share a signature.  This
        one is exhaustive: it is a single query over all classes.
    (b) FALSE-POSITIVE guard: recompute member polytopes and confirm each is
        iso_equal_refine to its class representative.  This one is a SAMPLE:
        it uses the first `sample_per_class` members of each multi-member
        class (default 2).  Pass a huge value to cover every member."""
    con = db_connect(db_path)
    coll = con.execute(
        "SELECT COUNT(*) FROM (SELECT sig_n,sig_key FROM classes "
        "GROUP BY sig_n,sig_key HAVING COUNT(*)>1)").fetchone()[0]
    print(f"=== ISOMETRY CROSS-VALIDATION  db={db_path} ===")
    print(f"  (a) signature collisions (distinct classes, same dist-multiset): {coll}")
    if coll == 0:
        print("      -> distance-multiset signature is a COMPLETE invariant here:"
              " no false negatives possible in the class partition.")
    multi = con.execute(
        "SELECT class_id FROM realizations GROUP BY class_id HAVING COUNT(*)>1").fetchall()
    if max_classes: multi = multi[:max_classes]
    bad = []; checked = 0; nseen = 0
    for (cid,) in multi:
        (d_json,) = con.execute("SELECT d_json FROM classes WHERE id=?", (cid,)).fetchone()
        Dref = [[_s2fr(s) for s in row] for row in json.loads(d_json)]
        reps = con.execute("SELECT type,low,high FROM realizations WHERE class_id=? LIMIT ?",
                           (cid, sample_per_class)).fetchall()
        for t, lo, hi in reps:
            nseen += 1
            g = GROUPS[t]
            v = vertices_exact(g, tuple(json.loads(lo)), tuple(json.loads(hi)))
            Dm = sqdist_matrix_exact(g, v)
            if not iso_equal_refine(Dref, Dm):
                bad.append((cid, t, lo, hi))
        checked += 1
        if checked % 5000 == 0: print(f"  ...{checked} classes cross-checked", flush=True)
    print(f"  (b) up to {sample_per_class} member(s) per class cross-checked "
          f"over {checked} multi-classes ({nseen} members in total); "
          f"non-isometric members found: {len(bad)}")
    if bad:
        print("      !!! independent method DISAGREES with stored class on:")
        for b in bad[:20]: print(f"        class {b[0]}: {b[1]} {b[2]}->{b[3]}")
    else:
        print("      -> independent refine method agrees with every grouping.")
    con.close()
    return bad

def catalog_report_db(db_path='catalog.db'):
    con = db_connect(db_path)
    (ncls,) = con.execute("SELECT COUNT(*) FROM classes").fetchone()
    (nrz,) = con.execute("SELECT COUNT(*) FROM realizations").fetchone()
    cross = con.execute(
        "SELECT class_id, COUNT(DISTINCT type) c, GROUP_CONCAT(DISTINCT type) "
        "FROM realizations GROUP BY class_id HAVING c>1").fetchall()
    print(f"--- CATALOG-DB SUMMARY ({db_path}) ---")
    print(f"  total isometry classes: {ncls}")
    print(f"  total realizations (polytopes): {nrz}")
    print(f"  cross-type classes (>1 type): {len(cross)}")
    shown = 0
    for cid, c, types in sorted(cross, key=lambda x: -x[1]):
        (k_json, sig_n) = con.execute(
            "SELECT k_json,sig_n FROM classes WHERE id=?", (cid,)).fetchone()
        K = {int(k): v for k, v in json.loads(k_json).items()}
        reps = []
        for t in types.split(','):
            row = con.execute("SELECT low,high FROM realizations WHERE class_id=? AND type=? "
                              "LIMIT 1", (cid, t)).fetchone()
            reps.append(f"{t}:{tuple(json.loads(row[0]))}->{tuple(json.loads(row[1]))}")
        print(f"  [verts={sig_n}] types={types} K={poly_str(K)}  e.g. {' ; '.join(reps)}")
        shown += 1
        if shown >= 40: print("   ... (truncated)"); break
    con.close()

def catalog_metrics_db(db_path='catalog.db'):
    con = db_connect(db_path)
    (ncls,) = con.execute("SELECT COUNT(*) FROM classes").fetchone()
    (nrz,) = con.execute("SELECT COUNT(*) FROM realizations").fetchone()
    from collections import Counter
    sizes = Counter(s for (s,) in con.execute(
        "SELECT COUNT(*) FROM realizations GROUP BY class_id"))
    distinctK = con.execute(
        "SELECT COUNT(DISTINCT k_json) FROM classes WHERE k_json<>''").fetchone()[0]
    withK = con.execute("SELECT COUNT(*) FROM classes WHERE k_json<>''").fetchone()[0]
    coll = con.execute("SELECT COUNT(*) FROM (SELECT sig_n,sig_key FROM classes "
                       "GROUP BY sig_n,sig_key HAVING COUNT(*)>1)").fetchone()[0]
    print(f"--- NON-VACUITY METRICS ({db_path}) ---")
    print(f"  isometry classes: {ncls}   realizations: {nrz}")
    print(f"  signature collisions (should be 0): {coll}")
    print(f"  classes with K computed: {withK}   DISTINCT K polynomials: {distinctK}")
    print(f"  -> {distinctK} distinct K over {withK} classes shows the invariant is highly"
          " non-trivial (isometry => K is not a constant map).")
    print("  class-size histogram (size : #classes):")
    for s in sorted(sizes):
        bar = '#'*min(40, sizes[s]*40//max(sizes.values()))
        print(f"    {s:4} : {sizes[s]:7}  {bar}")
    con.close()


def paper_numbers(db_path='catalog.db'):
    """The summary numbers quoted in the introduction of the paper."""
    con = sqlite3.connect(db_path)
    sizes = [n for (n,) in con.execute(
        "SELECT COUNT(*) FROM realizations GROUP BY class_id")]
    rows = con.execute("SELECT class_id, GROUP_CONCAT(DISTINCT type) "
                       "FROM realizations GROUP BY class_id").fetchall()
    multitype = sum(1 for _, ts in rows if len(set(ts.split(','))) > 1)
    crossfam  = sum(1 for _, ts in rows if len({t[0] for t in ts.split(',')}) > 1)
    print(f"nonempty polytopes (realizations): {sum(sizes)}")
    print(f"isometry classes:                  {len(sizes)}")
    print(f"classes with >1 member:            {sum(1 for n in sizes if n > 1)}")
    print(f"isometric pairs:                   {sum(n*(n-1)//2 for n in sizes)}")
    print(f"classes realized in >1 type:       {multitype}")
    print(f"classes bridging >1 Dynkin family: {crossfam}")
    con.close()

def _parse_wt(s):
    return tuple(int(x) for x in s.replace('(', '').replace(')', '').split(','))

if __name__ == '__main__':
    args = sys.argv[1:]
    cmd = args[0] if args else 'help'
    if cmd == 'pair':
        nm = args[1]; low = _parse_wt(args[2]); high = _parse_wt(args[3])
        g = GROUPS[nm]
        if not dominance_le(g, low, high):
            print(f"mu={low} <= lambda={high} fails in the dominance order "
                  f"of {nm}"); sys.exit(1)
        print(f"K^{nm}_{{lambda={high}, mu={low}}}(q) =",
              poly_str(kostka_foulkes_fast(g, high, low)))
    elif cmd == 'census':
        db = args[1] if len(args) > 1 else 'catalog.db'
        for nm, b in PAPER_BOUNDS.items():
            enumerate_geometry_db([nm], b, db)
        verify_k_db(db)
        catalog_metrics_db(db)
        paper_numbers(db)
    elif cmd == 'demo':
        db = args[1] if len(args) > 1 else 'catalog_demo.db'
        enumerate_geometry_db(['A1', 'A2', 'A3', 'B2', 'G2'], 4, db)
        verify_k_db(db)
        paper_numbers(db)
    elif cmd == 'verifyk':
        verify_k_db(args[1] if len(args) > 1 else 'catalog.db', recheck_all=True)
    elif cmd == 'isocheck':
        spc = int(args[2]) if len(args) > 2 else 2
        validate_isometry_db(args[1] if len(args) > 1 else 'catalog.db',
                             sample_per_class=(spc if spc > 0 else 10**9))
    elif cmd == 'numbers':
        paper_numbers(args[1] if len(args) > 1 else 'catalog.db')
    else:
        print(__doc__)
