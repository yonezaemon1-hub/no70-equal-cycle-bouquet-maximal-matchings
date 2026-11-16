#!/usr/bin/env python3
"""Exact rational-interval verification of No.70 v0.8 interlacing hypothesis.

No floating point root guesses. For every input q the program uses SymPy's
exact real-root isolator, tests multiplicity one, and checks the stipulated
root order using disjoint rational endpoints, allowing only a common zero
of a and b in the last slot. The mathematical interval lemma in
THEOREM_V08.md upgrades each successful q to every positive integer r.
"""
import argparse, hashlib, json, sys
from pathlib import Path
import sympy as S

z = S.Symbol('z')

def get_paths(limit):
    F = [S.Integer(1), S.Integer(1), z]
    for n in range(3, limit):
        F.append(S.expand(z*(F[n-2] + F[n-3])))
    return F

def norm_states(F, q):
    m=q//3
    return tuple(S.Poly(S.cancel(v/z**m), z) for v in
                 (z*z*F[q-5], 2*z*F[q-2], F[q-1]))

def intervals(poly):
    iv=poly.intervals(eps=S.Rational(1, 10**9))
    assert len(iv)==poly.degree(), ('real-rootedness failed',str(poly))
    assert all(mult==1 for _,mult in iv), ('multiple root', str(poly))
    return [tuple(pair) for pair,mult in iv]

def certify(q,F):
    a,b,c=norm_states(F,q)
    d=(q-1)//2-q//3
    e=1 if q%2==0 else 0
    assert (a.degree(),b.degree(),c.degree())==(d,d+e,d),('degrees',q)
    assert all(p.LC()>0 for p in (a,b,c))
    ra,rb,rc=[intervals(p) for p in (a,b,c)]
    less=lambda x,y:x[1]<y[0]
    eq0=lambda x,y:x==(S.Integer(0),S.Integer(0)) and y==x
    for i in range(d):
        if e==0:
            assert less(rc[i],ra[i]),('c<a',q,i)
            assert less(ra[i],rb[i]) or (i==d-1 and eq0(ra[i],rb[i])),('a<=b',q,i)
            if i+1<d:assert less(rb[i],rc[i+1]),('b<c',q,i)
        else:
            assert less(rb[i],rc[i]),('b<c',q,i)
            assert less(rc[i],ra[i]),('c<a',q,i)
            assert less(ra[i],rb[i+1]) or (i==d-1 and eq0(ra[i],rb[i+1])),('a<=b',q,i)
    assert all(upper <= 0 for roots in (ra,rb,rc) for lower,upper in roots)
    # Full intervals are stored verbatim as exact rational strings.
    as_str=lambda rs:[[str(lo),str(hi)] for lo,hi in rs]
    return {'q':q,'d':d,'e':e,'a_roots':as_str(ra),'b_roots':as_str(rb),'c_roots':as_str(rc),
            'trailing_common_zero':bool(d and eq0(ra[-1],rb[-1]))}

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--start',type=int,default=5)
    pa.add_argument('--end',type=int,default=200)
    pa.add_argument('--output-dir',type=Path,default=Path(__file__).parent)
    args=pa.parse_args()
    assert 5<=args.start<=args.end<=1000
    F=get_paths(args.end)
    records=[certify(q,F) for q in range(args.start,args.end+1)]
    out={'status':'PASS','scope':'EXACT Q RANGE; ALL R BY A SEPARATE ELEMENTARY CONDITIONAL THEOREM',
         'from_q':args.start,'to_q':args.end,'q_families':len(records),
         'sympy_version':S.__version__,
         'all_polynomial_roots_real_simple_and_ordered':True,
         'all_r_theorem':'Conditional real-root counting lemma with certified hypotheses for each q; not a global analytical interlacing proof for arbitrary q',
         'certificates':records}
    args.output_dir.mkdir(parents=True,exist_ok=True)
    fname=args.output_dir/f'INTERLACING_Q{args.start}_Q{args.end}.json'
    contents=json.dumps(out,ensure_ascii=False,indent=2)+'\n'
    fname.write_text(contents,encoding='utf-8')
    print('NO70 V08 INTERLACING PASS',args.start,args.end,'q_count',len(records),'sha256',hashlib.sha256(contents.encode('utf-8')).hexdigest())
    print('OUTPUT',fname)

if __name__=='__main__':
    main()
