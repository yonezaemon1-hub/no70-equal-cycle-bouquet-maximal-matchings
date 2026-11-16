#!/usr/bin/env python3
"""Supplemental finite exact interval audit of guidepost inequalities in v0.10.

This is NOT a proof oracle. It validates the interval placements newly made
explicit in Sections 2--3. Frozen v0.8/v0.9 auditors remain separate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as S

x=S.symbols('x')


def build(kmax):
    F=[S.Integer(1),S.Integer(1),x]
    for n in range(3,3*(kmax+2)+3):
        F.append(S.expand(x*(F[n-2]+F[n-3])))
    def series(k):
        return tuple(S.Poly(S.cancel(F[j]/x**power),x) for j,power in ((3*k,k),(3*k+1,k),(3*k+2,k+1)))
    return [series(k) for k in range(kmax+2)]


def realroots(poly):
    intervals=poly.intervals(eps=S.Rational(1,10**18))
    assert len(intervals)==poly.degree()
    assert all(multiplicity==1 for _,multiplicity in intervals)
    assert all(pair[1]<0 for pair,_ in intervals)
    return [pair for pair,_ in intervals]


def lt(u,v):
    return u[1]<v[0]


def chain(*roots):
    assert all(lt(u,v) for u,v in zip(roots,roots[1:])), roots


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--max-q',type=int,default=48)
    ap.add_argument('--output',type=Path,default=Path('BRIDGE_EVIDENCE_V010.json'))
    args=ap.parse_args()
    assert 7<=args.max_q<=130
    maxk=args.max_q//3+2
    polys=build(maxk)
    roots=[tuple(realroots(p) for p in triple) for triple in polys]
    count=0
    for k in range(1,maxk+1):
        old_a,old_b,old_c=roots[k]
        new_a,new_b,new_c=roots[k+1]
        if k%2==0:
            s=k//2
            for i in range(s):
                chain(new_b[i],old_c[i],new_c[i],new_a[i],old_b[i],new_b[i+1]); count+=1
                chain(old_c[i],old_a[i],new_a[i]); count+=1
        else:
            s=(k-1)//2
            for i in range(s):
                chain(new_c[i],new_a[i],old_b[i],new_b[i],old_c[i],new_c[i+1]); count+=1
                chain(old_c[i],old_a[i],new_a[i+1]); count+=1
            chain(new_c[s],new_a[s],old_b[s],new_b[s]); count+=1
    q_count=0
    for q in range(7,args.max_q+1):
        if q%3 !=1:continue
        t=(q-1)//3-1
        (ta,tb,tc)=roots[t-1]
        (ma,mb,mc)=roots[t]
        (na,nb,nc)=roots[t+1]
        if t%2==0:
            s=t//2
            for i in range(s-1):
                chain(mc[i],ma[i],na[i],mb[i],tc[i],mc[i+1]);count+=1
            chain(mc[s-1],ma[s-1],na[s-1],mb[s-1]);count+=1
        else:
            s=(t-1)//2
            chain(na[0],mb[0]);count+=1
            for i in range(s):
                chain(mb[i],tc[i],mc[i],ma[i],na[i+1],mb[i+1]);count+=1
        q_count+=1
    data={'status':'PASS','max_q':args.max_q,'braid_transitions_checked':maxk,
          'q_3k_plus_1_two_step_classes_checked':q_count,
          'strict_chain_assertions':count,
          'sympy':S.__version__,
          'scope':'finite exact rational-interval checks of the explicit old/new guidepost chains; not a substitute for the induction proof'}
    raw=json.dumps(data,indent=2,ensure_ascii=False)+'\n'
    args.output.write_text(raw,encoding='utf-8')
    print('NO70 V010 GUIDePOST CHECK PASS',data,'SHA256',hashlib.sha256(raw.encode()).hexdigest())


if __name__=='__main__':main()
