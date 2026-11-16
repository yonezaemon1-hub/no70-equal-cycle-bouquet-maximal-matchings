#!/usr/bin/env python3
"""No.70 v0.9 independent algebraic and exact root checks (not a proof oracle).
Requires SymPy. Checks two independent F construction methods, 3-residue
transfer, universal T(q) in sampled q, and exact Sturm root counts.
"""
import argparse, hashlib, json
from math import comb
from pathlib import Path
import sympy as S
x = S.symbols('x')

def recurrence(max_n):
    F=[S.Integer(1),S.Integer(1),x]
    for n in range(3,max_n+1): F.append(S.expand(x*(F[n-2]+F[n-3])))
    return F

def closed(n):
    return S.Add(*[S.Integer(comb(j+1,n-2*j))*x**j for j in range((n+1)//3,(n//2)+1) if 0<=n-2*j<=j+1]) if n else S.Integer(1)

def polys(F,q):
    m=q//3
    return [S.Poly(S.cancel(y/x**m),x) for y in (x*x*F[q-5],2*x*F[q-2],F[q-1])]

def roots(p):
    # exact rational intervals. May return point intervals for rational roots.
    iv=p.intervals(eps=S.Rational(1,10**12))
    assert len(iv)==p.degree() and all(mult==1 for _,mult in iv)
    assert all(pair[1]<=0 for pair,mult in iv)
    return [pair for pair,mult in iv]

def less(u,v): return u[1]<v[0]
def root_braid(k,A,B,C):
    ra,rb,rc=roots(A),roots(B),roots(C)
    if k%2==0:
        assert len(ra)==len(rb)==len(rc)==k//2
        for i in range(len(ra)):
            assert less(rc[i],ra[i]) and less(ra[i],rb[i])
            if i+1<len(ra):assert less(rb[i],rc[i+1])
    else:
        assert len(rb)==(k+1)//2 and len(ra)==len(rc)==(k-1)//2
        for i in range(len(ra)):
            assert less(rb[i],rc[i]) and less(rc[i],ra[i]) and less(ra[i],rb[i+1])

def hypothesis(q,F):
    a,b,c=polys(F,q); ra,rb,rc=map(roots,(a,b,c))
    d=(q-1)//2-q//3;e=int(q%2==0)
    assert (a.degree(),b.degree(),c.degree())==(d,d+e,d)
    assert len(ra)==len(rc)==d and len(rb)==d+e
    shared=lambda u,v:u==(S.Integer(0),S.Integer(0)) and u==v
    for i in range(d):
        if e:
            assert less(rb[i],rc[i]) and less(rc[i],ra[i])
            assert less(ra[i],rb[i+1]) or(i==d-1 and shared(ra[i],rb[i+1]))
        else:
            assert less(rc[i],ra[i])
            assert less(ra[i],rb[i]) or(i==d-1 and shared(ra[i],rb[i]))
            if i+1<d:assert less(rb[i],rc[i+1])
    return d,e

def count_nonreal(q,r,F):
    a,b,c=polys(F,q)
    P=a.as_expr()**r+r*b.as_expr()*c.as_expr()**(r-1)
    pp=S.Poly(S.expand(P),x)
    # algebraic real roots with multiplicity from real-root isolation
    real=sum(m for _,m in pp.intervals())
    return pp.degree()-real, pp.degree(), real

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--max-q',type=int,default=45)
    ap.add_argument('--max-r-sturm',type=int,default=6)
    ap.add_argument('--max-q-sturm',type=int,default=19)
    ap.add_argument('--output',type=Path,default=Path('v09_evidence.json'))
    args=ap.parse_args(); maxq=args.max_q
    assert 5<=args.max_q<=150 and 1<=args.max_r_sturm<=12 and 5<=args.max_q_sturm<=maxq
    F=recurrence(maxq+3)
    for n in range(maxq+1):
        assert S.expand(F[n]-closed(n))==0, ('binomial expression',n)
    for k in range(maxq//3+1):
        A=S.Poly(S.cancel(F[3*k]/x**k),x)
        B=S.Poly(S.cancel(F[3*k+1]/x**k),x)
        C=S.Poly(S.cancel(F[3*k+2]/x**(k+1)),x)
        if k==0:assert (A.as_expr(),B.as_expr(),C.as_expr())==(1,1,1)
        else:
            Ak=S.cancel(F[3*(k-1)]/x**(k-1))
            Bk=S.cancel(F[3*(k-1)+1]/x**(k-1))
            Ck=S.cancel(F[3*(k-1)+2]/x**k)
            assert S.expand(A.as_expr()-(Ak+Bk))==0
            assert S.expand(B.as_expr()-(Bk+x*Ck))==0
            assert S.expand(C.as_expr()-(A.as_expr()+Ck))==0
        root_braid(k,A,B,C)
    for q in range(5,maxq+1):
        a,b,c=polys(F,q);k=q//3
        A=lambda j:S.cancel(F[3*j]/x**j)
        B=lambda j:S.cancel(F[3*j+1]/x**j)
        C=lambda j:S.cancel(F[3*j+2]/x**(j+1))
        if q%3==0: expect=[B(k-2),2*B(k-1),C(k-1)]
        elif q%3==1:expect=[x*C(k-2),2*x*C(k-1),A(k)]
        else:expect=[x*A(k-1),2*x*A(k),B(k)]
        assert all(S.cancel(actual.as_expr()-want)==0 for actual,want in zip((a,b,c),expect)),('identities',q)
        hypothesis(q,F)
    sturm=0
    for q in range(5,args.max_q_sturm+1):
        for r in range(1,args.max_r_sturm+1):
            n,deg,real=count_nonreal(q,r,F)
            d=(q-1)//2-q//3
            assert n==2*d*((r-1)//2),(q,r,n)
            sturm+=1
    data={'status':'PASS','max_q':maxq,'F_binomial_vs_recurrence_n_count':maxq+1,'braid_k_count':maxq//3+1,'universal_T_q_sample_count':maxq-4,
          'sturm_exact_q_range':[5,args.max_q_sturm],'sturm_exact_r_range':[1,args.max_r_sturm],'sturm_case_count':sturm,
          'scope':'Finite spot checks only; the all-q/all-r claim is justified by the separate analytic induction proof, not by this verifier',
          'sympy':S.__version__}
    raw=json.dumps(data,ensure_ascii=False,indent=2)+'\n';args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(raw,encoding='utf8')
    print('NO70 V09 FINITE INDEPENDENT CHECK PASS',data,'SHA256',hashlib.sha256(raw.encode()).hexdigest())
if __name__=='__main__':main()
