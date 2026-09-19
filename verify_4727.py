#!/usr/bin/env python3
"""Self-contained exact verifier for the simplified B,L certificate in case (4,7,2,7).

A and B are Gaussian-integer matrices; L is Gaussian-rational and R=LL*.
All rank and positivity checks use exact integer arithmetic (Bareiss/Sylvester).
Standard library only.
"""
from __future__ import annotations
import argparse, base64, gzip, hashlib, json, time
from math import lcm
if not __debug__:
    raise RuntimeError("Run without python -O")

PAYLOAD = """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"""
DATA = json.loads(gzip.decompress(base64.b64decode(PAYLOAD)))

def transpose(x): return [list(t) for t in zip(*x)]
def product(ar,ai,br,bi):
    n,k,m=len(ar),len(br),len(br[0])
    cr=[[sum(ar[i][t]*br[t][j]-ai[i][t]*bi[t][j] for t in range(k)) for j in range(m)] for i in range(n)]
    ci=[[sum(ar[i][t]*bi[t][j]+ai[i][t]*br[t][j] for t in range(k)) for j in range(m)] for i in range(n)]
    return cr,ci
def adjoint(ar,ai): return transpose(ar),[[-v for v in row] for row in transpose(ai)]
def gram_rows(ar,ai): return product(ar,ai,*adjoint(ar,ai))
def gram_cols(ar,ai): return product(*adjoint(ar,ai),ar,ai)
def pt(x,n):
    N=4*n
    return [[x[(j//n)*n+i%n][(i//n)*n+j%n] for j in range(N)] for i in range(N)]
def realify(ar,ai):
    n=len(ar)
    for i in range(n):
        for j in range(n):
            assert ar[i][j]==ar[j][i] and ai[i][j]==-ai[j][i], "not Hermitian"
    return [ar[i]+[-v for v in ai[i]] for i in range(n)]+[ai[i]+ar[i] for i in range(n)]
def bareiss_pd(ar,ai):
    a=realify(ar,ai); n=len(a); prev=1; h=hashlib.sha256(); maxbits=0
    for k in range(n):
        p=a[k][k]
        if p<=0: raise AssertionError(f"nonpositive leading principal determinant at {k+1}")
        maxbits=max(maxbits,p.bit_length()); h.update((str(p)+"\n").encode())
        if k==n-1: break
        col=[a[i][k] for i in range(n)]
        for i in range(k+1,n):
            for j in range(i,n):
                val,rem=divmod(p*a[i][j]-col[i]*col[j],prev)
                assert rem==0, "Bareiss division not exact"
                a[i][j]=a[j][i]=val
        for i in range(k+1,n): a[i][k]=a[k][i]=0
        prev=p
    return {"real_size":n,"positive_leading_minors":n,"largest_minor_bits":maxbits,"minors_sha256":h.hexdigest()}
def check_matrix(re,im,rows,cols):
    assert len(re)==len(im)==rows and all(len(row)==cols for row in re+im)
    assert all(type(v)==int for row in re+im for v in row)
def verify_data(d):
    start=time.monotonic(); n,k,ell=d["n"],d["k"],d["ell"]; N=4*n
    ar,ai,br,bi=d["A_re"],d["A_im"],d["B_re"],d["B_im"]
    lr,li=d["L_re"],d["L_im"]; q=d["L_den"]; dn,dd=d["delta"]; r=len(lr[0])
    for re,im in zip(ar,ai): check_matrix(re,im,k,n)
    check_matrix(br,bi,N,ell); check_matrix(lr,li,N,r)
    aar=[[v for t in range(4) for v in ar[t][i]] for i in range(k)]
    aai=[[v for t in range(4) for v in ai[t][i]] for i in range(k)]
    agr,agi=gram_cols(aar,aai); bgr,bgi=gram_rows(br,bi); rr,ri=gram_rows(lr,li)
    fg,fi=pt(agr,n),pt(agi,n); rg,rgi=pt(rr,n),pt(ri,n); scale=lcm(q*q,dd)
    sr=[[scale*(fg[i][j]+bgr[i][j])-(scale//(q*q))*rg[i][j]-(dn*(scale//dd) if i==j else 0) for j in range(N)] for i in range(N)]
    si=[[scale*(fi[i][j]+bgi[i][j])-(scale//(q*q))*rgi[i][j] for j in range(N)] for i in range(N)]
    out={"case":[4,n,k,ell],"B_nnz":sum(br[i][j]!=0 or bi[i][j]!=0 for i in range(N) for j in range(ell)),
         "L_shape":[N,r],"L_nnz":sum(lr[i][j]!=0 or li[i][j]!=0 for i in range(N) for j in range(r)),
         "L_den":q,"L_numerator_height":max(abs(v) for row in lr+li for v in row),"delta":[dn,dd]}
    out["A_full_row_rank"]=bareiss_pd(*gram_rows(aar,aai))
    out["B_full_column_rank"]=bareiss_pd(*gram_cols(br,bi))
    out["R_positive_semidefinite"]="R = L L*, structural exact certificate"
    out["SOS_slack_positive_definite"]=bareiss_pd(sr,si)
    out["exact_verified"]=True; out["seconds"]=time.monotonic()-start
    return out

def main():
    result=verify_data(DATA["4727"])
    print(json.dumps(result,indent=2))
    print("EXACT PASS: 4727")
if __name__=="__main__": main()
