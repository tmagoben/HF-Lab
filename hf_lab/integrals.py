import math,numpy as np
BOHR_PER_ANGSTROM=1.8897261254578281
H_STO3G_EXP=np.array([3.42525091,0.62391373,0.16885540])
H_STO3G_COEF=np.array([0.15432897,0.53532814,0.44463454])
def primitive_norm(a): return (2*a/math.pi)**0.75
def boys0(T):
 T=float(T)
 if T<1e-10: return 1.0-T/3.0+T*T/10.0
 return 0.5*math.sqrt(math.pi/T)*math.erf(math.sqrt(T))
def overlap_primitive(a,A,b,B):
 A=np.asarray(A,float);B=np.asarray(B,float);p=a+b;mu=a*b/p;rab2=np.dot(A-B,A-B); return primitive_norm(a)*primitive_norm(b)*(math.pi/p)**1.5*math.exp(-mu*rab2)
def kinetic_primitive(a,A,b,B):
 A=np.asarray(A,float);B=np.asarray(B,float);p=a+b;mu=a*b/p;rab2=np.dot(A-B,A-B); return primitive_norm(a)*primitive_norm(b)*mu*(3-2*mu*rab2)*(math.pi/p)**1.5*math.exp(-mu*rab2)
def nuclear_attraction_primitive(a,A,b,B,C,Z):
 A=np.asarray(A,float);B=np.asarray(B,float);C=np.asarray(C,float);p=a+b;mu=a*b/p;rab2=np.dot(A-B,A-B);P=(a*A+b*B)/p
 return -Z*primitive_norm(a)*primitive_norm(b)*(2*math.pi/p)*math.exp(-mu*rab2)*boys0(p*np.dot(P-C,P-C))
def eri_primitive(a,A,b,B,c,C,d,D):
 A=np.asarray(A,float);B=np.asarray(B,float);C=np.asarray(C,float);D=np.asarray(D,float);p=a+b;q=c+d;mu=a*b/p;nu=c*d/q;P=(a*A+b*B)/p;Q=(c*C+d*D)/q
 pref=2*math.pi**2.5/(p*q*math.sqrt(p+q)); expo=math.exp(-mu*np.dot(A-B,A-B)-nu*np.dot(C-D,C-D)); return primitive_norm(a)*primitive_norm(b)*primitive_norm(c)*primitive_norm(d)*pref*expo*boys0((p*q/(p+q))*np.dot(P-Q,P-Q))
def contracted2(fn,A,B,*extra):
 return sum(H_STO3G_COEF[i]*H_STO3G_COEF[j]*fn(a,A,b,B,*extra) for i,a in enumerate(H_STO3G_EXP) for j,b in enumerate(H_STO3G_EXP))
def contracted_eri(A,B,C,D):
 s=0.0
 for i,a in enumerate(H_STO3G_EXP):
  for j,b in enumerate(H_STO3G_EXP):
   for k,c in enumerate(H_STO3G_EXP):
    for l,d in enumerate(H_STO3G_EXP): s+=H_STO3G_COEF[i]*H_STO3G_COEF[j]*H_STO3G_COEF[k]*H_STO3G_COEF[l]*eri_primitive(a,A,b,B,c,C,d,D)
 return s
