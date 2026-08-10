import numpy as np
from .integrals import BOHR_PER_ANGSTROM,contracted2,contracted_eri,overlap_primitive,kinetic_primitive,nuclear_attraction_primitive
def h2_sto3g(distance_angstrom=0.74):
 R=float(distance_angstrom)*BOHR_PER_ANGSTROM; centers=[np.array([0.,0.,-R/2]),np.array([0.,0.,R/2])]; S=np.zeros((2,2));T=np.zeros((2,2));V=np.zeros((2,2));eri=np.zeros((2,2,2,2))
 for p,A in enumerate(centers):
  for q,B in enumerate(centers):
   S[p,q]=contracted2(overlap_primitive,A,B); T[p,q]=contracted2(kinetic_primitive,A,B); V[p,q]=sum(contracted2(nuclear_attraction_primitive,A,B,C,1.0) for C in centers)
 for p,A in enumerate(centers):
  for q,B in enumerate(centers):
   for r,C in enumerate(centers):
    for s,D in enumerate(centers): eri[p,q,r,s]=contracted_eri(A,B,C,D)
 return {'S':S,'Hcore':T+V,'eri':eri,'nuclear_repulsion':1.0/R,'distance_bohr':R}
