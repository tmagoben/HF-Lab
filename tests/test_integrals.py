import numpy as np
from hf_lab.h2 import h2_sto3g
def test_overlap_and_eri_symmetry():
 d=h2_sto3g(.74); S=d['S']; eri=d['eri']; assert np.allclose(S,S.T); assert np.all(np.linalg.eigvalsh(S)>0)
 for p in range(2):
  for q in range(2):
   for r in range(2):
    for s in range(2):
     v=eri[p,q,r,s]; assert np.isclose(v,eri[q,p,r,s]); assert np.isclose(v,eri[p,q,s,r]); assert np.isclose(v,eri[r,s,p,q])
