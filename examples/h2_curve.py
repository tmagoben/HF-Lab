import numpy as np
from hf_lab import h2_sto3g,rhf_scf
for R in np.linspace(0.5,2.0,16):
 d=h2_sto3g(R); out=rhf_scf(d['Hcore'],d['S'],d['eri'],2,d['nuclear_repulsion']); print(f'{R:.4f} {out["total_energy"]:.12f}')
