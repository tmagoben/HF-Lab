try:
 from pyscf import gto,scf
except ImportError: raise SystemExit("install optional PySCF extra")
from hf_lab import h2_sto3g,rhf_scf
R=.74; d=h2_sto3g(R); ours=rhf_scf(d['Hcore'],d['S'],d['eri'],2,d['nuclear_repulsion'])['total_energy']; mol=gto.M(atom=f'H 0 0 {-R/2}; H 0 0 {R/2}',unit='Angstrom',basis='sto-3g',verbose=0); ref=scf.RHF(mol).run().e_tot; print('HF-Lab',ours); print('PySCF',ref); print('difference',ours-ref)
