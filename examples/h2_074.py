from hf_lab import h2_sto3g,rhf_scf
data=h2_sto3g(0.74); out=rhf_scf(data['Hcore'],data['S'],data['eri'],2,data['nuclear_repulsion']); print('E_total =',out['total_energy']); print('iterations =',out['iterations'])
