from hf_lab import h2_sto3g,rhf_scf
def test_h2_074_sto3g_reference():
 d=h2_sto3g(.74); out=rhf_scf(d['Hcore'],d['S'],d['eri'],2,d['nuclear_repulsion']); assert out['converged']; assert abs(out['total_energy']-(-1.1167593073781579))<2e-10
