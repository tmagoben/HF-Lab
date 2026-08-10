import numpy as np
def symmetric_orthogonalizer(S):
 e,U=np.linalg.eigh(S)
 if np.min(e)<=0: raise ValueError('overlap matrix must be positive definite')
 return U@np.diag(e**-0.5)@U.T
def build_fock(H,eri,P):
 J=np.einsum('pqrs,rs->pq',eri,P); K=np.einsum('prqs,rs->pq',eri,P); return H+J-0.5*K
def density_from_fock(F,S,n_electrons):
 X=symmetric_orthogonalizer(S); eps,Cp=np.linalg.eigh(X.T@F@X); C=X@Cp; Cocc=C[:,:n_electrons//2]; return eps,C,2*Cocc@Cocc.T
def electronic_energy(H,F,P): return 0.5*np.sum(P*(H+F))
def rhf_scf(H,S,eri,n_electrons,nuclear_repulsion=0.0,max_iter=100,energy_tol=1e-12,density_tol=1e-10):
 if n_electrons%2: raise ValueError('RHF requires an even number of electrons')
 _,_,P=density_from_fock(H,S,n_electrons); Eold=None
 for iteration in range(1,max_iter+1):
  F=build_fock(H,eri,P); eps,C,Pnew=density_from_fock(F,S,n_electrons); Fnew=build_fock(H,eri,Pnew); Eelec=electronic_energy(H,Fnew,Pnew); dP=np.linalg.norm(Pnew-P); dE=np.inf if Eold is None else abs(Eelec-Eold)
  if dE<energy_tol and dP<density_tol: return {'electronic_energy':float(Eelec),'total_energy':float(Eelec+nuclear_repulsion),'orbital_energies':eps,'C':C,'P':Pnew,'iterations':iteration,'converged':True}
  P=Pnew;Eold=Eelec
 raise RuntimeError('SCF did not converge')
