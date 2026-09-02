# Restricted Hartree-Fock equations

For a closed-shell density matrix $P$, the Fock matrix is

$$
F_{\mu\nu}
=
H_{\mu\nu}
+
\sum_{\lambda\sigma}
P_{\lambda\sigma}
\left[
(\mu\nu\mid\lambda\sigma)
-
\frac{1}{2}(\mu\lambda\mid\nu\sigma)
\right].
$$

Expanding the molecular orbitals in a finite basis and applying the variational
principle gives the Roothaan-Hall generalized eigenvalue problem,

$$
FC = SC\varepsilon.
$$

Here $F$ is the Fock matrix, $S$ is the overlap matrix, $C$ contains the molecular
orbital coefficients, and $\varepsilon$ is the diagonal matrix of orbital energies.
