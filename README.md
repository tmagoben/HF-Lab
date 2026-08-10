# HF-Lab

HF-Lab is a from-scratch restricted Hartree-Fock implementation for H$_2$ in the
STO-3G basis. Unlike the earlier toy-integral prototype, this version evaluates the
required one- and two-electron integrals analytically from contracted Gaussian basis
functions.

At an H-H distance of 0.74 Angstrom, the implementation gives approximately

$$
E_{\mathrm{RHF}}^{\mathrm{STO-3G}}=-1.1167593074\ E_h.
$$

## Implemented explicitly

- normalized primitive $s$ Gaussian functions;
- STO-3G contractions for hydrogen 1s;
- overlap, kinetic, nuclear-attraction, and electron-repulsion integrals;
- symmetric orthogonalization $X=S^{-1/2}$;
- closed-shell Coulomb and exchange matrices;
- Roothaan-Hall SCF iteration;
- electronic and nuclear-repulsion energies;
- regression tests for integral symmetries and the H$_2$ reference energy.

```bash
pip install -e ".[dev]"
python examples/h2_074.py
pytest -q
```

Atomic units are used internally; input bond lengths in the examples are Angstrom.
