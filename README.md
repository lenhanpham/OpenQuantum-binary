<p align="center">
  <a href="https://github.com/lenhanpham/OpenQuantum-binary">
    <picture>
      <img src="resources/oquantum-logo.svg" alt="OpenQuantum" style="width: 50%;">
    </picture>
  </a>
</p>
<p align="center">
  <a href="https://www.rust-lang.org/">
    <img src="https://img.shields.io/badge/Rust-1.70+-orange.svg" alt="Rust">
  </a>
</p

OpenQuantum is an ab initio molecular orbital (MO) calculation program written in Rust.
It performs Hartree-Fock and post-Hartree-Fock calculations within the Linear Combination
of Atomic Orbitals (LCAO) framework. Supports RHF, UHF, and ROHF for closed- and open-shell
molecules, MP2 and CCSD(T) correlation energies, geometry optimization (BFGS, Berny RFO,
primitive/DLC/TRIC internal coordinates) with analytic gradients, IRC and NEB pathways,
harmonic frequency analysis with analytical Hessians and thermochemistry, and effective
core potentials (ECPs).

> **Note**: The project is under active development and not yet ready for production use.

## Full Features, examples, and instructions
[https://lenhanpham.github.io/OpenQuantum](https://lenhanpham.github.io/OpenQuantum)


## Features

- **Molecular Input**: XYZ Cartesian and Z-matrix geometry formats, Angstrom or Bohr units
- **Basis Sets**: STO-3G, 6-31G, 6-31G\*, cc-pVDZ, cc-pVTZ, def2-SVP, def2-TZVP, and 500+ more (H–Xe and beyond)
- **One-Electron Integrals**: Overlap (S), kinetic energy (T), and nuclear attraction (V) via Obara-Saika recurrences
- **Two-Electron Integrals**: McMurchie-Davidson algorithm with 8-fold permutational symmetry; SP analytical fast paths for S/P shells (up to 10× speedup); Rys quadrature for D/F+ shells; direct SCF mode for large systems with symmetry-accelerated cache lookup; engine selection via `INT section` — applies uniformly to in-core ERI build, direct SCF, analytical gradients/Hessians, and semi-analytical (FD) Hessians
- **SCF Methods**: RHF, UHF, and ROHF with DIIS, level shifting, damping, Fermi broadening, multiple initial guesses (core Hamiltonian, Hückel, SAD), and Quadratic Convergence SCF (QC-SCF) with Newton-Raphson orbital optimization
- **Post-HF**: MP2 and CCSD(T) correlation energies
- **Analytical Gradients**: RHF and UHF analytic nuclear gradients with symmetry-accelerated ERI derivatives (skips symmetry-equivalent shell quartets); used by BFGS optimizer
- **Analytical Hessians**: Fully analytical RHF and UHF Hessians including CPHF response, occupied-occupied reorthonormalization, and analytic d²ERI integrals with symmetry acceleration (skips symmetry-equivalent quartets); semi-analytical (FD) path also propagates the selected ERI engine to all displaced SCF evaluations
- **Geometry Optimization**: BFGS optimizer with analytic nuclear gradients; Berny RFO algorithm with trust-radius step control via `opt=(algorithm=berny)`; **Rust Berny backend** selectable with `r`-prefixed route keywords (`rberny`, `opt=(algorithm=rberny)`, `ropt=(...)`) and backend-specific controls (`dihedral`, `superweakdih`, `energynoise`); Primitive, DLC, and TRIC internal-coordinate back-ends with topology-inferred primitive sets (bonds, angles, dihedrals, out-of-plane, linear-angle supplementary) and IC back-transform; Lindh model Hessian guess with BFGS/MSP/SR1/DFP/Bofill/PSB/MS update; 4-criterion style convergence; Kabsch-aligned RMSD displacement tracking; IRC path following (Gonzalez-Schlegel, mass-weighted, predictor-corrector, bidirectional); NEB pathway optimization with Henkelman energy-weighted tangent and climbing-image CI-NEB; **Transition-state search** with P-RFO step, Bofill Hessian update, initial TS Hessian initialization, periodic eigenvalue correction, and TS-specific trust radius (0.01 Å); **saddle-point-aware IC optimizer** with TS-BFGS Hessian update, P-RFO or Minimum Mode Following (MMF) step, sigma-based trust-radius schedule, and Davidson partial eigensolver (`sella` keyword); **GDIIS/GEDIIS** geometry-space DIIS extrapolation in IC optimizers (`opt=(diis)`); **Bofill Hessian update** for transition-state optimization (`opt=(update=bofill)`)
- **Frequency Analysis**: Harmonic vibrational frequencies via semi-analytical (finite-difference of gradients) or fully analytical Hessian; IR intensities, thermochemistry
- **Thermochemistry**: Zero-point energy, thermal corrections (U, H, G), entropy (translational + rotational + vibrational) via RRHO model with symmetry number
- **ECP Support**: Effective core potentials (LANL2DZ, Stuttgart, etc.)
- **Symmetry**: Point group detection, character tables, irrep assignment; symmetry acceleration for all ERI-heavy computations — energy (in-core + direct SCF), analytical gradient, and analytical Hessian (matching GRAD2E SymShl — skips symmetry-equivalent shell quartets before expensive integral evaluation)
- **CPHF Solver**: Coupled-perturbed Hartree-Fock equations for response properties; RHF and coupled 2×2 UHF spins
- **Analysis**: Mulliken population analysis, orbital energies, spin contamination ⟨S²⟩, dipole moments, EFG
- **Checkpointing**: Save and restart SCF from binary checkpoint files; checkpoint and other temporary files are written to a configurable scratch directory (see [Environment Variables](#environment-variables))

## Installation

### Download the binary file that fits your system
OpenQuantum binaries: 

| Platform | Target | Recommendation |
| :--- | :--- | :--- |
| Windows |[x86_64-pc-windows-msvc.zip](https://github.com/user-attachments/files/28456251/x86_64-pc-windows-msvc.zip) | Best compatibility |
| Ubuntu / Debian | [x86_64-unknown-linux-gnu.zip](https://github.com/user-attachments/files/28456254/x86_64-unknown-linux-gnu.zip) | Good |
| Rocky / RHEL | [x86_64-unknown-linux-gnu.zip](https://github.com/user-attachments/files/28456254/x86_64-unknown-linux-gnu.zip)| `musl` better for distribution |
| Rocky / RHEL | [x86_64-unknown-linux-musl.zip](https://github.com/user-attachments/files/28456257/x86_64-unknown-linux-musl.zip) | `musl` better for distribution |
| General Linux | [x86_64-unknown-linux-musl.zip](https://github.com/user-attachments/files/28456257/x86_64-unknown-linux-musl.zip)| Best (static, very portable) |
| macOS Intel | [x86_64-apple-darwin.zip](https://github.com/user-attachments/files/28456248/x86_64-apple-darwin.zip) | Good |
| macOS Apple Silicon | [aarch64-apple-darwin.zip](https://github.com/user-attachments/files/28456242/aarch64-apple-darwin.zip) | Good |




## License

NOT DECIDED YET

## Contributing

Le Nhan Pham
