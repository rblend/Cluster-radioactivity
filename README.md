# Cluster Decay Calculator

A Python implementation of the semi-classical model for nuclear cluster decay
half-lives and related quantities.

## Physical Model

The code implements the WKB (Gamow) tunnelling model described in:

- R. Blendowske, H. Walliser,  
  *Systematics of cluster-radioactivity-decay constants as suggested by  
  microscopic calculations*,  
  **Phys. Rev. Lett. 61 (1988) 1930**

- R. Blendowske, T. Fliessbach, H. Walliser,  
  *Microscopic calculation of the 14C decay of Ra isotopes*,  
  **Z. Phys. A339 (1991) 121**



### Key ingredients

| Quantity | Description |
|---|---|
| Nuclear potential | Elastic scattering potential (Christensen & Winter, Phys. Lett. 65B, 1976) |
| Coulomb potential | Point-charge approximation |
| Electronic correction | Lunney, Pearson & Thibault, Rev. Mod. Phys. 75 (2003) 1021 |
| Cluster binding energy | Wapstra et al., Nucl. Phys. A432 (1985) 1; semi-empirical fallback |
| Spectroscopic factor | S = (6.3×10⁻³)^((Ac−1)/3) for even, (3.2×10⁻³)^((Ac−1)/3) for odd |

## Installation

No special installation required. Dependencies:

```
numpy
scipy
```

Install with:

```bash
pip install numpy scipy
```

## Usage

### Single decay

```python
from cluster_decays import cluster_decays

# Ra-226 --> C-14 + Pb-212
TAUlog, E, PEN, S, WKB, RI, RA, ES = cluster_decays(
    AP=226, ZP=88,   # parent nucleus
    AC=14,  ZC=6,    # emitted cluster
)

print(f"log10(T½ / s) = {TAUlog:.2f}")   # experimental: ~21.2
print(f"Q-value       = {E:.3f} MeV")
print(f"Penetrability = {PEN:.3e}")
```

### With explicit Q-value

```python
result = cluster_decays(226, 88, 14, 6, ENC=28.20)
```

### Odd-A parent (hindrance)

```python
result = cluster_decays(223, 88, 14, 6, odd=1)
```

### Full parameter signature

```python
cluster_decays(AP, ZP, AC, ZC, odd=0, L=0, ENC=None)
```

| Parameter | Description |
|---|---|
| `AP`, `ZP` | Mass number and atomic number of parent |
| `AC`, `ZC` | Mass number and atomic number of emitted cluster |
| `odd` | `0` (default): even-even parent; `1`: odd parent (smaller S) |
| `L` | Angular momentum of the emitted cluster (default `0`) |
| `ENC` | Q-value in MeV; computed from `atomic_mass_data.txt` if omitted |

### Return values

| Name | Description |
|---|---|
| `TAUlog` | log₁₀(T½ / s) |
| `E` | Decay energy including electronic correction (MeV) |
| `PEN` | Barrier penetrability |
| `S` | Spectroscopic (preformation) factor |
| `WKB` | Gamow decay constant λ (s⁻¹) |
| `RI` | Inner classical turning point (fm) |
| `RA` | Outer classical turning point (fm) |
| `ES` | Data-quality flags `[parent, daughter, cluster]` (1=exact, 0=estimated) |

### Systematics table

```bash
python cluster_decay_systematics.py
```

Prints a comprehensive table covering Ba, Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am,
and Cm parent nuclei with various cluster species.

## Atomic mass data

Q-values are read from `atomic_mass_data.txt` (whitespace-separated columns:
`Z  A  Ebind/A [keV]  status`). The file should be based on the Atomic Mass
Evaluation (AME); the status flag marks whether a value is experimentally known
(`1`) or extrapolated (`0`).

## Validation

```bash
python validate.py
```

Compares computed half-lives against experimental values for a set of
well-measured cluster decays and prints residuals. 

## Repository layout

```
cluster_decays.py              Core physics module
cluster_decay_systematics.py   Full systematics table
validate.py                    Validation against experiment
atomic_mass_data.txt           Atomic mass input data (not included, see below)
README.md                      This file
LICENSE                        MIT License
```

## License

MIT License. See `LICENSE`.

The underlying physical model and original MATLAB code are the intellectual
property of R. Blendowske, T. Fliessbach, and H. Walliser. This Python
translation is published with their explicit permission.

## Citation

R. Blendowske, H. Walliser,
Phys. Rev. Lett. 61, 1930 (1988).
https://doi.org/10.1103/PhysRevLett.61.1930

R. Blendowske, T. Fliessbach, H. Walliser,
Z. Phys. A 339, 121 (1991).
https://doi.org/10.1007/BF01282943


