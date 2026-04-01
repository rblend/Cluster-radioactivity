"""
validate.py  –  Comparison of computed half-lives against experimental values
               for a representative set of well-measured cluster decays.

Experimental data compiled from:
  - B.A. Brown et al., Phys. Rev. C 45 (1992) R1
  - A. Sandulescu et al., Phys. Rev. C 54 (1996) 258
  - NUBASE / ENSDF evaluated nuclear data

Run:
    python validate.py

All cases use experimental Q-values (ENC) so the test is purely on the
WKB tunnelling model, independent of the mass table.
"""

from cluster_decays import cluster_decays

# ---------------------------------------------------------------------------
# Experimental reference data
# Each entry: (AP, ZP, AC, ZC, odd, L, Q_MeV, log10_T_exp, label)
# ---------------------------------------------------------------------------
EXPERIMENTS = [
    # Ra isotopes -> C-14  (best-measured series)
    (221, 88, 14, 6, 0, 0, 32.394, 14.6,  "Ra-221 -> C-14"),
    (222, 88, 14, 6, 0, 0, 33.049, 11.0,  "Ra-222 -> C-14"),
    (223, 88, 14, 6, 1, 0, 31.829, 15.1,  "Ra-223 -> C-14  (odd)"),
    (224, 88, 14, 6, 0, 0, 30.535, 15.9,  "Ra-224 -> C-14"),
    (226, 88, 14, 6, 0, 0, 28.196, 21.3,  "Ra-226 -> C-14"),

    # Fr -> C-14
    (221, 87, 14, 6, 1, 0, 31.290, 14.6,  "Fr-221 -> C-14  (odd)"),

    # Ac -> C-14 / N-15
    (225, 89, 14, 6, 1, 0, 30.476, 17.3,  "Ac-225 -> C-14  (odd)"),
    (223, 89, 15, 7, 1, 0, 39.470, 14.7,  "Ac-223 -> N-15  (odd)"),

    # Th -> O-20 / Ne-24
    (228, 90, 20, 8, 0, 0, 44.723, 20.7,  "Th-228 -> O-20"),
    (230, 90, 24,10, 0, 0, 57.758, 24.6,  "Th-230 -> Ne-24"),
    (232, 90, 26,10, 0, 0, 55.972, 29.2,  "Th-232 -> Ne-26"),

    # U -> Ne / Mg
    (232, 92, 24,10, 0, 0, 62.309, 20.5,  "U-232  -> Ne-24"),
    (234, 92, 24,10, 0, 0, 58.826, 25.9,  "U-234  -> Ne-24"),
    (234, 92, 28,12, 0, 0, 74.108, 25.7,  "U-234  -> Mg-28"),
    (236, 92, 28,12, 0, 0, 70.560, 27.6,  "U-236  -> Mg-28"),

    # Pu -> Mg / Si
    (238, 94, 28,12, 0, 0, 75.910, 25.7,  "Pu-238 -> Mg-28"),
    (238, 94, 32,14, 0, 0, 91.188, 25.3,  "Pu-238 -> Si-32"),
    (240, 94, 32,14, 0, 0, 88.490, 27.4,  "Pu-240 -> Si-32"),
]

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

header = (f"{'Decay':<26}  {'log(T) calc':>11}  {'log(T) exp':>10}  "
          f"{'Δ log(T)':>9}  {'E calc':>8}")
print(header)
print("-" * len(header))

max_dev = 0.0
for AP, ZP, AC, ZC, odd, L, Q, T_exp, label in EXPERIMENTS:
    TAUlog, E, *_ = cluster_decays(AP, ZP, AC, ZC, odd, L, ENC=Q)
    delta = TAUlog - T_exp
    max_dev = max(max_dev, abs(delta))
    flag = "  <--" if abs(delta) > 3.0 else ""
    print(f"{label:<26}  {TAUlog:11.2f}  {T_exp:10.1f}  {delta:+9.2f}  {E:8.3f} MeV{flag}")

print("-" * len(header))
print(f"  Maximum deviation: {max_dev:.2f} in log10(T½)")
print()
print("Note: deviations > ~2 in log10 are expected for transitions to excited")
print("states or strongly hindered (odd-A) decays not individually tuned.")
