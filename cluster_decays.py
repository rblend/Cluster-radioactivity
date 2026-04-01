"""
CALCULATION OF DECAY CONSTANTS FOR CLUSTER DECAYS
R.Blendowske, H.Walliser  Phys.Rev.Lett. 61 (1988) 1930
R.Blendowske, T.Fliessbach, H.Walliser  Z.Phys. A339 (1991) 121

Python translation of the original MATLAB code.
"""

import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad


# ---------------------------------------------------------------------------
# Internal helper functions
# ---------------------------------------------------------------------------

def _fpent(x, REFF, R1, R2, HM, L, GA, E):
    """Potential function – its zeros are the classical turning points."""
    VN = -REFF * np.exp(-(x - R1 - R2) / 0.63)   # nuclear (Christensen & Winter)
    VZ = HM * L * (L + 1) / (x * x)               # centrifugal
    VC = GA / x                                     # Coulomb
    return (VN + VC + VZ - E) / HM


def _fpent2(x, REFF, R1, R2, HM, L, GA, E):
    """Integrand for the WKB barrier integral."""
    VN = -REFF * np.exp(-(x - R1 - R2) / 0.63)
    VZ = HM * L * (L + 1) / (x * x)
    VC = GA / x
    return np.sqrt(np.abs((VN + VC + VZ - E) / HM))


def _bracket_and_solve(f, x0, args, dx0=None, factor=1.05, max_steps=2000):
    """
    Find a root of f near x0 by scanning outward in both directions until a
    sign change is bracketed, then refine with brentq.
    Mimics MATLAB fzero's automatic bracketing strategy.
    """
    if dx0 is None:
        dx0 = abs(x0) * 0.005 if x0 != 0 else 0.01

    f0 = f(x0, *args)
    if f0 == 0.0:
        return x0

    # Scan symmetrically in both directions
    step = dx0
    for _ in range(max_steps):
        for sign in (+1, -1):
            xb = x0 + sign * step
            if xb <= 0:
                continue
            fb = f(xb, *args)
            if f0 * fb < 0:
                a, b = (x0, xb) if xb > x0 else (xb, x0)
                return brentq(f, a, b, args=args, xtol=1e-10, rtol=1e-10)
        step *= factor

    raise RuntimeError(
        f"_bracket_and_solve: no bracket found near x0={x0:.6g} "
        f"after {max_steps} steps."
    )


def _ecorr(E, Z1, Z2):
    """
    Electronic correction to the Q-value (MeV).
    D. Lunney, J.M. Pearson and C. Thibault, Rev. Mod. Phys. 75 (2003) 1021
    """
    Z  = Z1 + Z2
    EM = 14.4381 * Z**2.39  + 1.55468e-6 * Z**5.35
    E1 = 14.4381 * Z1**2.39 + 1.55468e-6 * Z1**5.35
    E2 = 14.4381 * Z2**2.39 + 1.55468e-6 * Z2**5.35
    return E + (EM - E1 - E2) / 1e6             # eV -> MeV


def _ebind(A, Z):
    """
    Binding energies for the emitted cluster (MeV).
    A.H. Wapstra et al. Nucl. Phys. A432 (1985) 1
    Falls back to a semi-empirical formula when not tabulated.
    """
    table = {
        ( 4,  2):  28.296,
        (10,  4):  64.977,
        (12,  6):  92.162,
        (14,  6): 105.285,
        (15,  7): 115.492,
        (16,  8): 127.620,
        (18,  8): 139.808,
        (20,  8): 151.372,
        (22,  8): 161.870,
        (22, 10): 177.773,
        (23,  9): 175.290,
        (24, 10): 191.839,
        (25, 10): 196.120,
        (26, 10): 201.590,
        (28, 12): 231.629,
        (29, 12): 235.300,
        (30, 12): 241.850,
        (32, 14): 271.412,
        (34, 12): 258.100,
        (34, 14): 283.330,
        (36, 16): 308.716,
    }
    if (A, Z) in table:
        return table[(A, Z)]

    print(f"  WARNING: Binding energy of cluster (A={A}, Z={Z}) not stored – "
          "using semi-empirical estimate.")
    return (14.1  * A
            - 13.0 * A**0.6666
            - 0.595 * Z * Z / A**0.3333
            - 19.0 * (A - 2 * Z)**2 / A
            - 33.5 / A**0.75)


def load_atomic_mass_tables(AP, ZP, AC, ZC, filename="atomic_mass_data.txt"):
    """
    Load Q-value from an atomic mass data file.

    File format (whitespace-separated): Z  A  Ebind/nucleon[keV]  status

    Returns
    -------
    Q : float        Q-value in MeV
    ES : list[int]   status flags [parent, daughter, cluster]
    """
    data = np.loadtxt(filename)
    Z_col      = data[:, 0].astype(int)
    A_col      = data[:, 1].astype(int)
    E_col      = data[:, 2]               # keV per nucleon
    status_col = data[:, 3].astype(int)

    AD = AP - AC
    ZD = ZP - ZC

    def _lookup(Z, A):
        idx = np.where((Z_col == Z) & (A_col == A))[0]
        if len(idx) == 0:
            raise ValueError(f"Nucleus Z={Z}, A={A} not found in mass table.")
        return A * E_col[idx[0]] / 1000.0, int(status_col[idx[0]])  # -> MeV

    EP, es_p = _lookup(ZP, AP)
    ED, es_d = _lookup(ZD, AD)
    EC, es_c = _lookup(ZC, AC)

    return ED + EC - EP, [es_p, es_d, es_c]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def cluster_decays(AP, ZP, AC, ZC, odd=0, L=0, ENC=None):
    """
    Calculate decay constants for cluster decays.

    Parameters
    ----------
    AP, ZP : int   Mass and atomic number of the parent nucleus.
    AC, ZC : int   Mass and atomic number of the emitted cluster.
    odd    : int   1 -> S = 3.2e-3 (odd nucleus);  0 (default) -> S = 6.3e-3.
    L      : int   Angular momentum (default 0).
    ENC    : float or None
                   Q-value in MeV.  Computed from mass tables when None.

    Returns
    -------
    TAUlog : float   log10 of the half-life [s]
    E      : float   Decay energy with electronic correction [MeV]
    PEN    : float   Penetrability
    S      : float   Preformation probability (spectroscopic factor)
    WKB    : float   Gamow decay constant [s^-1]
    RI     : float   Inner classical turning point [fm]
    RA     : float   Outer classical turning point [fm]
    ES     : list    Status flags [parent, daughter, cluster] (0 or 1 each)
    """

    AD = AP - AC
    ZD = ZP - ZC

    ES = [0, 0, 0]
    if ENC is None:
        ENC, ES = load_atomic_mass_tables(AP, ZP, AC, ZC)

    # Electronic correction
    E = _ecorr(ENC, ZD, ZC)

    # Reduced mass  HM = hbar^2/(2*mu)  [MeV*fm^2], with cluster binding correction
    XMASS = (AC - ZC) * 939.5732 + ZC * 938.2796 - _ebind(AC, ZC)
    HM = 0.5 * 197.32858**2 / XMASS
    HM = HM * (AD + AC) / AD

    # Coulomb parameter
    GA = 1.4399 * ZD * ZC

    # Radii (Phys. Lett. 65B (1976) 19)
    T    = 1.0 / 3.0
    R1   = 1.233 * AD**T - 0.978 / AD**T
    R2   = 1.233 * AC**T - 0.978 / AC**T
    REFF = 50.0 * R1 * R2 / (R1 + R2)

    args = (REFF, R1, R2, HM, L, GA, E)

    # Outer turning point: pure Coulomb gives GA/E, search bidirectionally
    RA = _bracket_and_solve(_fpent, GA / E, args)

    # Inner turning point: near nuclear surface, search bidirectionally
    RI = _bracket_and_solve(_fpent, R1 + R2 + 1.0, args)

    # Ensure RI < RA
    if RI > RA:
        RI, RA = RA, RI

    # WKB barrier integral
    SUM, _ = quad(_fpent2, RI, RA, args=args,
                  limit=200, epsabs=1e-7, epsrel=1e-7)

    EU  = np.log10(np.e)
    EXX = -2.0 * EU * SUM
    PEN = 10.0**EXX

    # Knocking frequency (25 MeV/nucleon)
    VIN = 25.0 * AC
    WKB = np.sqrt(HM * VIN) * 10.0**(EXX + 22) / (RI * 6.582173)

    # Spectroscopic factor
    S = 6.3e-3 if odd == 0 else 3.2e-3
    S = S**((AC - 1) / 3.0)

    # Half-life
    XLA    = WKB * S
    TAUlog = np.log10(np.log(2.0) / XLA)

    return TAUlog, E, PEN, S, WKB, RI, RA, ES
