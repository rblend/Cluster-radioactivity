"""
Cluster decay systematics – Python translation of cluster_decay_systematics.m

Prints a table of calculated half-lives and related quantities for a
comprehensive set of cluster-decay reactions.

Usage
-----
    python cluster_decay_systematics.py

Requires cluster_decays.py (and atomic_mass_data.txt for automatic Q-values).
"""

from cluster_decays import cluster_decays


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def cluster_out(A, Z, AC, ZC, L=0, odd=0, ENC=None):
    """Calculate and print one row of the systematics table."""
    if ENC is not None:
        TAUlog, E, PEN, S, WKB, RI, RA, ES = cluster_decays(
            A, Z, AC, ZC, odd, L, ENC)
    else:
        TAUlog, E, PEN, S, WKB, RI, RA, ES = cluster_decays(
            A, Z, AC, ZC, odd, L)

    print(
        f"{A:5g} {Z:5g} {AC:5g} {ZC:5g}   "
        f"{TAUlog:9.3f} {E:9.3f}   "
        f"{PEN:12.4g} {S:12.4g} {WKB:12.4e}  "
        f"{RI:9.3f} {RA:9.3f}   "
        f"{L:3g}  {sum(ES):3g}  {odd:10g}"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    header = (
        "   A     Z     AC     ZC     log(t)      E          PEN            S"
        "       Lambda-Gamow     RI        RA        L    E:status   odd"
    )
    print(header)
    print()

    # --- Ba parent ---
    cluster_out(114, 56, 12, 6)
    print()

    # --- Fr parent ---
    cluster_out(221, 87, 14, 6, odd=0)
    cluster_out(221, 87, 14, 6, odd=1)
    print()

    # --- Ra parent ---
    cluster_out(221, 88, 14, 6, odd=0)
    cluster_out(222, 88, 14, 6, odd=0)
    cluster_out(223, 88, 14, 6, odd=0)
    cluster_out(223, 88, 14, 6, odd=1)          # odd
    cluster_out(223, 88, 14, 6, odd=0, ENC=31)  # favoured, excited-state transition
    cluster_out(224, 88, 14, 6, odd=0)
    cluster_out(226, 88, 14, 6, odd=0)
    print()

    # --- Ac parent ---
    cluster_out(223, 89, 14, 6,  odd=0)
    cluster_out(223, 89, 15, 7,  odd=0)
    cluster_out(223, 89, 15, 7,  odd=1)         # odd
    cluster_out(225, 89, 14, 6,  odd=0)
    print()

    # --- Th parent ---
    cluster_out(228, 90, 20,  8, odd=0)
    cluster_out(230, 90, 24, 10, odd=0)
    cluster_out(232, 90, 24, 10, odd=0)
    cluster_out(232, 90, 26, 10, odd=0)
    print()

    # --- Pa parent ---
    cluster_out(231, 91, 24, 10, odd=0)
    cluster_out(231, 91, 24, 10, odd=1)         # odd
    print()

    # --- U parent (Ne clusters) ---
    cluster_out(230, 92, 22, 10, odd=0)
    cluster_out(232, 92, 24, 10, odd=0)
    cluster_out(233, 92, 24, 10, odd=0)
    cluster_out(233, 92, 24, 10, odd=1)         # odd
    cluster_out(233, 92, 25, 10, odd=0)
    cluster_out(233, 92, 25, 10, odd=1)         # odd
    cluster_out(234, 92, 24, 10, odd=0)
    cluster_out(234, 92, 26, 10, odd=0)
    cluster_out(235, 92, 24, 10, odd=0)
    cluster_out(235, 92, 25, 10, odd=0)
    cluster_out(236, 92, 24, 10, odd=0)
    cluster_out(236, 92, 26, 10, odd=0)
    print()

    # --- U parent (Mg clusters) ---
    cluster_out(232, 92, 28, 12, odd=0)
    cluster_out(233, 92, 28, 12, odd=0)
    cluster_out(233, 92, 28, 12, odd=1)         # odd
    cluster_out(234, 92, 28, 12, odd=0)
    cluster_out(235, 92, 28, 12, odd=0)
    cluster_out(235, 92, 28, 12, odd=1)         # odd
    cluster_out(235, 92, 29, 12, odd=0)
    cluster_out(235, 92, 29, 12, odd=1)         # odd
    cluster_out(236, 92, 28, 12, odd=0)
    cluster_out(236, 92, 28, 12, odd=0, ENC=71.69)
    cluster_out(236, 92, 30, 12, odd=0)
    print()

    # --- Np parent ---
    cluster_out(237, 93, 30, 12, odd=0)
    cluster_out(237, 93, 30, 12, odd=1)         # odd
    print()

    # --- Pu parent ---
    cluster_out(236, 94, 28, 12, odd=0)
    cluster_out(238, 94, 28, 12, odd=0)
    cluster_out(238, 94, 30, 12, odd=0)
    cluster_out(238, 94, 32, 14, odd=0)
    cluster_out(240, 94, 32, 14, odd=0)
    print()

    # --- Am parent ---
    cluster_out(241, 95, 34, 14, odd=0)
    cluster_out(241, 95, 34, 14, odd=1)         # odd
    print()

    # --- Cm parent ---
    cluster_out(242, 96, 34, 14, odd=0)


if __name__ == "__main__":
    main()
