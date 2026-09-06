#!/usr/bin/env python3
"""Generate bbn/rates_reaclib.dat from the JINA REACLIB database.

The BBN solver (tools/bbn_rates.cpp) reads FORWARD reaction rates only, in the
REACLIB seven-parameter form. Reverse rates are derived in C++ by detailed
balance and have no representation in the data file, so this script must never
emit a REACLIB "reverse" entry.

Provenance is the point of this script: the data file it writes is generated,
not hand-transcribed, and the source label of every rate is carried into the
file so a reader can trace each number back to its publication.

Usage:
    pip install pynucastro
    python3 scripts/gen_bbn_rates.py [-o bbn/rates_reaclib.dat]
"""

import argparse
import math
import pathlib
import sys

# The Kawano key network: nine nuclides, fourteen strong reactions. Free neutron
# decay and the two-directional weak n <-> p rates are NOT here -- they are
# computed from phase-space integrals in tools/bbn_weak.cpp, normalized to the
# tau_n input, and so are not rate-file data.
REACTIONS = [
    "n(p,g)d",
    "d(p,g)he3",
    "d(d,n)he3",
    "d(d,p)t",
    "he3(n,p)t",
    "t(d,n)he4",
    "he3(d,p)he4",
    "he3(he4,g)be7",
    "t(he4,g)li7",
    "be7(n,p)li7",
    "li7(p,he4)he4",
    "he4(d,g)li6",
    "li6(p,he4)he3",
    "t(p,g)he4",
]

# Expanded source labels, so the emitted file is readable without a REACLIB manual.
LABEL_SOURCES = {
    "an06": "Ando et al. 2006, PRC 74, 025809",
    "de04": "Descouvemont et al. 2004, ADNDT 88, 203",
    "gi17": "Gomez Inesta et al. 2017, ApJ 849, 134",
    "go17": "Gomez Inesta et al. 2017, ApJ 849, 134",
    "cd08": "Cyburt & Davids 2008, PRC 78, 064614",
    "db18": "Damone et al. 2018 (n_TOF), PRL 121, 042701",
    "tu19": "Tumino et al. 2019",
    "pt05": "Pizzone et al. 2005, A&A 438, 779",
    "cf88": "Caughlan & Fowler 1988, ADNDT 40, 283",
}


def reaclib_eval(sets, t9):
    """Evaluate the seven-parameter REACLIB sum -- the same expression that
    tools/bbn_rates.cpp implements. Used below to verify that our reading of the
    coefficient ordering agrees with pynucastro's own evaluator."""
    total = 0.0
    for a in sets:
        total += math.exp(
            a[0]
            + a[1] / t9
            + a[2] / t9 ** (1.0 / 3.0)
            + a[3] * t9 ** (1.0 / 3.0)
            + a[4] * t9
            + a[5] * t9 ** (5.0 / 3.0)
            + a[6] * math.log(t9)
        )
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--output", type=pathlib.Path,
                    default=pathlib.Path(__file__).parent.parent / "bbn" / "rates_reaclib.dat")
    args = ap.parse_args()

    try:
        import pynucastro as pyna
    except ImportError:
        sys.exit("pynucastro is required: pip install pynucastro")

    lib = pyna.ReacLibLibrary()
    rows, labels_seen = [], set()

    for name in REACTIONS:
        rate = lib.get_rate_by_name(name)
        sets = [list(s.a) for s in rate.sets]

        # Guard: the coefficient ordering we write must reproduce pynucastro's
        # evaluator across the whole BBN temperature range. A silent transposition
        # here would be a wrong rate that still looks like a plausible number.
        for t9 in (0.01, 0.1, 1.0, 10.0):
            ours, theirs = reaclib_eval(sets, t9), rate.eval(t9 * 1e9)
            if not math.isclose(ours, theirs, rel_tol=1e-10):
                sys.exit(f"{name}: fit form disagrees at T9={t9}: {ours} vs {theirs}")

        for s in rate.sets:
            label = s.labelprops.strip()
            labels_seen.add(label.rstrip("nr"))
            rows.append((name, label, list(s.a)))

    with open(args.output, "w") as f:
        f.write("# CLASS++ BBN forward reaction rates, JINA REACLIB seven-parameter form.\n")
        f.write("#\n")
        f.write("#   lambda = sum_sets exp( a0 + a1/T9 + a2/T9^(1/3) + a3*T9^(1/3)\n")
        f.write("#                            + a4*T9 + a5*T9^(5/3) + a6*ln(T9) )\n")
        f.write("#\n")
        f.write("# Units: two-body reactions  N_A<sigma v>  [cm^3 mol^-1 s^-1]\n")
        f.write("#        one-body reactions  lambda        [s^-1]\n")
        f.write("#\n")
        f.write("# REVERSE RATES ARE NOT IN THIS FILE and have no syntax here. They are\n")
        f.write("# derived by detailed balance in tools/bbn_rates.cpp from the nuclide\n")
        f.write("# degeneracies and masses in tools/bbn_nuclides.h, so that editing this\n")
        f.write("# file cannot break thermodynamic consistency.\n")
        f.write("#\n")
        f.write("# Lines with the same reaction name are summed. A reaction the code knows\n")
        f.write("# about but which is missing here is a fatal error, as is a reaction named\n")
        f.write("# here that the code does not know about.\n")
        f.write("#\n")
        f.write("# Weak n <-> p rates and free neutron decay are deliberately absent: they\n")
        f.write("# come from phase-space integrals normalized to the tau_n input.\n")
        f.write("#\n")
        f.write("# GENERATED by scripts/gen_bbn_rates.py -- do not hand-edit; regenerate.\n")
        f.write(f"# Source: JINA REACLIB via pynucastro {pyna.__version__}\n")
        f.write("#\n")
        f.write("# Rate sources by label:\n")
        for lab in sorted(labels_seen):
            f.write(f"#   {lab:6s} {LABEL_SOURCES.get(lab, 'see JINA REACLIB')}\n")
        f.write("#\n")
        f.write(f"# {'reaction':<16s} {'label':<7s}"
                + "".join(f"{'a'+str(i):>16s}" for i in range(7)) + "\n")
        for name, label, a in rows:
            f.write(f"  {name:<16s} {label:<7s}"
                    + "".join(f"{v:16.8e}" for v in a) + "\n")

    print(f"wrote {args.output} ({len(rows)} sets over {len(REACTIONS)} reactions)")


if __name__ == "__main__":
    main()
