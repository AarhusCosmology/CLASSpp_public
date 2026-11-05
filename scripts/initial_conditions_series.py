"""Superhorizon initial conditions as power series in tau, checked term by term against CLASS++.

Solves the synchronous-gauge perturbation equations for one regular superhorizon mode as power series in
conformal time, and compares every coefficient with the initial conditions that perturb_initial_conditions
(source/perturbations_module.cpp) and the species' ApplyInitialConditions set. Written in 2026-11 to recheck
the non-adiabatic modes (the niv shear, l=3 and eta values were wrong, and the adiabatic l=3 value was half
of what the shear implies).

Equations (Ma & Bertschinger 1995 in CLASS's variables; flat; photons and baryons tightly coupled at zeroth
order; massless neutrinos with the hierarchy truncated at LMAX; radiation + matter background):

    a ~ tau (1 + om tau / 4),   y = rho_m / rho_r = om tau (1 + om tau / 4),
    aH = (1 + om tau / 2) / (tau (1 + om tau / 4)),   om = a rho_m / sqrt(rho_r)   (8 pi G / 3 = 1),
    f_nu = rho_nu / rho_r,  f_g = 1 - f_nu,  f_b = rho_b / rho_m,  f_c = 1 - f_b,  R = 3/4 f_b y / f_g,

    aH h'      = 2 k^2 eta + 3 (aH)^2 delta_tot                           (00 Einstein equation)
    eta'       = 3/2 (aH)^2 / k^2 (rho + p) theta_tot / rho_tot           (0i Einstein equation)
    delta_c'   = -h'/2
    delta_b'   = -theta_g - h'/2                                          (theta_b = theta_g)
    delta_g'   = -4/3 theta_g - 2/3 h'
    (1 + R) theta_g' + R aH theta_g = k^2 delta_g / 4
    delta_nu'  = -4/3 theta_nu - 2/3 h'
    theta_nu'  = k^2 (delta_nu / 4 - sigma_nu)
    sigma_nu'  = 1/2 (8/15 theta_nu - 3/5 k F_3 + 4/15 h' + 8/5 eta')
    F_l'       = k / (2l + 1) (l F_{l-1} - (l + 1) F_{l+1}),   F_2 = 2 sigma_nu,  F_{LMAX+1} = 0.

Each variable is sum_j c_j tau^j. The equations, multiplied by their denominators, give linear equations
for the c_j order by order. The mode is fixed by its tau^0 data:

    ad: eta = S,  cdi: delta_c = S,  bi: delta_b = S,  nid: delta_nu = S,  niv: theta_nu = 3/4 k S

and regularity (sigma_nu and F_l start at tau^1 and tau^(l-1)); the other tau^0 values follow, e.g.
theta_g = -3/4 k S f_nu / f_g for niv. With k = 1, a term om^n tau^j is (k tau)^(j-n) (om tau)^n. Truncating
at order N leaves the top coefficients free; a coefficient is reported as exact while it does not depend
on them.

Output: per variable, every exact coefficient; CLASS++'s value is printed next to it where it differs,
marked "not included" when CLASS++ has no such term at a higher order (the series is truncated) and "WRONG"
when CLASS++ has a different value or misses a term of the lowest order. The exit status is the number of WRONG coefficients.

usage: python scripts/initial_conditions_series.py MODE [N [F_NU F_B]]
       MODE = ad | cdi | bi | nid | niv,  N = order in tau (default 4),
       F_NU, F_B rationals such as 41/100 4/25 (default: symbolic, much slower).   Needs sympy.
"""
import sys

import sympy as sp

LMAX = 8

MODE = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 4
tau, om, S = sp.symbols("tau om S", positive=True)
if len(sys.argv) > 4:
    fn, fb = sp.Rational(sys.argv[3]), sp.Rational(sys.argv[4])
else:
    fn, fb = sp.symbols("f_nu f_b", positive=True)
fg, fc = 1 - fn, 1 - fb
k = 1

# --- the series -------------------------------------------------------------------------------------------
names = ["dc", "db", "dg", "thg", "dn", "thn", "sn"] + [f"F{l}" for l in range(3, LMAX + 1)] + ["eta"]
c = {v: sp.symbols(f"{v}_0:{N + 1}") for v in names}
Y = {v: sum(c[v][j] * tau**j for j in range(N + 1)) for v in names}
D = {v: sp.diff(Y[v], tau) for v in names}

q = 1 + om * tau / 4
aH = (1 + om * tau / 2) / (tau * q)
y = om * tau * q
R = sp.Rational(3, 4) * fb * y / fg
delta_tot = (fg * Y["dg"] + fn * Y["dn"] + y * (fb * Y["db"] + fc * Y["dc"])) / (1 + y)
theta_tot = (sp.Rational(4, 3) * (fg * Y["thg"] + fn * Y["thn"]) + y * fb * Y["thg"]) / (1 + y)
hp = 2 * k**2 * Y["eta"] / aH + 3 * aH * delta_tot


def F(l):
    return 2 * Y["sn"] if l == 2 else (Y[f"F{l}"] if 3 <= l <= LMAX else 0)


eqs = {
    "eta": D["eta"] - sp.Rational(3, 2) * aH**2 / k**2 * theta_tot,
    "dc": D["dc"] + hp / 2,
    "db": D["db"] + Y["thg"] + hp / 2,
    "dg": D["dg"] + sp.Rational(4, 3) * Y["thg"] + sp.Rational(2, 3) * hp,
    "thg": (1 + R) * D["thg"] + R * aH * Y["thg"] - k**2 * Y["dg"] / 4,
    "dn": D["dn"] + sp.Rational(4, 3) * Y["thn"] + sp.Rational(2, 3) * hp,
    "thn": D["thn"] - k**2 * (Y["dn"] / 4 - Y["sn"]),
    "sn": D["sn"] - (sp.Rational(8, 15) * Y["thn"] - sp.Rational(3, 5) * k * F(3) + sp.Rational(4, 15) * hp
                     + sp.Rational(8, 5) * D["eta"]) / 2,
}
for l in range(3, LMAX + 1):
    eqs[f"F{l}"] = D[f"F{l}"] - sp.Rational(k, 2 * l + 1) * (l * F(l - 1) - (l + 1) * F(l + 1))

lin = []
for v, e in eqs.items():
    num, den = sp.fraction(sp.together(e * tau**2 * q**2 * (1 + om * tau / 2)**2 * (1 + y)))
    assert not den.has(tau), (v, den)
    num = sp.Poly(sp.expand(num), tau)
    lin += [co for j in range(N + 1) if (co := num.coeff_monomial(tau**j)) != 0]
# regularity: sigma_nu ~ tau, F_l ~ tau^(l-1)
lin += [c["sn"][0]] + [c[f"F{l}"][j] for l in range(3, LMAX + 1) for j in range(min(l - 1, N + 1))]
FIX = {"ad": dict(eta=S, dc=0, db=0, dn=0, thn=0),
       "cdi": dict(eta=0, dc=S, db=0, dn=0, thn=0),
       "bi": dict(eta=0, dc=0, db=S, dn=0, thn=0),
       "nid": dict(eta=0, dc=0, db=0, dn=S, thn=0),
       "niv": dict(eta=0, dc=0, db=0, dn=0, thn=sp.Rational(3, 4) * S * k)}[MODE]
lin += [c[v][0] - val for v, val in FIX.items()]
unknowns = [s for v in names for s in c[v]]
sol = sp.solve(lin, unknowns, dict=True, rational=True)
assert len(sol) == 1, "the mode is not unique"
sol = sol[0]
free = [s for s in unknowns if s not in sol]

# --- CLASS++'s initial conditions (flat: s2_squared = s_l = 1; curvature_ini = entropy_ini = S) -----------
kt, kt2, kt3 = k * tau, (k * tau)**2, (k * tau)**3
CL = {}
if MODE == "ad":
    CL = dict(
        dg=-kt2 / 3 * (1 - om * tau / 5) * S,
        thg=-k * kt3 / 36 * (1 - 3 * (1 + 5 * fb - fn) / 20 / (1 - fn) * om * tau) * S,
        dn=-kt2 / 3 * (1 - om * tau / 5) * S,
        thn=-k * kt3 / 36 / (4 * fn + 15) * (4 * fn + 23 - 3 * (8 * fn**2 + 50 * fn + 275) / 20 / (2 * fn + 15)
                                             * tau * om) * S,
        sn=2 * kt2 / (45 + 12 * fn) * (1 + (4 * fn - 5) / 4 / (2 * fn + 15) * tau * om) * S,
        F3=4 * kt3 / 7 / (12 * fn + 45) * S,
        eta=(1 - kt2 / 12 / (15 + 4 * fn) * (5 + 4 * fn - (16 * fn**2 + 280 * fn + 325) / 10 / (2 * fn + 15)
                                              * tau * om)) * S)
    CL["db"] = CL["dc"] = sp.Rational(3, 4) * CL["dg"]
elif MODE in ("cdi", "bi"):
    f = fc if MODE == "cdi" else fb
    dg = S * f * om * tau * (-sp.Rational(2, 3) + om * tau / 4)
    CL = dict(dg=dg, thg=-S * f * om * kt2 / 12, dn=dg, thn=-S * f * om * kt2 / 12,
              sn=-S * f * kt2 * tau * om / 6 / (2 * fn + 15), eta=-S * f * om * tau * (sp.Rational(1, 6) - om * tau / 16))
    CL["dc"] = (S if MODE == "cdi" else 0) + sp.Rational(3, 4) * dg
    CL["db"] = (S if MODE == "bi" else 0) + sp.Rational(3, 4) * dg
elif MODE == "nid":
    CL = dict(
        dg=S * fn / fg * (-1 + kt2 / 6),
        thg=-S * fn / fg * k * kt * (sp.Rational(1, 4) - fb / fg * 3 * om * tau / 16),
        db=S * fn / fg / 8 * kt2,
        dc=-S * fn * fb / fg / 80 * kt2 * om * tau,
        dn=S * (1 - kt2 / 6), thn=S * k * kt / 4, sn=S * kt2 / (4 * fn + 15) / 2,
        F3=S * kt3 / 7 / (4 * fn + 15),
        eta=-S * fn / (4 * fn + 15) / 6 * kt2)
elif MODE == "niv":
    CL = dict(
        dg=S * kt * fn / fg * (1 - sp.Rational(3, 16) * fb * (2 + fg) / fg * om * tau),
        thg=S * fn / fg * sp.Rational(3, 4) * k * (-1 + sp.Rational(3, 4) * fb / fg * om * tau
                                                   + sp.Rational(3, 16) * (om * tau)**2 * fb / fg**2 * (fg - 3 * fb)
                                                   + kt2 / 6),
        dc=-S * sp.Rational(9, 64) * fn * fb / fg * kt * om * tau,
        dn=-S * kt * (1 + sp.Rational(3, 16) * fb * fn / fg * om * tau),
        thn=S * sp.Rational(3, 4) * k * (1 - kt2 / 6 * (4 * fn + 9) / (4 * fn + 5)),
        sn=S / (4 * fn + 5) * kt * (1 + 3 * om * tau * fn / (4 * fn + 15)),
        F3=S * sp.Rational(3, 7) / (4 * fn + 5) * kt2,
        eta=S * fn * kt * (-1 / (4 * fn + 5) + (-sp.Rational(3, 64) * fb / fg
                                                + sp.Rational(15, 4) / (4 * fn + 15) / (4 * fn + 5)) * om * tau))
    CL["db"] = sp.Rational(3, 4) * CL["dg"]


# --- compare ----------------------------------------------------------------------------------------------
def by_order(e):
    """{(a, b): coefficient of (k tau)^a (om tau)^b}"""
    e, out = sp.expand(e), {}
    for j in range(N + 1):
        cj = sp.expand(e.coeff(tau, j))
        for n in range(j + 1):
            cn = sp.factor(cj.coeff(om, n))
            if cn != 0:
                out[(j - n, n)] = cn
    return out


CLASS_NAME = dict(dc="delta_cdm", db="delta_b", dg="delta_g", thg="theta_g", dn="delta_ur", thn="theta_ur",
                  sn="shear_ur", F3="l3_ur", eta="eta")
print(f"# {MODE}, series to tau^{N}, f_nu = {fn}, f_b = {fb}, l_max = {LMAX}")
wrong = 0
for v in CLASS_NAME:
    exact = {}
    for j, cj in enumerate(c[v]):
        val = sp.expand(sol.get(cj, cj))
        if any(val.has(f) for f in free):
            break
        exact[j] = val
    jmax = max(exact) if exact else -1
    ex = by_order(sum(val * tau**j for j, val in exact.items()))
    cl = {key: val for key, val in by_order(CL.get(v, 0)).items() if sum(key) <= jmax}
    print(f"\n{CLASS_NAME[v]}  (exact through tau^{jmax})")
    lowest = min((sum(key) for key in ex), default=None)  # a term missing at the lowest order is an error too
    for key in sorted(set(ex) | set(cl), key=lambda kk: (sum(kk), kk[1])):
        a, b = ex.get(key, 0), cl.get(key, 0)
        flag = ""
        if sp.simplify(a - b) != 0:
            bad = b != 0 or sum(key) == lowest
            flag = f"   <-- CLASS++: {b}  WRONG" if bad else "   <-- CLASS++: not included"
            wrong += bad
        print(f"   (k tau)^{key[0]} (om tau)^{key[1]}:  {a}{flag}")
print(f"\n{wrong} wrong coefficient(s)")
sys.exit(wrong)
