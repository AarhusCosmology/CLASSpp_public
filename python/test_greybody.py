"""Grey-body NCDM: FD limit and moments sanity via classy."""
from classy import Class

COMMON = {
    'output': 'mPk',
    'P_k_max_1/Mpc': 1.0,
    'N_ur': 2.0308,
    'h': 0.67556,
    'omega_b': 0.022032,
    'omega_cdm': 0.12038,
}


def _run(extra):
    cosmo = Class()
    cosmo.set({**COMMON, **extra})
    cosmo.compute()
    return cosmo


def test_greybody_fd_limit():
    # Standard NCDM neutrino.
    std = _run({
        'nu.type': 'ncdm_standard', 'nu.m': 0.06, 'nu.deg': 1.0,
        'nu.momenta_bins': 7, 'nu.quadrature_strategy': 2, 'nu.max_q': 18.0,
    })
    # Grey-body in the FD limit. At alpha=1, x=1, q0=1 the grey-body PSD reduces
    # to 1/(e^q + 1) (times the same 2/(2pi)^3 normalization the FD path carries
    # via the ksi=0 particle+antiparticle sum), so with identical m and deg the
    # two species are physically identical -> SAME deg (1.0).
    gb = _run({
        'gb.type': 'ncdm_greybody', 'gb.parameterization': 'direct',
        'gb.alpha': 1.0, 'gb.x': 1.0, 'gb.q0': 1.0,
        'gb.m': 0.06, 'gb.deg': 1.0,
        'gb.momenta_bins': 7, 'gb.quadrature_strategy': 2, 'gb.max_q': 18.0,
    })
    om_std = std.get_current_derived_parameters(['Omega_m'])['Omega_m']
    om_gb = gb.get_current_derived_parameters(['Omega_m'])['Omega_m']
    assert abs(om_gb - om_std) / om_std < 1e-3, (om_gb, om_std)
    std.struct_cleanup()
    gb.struct_cleanup()


def test_greybody_moments_sane():
    # Direct params; read back the moments derived parameters; confirm finite,
    # positive, and that M2*M4/M3^2 > 1 (Cauchy-Schwarz).
    gb = _run({
        'gb.type': 'ncdm_greybody', 'gb.parameterization': 'direct',
        'gb.alpha': 2.5, 'gb.x': 0.8, 'gb.q0': 1.3, 'gb.m': 0.1,
        'gb.momenta_bins': 15,
    })
    d = gb.get_current_derived_parameters(
        ['gb.M_2', 'gb.M_3', 'gb.M_4', 'gb.alpha', 'gb.q0', 'gb.x'])
    M2, M3, M4 = d['gb.M_2'], d['gb.M_3'], d['gb.M_4']
    assert M2 > 0 and M3 > 0 and M4 > 0, d
    assert M2 * M4 / (M3 * M3) > 1.0, d
    # Direct params must round-trip exactly through the per-species GetParam path.
    assert abs(d['gb.alpha'] - 2.5) < 1e-10, d
    assert abs(d['gb.q0'] - 1.3) < 1e-10, d
    assert abs(d['gb.x'] - 0.8) < 1e-10, d
    gb.struct_cleanup()


if __name__ == '__main__':
    test_greybody_fd_limit()
    test_greybody_moments_sane()
    print('grey-body classy tests passed')
