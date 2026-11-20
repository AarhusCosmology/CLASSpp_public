"""Class: a cosmology computed by CLASS++, as Python sees it.

Class subclasses the compiled core, so the core's methods (the per-redshift ones that Cobaya
calls thousands of times per point among them) are called directly, with no Python layer in
between. What is defined here is composition over those methods; the numbers themselves come
from C++ (docs/superpowers/specs/2026-11-20-classy-package-design.md).
"""
import numpy as np

from ._core import PyCosmology
from ._errors import CosmoSevereError


class Class(PyCosmology):
    """A cosmology: set its parameters with set(), then ask it for anything.

    The modules are computed on first use, so a question about the background never pays for
    the perturbations.
    """

    # ---------------------------------------------------------- derived parameters --

    def get_current_derived_parameters(self, names):
        """
        get_current_derived_parameters(names)

        Return a dictionary containing an entry for all the names defined in the
        input list.

        Parameters
        ----------
        names : list
                Derived parameters that can be asked from Monte Python, or
                elsewhere. The modules answer their own names (h, z_rec,
                sigma8, ...), and a species answers the names it knows, such
                as a_c of a pheno-axion fluid. "<species key>.<field>" asks one
                species instance alone: "Fluid.f_peak" is the peak EDE fraction
                f_EDE of a pheno-axion fluid, "Fluid.z_peak" the redshift where
                it peaks, "NEDE.z_star" the NEDE transition (a plain "z_star" is
                recombination's). Two species answering one plain name is a
                CosmoSevereError; ask for one of them with the dotted name.

        Returns
        -------
        derived : dict

        .. warning::

            This method used to take as an argument directly the data class from
            Monte Python. To maintain compatibility with this old feature, a
            check is performed to verify that names is indeed a list. If not, it
            returns a TypeError. The old version of this function, when asked
            with the new argument, will raise an AttributeError.

        """
        if not isinstance(names, list):
            raise TypeError("Deprecated")
        return {name: self._derived_parameter(name) for name in names}

    # One derived parameter each.
    def age(self):
        """Age of the universe today, in Gyr."""
        return self._derived_parameter('age')

    def h(self):
        """Reduced Hubble parameter today."""
        return self._derived_parameter('h')

    def n_s(self):
        """Scalar spectral index at the pivot scale."""
        return self._derived_parameter('n_s')

    def A_s(self):
        """Scalar amplitude at the pivot scale."""
        return self._derived_parameter('A_s')

    def tau_reio(self):
        """Optical depth to reionization."""
        return self._derived_parameter('tau_reio')

    def Omega_m(self):
        """Total non-relativistic matter today."""
        return self._derived_parameter('Omega_m')

    def Omega_r(self):
        """Total ultra-relativistic radiation today."""
        return self._derived_parameter('Omega_r')

    def theta_s_100(self):
        """100 times the angular sound horizon at recombination."""
        return self._derived_parameter('100*theta_s')

    def theta_star_100(self):
        """100 times the angular sound horizon at z_star."""
        return self._derived_parameter('100*theta_star')

    def theta_d_100(self):
        """100 times the angular photon damping scale at recombination; needs
        'compute damping scale': 'yes'."""
        return self._derived_parameter('100*theta_d')

    def Omega_Lambda(self):
        """Omega of the cosmological constant today (0 without one)."""
        return self._derived_parameter('Omega_Lambda')

    def Omega_g(self):
        """Omega of the photons today."""
        return self._derived_parameter('Omega_g')

    def Omega_b(self):
        """Omega of the baryons today."""
        return self._derived_parameter('Omega_b')

    def omega_b(self):
        """Omega_b h^2."""
        return self._derived_parameter('omega_b')

    def Neff(self):
        """Effective number of neutrino species."""
        return self._derived_parameter('Neff')

    def f_idr_drmd(self):
        """Dark-radiation fraction of the DRMD model (0 without it)."""
        return self._derived_parameter('f_idr_drmd')

    def z_dec_drmd(self):
        """Decoupling redshift of the DRMD model (-1 without it or without decoupling)."""
        return self._derived_parameter('z_dec_drmd')

    def k_eq(self):
        """Wavenumber of matter-radiation equality, a_eq H_eq, in 1/Mpc."""
        return self._derived_parameter('k_eq')

    def sigma8(self):
        """sigma8 of all matter."""
        return self._derived_parameter('sigma8')

    def sigma8_cb(self):
        """sigma8 of cdm and baryons."""
        return self._derived_parameter('sigma8_cb')

    def rs_drag(self):
        """Comoving sound horizon at baryon drag, in Mpc."""
        return self._derived_parameter('rs_d')

    def z_reio(self):
        """Redshift of reionization."""
        return self._derived_parameter('z_reio')

    def T_cmb(self):
        """CMB temperature today, in K."""
        return self._derived_parameter('T_cmb')

    def Omega0_m(self):
        """The sum of Omega0 of all non-relativistic components (as Omega_m)."""
        return self._derived_parameter('Omega_m')

    def Omega0_k(self):
        """Curvature contribution today."""
        return self._derived_parameter('Omega0_k')

    def Omega0_cdm(self):
        """Omega of cold dark matter today (0 without it)."""
        return self._derived_parameter('Omega0_cdm')

    @property
    def Omega_nu(self):
        """Omega of all species of the NCDM family (massive neutrinos and the like) today."""
        return self._derived_parameter('Omega_nu')

    # ------------------------------------------------------------------ the grids --

    def get_pk_and_k_and_z(self, nonlinear=True, only_clustering_species=False, h_units=False):
        """
        Return the stored grid of matter power spectrum values with the k and z it
        was computed on. Useful for building interpolators.

        Parameters
        ----------
        nonlinear : bool
                Return the non-linear spectrum (default) or the linear one
        only_clustering_species : bool
                Exclude massive neutrinos, as galaxy clustering does, rather than
                including everything (default)
        h_units : bool
                Return k in h/Mpc rather than 1/Mpc (default)

        Returns
        -------
        pk : grid of power spectrum values, pk[index_k, index_z], in Mpc**3
        k : vector of k values, k[index_k]
        z : vector of z values, z[index_z]
        """
        pk, k = self._pk_grid(nonlinear, only_clustering_species)
        z = self._z_grid_for_output()
        # The grid comes back in its stored order, z-major.
        pk = pk.reshape(z.size, k.size).T.copy()
        if h_units:
            k = k/self.h()
        return pk, k, z

    def get_transfer_and_k_and_z(self, output_format='class', h_units=False):
        """
        Return a dictionary of density and/or velocity transfer function grids with
        the k and z they were computed on. Useful for building interpolators.

        Include at least one of 'dTk' or 'vTk' in 'output'. With the default
        output_format='class' the transfer functions are normalised to curvature
        R=1 at initial time, and the newtonian metric fluctuations phi and psi are
        included; output_format='camb' switches to the CAMB convention, which
        outputs density transfer functions only.

        Parameters
        ----------
        output_format : 'class' or 'camb'
        h_units : bool
                Return k in h/Mpc rather than 1/Mpc (default)

        Returns
        -------
        tk : dict of grids, tk[name][index_k, index_z]
        k : vector of k values, k[index_k]
        z : vector of z values, z[index_z]
        """
        k = self._scalar_transfer_k()
        z = self._z_grid_for_output()
        k_size = k.size
        if h_units:
            k = k/self.h()

        tk = {}
        for index_z in range(z.size):
            # By node, not by redshift: the stored table is what these grids are,
            # and its topmost node does not survive the round trip through z.
            transfers = self._transfer_at_node(index_z, output_format)
            # get_transfer leads with the k column; k is returned separately here.
            transfers.pop('k (h/Mpc)', None)
            for name, values in transfers.items():
                if name not in tk:
                    tk[name] = np.empty((k_size, z.size), np.double)
                tk[name][:, index_z] = values

        return tk, k, z

    def get_Weyl_pk_and_k_and_z(self, nonlinear=False, h_units=False):
        """
        Return the grid of k**4 times the power spectrum of the Weyl potential
        (phi+psi)/2, with the k and z it was computed on.

        This is the matter power spectrum rescaled by [(phi+psi)/2/d_m]**2. The
        factor k**4 is a convention: the Poisson equation carries a k**2, so the
        rescaled spectrum has a shape similar to the matter one.

        Parameters
        ----------
        nonlinear : bool
                Rescale the non-linear matter spectrum; the default is to rescale
                the linear one
        h_units : bool
                Return k in h/Mpc rather than 1/Mpc (default)

        Returns
        -------
        Weyl_pk : grid of values, Weyl_pk[index_k, index_z]
        k : vector of k values, k[index_k]
        z : vector of z values, z[index_z]
        """
        pk, k, z = self.get_pk_and_k_and_z(nonlinear=nonlinear,
                                           only_clustering_species=False,
                                           h_units=h_units)
        tk, _, _ = self.get_transfer_and_k_and_z(output_format='class', h_units=h_units)
        for name in ('phi', 'psi', 'd_m'):
            if name not in tk:
                raise CosmoSevereError(
                    "The Weyl spectrum is built from the '%s' transfer function, which was not "
                    "computed. It needs 'dTk' in 'output' (and, for d_m, 'mPk')." % name)
        # k enters in whichever units the caller asked for, so that the convention
        # is the same on both sides of the rescaling.
        return pk*((tk['phi'] + tk['psi'])/2./tk['d_m'])**2*k[:, None]**4, k, z

    # ---------------------------------------------------------------------- sigma --

    def sigma(self, R, z, h_units=False):
        """
        Gives sigma (total matter) for a given R and z

        R is the radius in units of Mpc, so R=8/h gives the usual sigma8(z) --
        unless h_units is set, in which case R is in Mpc/h and R=8 gives it.

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        if h_units:
            R = R/self.h()
        return self._sigma(R, z, False)

    def sigma_cb(self, R, z, h_units=False):
        """
        Gives sigma (cdm+b) for a given R and z

        R is the radius in units of Mpc, so R=8/h gives the usual sigma8(z) --
        unless h_units is set, in which case R is in Mpc/h and R=8 gives it.

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        if h_units:
            R = R/self.h()
        return self._sigma(R, z, True)

    def effective_f_sigma8(self, z, z_step=0.1):
        """
        Estimate f(z)*sigma8(z) as (d sigma8 / d ln a)(z), by finite difference.

        Parameters
        ----------
        z : float
                Desired redshift
        z_step : float
                Step of the two-sided derivative. For z below z_step the step
                shrinks to z, so the derivative stays two-sided; below z_step/10
                it goes one-sided instead.

        Returns
        -------
        (d sigma8/d ln a)(z) (dimensionless)
        """
        # d sigma8/d ln a = -(d sigma8/dz)*(1+z)
        if z < z_step/10.:
            z_step = z_step/10.
            return (self.sigma(8, z, h_units=True)
                    - self.sigma(8, z + z_step, h_units=True))/z_step*(1 + z)

        if z < z_step:
            z_step = z

        return (self.sigma(8, z - z_step, h_units=True)
                - self.sigma(8, z + z_step, h_units=True))/(2.*z_step)*(1 + z)

    # ---------------------------------------------------------------- CosmoHammer --

    def __call__(self, ctx):
        """
        Function to interface with CosmoHammer

        Parameters
        ----------
        ctx : context
                Contains several dictionaries storing data and cosmological
                information

        """
        data = ctx.get('data')  # recover data from the context

        # Set the module to the current values
        self.empty()
        self.set(data.cosmo_arguments)
        self.compute()

        # Compute the derived paramter value and store them
        params = ctx.getData()
        self.get_current_derived_parameters(
            data.get_mcmc_parameters(['derived']))
        for elem in data.get_mcmc_parameters(['derived']):
            data.mcmc_parameters[elem]['current'] /= \
                data.mcmc_parameters[elem]['scale']
            params[elem] = data.mcmc_parameters[elem]['current']

        ctx.add('boundary', True)
        # Store itself into the context, to be accessed by the likelihoods
        ctx.add('cosmo', self)
