
from libcpp cimport bool
from libcpp.map cimport map
from libcpp.memory cimport shared_ptr
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "class.h":
    pair[string, string] get_my_py_error_message()

    ctypedef char FileArg[1024]
    ctypedef char ErrorMsg[2048]
    ctypedef char FileName[256]

    ctypedef enum file_format :class_format,camb_format
    ctypedef enum equation_of_state :CLP,EDE
    ctypedef enum selection_type :gaussian,tophat,dirac
    ctypedef enum primordial_spectrum_type :
       analytic_Pk,
       two_scales,
       inflation_V,
       inflation_H,
       inflation_V_end,
       external_Pk
    ctypedef enum linear_or_logarithmic :
       linear,
       logarithmic
    ctypedef enum non_linear_method :nl_none,nl_halofit,nl_HMcode
    ctypedef enum pk_outputs :pk_linear,pk_nonlinear
    ctypedef enum source_extrapolation :extrap_zero,extrap_only_max,extrap_only_max_units,extrap_max_scaled,extrap_hmcode,extrap_user_defined
    ctypedef enum halofit_integral_type :halofit_integral_one, halofit_integral_two, halofit_integral_three
    ctypedef enum hmcode_baryonic_feedback_model :nl_emu_dmonly, nl_owls_dmonly, nl_owls_ref, nl_owls_agn, nl_owls_dblim, nl_user_defined
    ctypedef enum out_sigmas :out_sigma,out_sigma_prime,out_sigma_disp

    cdef struct precision:
        FileArg class_dir
        double smallest_allowed_variation

    cdef struct background:
        double H0
        double Omega0_g
        double T_cmb
        double Omega0_b
        double Omega0_cdm
        double Omega0_lambda
        double Omega0_fld
        equation_of_state fluid_equation_of_state
        double w0_fld
        double wa_fld
        double Omega_EDE
        double cs2_fld
        short use_ppf
        double c_gamma_over_c_fld
        double Omega0_ur
        double Omega0_idr
        double T_idr
        double Omega0_idm_dr
        double Omega0_dcdmdr
        double Gamma_dcdm
        double Omega_ini_dcdm
        double Omega0_scf
        short attractor_ic_scf
        double phi_ini_scf
        double phi_prime_ini_scf
        int scf_tuning_index
        double Omega0_k
        int N_ncdm
        double Omega0_ncdm_tot
        int N_decay_dr
        double h
        double K
        int sgnK
        double a_today
        short has_cdm
        short has_dcdm
        short has_ncdm_decay_dr
        short has_dr
        short has_scf
        short has_ncdm
        short has_lambda
        short has_fld
        short has_ur
        short has_idr
        short has_idm_dr
        short has_curvature
        short has_idm_drmd
        short has_idr_drmd
        short has_idm
        double Omega0_idr_drmd
        double Omega0_idm_drmd
        double f_idm_drmd
        double G_over_aH_drmd
        double delta_Neff_drmd
        double z_stop
        short short_info
        short normal_info
        short long_info
        short inter_normal
        short inter_closeby
        short background_verbose

    cdef struct lensing:
        short has_lensed_cls
        short lensing_verbose

    cdef struct perturbs:
        short has_perturbations
        short has_cls
        short has_scalars
        short has_vectors
        short has_tensors
        short has_ad
        short has_bi
        short has_cdi
        short has_nid
        short has_niv
        short has_perturbed_recombination
        short has_cl_cmb_temperature
        short has_cl_cmb_polarization
        short has_cl_cmb_lensing_potential
        short has_cl_lensing_potential
        short has_cl_number_count
        short has_pk_matter
        short has_density_transfers
        short has_velocity_transfers
        short has_Nbody_gauge_transfers
        short has_nl_corrections_based_on_delta_m
        short has_nc_density
        short has_nc_rsd
        short has_nc_lens
        short has_nc_gr
        int l_scalar_max
        int l_vector_max
        int l_tensor_max
        int l_lss_max
        double k_max_for_pk
        int selection_num
        selection_type selection
        double selection_mean[100]
        double selection_width[100]
        int switch_sw
        int switch_eisw
        int switch_lisw
        int switch_dop
        int switch_pol
        double eisw_lisw_split_z
        int store_perturbations
        int k_output_values_num
        double k_output_values[30]
        double G_eff_ur
        double z_max_pk
        int idr_nature
        double selection_min_of_tau_min
        double selection_max_of_tau_max
        double selection_delta_tau
        short perturbations_verbose

    cdef struct thermo:
        double YHe
        double tau_reio
        double z_reio
        short compute_cb2_derivatives
        short compute_damping_scale
        double reionization_width
        double reionization_exponent
        double helium_fullreio_redshift
        double helium_fullreio_width
        int binned_reio_num
        double binned_reio_step_sharpness
        int many_tanh_num
        double many_tanh_width
        int reio_inter_num
        double annihilation
        short has_on_the_spot
        double decay
        double annihilation_variation
        double annihilation_z
        double annihilation_zmax
        double annihilation_zmin
        double annihilation_f_halo
        double annihilation_z_halo
        double a_idm_dr
        double b_idr
        double nindex_idm_dr
        double m_idm
        short thermodynamics_verbose

    cdef struct primordial:
        double k_pivot
        primordial_spectrum_type primordial_spec_type
        double A_s
        double sigma8
        double n_s
        double alpha_s
        double r
        double n_t
        double alpha_t
        double f_bi
        double n_bi
        double alpha_bi
        double f_cdi
        double n_cdi
        double alpha_cdi
        double f_nid
        double n_nid
        double alpha_nid
        double f_niv
        double n_niv
        double alpha_niv
        double c_ad_bi
        double n_ad_bi
        double alpha_ad_bi
        double c_ad_cdi
        double n_ad_cdi
        double alpha_ad_cdi
        double c_ad_nid
        double n_ad_nid
        double alpha_ad_nid
        double c_ad_niv
        double n_ad_niv
        double alpha_ad_niv
        double c_bi_cdi
        double n_bi_cdi
        double alpha_bi_cdi
        double c_bi_nid
        double n_bi_nid
        double alpha_bi_nid
        double c_bi_niv
        double n_bi_niv
        double alpha_bi_niv
        double c_cdi_nid
        double n_cdi_nid
        double alpha_cdi_nid
        double c_cdi_niv
        double n_cdi_niv
        double alpha_cdi_niv
        double c_nid_niv
        double n_nid_niv
        double alpha_nid_niv
        double V0
        double V1
        double V2
        double V3
        double V4
        double H0
        double H1
        double H2
        double H3
        double H4
        double phi_end
        double phi_pivot_target
        double custom1
        double custom2
        double custom3
        double custom4
        double custom5
        double custom6
        double custom7
        double custom8
        double custom9
        double custom10
        short primordial_verbose

    cdef struct spectra:
        double z_max_pk
        int non_diag
        short spectra_verbose

    cdef struct transfers:
        double lcmb_rescale
        double lcmb_tilt
        double lcmb_pivot
        double selection_bias[100]
        double selection_magnification_bias[100]
        short has_nz_file
        short has_nz_analytic
        short has_nz_evo_file
        short has_nz_evo_analytic
        short initialise_HIS_cache
        short transfer_verbose

    cdef struct output:
        int z_pk_num
        double z_pk[100]
        short write_header
        file_format output_format
        short write_background
        short write_thermodynamics
        short write_perturbations
        short write_primordial
        short output_verbose

    cdef struct nonlinear:
        non_linear_method method
        source_extrapolation extrapolation_method
        hmcode_baryonic_feedback_model feedback
        double c_min
        double eta_0
        double z_infinity
        short has_pk_eq
        short nonlinear_verbose


cdef extern from "parser.h":
    cdef cppclass FileContent:
        bool is_shooting
        void set(const string& name, const string& value) except +
        int size() except +
        bool read_int(const string& name, int& value) except +
        bool read_double(const string& name, double& value) except +
        bool read_string(const string& name, string& value) except +
        bool read_list_of_doubles(const string& name, vector[double]& values) except +
        bool read_list_of_integers(const string& name, vector[int]& values) except +
        bool read_list_of_strings(const string& name, vector[string]& values) except +
        void mark_all_unread() except +
        bool was_read(const string& name) except +
        vector[string] unread_parameters() except +

cdef extern from "thermodynamics_module.h":
    cdef cppclass ThermodynamicsModule:
        int thermodynamics_output_titles(char titles[8000]) except +
        int thermodynamics_output_data(int number_of_titles, double* data) except +
        int thermodynamics_at_z(double z, short inter_mode, int* last_index, double* pvecback, double* pvecthermo) except +
        double tau_ini_
        double YHe_
        short inter_normal_
        short inter_closeby_
        int tt_size_
        double z_rec_
        int th_size_
        double tau_rec_
        double angular_rescaling_
        double tau_free_streaming_
        double tau_idr_free_streaming_
        double tau_cut_
        double tau_reionization_
        double z_reionization_
        double n_e_
        double ds_rec_
        double da_rec_
        double rd_rec_
        double rs_rec_
        double ra_rec_
        double rs_star_
        double ra_star_
        double z_star_
        double tau_star_
        double ds_star_
        double da_star_
        double rd_star_
        double z_d_
        double tau_d_
        double rs_d_
        double ds_d_
        int index_th_xe_
        int index_th_rate_
        int index_th_tau_d_
        int index_th_dkappa_
        int index_th_ddkappa_
        int index_th_dddkappa_
        int index_th_exp_m_kappa_
        int index_th_g_
        int index_th_dg_
        int index_th_ddg_
        int index_th_dmu_idm_dr_
        int index_th_ddmu_idm_dr_
        int index_th_dddmu_idm_dr_
        int index_th_dmu_idr_
        int index_th_tau_idm_dr_
        int index_th_tau_idr_
        int index_th_g_idm_dr_
        int index_th_cidm_dr2_
        int index_th_Tidm_dr_
        int index_th_Tb_
        int index_th_wb_
        int index_th_cb2_
        int index_th_dcb2_
        int index_th_ddcb2_
        int index_th_r_d_
        ErrorMsg error_message_

cdef extern from "primordial_module.h":
    cdef cppclass PrimordialModule:
        int primordial_spectrum_at_k(int index_md, linear_or_logarithmic mode, double k, double* pk) except +
        int primordial_output_titles(char titles[8000]) except +
        int primordial_output_data(int number_of_titles, double* data) except +
        int* ic_size_
        int* ic_ic_size_
        short** is_non_zero_
        int lnk_size_
        double phi_pivot_
        double phi_min_
        double phi_max_
        double phi_stop_
        double A_s_
        double n_s_
        double alpha_s_
        double beta_s_
        double r_
        double n_t_
        double alpha_t_
        ErrorMsg error_message_

cdef extern from "spectra_module.h":
    cdef cppclass SpectraModule:
        int spectra_cl_at_l(double l, double * cl, double ** cl_md, double ** cl_md_ic) except +
        map[string, int] cl_output_index_map() except +
        map[string, vector[double]] cl_output(int lmax) except +
        void cl_output_no_copy(int lmax, vector[double*]& output_pointers) except +
        int md_size_
        int* ic_size_
        int* ic_ic_size_
        short** is_non_zero_
        int has_tt_
        int has_ee_
        int has_te_
        int has_bb_
        int has_pp_
        int has_tp_
        int has_ep_
        int has_dd_
        int has_td_
        int has_pd_
        int has_ll_
        int has_tl_
        int has_dl_
        int index_ct_tt_
        int index_ct_ee_
        int index_ct_te_
        int index_ct_bb_
        int index_ct_pp_
        int index_ct_tp_
        int index_ct_ep_
        int index_ct_dd_
        int index_ct_td_
        int index_ct_pd_
        int index_ct_ll_
        int index_ct_tl_
        int index_ct_dl_
        int ct_size_
        int d_size_
        int l_size_max_
        double* l_
        int** l_max_ct_
        int* l_max_
        int l_max_tot_
        ErrorMsg error_message_

cdef extern from "lensing_module.h":
    cdef cppclass LensingModule:
        map[string, vector[double]] cl_output(int lmax) except +
        map[string, vector[double]] cl_output_computed() except +
        int lensing_cl_at_l(int l, double * cl_lensed) except +
        int l_unlensed_max_
        int l_lensed_max_
        ErrorMsg error_message_

cdef extern from "input_module.h":
    cdef cppclass InputModule:
        int file_content_from_arguments(int argc, char** argv, FileContent& fc, ErrorMsg errmsg) except +
        precision precision_
        background background_
        thermo thermodynamics_
        perturbs perturbations_
        transfers transfers_
        primordial primordial_
        spectra spectra_
        nonlinear nonlinear_
        lensing lensing_
        output output_
        shared_ptr[NonColdDarkMatter] ncdm_
        ErrorMsg error_message_

cdef extern from "perturbations_module.h":
    cdef cppclass PerturbationsModule:
        int perturb_output_data(file_format output_format, double z, int number_of_titles, double* data) except +
        int perturb_output_titles(file_format output_format, char titles[8000]) except +
        int perturb_output_firstline_and_ic_suffix(int index_ic, char first_line[1024], FileName ic_suffix) except +
        int index_md_scalars_
        int index_md_tensors_
        int index_md_vectors_
        int md_size_
        int index_ic_ad_
        int index_ic_cdi_
        int index_ic_bi_
        int index_ic_nid_
        int index_ic_niv_
        int index_ic_ten_
        int* ic_size_
        int index_tp_t0_
        int index_tp_t1_
        int index_tp_t2_
        int index_tp_p_
        int index_tp_delta_tot_
        int index_tp_delta_g_
        int index_tp_delta_b_
        int index_tp_delta_cdm_
        int index_tp_delta_dcdm_
        int index_tp_delta_fld_
        int index_tp_delta_scf_
        int index_tp_delta_dr_
        int index_tp_delta_ur_
        int index_tp_delta_idr_
        int index_tp_delta_idm_dr_
        int index_tp_delta_idr_drmd_
        int index_tp_delta_idm_drmd_
        int index_tp_delta_ncdm1_
        int index_tp_perturbed_recombination_delta_temp_
        int index_tp_perturbed_recombination_delta_chi_
        int index_tp_theta_m_
        int index_tp_theta_cb_
        int index_tp_theta_tot_
        int index_tp_theta_g_
        int index_tp_theta_b_
        int index_tp_theta_cdm_
        int index_tp_theta_dcdm_
        int index_tp_theta_fld_
        int index_tp_theta_scf_
        int index_tp_theta_ur_
        int index_tp_theta_idr_
        int index_tp_theta_idm_dr_
        int index_tp_theta_idr_drmd_
        int index_tp_theta_idm_drmd_
        int index_tp_theta_dr_
        int index_tp_theta_ncdm1_
        int index_tp_phi_
        int index_tp_phi_prime_
        int index_tp_phi_plus_psi_
        int index_tp_psi_
        int index_tp_h_
        int index_tp_h_prime_
        int index_tp_eta_
        int index_tp_eta_prime_
        int index_tp_H_T_Nb_prime_
        int index_tp_k2gamma_Nb_
        int index_tp_delta_m_
        int index_tp_delta_cb_
        int* tp_size_
        short has_source_t_
        short has_source_p_
        short has_source_delta_m_
        short has_source_delta_cb_
        short has_source_delta_tot_
        short has_source_delta_g_
        short has_source_delta_b_
        short has_source_delta_cdm_
        short has_source_delta_dcdm_
        short has_source_delta_fld_
        short has_source_delta_scf_
        short has_source_delta_dr_
        short has_source_delta_ur_
        short has_source_delta_idr_
        short has_source_delta_idm_dr_
        short has_source_delta_idr_drmd_
        short has_source_delta_idm_drmd_
        short has_source_delta_ncdm_
        short has_source_theta_m_
        short has_source_theta_cb_
        short has_source_theta_tot_
        short has_source_theta_g_
        short has_source_theta_b_
        short has_source_theta_cdm_
        short has_source_theta_dcdm_
        short has_source_theta_fld_
        short has_source_theta_scf_
        short has_source_theta_dr_
        short has_source_theta_ur_
        short has_source_theta_idr_
        short has_source_theta_idm_dr_
        short has_source_theta_idr_drmd_
        short has_source_theta_idm_drmd_
        short has_source_theta_ncdm_
        short has_source_phi_
        short has_source_phi_prime_
        short has_source_phi_plus_psi_
        short has_source_psi_
        short has_source_h_
        short has_source_h_prime_
        short has_source_eta_
        short has_source_eta_prime_
        short has_source_H_T_Nb_prime_
        short has_source_k2gamma_Nb_
        int* index_k_output_values_
        char scalar_titles_[8000]
        char vector_titles_[8000]
        char tensor_titles_[8000]
        double* scalar_perturbations_data_[30]
        double* vector_perturbations_data_[30]
        double* tensor_perturbations_data_[30]
        int size_scalar_perturbation_data_[30]
        int size_vector_perturbation_data_[30]
        int size_tensor_perturbation_data_[30]
        double*** sources_
        double* ln_tau_
        int ln_tau_size_
        double* tau_sampling_
        int tau_size_
        int* k_size_cl_
        int* k_size_
        double** k_
        double k_min_
        double k_max_
        ErrorMsg error_message_

cdef extern from "nonlinear_module.h":
    cdef cppclass NonlinearModule:
        int nonlinear_pk_at_z(linear_or_logarithmic mode, pk_outputs pk_output, double z, int index_pk, double* out_pk, double* out_pk_ic) except +
        int nonlinear_pks_at_z(linear_or_logarithmic mode, pk_outputs pk_output, double z, double* out_pk, double* out_pk_ic, double* out_pk_cb, double* out_pk_cb_ic) except +
        int nonlinear_pk_at_k_and_z(pk_outputs pk_output, double k, double z, int index_pk, double* out_pk, double* out_pk_ic) except +
        int nonlinear_pks_at_k_and_z(pk_outputs pk_output, double k, double z, double* out_pk, double* out_pk_ic, double* out_pk_cb, double* out_pk_cb_ic) except +
        int nonlinear_pks_at_kvec_and_zvec(pk_outputs pk_output, double* kvec, int kvec_size, double* zvec, int zvec_size, double* out_pk, double* out_pk_cb) except +
        int nonlinear_sigmas_at_z(double R, double z, int index_pk, out_sigmas sigma_output, double* result) except +
        int nonlinear_pk_tilt_at_k_and_z(pk_outputs pk_output, double k, double z, int index_pk, double* pk_tilt) except +
        int nonlinear_k_nl_at_z(double z, double* k_nl, double* k_nl_cb) except +
        int nonlinear_sigma_at_z(double R, double z, int index_pk, double k_per_decade, double* result) except +
        int k_size_
        double* ln_k_
        double** nl_corr_density_
        short* is_non_zero_
        int ic_size_
        int ic_ic_size_
        short has_pk_m_
        short has_pk_cb_
        int index_pk_m_
        int index_pk_cb_
        int pk_size_
        double* sigma8_
        ErrorMsg error_message_

cdef extern from "background_module.h":
    cdef cppclass BackgroundModule:
        int background_output_titles(char titles[8000]) except +
        int background_output_data(int number_of_titles, double* data) except +
        int background_at_tau(double tau, short return_format, short inter_mode, int* last_index, double* pvecback) except +
        int background_tau_of_z(double z, double* tau) except +
        int background_w_fld(double a, double* w_fld, double* dw_over_da_fld, double* integral_fld) except +
        int background_free_noinput() except +
        double dV_scf(double phi) except +
        int background_idm_drmd(double a, double rho_idm_over_rho_idr, double *Rint, double *csp2, double *Gint) except +
        int index_bg_a_
        int index_bg_H_
        int index_bg_H_prime_
        int index_bg_rho_g_
        int index_bg_rho_b_
        int index_bg_rho_cdm_
        int index_bg_rho_lambda_
        int index_bg_rho_fld_
        int index_bg_w_fld_
        int index_bg_rho_ur_
        int index_bg_rho_idm_dr_
        int index_bg_rho_idr_
        int index_bg_rho_dcdm_
        int index_bg_lnf_ncdm_decay_dr1_
        int index_bg_dlnfdlnq_ncdm_decay_dr1_
        int index_bg_dlnfdlnq_separate_ncdm_decay_dr1_
        int index_bg_rho_idm_drmd_
        int index_bg_G_over_aH_drmd_
        int index_bg_rho_idr_drmd_
        int index_bg_Gamma0_drmd_
        int index_bg_rho_dr_species_
        int index_bg_rho_dr_
        int index_bg_phi_scf_
        int index_bg_phi_prime_scf_
        int index_bg_V_scf_
        int index_bg_dV_scf_
        int index_bg_ddV_scf_
        int index_bg_rho_scf_
        int index_bg_p_scf_
        int index_bg_p_prime_scf_
        int index_bg_number_ncdm1_
        int index_bg_rho_ncdm1_
        int index_bg_p_ncdm1_
        int index_bg_pseudo_p_ncdm1_
        int index_bg_rho_tot_
        int index_bg_p_tot_
        int index_bg_p_tot_prime_
        int index_bg_Omega_r_
        int index_bg_rho_crit_
        int index_bg_Omega_m_
        int index_bg_conf_distance_
        int index_bg_ang_distance_
        int index_bg_lum_distance_
        int index_bg_time_
        int index_bg_rs_
        int index_bg_D_
        int index_bg_f_
        int bg_size_short_
        int bg_size_normal_
        int bg_size_
        int bt_size_
        double* tau_table_
        double* background_table_
        double conformal_age_
        double Neff_
        double a_eq_
        double H_eq_
        double Omega0_m_
        double Omega0_de_
        double age_
        double Omega0_r_
        double z_eq_
        double tau_eq_
        double Omega0_dcdm_
        double Omega0_dr_
        double Omega0_idr_drmd_
        double Omega0_idm_drmd_
        double Omega0_idm_
        double G_over_aH_tmp_
        double Gamma0_drmd_
        double f_idr_drmd_
        double z_dec_drmd_
        ErrorMsg error_message_

cdef extern from "transfer_module.h":
    cdef cppclass TransferModule:
        int index_tt_t0_
        int index_tt_t1_
        int index_tt_t2_
        int index_tt_e_
        int index_tt_b_
        int index_tt_lcmb_
        int index_tt_density_
        int index_tt_lensing_
        int index_tt_rsd_
        int index_tt_d0_
        int index_tt_d1_
        int index_tt_nc_lens_
        int index_tt_nc_g1_
        int index_tt_nc_g2_
        int index_tt_nc_g3_
        int index_tt_nc_g4_
        int index_tt_nc_g5_
        int* tt_size_
        int l_size_max_
        int** l_size_tt_
        int* l_size_
        int* l_
        int q_size_
        double* q_
        double** k_
        int index_q_flat_approximation_
        double** transfer_
        ErrorMsg error_message_

cdef extern from "non_cold_dark_matter.h":
    cdef cppclass NcdmSettings:
        double h
        double T_cmb
        double tol_ncdm
        double tol_ncdm_bg
        double tol_M_ncdm

cdef extern from "non_cold_dark_matter.h":
    cdef cppclass NonColdDarkMatter:
        shared_ptr[NonColdDarkMatter] Create(FileContent* pfc, const NcdmSettings&) except +
        int background_ncdm_momenta(int n_ncdm, double z, double* n, double* rho, double* p, double* drho_dM, double* pseudo_p) except +
        int background_ncdm_momenta_deg(int n_ncdm, double deg, double z, double T_cmb, double* n, double* rho, double* p, double* drho_ddeg, double* pseudo_p) except +
        double GetOmega0() except +
        double GetNeff(double z) except +
        double GetMassInElectronvolt(int n_ncdm) except +
        void PrintNeffInfo() except +
        void PrintMassInfo() except +
        void PrintOmegaInfo() except +
        void SetBackgroundWeight(int n_ncdm, int q_index, double weight) except +
        void SetOmega0(int n_ncdm, double Omega0, double h) except +
        void SetDegAndFactor(int n_ncdm, double deg, double T_cmb) except +
        void SetDeg_from_Omega_ini(int n_ncdm, double z_ini, double H0, double Omega_ini, double T_cmb) except +
        double GetIni(double a, double a_today, double tol_ncdm_initial_w) except +
        double GetDeg(int n_ncdm) except +
        double GetRescalingFactor(int n_ncdm, double* pvecback_begin) except +
        int N_ncdm_
        int N_ncdm_standard_
        int N_ncdm_decay_dr_
        double* m_ncdm_in_eV_
        double* M_ncdm_
        int* q_size_ncdm_
        double* factor_ncdm_
        double** q_ncdm_
        double** w_ncdm_
        double** dlnf0_dlnq_ncdm_
        double* Omega_dncdmdr_
        int q_total_size_dncdm_

cdef extern from "class.h":
    cdef cppclass ClassConstants:
        int sMAXTITLESTRINGLENGTH
        int sFALSE
        int sARGUMENT_LENGTH_MAX
        int sFAILURE


