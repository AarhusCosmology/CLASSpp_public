/** @file input.c Documented input module.
 *
 * Julien Lesgourgues, 27.08.2010
 */

#include "input_module.h"

#include <algorithm>
#include <filesystem>

#include "background_module.h"
#include "cosmology.h"
#include "lensing_module.h"
#include "nonlinear_module.h"
#include "perturbations_module.h"
#include "primordial_module.h"
#include "spectra_module.h"
#include "thermodynamics_module.h"

// All species — single consolidated header
#include <thread>

#include "../species/all_species.h"
#include "../species/dcdm_dr_species.h"
#include "../species/species_input.h"

namespace {

int readDoubleList(
    FileContent* pfc, const char* name, std::vector<double>& values, int* found, ErrorMsg errmsg) {
  try {
    *found = pfc->read_list_of_doubles(name, values) ? _TRUE_ : _FALSE_;
  }
  catch (const std::exception& e) {
    class_stop(errmsg, "%s", e.what());
  }
  return _SUCCESS_;
}

}  // namespace
/**
 * Use this routine to extract initial parameters from files 'xxx.ini'
 * and/or 'xxx.pre'. They can be the arguments of the main() routine.
 *
 * If class is embedded into another code, you will probably prefer to
 * call directly input_init() in order to pass input parameters
 * through a 'file_content' structure.
 */

int InputModule::file_content_from_arguments(int argc,
                                             char** argv,
                                             FileContent& fc,
                                             ErrorMsg errmsg) {
  /** Summary: */

  /** - define local variables */
  FileContent fc_input;     /** - --> a temporary structure with all input parameters */
  FileContent fc_precision; /** - --> a temporary structure with all precision parameters */
  FileContent fc_root;      /** - --> a temporary structure with only the root name */
  FileContent fc_inputroot; /** - --> sum of fc_inoput and fc_root */
  FileContent* pfc_input;   /** - --> a pointer to either fc_root or fc_inputroot */

  char input_file[_ARGUMENT_LENGTH_MAX_];
  char precision_file[_ARGUMENT_LENGTH_MAX_];
  const size_t tmp_file_size =
      _ARGUMENT_LENGTH_MAX_ +
      26;  // 26 is enough to extend the file name [...] with the characters "output/[...]%02d_parameters.ini" (as done below)
  char tmp_file[tmp_file_size];
  char extension[5];
  FileArg stringoutput, inifilename;
  int flag1, filenum;

  pfc_input = &fc_input;

  /** - Initialize the two file_content structures (for input
      parameters and precision parameters) to some null content. If no
      arguments are passed, they will remain null and inform
      init_params() that all parameters take default values. */

  fc                = FileContent();
  fc_input          = FileContent();
  fc_precision      = FileContent();
  input_file[0]     = '\0';
  precision_file[0] = '\0';

  /** - If some arguments are passed, identify eventually some 'xxx.ini'
      and 'xxx.pre' files, and store their name. */

  if (argc > 1) {
    for (int i = 1; i < argc; i++) {
      strncpy(extension, (argv[i] + strlen(argv[i]) - 4), 4);
      extension[4] = '\0';
      if (strcmp(extension, ".ini") == 0) {
        class_test(input_file[0] != '\0',
                   errmsg,
                   "You have passed more than one input file with extension '.ini', choose one.");
        strcpy(input_file, argv[i]);
      }
      else if (strcmp(extension, ".pre") == 0) {
        class_test(precision_file[0] != '\0',
                   errmsg,
                   "You have passed more than one precision with extension '.pre', choose one.");
        strcpy(precision_file, argv[i]);
      }
      else {
        fprintf(stdout,
                "Warning: the file '%s' has an extension different from .ini and .pre, so it has "
                "been ignored\n",
                argv[i]);
      }
    }
  }

  /** - if there is an 'xxx.ini' file, read it and store its content. */

  if (input_file[0] != '\0') {
    class_call(parser_read_file(input_file, &fc_input, errmsg), errmsg, errmsg);

    /** - check whether a root name has been set */

    class_call(parser_read_string(&fc_input, "root", &stringoutput, &flag1, errmsg),
               errmsg,
               errmsg);

    /** - if root has not been set, use root=output/inputfilennameN_ */

    if (flag1 == _FALSE_) {
      //printf("strlen-4 = %zu\n",strlen(input_file)-4);
      strncpy(inifilename, input_file, strlen(input_file) - 4);
      inifilename[strlen(input_file) - 4] = '\0';
      for (filenum = 0; filenum < 100; filenum++) {
        std::error_code ec;
        snprintf(tmp_file, tmp_file_size, "output/%s%02d_cl.dat", inifilename, filenum);
        if (std::filesystem::exists(tmp_file, ec))
          continue;
        snprintf(tmp_file, tmp_file_size, "output/%s%02d_pk.dat", inifilename, filenum);
        if (std::filesystem::exists(tmp_file, ec))
          continue;
        snprintf(tmp_file, tmp_file_size, "output/%s%02d_tk.dat", inifilename, filenum);
        if (std::filesystem::exists(tmp_file, ec))
          continue;
        snprintf(tmp_file, tmp_file_size, "output/%s%02d_parameters.ini", inifilename, filenum);
        if (std::filesystem::exists(tmp_file, ec))
          continue;
        break;
      }
      snprintf(tmp_file, tmp_file_size, "output/%s%02d_", inifilename, filenum);
      fc_root.set("root", tmp_file);
      class_call(parser_cat(&fc_input, &fc_root, &fc_inputroot, errmsg), errmsg, errmsg);
      pfc_input = &fc_inputroot;
    }
  }

  /** - if there is an 'xxx.pre' file, read it and store its content. */

  if (precision_file[0] != '\0')

    class_call(parser_read_file(precision_file, &fc_precision, errmsg), errmsg, errmsg);

  /** - if one or two files were read, merge their contents in a
      single 'file_content' structure. */

  if ((input_file[0] != '\0') || (precision_file[0] != '\0'))

    class_call(parser_cat(pfc_input, &fc_precision, &fc, errmsg), errmsg, errmsg);

  return _SUCCESS_;
}

InputModule::InputModule(FileContent& fc) : file_content_(fc) {
  file_content_.mark_all_unread();
  try {
    input_init();
  }
  catch (const std::runtime_error& e) {
    throw std::invalid_argument(e.what());
  }
  ConstructSpecies();
}

void InputModule::ConstructSpecies() {
  background* pba = &background_;  // need non-const to write closure Omega

  NcdmSettings ncdm_settings;
  ncdm_settings.h           = pba->h;
  ncdm_settings.T_cmb       = pba->T_cmb;
  ncdm_settings.tol_ncdm    = precision_.tol_ncdm;
  ncdm_settings.tol_ncdm_bg = precision_.tol_ncdm_bg;
  ncdm_settings.tol_M_ncdm  = precision_.tol_M_ncdm;

  SpeciesBuildContext ctx{
      /*pfc=*/&file_content_,
      /*pba=*/pba,
      /*ppr=*/&precision_,
      /*ncdm_settings=*/&ncdm_settings,
      /*bgm=*/nullptr,
      /*all_species=*/&all_species_,
      /*omega_budget=*/&omega_budget_,
  };

  // Read input_verbose for the closure verbose message.
  int input_verbose = 0;
  {
    int flag = _FALSE_;
    int val  = 0;
    class_call(parser_read_int(&file_content_, "input_verbose", &val, &flag, error_message_),
               error_message_,
               error_message_);
    if (flag == _TRUE_)
      input_verbose = val;
  }

  const std::string_view closure_name = ClosureSpeciesName(pba->closure_species);

  // Pass 1: build every non-closure species, summing Omega0 contributions.
  // Also tally NCDM-family counters as we go: we know which factories are
  // NCDM-family by their entry name, so no downcast is needed.
  double omega0_sum = 0.0;
  int n_ncdm        = 0;
  for (const auto& entry : kAllSpeciesFactories) {
    if (entry.name == closure_name)
      continue;
    auto produced             = entry.create_all(ctx);
    const bool is_ncdm_family = (entry.name == "NCDM" || entry.name == "DNCDM_DR" ||
                                 entry.name == "NCDMInt");
    for (auto& e : produced) {
      omega0_sum += e.species->GetOmega0();
      if (is_ncdm_family) {
        n_ncdm += 1;
      }
      all_species_.insert(e.key, std::move(e.species));
    }
  }

  // Compute closure value and hand it to the closure species via the context.
  // Closure species' factories pick it up from ctx.omega0_closure_override
  // rather than reading pba->Omega0_X (which is no longer written here).
  if (pba->closure_species != ClosureSpecies::None) {
    const double closure_value = 1.0 - pba->Omega0_k - omega0_sum;
    if (input_verbose > 0) {
      const char* name_for_msg = (pba->closure_species == ClosureSpecies::Lambda)  ? "Omega_Lambda"
                                 : (pba->closure_species == ClosureSpecies::Fluid) ? "Omega_fld"
                                 : (pba->closure_species == ClosureSpecies::ScalarField)
                                     ? "Omega_scf"
                                     : "";
      printf(" -> matched budget equations by adjusting %s = %e\n", name_for_msg, closure_value);
    }

    // Pass 2: build the closure species. Override is consumed by its CreateAll
    // and cleared afterward so subsequent callers (shooting guess, etc.) don't see it.
    ctx.omega0_closure_override = closure_value;
    for (const auto& entry : kAllSpeciesFactories) {
      if (entry.name != closure_name)
        continue;
      auto produced = entry.create_all(ctx);
      for (auto& e : produced) {
        all_species_.insert(e.key, std::move(e.species));
      }
      break;
    }
    ctx.omega0_closure_override.reset();
  }

  all_species_.freeze();

  // Precision-parameter consistency check that fires when any NCDM-family
  // species is present (NCDMSpecies + DNCDMSpecies + NCDMInteractingSpecies).
  const precision* ppr = &precision_;
  if (n_ncdm > 0) {
    if (ppr->ncdm_fluid_trigger_tau_over_tau_k == ppr->radiation_streaming_trigger_tau_over_tau_k) {
      throw std::invalid_argument(
          "please choose different values for precision parameters "
          "ncdm_fluid_trigger_tau_over_tau_k and radiation_streaming_trigger_tau_over_tau_k, "
          "in order to avoid switching two approximation schemes at the same time");
    }
    if (ppr->ncdm_fluid_trigger_tau_over_tau_k == ppr->ur_fluid_trigger_tau_over_tau_k) {
      throw std::invalid_argument(
          "please choose different values for precision parameters "
          "ncdm_fluid_trigger_tau_over_tau_k and ur_fluid_trigger_tau_over_tau_k, in order to "
          "avoid switching two approximation schemes at the same time");
    }
  }
}

int InputModule::ReadCoupledOmegaBudget() {
  char* errmsg        = error_message_;
  FileContent* pfc    = &file_content_;
  precision* ppr      = &precision_;
  background* pba     = &background_;
  thermo* /*pth*/ pth = &thermodynamics_;
  perturbs* ppt       = &perturbations_;

  (void) pth;  // future-proof; the budget computation itself doesn't touch pth.

  int flag1, flag2, flag3, flag4;
  double param1, param2, param3, param4;
  int int1;
  int input_verbose = 0;
  class_read_int("input_verbose", input_verbose);

  // ── IDR: stat_f_idr * (T_idr / T_cmb)^4 * Omega0_g ───────────────────────────
  // T_idr is a physics param owned by IDRSpecies; compute a local copy here
  // solely to derive the omega budget entry for IDR.
  double stat_f_idr = 7. / 8.;
  class_read_double("stat_f_idr", stat_f_idr);

  class_call(parser_read_double(pfc, "N_idr", &param1, &flag1, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "N_dg", &param2, &flag2, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "xi_idr", &param3, &flag3, errmsg), errmsg, errmsg);
  class_test(class_at_least_two_of_three(flag1, flag2, flag3),
             errmsg,
             "In input file, you can only enter one of N_idr, N_dg or xi_idr, choose one");

  double T_idr_local = 0.;
  if (flag1 == _TRUE_) {
    T_idr_local = pow(param1 / stat_f_idr * (7. / 8.) / pow(11. / 4., (4. / 3.)), (1. / 4.)) *
                  pba->T_cmb;
    if (input_verbose > 1)
      printf(
          "You passed N_idr = N_dg = %e, this is equivalent to xi_idr = %e in the ETHOS notation. "
          "\n",
          param1,
          T_idr_local / pba->T_cmb);
  }
  else if (flag2 == _TRUE_) {
    T_idr_local = pow(param2 / stat_f_idr * (7. / 8.) / pow(11. / 4., (4. / 3.)), (1. / 4.)) *
                  pba->T_cmb;
    if (input_verbose > 2)
      printf(
          "You passed N_dg = N_idr = %e, this is equivalent to xi_idr = %e in the ETHOS notation. "
          "\n",
          param2,
          T_idr_local / pba->T_cmb);
  }
  else if (flag3 == _TRUE_) {
    T_idr_local = param3 * pba->T_cmb;
    if (input_verbose > 1)
      printf(
          "You passed xi_idr = %e, this is equivalent to N_idr = N_dg = %e in the NADM notation. "
          "\n",
          param3,
          stat_f_idr * pow(param3, 4.) / (7. / 8.) * pow(11. / 4., (4. / 3.)));
  }

  // Mark the budget slot present iff any of the three IDR-temperature inputs was given.
  // (Matches the legacy semantics: pba->Omega0_idr was always written, but only became
  // nonzero when T_idr was set by one of these inputs.)
  if (flag1 == _TRUE_ || flag2 == _TRUE_ || flag3 == _TRUE_) {
    omega_budget_.idr = stat_f_idr * pow(T_idr_local / pba->T_cmb, 4.) * pba->Omega0_g;
  }

  // ── CDM: parser value (or default), then synchronous-gauge minimum ───────────
  // input_default_params() already set pba->Omega0_cdm = 0.12038/h^2 as the fallback
  // for the closure-Omega budget computation; honour that as the budget default too.
  double omega0_cdm = 0.12038 / (pba->h * pba->h);
  bool cdm_user_set = false;

  class_call(parser_read_double(pfc, "Omega_cdm", &param1, &flag1, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "omega_cdm", &param2, &flag2, errmsg), errmsg, errmsg);
  class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
             errmsg,
             "In input file, you can only enter one of Omega_cdm or omega_cdm, choose one");
  if (flag1 == _TRUE_) {
    omega0_cdm   = param1;
    cdm_user_set = true;
  }
  if (flag2 == _TRUE_) {
    omega0_cdm   = param2 / pba->h / pba->h;
    cdm_user_set = true;
  }

  // Synchronous-gauge minimum: if CDM ends up zero in synchronous gauge, bump it
  // to Omega0_cdm_min_synchronous. Track presence: the gauge minimum kicks in
  // even if the user supplied nothing, so the budget slot becomes present.
  bool cdm_present = cdm_user_set || omega0_cdm != 0.;
  if ((ppt->gauge == synchronous) && (omega0_cdm == 0.)) {
    omega0_cdm  = ppr->Omega0_cdm_min_synchronous;
    cdm_present = true;
  }
  if (cdm_present)
    omega_budget_.cdm = omega0_cdm;

  // ── IDM_DR: Omega_idm_dr / omega_idm_dr / f_idm_dr ───────────────────────────
  class_call(parser_read_double(pfc, "Omega_idm_dr", &param1, &flag1, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "omega_idm_dr", &param2, &flag2, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "f_idm_dr", &param3, &flag3, errmsg), errmsg, errmsg);
  class_test(class_at_least_two_of_three(flag1, flag2, flag3),
             errmsg,
             "In input file, you can only enter one of Omega_idm_dr, omega_idm_dr or f_idm_dr, "
             "choose one");

  if (flag1 == _TRUE_)
    omega_budget_.idm_dr = param1;
  if (flag2 == _TRUE_)
    omega_budget_.idm_dr = param2 / pba->h / pba->h;

  if (flag3 == _TRUE_) {
    class_test((param3 < 0.) || (param3 > 1.),
               errmsg,
               "The fraction of interacting DM with DR must be between 0 and 1, you asked for "
               "f_idm_dr=%e",
               param3);
    const double cdm_for_frac = omega_budget_.cdm.value_or(0.);
    class_test((param3 > 0.) && (cdm_for_frac == 0.),
               errmsg,
               "If you want a fraction of interacting DM with DR, to be consistent, you should not "
               "set the fraction of CDM to zero");

    double new_idm_dr = param3 * cdm_for_frac;
    double new_cdm    = cdm_for_frac - new_idm_dr;
    // avoid CDM=0 in synchronous gauge after the subtraction
    if ((ppt->gauge == synchronous) && (new_cdm == 0.)) {
      new_cdm    += ppr->Omega0_cdm_min_synchronous;
      new_idm_dr -= ppr->Omega0_cdm_min_synchronous;
    }
    omega_budget_.idm_dr = new_idm_dr;
    omega_budget_.cdm    = new_cdm;
  }

  // ── DCDM_DR: Omega_dcdmdr / omega_dcdmdr ────────────────────────────────────
  // Gamma_dcdm and Omega_ini_dcdm are physics params owned by DCDMSpecies;
  // they are parsed in DCDM_DR_Species::CreateAll, not here.
  class_call(parser_read_double(pfc, "Omega_dcdmdr", &param1, &flag1, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "omega_dcdmdr", &param2, &flag2, errmsg), errmsg, errmsg);
  class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
             errmsg,
             "In input file, you can only enter one of Omega_dcdmdr or omega_dcdmdr, choose one");
  if (flag1 == _TRUE_)
    omega_budget_.dcdmdr = param1;
  if (flag2 == _TRUE_)
    omega_budget_.dcdmdr = param2 / pba->h / pba->h;

  // ── DRMD: z_stop, G_over_aH_drmd_ini, f_idm_drmd, delta_Neff_drmd ────────────
  // These fields now live on IDM_DRMD_IDR_DRMD_Species; we parse them locally here
  // only for the budget math (delta_Neff_drmd → omega_budget_.idr_drmd;
  // f_idm_drmd → omega_budget_.idm_drmd with the CDM subtraction).
  class_call(parser_read_double(pfc, "z_stop", &param1, &flag1, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "G_over_aH_drmd_ini", &param2, &flag2, errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc, "f_idm_drmd", &param3, &flag3, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "delta_Neff_drmd", &param4, &flag4, errmsg), errmsg, errmsg);

  const double z_stop_drmd           = (flag1 == _TRUE_) ? param1 : 0.;
  const double f_idm_drmd_local      = (flag3 == _TRUE_) ? param3 : 0.;
  const double delta_Neff_drmd_local = (flag4 == _TRUE_) ? param4 : 0.;

  const int any_drmd = (flag1 == _TRUE_) || (flag2 == _TRUE_) || (flag3 == _TRUE_) ||
                       (flag4 == _TRUE_);
  const int all_drmd = (flag1 == _TRUE_) && (flag2 == _TRUE_) && (flag3 == _TRUE_) &&
                       (flag4 == _TRUE_);
  class_test(any_drmd && !all_drmd,
             errmsg,
             "If any DRMD parameter is set, all of them must be non-zero.\nDRMD parameters are "
             "'z_stop', 'G_over_aH_drmd_ini', 'f_idm_drmd' and 'delta_Neff_drmd'.");

  if (delta_Neff_drmd_local > 0.) {
    omega_budget_.idr_drmd = delta_Neff_drmd_local * 7. / 8. * pow(4. / 11., 4. / 3.) *
                             pba->Omega0_g;
    if (f_idm_drmd_local > 0) {
      class_test((z_stop_drmd > 200000.),
                 errmsg,
                 "z_stop is chosen too large. If you want to probe z_stop > 1000000 you need to "
                 "start evolving perturbations earlier in CLASS by changing the precision "
                 "settings. Also you should check that the exponential suppression factor does not "
                 "lead to numerical problems.");
    }
  }

  if (f_idm_drmd_local > 0) {
    class_test((f_idm_drmd_local > 1.),
               errmsg,
               "The fraction of interacting DM with DR must be between 0 and 1, you asked for "
               "f_idm_drmd=%e",
               f_idm_drmd_local);
    const double cdm_for_drmd = omega_budget_.cdm.value_or(0.);
    class_test((cdm_for_drmd == 0.),
               errmsg,
               "If you want a fraction of interacting DM with DRMD, to be consistent, you should "
               "not set the fraction of CDM to zero");

    double new_idm_drmd = f_idm_drmd_local * cdm_for_drmd;
    double new_cdm      = cdm_for_drmd - new_idm_drmd;
    if ((ppt->gauge == synchronous) && (new_cdm == 0.)) {
      new_cdm      += ppr->Omega0_cdm_min_synchronous;
      new_idm_drmd -= ppr->Omega0_cdm_min_synchronous;
    }
    omega_budget_.idm_drmd = new_idm_drmd;
    omega_budget_.cdm      = new_cdm;
  }

  return _SUCCESS_;
}

/**
 * Initialize each parameter, first to its default values, and then
 * from what can be interpreted from the values passed in the input
 * 'file_content' structure. If its size is null, all parameters keep
 * their default values.
 *
 */

int InputModule::input_init() {
  char* errmsg     = error_message_;
  FileContent* pfc = &file_content_;

  int flag1;

  char string1[_ARGUMENT_LENGTH_MAX_];

  /**
   * Before getting into the assignment of parameters,
   * and before the shooting, we want to already fix our precision parameters.
   *
   * No precision parameter should depend on any input parameter
   *
   */

  class_call(input_read_precisions(), error_message_, error_message_);

  /**
   * 'Shooting' resolves inputs that can't be set directly — a condition (e.g. the
   *  angular sound-horizon scale 100*theta_s, or a species' today density) is satisfied
   *  by root-finding an unknown (e.g. h, or Omega_ini_dcdm) through repeated CLASS runs.
   *  This no longer happens in the constructor: each shooting-capable species guesses its
   *  own unknown during ConstructSpecies and reports its target, and Cosmology::GetInputModule
   *  lazily calls InputModule::DoShooting to solve the coupled system. See DoShooting.
   */

  int input_verbose = 0, int1;
  class_read_int("input_verbose", input_verbose);
  if (input_verbose > 0)
    printf("Reading input parameters\n");

  /** - -->  read all parameters from input pfc: */
  class_call(input_read_parameters(), errmsg, errmsg);

  /** - eventually write all the read parameters in a file, unread parameters in another file, and warnings about unread parameters */

  class_call(parser_read_string(pfc, "write parameters", &string1, &flag1, errmsg), errmsg, errmsg);

  if ((flag1 == _TRUE_) &&
      ((strstr(string1, "y") != nullptr) || (strstr(string1, "Y") != nullptr))) {
    output* pop = &output_;
    char param_output_name[_LINE_LENGTH_MAX_];
    char param_unused_name[_LINE_LENGTH_MAX_];
    snprintf(param_output_name, _LINE_LENGTH_MAX_, "%s%s", pop->root, "parameters.ini");
    snprintf(param_unused_name, _LINE_LENGTH_MAX_, "%s%s", pop->root, "unused_parameters");

    FILE* param_output;
    FILE* param_unused;
    class_open(param_output, param_output_name, "w", errmsg);
    class_open(param_unused, param_unused_name, "w", errmsg);

    fprintf(param_output, "# List of input/precision parameters actually read\n");
    fprintf(param_output, "# (all other parameters set to default values)\n");
    fprintf(param_output,
            "# Obtained with CLASS %s (for developers: svn version %s)\n",
            _VERSION_,
            _SVN_VERSION_);
    fprintf(param_output, "#\n");
    fprintf(param_output, "# This file can be used as the input file of another run\n");
    fprintf(param_output, "#\n");

    fprintf(param_unused, "# List of input/precision parameters passed\n");
    fprintf(param_unused, "# but not used (just for info)\n");
    fprintf(param_unused, "#\n");

    pfc->for_each([&](const std::string& name, const std::string& value, bool read) {
      if (read)
        fprintf(param_output, "%s = %s\n", name.c_str(), value.c_str());
      else
        fprintf(param_unused, "%s = %s\n", name.c_str(), value.c_str());
    });
    fprintf(param_output, "#\n");

    fclose(param_output);
    fclose(param_unused);
  }

  class_call(parser_read_string(pfc, "write warnings", &string1, &flag1, errmsg), errmsg, errmsg);

  if ((flag1 == _TRUE_) &&
      ((strstr(string1, "y") != nullptr) || (strstr(string1, "Y") != nullptr))) {
    pfc->for_each([](const std::string& name, const std::string& value, bool read) {
      if (!read)
        fprintf(stdout,
                "[WARNING: input line not recognized and not taken into account: '%s=%s']\n",
                name.c_str(),
                value.c_str());
    });
  }

  return _SUCCESS_;
}
int InputModule::input_read_precisions() {
  precision* ppr = &precision_;
  int flag1;

  class_call(parser_read_string(&file_content_,
                                "class_dir",
                                &(ppr->class_dir),
                                &flag1,
                                error_message_),
             error_message_,
             error_message_);
  if (flag1 == _FALSE_) {
    strncpy(ppr->class_dir, __CLASSDIR__, _ARGUMENT_LENGTH_MAX_);
  }

  /** - set string parameter defaults (require runtime path concatenation) */
  class_call(input_default_precision(), error_message_, error_message_);

  const auto standard_ncdm_instances = file_content_.instances_with("type", "ncdm_standard");
  SynthesiseIdenticalScalarField(&file_content_,
                                 standard_ncdm_instances,
                                 "fluid_approximation",
                                 "ncdm_fluid_approximation",
                                 "dot-syntax standard NCDM species");

  /** - parse all precision parameters from config file */
  ppr->parse(file_content_);

  return _SUCCESS_;
}
int InputModule::input_read_parameters() {
  /** Summary: */

  /** - define local variables */
  char* errmsg     = error_message_;
  FileContent* pfc = &file_content_;
  precision* ppr   = &precision_;      /* for precision parameters */
  background* pba  = &background_;     /* for cosmological background */
  thermo* pth      = &thermodynamics_; /* for thermodynamics */
  perturbs* ppt    = &perturbations_;  /* for source functions */
  primordial* ppm  = &primordial_;     /* for primordial spectra */
  nonlinear* pnl   = &nonlinear_;      /* for non-linear spectra */
  transfers* ptr   = &transfers_;      /* for transfer functions */
  spectra* psp     = &spectra_;        /* for output spectra */
  lensing* ple     = &lensing_;        /* for lensed spectra */
  output* pop      = &output_;

  int flag1, flag2, flag3;
  double param1, param2, param3;
  int entries_read;
  int int1;
  double* pointer1;
  char string1[_ARGUMENT_LENGTH_MAX_];
  char string2[_ARGUMENT_LENGTH_MAX_];
  double k1         = 0.;
  double k2         = 0.;
  double prr1       = 0.;
  double prr2       = 0.;
  double pii1       = 0.;
  double pii2       = 0.;
  double pri1       = 0.;
  double pri2       = 0.;
  double n_iso      = 0.;
  double f_iso      = 0.;
  double n_cor      = 0.;
  double c_cor      = 0.;
  double stat_f_idr = 7. / 8.;

  int i;

  double sigma_B; /* Stefan-Boltzmann constant in \f$ W/m^2/K^4 = Kg/K^4/s^3 \f$*/

  double z_max      = 0.;
  int input_verbose = 0;

  sigma_B = 2. * pow(_PI_, 5) * pow(_k_B_, 4) / 15. / pow(_h_P_, 3) / pow(_c_, 2);

  /** - set all input parameters to default values */

  class_call(input_default_params(), error_message_, error_message_);

  /** - if entries passed in file_content structure, carefully read
      and interpret each of them, and tune the relevant input
      parameters accordingly*/

  class_read_int("input_verbose", input_verbose);

  class_call(parser_read_int(pfc, "threads", &int1, &flag1, errmsg), errmsg, errmsg);
  if (flag1 == _TRUE_) {
    pba->number_of_threads = int1;
  }
  else {
    // If threads was not specified an input, use std::thread::hardware_concurrency() as
    // default. Then loop through a list of environment variable names, and if one of
    // them is set to value 0 < threads <= 8192, use that one instead.
    pba->number_of_threads = std::thread::hardware_concurrency();
    for (const std::string& env_var_name : {"SLURM_CPUS_PER_TASK", "OMP_NUM_THREADS"}) {
      if (char* s = std::getenv(env_var_name.c_str())) {
        int threads = std::atoi(s);
        if ((threads > 0) && (threads <= 8192)) {
          pba->number_of_threads = threads;
          break;
        }
      }
    }
  }

  /** Knowing the gauge from the very beginning is useful (even if
      this could be a run not requiring perturbations at all: even in
      that case, knowing the gauge is important e.g. for fixing the
      sampling in momentum space for non-cold dark matter) */

  class_call(parser_read_string(pfc, "gauge", &string1, &flag1, errmsg), errmsg, errmsg);

  if (flag1 == _TRUE_) {
    if ((strstr(string1, "newtonian") != nullptr) || (strstr(string1, "Newtonian") != nullptr) ||
        (strstr(string1, "new") != nullptr)) {
      ppt->gauge = newtonian;
    }

    if ((strstr(string1, "synchronous") != nullptr) || (strstr(string1, "sync") != nullptr) ||
        (strstr(string1, "Synchronous") != nullptr)) {
      ppt->gauge = synchronous;
    }
  }

  /** (a) background parameters */

  /** - scale factor today (arbitrary) */
  class_read_double("a_today", pba->a_today);

  /** - h (dimensionless) and [\f$ H_0/c\f$] in \f$ Mpc^{-1} = h / 2997.9... = h * 10^5 / c \f$ */
  class_call(parser_read_double(pfc, "H0", &param1, &flag1, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "h", &param2, &flag2, errmsg), errmsg, errmsg);
  class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
             errmsg,
             "In input file, you cannot enter both h and H0, choose one");
  if (flag1 == _TRUE_) {
    pba->H0 = param1 * 1.e3 / _c_;
    pba->h  = param1 / 100.;
  }
  if (flag2 == _TRUE_) {
    pba->H0 = param2 * 1.e5 / _c_;
    pba->h  = param2;
  }

  /** - 100*theta_s is a module-level shooting target (varies h; resolved later by
   *  DoShooting). Consume it here so the unread-parameter check passes, and reject
   *  combining it with a direct h/H0 — except inside a shooting build, where DoShooting
   *  has itself set the trial h. h keeps its default; DoShooting seeds and solves it. */
  class_call(parser_read_double(pfc, "100*theta_s", &param3, &flag3, errmsg), errmsg, errmsg);
  class_test((flag3 == _TRUE_) && (flag1 == _TRUE_ || flag2 == _TRUE_) && !pfc->is_shooting,
             errmsg,
             "In input file, you cannot enter both 100*theta_s and h (or H0), choose one");

  /** - Omega_0_g (photons) and T_cmb */
  class_call(parser_read_double(pfc, "T_cmb", &param1, &flag1, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "Omega_g", &param2, &flag2, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "omega_g", &param3, &flag3, errmsg), errmsg, errmsg);
  class_test(class_at_least_two_of_three(flag1, flag2, flag3),
             errmsg,
             "In input file, you can only enter one of T_cmb, Omega_g or omega_g, choose one");

  if (class_none_of_three(flag1, flag2, flag3)) {
    pba->Omega0_g = (4. * sigma_B / _c_ * pow(pba->T_cmb, 4.)) /
                    (3. * _c_ * _c_ * 1.e10 * pba->h * pba->h / _Mpc_over_m_ / _Mpc_over_m_ / 8. /
                     _PI_ / _G_);
  }
  else {
    if (flag1 == _TRUE_) {
      /** - Omega0_g = rho_g / rho_c0, each of them expressed in \f$ Kg/m/s^2 \f$*/
      /** - rho_g = (4 sigma_B / c) \f$ T^4 \f$*/
      /** - rho_c0 \f$ = 3 c^2 H_0^2 / (8 \pi G) \f$*/
      pba->Omega0_g = (4. * sigma_B / _c_ * pow(param1, 4.)) /
                      (3. * _c_ * _c_ * 1.e10 * pba->h * pba->h / _Mpc_over_m_ / _Mpc_over_m_ / 8. /
                       _PI_ / _G_);
      pba->T_cmb    = param1;
    }

    if (flag2 == _TRUE_) {
      pba->Omega0_g = param2;
      pba->T_cmb = pow(pba->Omega0_g *
                           (3. * _c_ * _c_ * 1.e10 * pba->h * pba->h / _Mpc_over_m_ / _Mpc_over_m_ /
                            8. / _PI_ / _G_) /
                           (4. * sigma_B / _c_),
                       0.25);
    }

    if (flag3 == _TRUE_) {
      pba->Omega0_g = param3 / pba->h / pba->h;
      pba->T_cmb = pow(pba->Omega0_g *
                           (3. * _c_ * _c_ * 1.e10 * pba->h * pba->h / _Mpc_over_m_ / _Mpc_over_m_ /
                            8. / _PI_ / _G_) /
                           (4. * sigma_B / _c_),
                       0.25);
    }
  }

  /** - Omega_0_b (baryons) */
  class_call(parser_read_double(pfc, "Omega_b", &param1, &flag1, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "omega_b", &param2, &flag2, errmsg), errmsg, errmsg);
  class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
             errmsg,
             "In input file, you can only enter one of Omega_b or omega_b, choose one");
  if (flag1 == _TRUE_)
    pba->Omega0_b = param1;
  if (flag2 == _TRUE_)
    pba->Omega0_b = param2 / pba->h / pba->h;

  // NOTE: N_ur / N_eff / Omega_ur / omega_ur parsing has been moved to
  // UltraRelativisticSpecies::CreateAll, which reads pfc directly.

  class_call(parser_read_double(pfc, "ceff2_ur", &param1, &flag1, errmsg), errmsg, errmsg);
  if (flag1 == _TRUE_)
    ppt->three_ceff2_ur = 3. * param1;

  class_call(parser_read_double(pfc, "cvis2_ur", &param1, &flag1, errmsg), errmsg, errmsg);
  if (flag1 == _TRUE_)
    ppt->three_cvis2_ur = 3. * param1;

  class_call(parser_read_double(pfc, "G_eff_ur", &ppt->G_eff_ur, &flag1, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "log10_G_eff_ur", &ppt->G_eff_ur, &flag2, errmsg),
             errmsg,
             errmsg);
  if (flag2 == _TRUE_) {
    class_test(flag1 == _TRUE_,
               errmsg,
               "In input file, you cannot enter both log10_G_eff_ur and G_eff_ur, choose one");
    ppt->G_eff_ur = pow(10.0, ppt->G_eff_ur);
  }

  // Coupled-species cluster (CDM, IDR, IDM_DR, DCDM_DR, IDM_DRMD, IDR_DRMD)
  // Omega-budget parsing has been moved to ReadCoupledOmegaBudget so each
  // species's CreateAll can read its slot from the budget instead of from
  // pba->Omega0_X. Sub-parameters that are not Omega0 (pth->a_idm_dr,
  // thermo block, etc.) remain below and consume the budget.
  class_call(ReadCoupledOmegaBudget(), errmsg, errmsg);

  if (omega_budget_.idm_dr.value_or(0.) > 0.) {
    const double omega0_idr_budget = omega_budget_.idr.value_or(0.);
    class_test(omega0_idr_budget == 0.0,
               errmsg,
               "You have requested interacting DM ith DR, this requires a non-zero density of "
               "interacting DR. Please set either N_idr or xi_idr");

    class_call(parser_read_double(pfc, "a_idm_dr", &param1, &flag1, errmsg), errmsg, errmsg);
    class_call(parser_read_double(pfc, "a_dark", &param2, &flag2, errmsg), errmsg, errmsg);
    class_call(parser_read_double(pfc, "Gamma_0_nadm", &param3, &flag3, errmsg), errmsg, errmsg);
    class_test(class_at_least_two_of_three(flag1, flag2, flag3),
               errmsg,
               "In input file, you can only enter one of a_idm_dr, a_dark or Gamma_0_nadm, choose "
               "one");

    if (flag1 == _TRUE_) {
      pth->a_idm_dr = param1;
      if (input_verbose > 1)
        printf(
            "You passed a_idm_dr = a_dark = %e, this is equivalent to Gamma_0_nadm = %e in the "
            "NADM notation. \n",
            param1,
            param1 * (4. / 3.) * (pba->h * pba->h * omega0_idr_budget));
    }
    else if (flag2 == _TRUE_) {
      pth->a_idm_dr = param2;
      if (input_verbose > 1)
        printf(
            "You passed a_dark = a_idm_dr = %e, this is equivalent to Gamma_0_nadm = %e in the "
            "NADM notation. \n",
            param2,
            param2 * (4. / 3.) * (pba->h * pba->h * omega0_idr_budget));
    }
    else if (flag3 == _TRUE_) {
      pth->a_idm_dr = param3 * (3. / 4.) / (pba->h * pba->h * omega0_idr_budget);
      if (input_verbose > 1)
        printf(
            "You passed Gamma_0_nadm = %e, this is equivalent to a_idm_dr = a_dark = %e in the "
            "ETHOS notation. \n",
            param3,
            pth->a_idm_dr);
    }

    /** - Load the rest of the parameters for idm and idr */

    if (flag3 ==
        _TRUE_) { /* If the user passed Gamma_0_nadm, assume they want nadm parameterisation*/
      pth->nindex_idm_dr = 0;
      ppt->idr_nature    = idr_fluid;
      if (input_verbose > 1)
        printf("NADM requested. Defaulting on nindex_idm_dr = %e and idr_nature = fluid \n",
               pth->nindex_idm_dr);
    }

    else {
      class_read_double_one_of_two("nindex_dark", "nindex_idm_dr", pth->nindex_idm_dr);

      class_call(parser_read_string(pfc, "idr_nature", &string1, &flag1, errmsg), errmsg, errmsg);

      if (flag1 == _TRUE_) {
        if ((strstr(string1, "free_streaming") != nullptr) ||
            (strstr(string1, "Free_Streaming") != nullptr) ||
            (strstr(string1, "Free_streaming") != nullptr) ||
            (strstr(string1, "FREE_STREAMING") != nullptr)) {
          ppt->idr_nature = idr_free_streaming;
        }
        if ((strstr(string1, "fluid") != nullptr) || (strstr(string1, "Fluid") != nullptr) ||
            (strstr(string1, "FLUID") != nullptr)) {
          ppt->idr_nature = idr_fluid;
        }
      }
    }

    class_read_double_one_of_two("m_idm", "m_dm", pth->m_idm);

    class_read_double_one_of_two("b_dark", "b_idr", pth->b_idr);

    /* Read alpha_idm_dr or alpha_dark */

    std::vector<double> alpha_values;
    class_call(readDoubleList(pfc, "alpha_idm_dr", alpha_values, &flag1, errmsg), errmsg, errmsg);

    /* try with the other syntax */
    if (flag1 == _FALSE_) {
      class_call(readDoubleList(pfc, "alpha_dark", alpha_values, &flag1, errmsg), errmsg, errmsg);
    }

    if (flag1 == _TRUE_) {
      entries_read              = static_cast<int>(alpha_values.size());
      ppt->alpha_idm_dr_storage = alpha_values;
      if (entries_read != (ppr->l_max_idr - 1)) {
        ppt->alpha_idm_dr_storage.resize(ppr->l_max_idr - 1,
                                         ppt->alpha_idm_dr_storage[entries_read - 1]);
      }
    }
    else {
      ppt->alpha_idm_dr_storage.assign(ppr->l_max_idr - 1, 1.5);
    }
    ppt->alpha_idm_dr = ppt->alpha_idm_dr_storage.data();

    /* Read alpha_idm_dr or alpha_dark */

    std::vector<double> beta_values;
    class_call(readDoubleList(pfc, "beta_idr", beta_values, &flag1, errmsg), errmsg, errmsg);

    /* try with the other syntax */
    if (flag1 == _FALSE_) {
      class_call(readDoubleList(pfc, "beta_dark", beta_values, &flag1, errmsg), errmsg, errmsg);
    }

    if (flag1 == _TRUE_) {
      entries_read          = static_cast<int>(beta_values.size());
      ppt->beta_idr_storage = beta_values;
      if (entries_read != (ppr->l_max_idr - 1)) {
        ppt->beta_idr_storage.resize(ppr->l_max_idr - 1, ppt->beta_idr_storage[entries_read - 1]);
      }
    }
    else {
      ppt->beta_idr_storage.assign(ppr->l_max_idr - 1, 1.5);
    }
    ppt->beta_idr = ppt->beta_idr_storage.data();
  }

  // Omega_dcdmdr parsing has been moved to ReadCoupledOmegaBudget
  // (Omega_dcdmdr → omega_budget_.dcdmdr); Gamma_dcdm and Omega_ini_dcdm are
  // owned by DCDMSpecies and parsed in DCDM_DR_Species::CreateAll.
  // T_idr and l_max_idr are owned by IDRSpecies and parsed in
  // IDM_DR_IDR_Species::CreateAll.

  /** - Omega_0_k (effective fractional density of curvature) */
  class_read_double("Omega_k", pba->Omega0_k);
  /** - Set curvature parameter K */
  pba->K = -pba->Omega0_k * pow(pba->a_today * pba->H0, 2);
  /** - Set curvature sign */
  if (pba->K > 0.)
    pba->sgnK = 1;
  else if (pba->K < 0.)
    pba->sgnK = -1;

  // DRMD parameter block (z_stop, G_over_aH_drmd_ini, f_idm_drmd, delta_Neff_drmd,
  // and the resulting Omega0_idr_drmd / Omega0_idm_drmd contributions) has been
  // moved to ReadCoupledOmegaBudget.

  /** - Omega_0_lambda (cosmological constant), Omega0_fld (dark energy fluid), Omega0_scf (scalar field) */

  class_call(parser_read_double(pfc, "Omega_Lambda", &param1, &flag1, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "Omega_fld", &param2, &flag2, errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "Omega_scf", &param3, &flag3, errmsg), errmsg, errmsg);

  class_test((flag1 == _TRUE_) && (flag2 == _TRUE_) && ((flag3 == _FALSE_) || (param3 >= 0.)),
             errmsg,
             "In input file, either Omega_Lambda or Omega_fld must be left unspecified, except if "
             "Omega_scf is set and <0.0, in which case the contribution from the scalar field will "
             "be the free parameter.");

  /** - --> (flag3 == _FALSE_) || (param3 >= 0.) explained:
   *  it means that either we have not read Omega_scf so we are ignoring it
   *  (unlike lambda and fld!) OR we have read it, but it had a
   *  positive value and should not be used for filling.
   *  We only record which component closes the budget here; the actual
   *  Omega0_X value is computed in ConstructSpecies after all non-closure
   *  species have been built and their GetOmega0() contributions summed.
   *  The non-closure species' factories read their Omega from pfc directly.
   */

  // Record which species (if any) is the closure species.
  // NOTE: pba->Omega0_lambda/fld/scf are no longer written from this block;
  // FluidSpecies::CreateAll and ScalarFieldSpecies::CreateAll read pfc directly,
  // and Pass 2 of ConstructSpecies hands the closure species its closure value via
  // SpeciesBuildContext::omega0_closure_override.
  if (flag1 == _FALSE_) {
    pba->closure_species = ClosureSpecies::Lambda;
  }
  else if (flag2 == _FALSE_) {
    pba->closure_species = ClosureSpecies::Fluid;
  }
  else if ((flag3 == _TRUE_) && (param3 < 0.)) {
    pba->closure_species = ClosureSpecies::ScalarField;
  }
  else {
    pba->closure_species = ClosureSpecies::None;
  }

  /** - Test that the user have not specified Omega_scf = -1 but left either
      Omega_lambda or Omega_fld unspecified:*/
  class_test(((flag1 == _FALSE_) || (flag2 == _FALSE_)) && ((flag3 == _TRUE_) && (param3 < 0.)),
             errmsg,
             "It looks like you want to fulfil the closure relation sum Omega = 1 using the scalar "
             "field, so you have to specify both Omega_lambda and Omega_fld in the .ini file");

  // Snapshot fluid presence here: later parser blocks reuse flag1/flag2/flag3
  // and param1/param2/param3, and pba->Omega0_fld is no longer written.
  // This flag decides whether to apply the halofit gate (see below).
  const bool fluid_present_pfc = (flag2 == _TRUE_) ||
                                 (pba->closure_species == ClosureSpecies::Fluid);

  /* Fluid physics params (use_ppf, fluid_equation_of_state, w0_fld, wa_fld,
     cs2_fld, Omega_EDE, c_gamma_over_c_fld) are parsed inside
     FluidSpecies::CreateAll directly from pfc; no per-key writes to pba here.

     Scalar field physics params (scf_parameters, scf_tuning_index,
     scf_shooting_parameter, attractor_ic_scf, phi_ini_scf, phi_prime_ini_scf)
     are parsed inside ScalarFieldSpecies::CreateAll directly from pfc; no
     per-key writes to pba here. */

  /** (b) assign values to thermodynamics cosmological parameters */

  /** - primordial helium fraction */
  class_call(parser_read_string(pfc, "YHe", &string1, &flag1, errmsg), errmsg, errmsg);

  if (flag1 == _TRUE_) {
    if ((strstr(string1, "BBN") != nullptr) || (strstr(string1, "bbn") != nullptr)) {
      pth->YHe = _BBN_;
    }
    else {
      class_read_double("YHe", pth->YHe);
    }
  }

  /** - recombination parameters */
  class_call(parser_read_string(pfc, "recombination", &string1, &flag1, errmsg), errmsg, errmsg);

  if (flag1 == _TRUE_) {
    if ((strstr(string1, "HYREC") != nullptr) || (strstr(string1, "hyrec") != nullptr) ||
        (strstr(string1, "HyRec") != nullptr)) {
      pth->recombination = hyrec;
    }
  }

  /** - reionization parametrization */
  class_call(parser_read_string(pfc, "reio_parametrization", &string1, &flag1, errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {
    flag2 = _FALSE_;
    if (strcmp(string1, "reio_none") == 0) {
      pth->reio_parametrization = reio_none;
      flag2                     = _TRUE_;
    }
    if (strcmp(string1, "reio_camb") == 0) {
      pth->reio_parametrization = reio_camb;
      flag2                     = _TRUE_;
    }
    if (strcmp(string1, "reio_bins_tanh") == 0) {
      pth->reio_parametrization = reio_bins_tanh;
      flag2                     = _TRUE_;
    }
    if (strcmp(string1, "reio_half_tanh") == 0) {
      pth->reio_parametrization = reio_half_tanh;
      flag2                     = _TRUE_;
    }
    if (strcmp(string1, "reio_many_tanh") == 0) {
      pth->reio_parametrization = reio_many_tanh;
      flag2                     = _TRUE_;
    }
    if (strcmp(string1, "reio_inter") == 0) {
      pth->reio_parametrization = reio_inter;
      flag2                     = _TRUE_;
    }

    class_test(flag2 == _FALSE_,
               errmsg,
               "could not identify reionization_parametrization value, check that it is one of "
               "'reio_none', 'reio_camb', 'reio_bins_tanh', 'reio_half_tanh', 'reio_many_tanh', "
               "'reio_inter'...");
  }

  /** - reionization parameters if reio_parametrization=reio_camb */
  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {
    class_call(parser_read_double(pfc, "z_reio", &param1, &flag1, errmsg), errmsg, errmsg);
    class_call(parser_read_double(pfc, "tau_reio", &param2, &flag2, errmsg), errmsg, errmsg);
    class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
               errmsg,
               "In input file, you can only enter one of z_reio or tau_reio, choose one");
    if (flag1 == _TRUE_) {
      pth->z_reio        = param1;
      pth->reio_z_or_tau = reio_z;
    }
    if (flag2 == _TRUE_) {
      pth->tau_reio      = param2;
      pth->reio_z_or_tau = reio_tau;
    }

    class_read_double("reionization_exponent", pth->reionization_exponent);
    class_read_double("reionization_width", pth->reionization_width);
    class_read_double("helium_fullreio_redshift", pth->helium_fullreio_redshift);
    class_read_double("helium_fullreio_width", pth->helium_fullreio_width);
  }

  /** - reionization parameters if reio_parametrization=reio_bins_tanh */
  if (pth->reio_parametrization == reio_bins_tanh) {
    class_read_int("binned_reio_num", pth->binned_reio_num);
    class_call(readDoubleList(pfc, "binned_reio_z", pth->binned_reio_z_storage, &flag1, errmsg),
               errmsg,
               errmsg);
    class_test(flag1 == _FALSE_ ||
                   static_cast<int>(pth->binned_reio_z_storage.size()) != pth->binned_reio_num,
               errmsg,
               "Number of entries in binned_reio_z does not match expected number, %d.",
               pth->binned_reio_num);
    pth->binned_reio_z = pth->binned_reio_z_storage.data();
    class_call(readDoubleList(pfc, "binned_reio_xe", pth->binned_reio_xe_storage, &flag1, errmsg),
               errmsg,
               errmsg);
    class_test(flag1 == _FALSE_ ||
                   static_cast<int>(pth->binned_reio_xe_storage.size()) != pth->binned_reio_num,
               errmsg,
               "Number of entries in binned_reio_xe does not match expected number, %d.",
               pth->binned_reio_num);
    pth->binned_reio_xe = pth->binned_reio_xe_storage.data();
    class_read_double("binned_reio_step_sharpness", pth->binned_reio_step_sharpness);
  }

  /** - reionization parameters if reio_parametrization=reio_many_tanh */
  if (pth->reio_parametrization == reio_many_tanh) {
    class_read_int("many_tanh_num", pth->many_tanh_num);
    class_call(readDoubleList(pfc, "many_tanh_z", pth->many_tanh_z_storage, &flag1, errmsg),
               errmsg,
               errmsg);
    class_test(flag1 == _FALSE_ ||
                   static_cast<int>(pth->many_tanh_z_storage.size()) != pth->many_tanh_num,
               errmsg,
               "Number of entries in many_tanh_z does not match expected number, %d.",
               pth->many_tanh_num);
    pth->many_tanh_z = pth->many_tanh_z_storage.data();
    class_call(readDoubleList(pfc, "many_tanh_xe", pth->many_tanh_xe_storage, &flag1, errmsg),
               errmsg,
               errmsg);
    class_test(flag1 == _FALSE_ ||
                   static_cast<int>(pth->many_tanh_xe_storage.size()) != pth->many_tanh_num,
               errmsg,
               "Number of entries in many_tanh_xe does not match expected number, %d.",
               pth->many_tanh_num);
    pth->many_tanh_xe = pth->many_tanh_xe_storage.data();
    class_read_double("many_tanh_width", pth->many_tanh_width);
  }

  /** - reionization parameters if reio_parametrization=reio_many_tanh */
  if (pth->reio_parametrization == reio_inter) {
    class_read_int("reio_inter_num", pth->reio_inter_num);
    class_call(readDoubleList(pfc, "reio_inter_z", pth->reio_inter_z_storage, &flag1, errmsg),
               errmsg,
               errmsg);
    class_test(flag1 == _FALSE_ ||
                   static_cast<int>(pth->reio_inter_z_storage.size()) != pth->reio_inter_num,
               errmsg,
               "Number of entries in reio_inter_z does not match expected number, %d.",
               pth->reio_inter_num);
    pth->reio_inter_z = pth->reio_inter_z_storage.data();
    class_call(readDoubleList(pfc, "reio_inter_xe", pth->reio_inter_xe_storage, &flag1, errmsg),
               errmsg,
               errmsg);
    class_test(flag1 == _FALSE_ ||
                   static_cast<int>(pth->reio_inter_xe_storage.size()) != pth->reio_inter_num,
               errmsg,
               "Number of entries in reio_inter_xe does not match expected number, %d.",
               pth->reio_inter_num);
    pth->reio_inter_xe = pth->reio_inter_xe_storage.data();
  }

  /** - energy injection parameters from CDM annihilation/decay */

  class_read_double("annihilation", pth->annihilation);

  if (pth->annihilation > 0.) {
    class_read_double("annihilation_variation", pth->annihilation_variation);
    class_read_double("annihilation_z", pth->annihilation_z);
    class_read_double("annihilation_zmax", pth->annihilation_zmax);
    class_read_double("annihilation_zmin", pth->annihilation_zmin);
    class_read_double("annihilation_f_halo", pth->annihilation_f_halo);
    class_read_double("annihilation_z_halo", pth->annihilation_z_halo);

    class_call(parser_read_string(pfc, "on the spot", &(string1), &(flag1), errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {
      if ((strstr(string1, "y") != nullptr) || (strstr(string1, "Y") != nullptr)) {
        pth->has_on_the_spot = _TRUE_;
      }
      else {
        if ((strstr(string1, "n") != nullptr) || (strstr(string1, "N") != nullptr)) {
          pth->has_on_the_spot = _FALSE_;
        }
        else {
          class_stop(errmsg, "incomprehensible input '%s' for the field 'on the spot'", string1);
        }
      }
    }
  }

  class_read_double("decay", pth->decay);

  class_call(parser_read_string(pfc, "compute damping scale", &(string1), &(flag1), errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {
    if ((strstr(string1, "y") != nullptr) || (strstr(string1, "Y") != nullptr)) {
      pth->compute_damping_scale = _TRUE_;
    }
    else {
      if ((strstr(string1, "n") != nullptr) || (strstr(string1, "N") != nullptr)) {
        pth->compute_damping_scale = _FALSE_;
      }
      else {
        class_stop(errmsg,
                   "incomprehensible input '%s' for the field 'compute damping scale'",
                   string1);
      }
    }
  }

  /** (c) define which perturbations and sources should be computed, and down to which scale */

  ppt->has_perturbations = _FALSE_;
  ppt->has_cls           = _FALSE_;

  class_call(parser_read_string(pfc, "output", &string1, &flag1, errmsg), errmsg, errmsg);

  if (flag1 == _TRUE_) {
    if ((strstr(string1, "tCl") != nullptr) || (strstr(string1, "TCl") != nullptr) ||
        (strstr(string1, "TCL") != nullptr)) {
      ppt->has_cl_cmb_temperature = _TRUE_;
      ppt->has_perturbations      = _TRUE_;
      ppt->has_cls                = _TRUE_;
    }

    if ((strstr(string1, "pCl") != nullptr) || (strstr(string1, "PCl") != nullptr) ||
        (strstr(string1, "PCL") != nullptr)) {
      ppt->has_cl_cmb_polarization = _TRUE_;
      ppt->has_perturbations       = _TRUE_;
      ppt->has_cls                 = _TRUE_;
    }

    if ((strstr(string1, "lCl") != nullptr) || (strstr(string1, "LCl") != nullptr) ||
        (strstr(string1, "LCL") != nullptr)) {
      ppt->has_cl_cmb_lensing_potential = _TRUE_;
      ppt->has_perturbations            = _TRUE_;
      ppt->has_cls                      = _TRUE_;
    }

    if ((strstr(string1, "nCl") != nullptr) || (strstr(string1, "NCl") != nullptr) ||
        (strstr(string1, "NCL") != nullptr) || (strstr(string1, "dCl") != nullptr) ||
        (strstr(string1, "DCl") != nullptr) || (strstr(string1, "DCL") != nullptr)) {
      ppt->has_cl_number_count = _TRUE_;
      ppt->has_perturbations   = _TRUE_;
      ppt->has_cls             = _TRUE_;
    }

    if ((strstr(string1, "sCl") != nullptr) || (strstr(string1, "SCl") != nullptr) ||
        (strstr(string1, "SCL") != nullptr)) {
      ppt->has_cl_lensing_potential = _TRUE_;
      ppt->has_perturbations        = _TRUE_;
      ppt->has_cls                  = _TRUE_;
    }

    if ((strstr(string1, "mPk") != nullptr) || (strstr(string1, "MPk") != nullptr) ||
        (strstr(string1, "MPK") != nullptr)) {
      ppt->has_pk_matter     = _TRUE_;
      ppt->has_perturbations = _TRUE_;

      /*if (pba->Omega0_ncdm_tot != 0.0){
        class_call(parser_read_string(pfc,"pk_only_cdm_bar",&string1,&flag1,errmsg),
        errmsg,
        errmsg);
        if (flag1 == _TRUE_){
        if((strstr(string1,"y") != nullptr) || (strstr(string1,"Y") != nullptr)){
        ppt->pk_only_cdm_bar = _TRUE_;
        }
        else {
        ppt->pk_only_cdm_bar = _FALSE_;
        }
        }
        }*/
    }

    if ((strstr(string1, "mTk") != nullptr) || (strstr(string1, "MTk") != nullptr) ||
        (strstr(string1, "MTK") != nullptr) || (strstr(string1, "dTk") != nullptr) ||
        (strstr(string1, "DTk") != nullptr) || (strstr(string1, "DTK") != nullptr)) {
      ppt->has_density_transfers = _TRUE_;
      ppt->has_perturbations     = _TRUE_;
    }

    if ((strstr(string1, "vTk") != nullptr) || (strstr(string1, "VTk") != nullptr) ||
        (strstr(string1, "VTK") != nullptr)) {
      ppt->has_velocity_transfers = _TRUE_;
      ppt->has_perturbations      = _TRUE_;
    }
  }

  if (ppt->has_density_transfers == _TRUE_) {
    class_call(parser_read_string(pfc, "extra metric transfer functions", &string1, &flag1, errmsg),
               errmsg,
               errmsg);

    if ((flag1 == _TRUE_) &&
        ((strstr(string1, "y") != nullptr) || (strstr(string1, "y") != nullptr))) {
      ppt->has_metricpotential_transfers = _TRUE_;
    }
  }

  if (ppt->has_cl_cmb_temperature == _TRUE_) {
    class_call(parser_read_string(pfc, "temperature contributions", &string1, &flag1, errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {
      ppt->switch_sw   = 0;
      ppt->switch_eisw = 0;
      ppt->switch_lisw = 0;
      ppt->switch_dop  = 0;
      ppt->switch_pol  = 0;

      if ((strstr(string1, "tsw") != nullptr) || (strstr(string1, "TSW") != nullptr))
        ppt->switch_sw = 1;
      if ((strstr(string1, "eisw") != nullptr) || (strstr(string1, "EISW") != nullptr))
        ppt->switch_eisw = 1;
      if ((strstr(string1, "lisw") != nullptr) || (strstr(string1, "LISW") != nullptr))
        ppt->switch_lisw = 1;
      if ((strstr(string1, "dop") != nullptr) || (strstr(string1, "Dop") != nullptr))
        ppt->switch_dop = 1;
      if ((strstr(string1, "pol") != nullptr) || (strstr(string1, "Pol") != nullptr))
        ppt->switch_pol = 1;

      class_test((ppt->switch_sw == 0) && (ppt->switch_eisw == 0) && (ppt->switch_lisw == 0) &&
                     (ppt->switch_dop == 0) && (ppt->switch_pol == 0),
                 errmsg,
                 "In the field 'output', you selected CMB temperature, but in the field "
                 "'temperature contributions', you removed all contributions");

      class_read_double("early/late isw redshift", ppt->eisw_lisw_split_z);
    }
  }

  if (ppt->has_cl_number_count == _TRUE_) {
    class_call(parser_read_string(pfc, "number count contributions", &string1, &flag1, errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {
      if (strstr(string1, "density") != nullptr)
        ppt->has_nc_density = _TRUE_;
      if (strstr(string1, "rsd") != nullptr)
        ppt->has_nc_rsd = _TRUE_;
      if (strstr(string1, "lensing") != nullptr)
        ppt->has_nc_lens = _TRUE_;
      if (strstr(string1, "gr") != nullptr)
        ppt->has_nc_gr = _TRUE_;

      class_test((ppt->has_nc_density == _FALSE_) && (ppt->has_nc_rsd == _FALSE_) &&
                     (ppt->has_nc_lens == _FALSE_) && (ppt->has_nc_gr == _FALSE_),
                 errmsg,
                 "In the field 'output', you selected number count Cl's, but in the field 'number "
                 "count contributions', you removed all contributions");
    }

    else {
      /* default: only the density contribution */
      ppt->has_nc_density = _TRUE_;
    }
  }

  if (ppt->has_perturbations == _TRUE_) {
    /* perturbed recombination */
    class_call(parser_read_string(pfc, "perturbed recombination", &(string1), &(flag1), errmsg),
               errmsg,
               errmsg);

    if ((flag1 == _TRUE_) &&
        ((strstr(string1, "y") != nullptr) || (strstr(string1, "Y") != nullptr))) {
      ppt->has_perturbed_recombination = _TRUE_;
    }

    /* modes */
    class_call(parser_read_string(pfc, "modes", &string1, &flag1, errmsg), errmsg, errmsg);

    if (flag1 == _TRUE_) {
      /* if no modes are specified, the default is has_scalars=_TRUE_;
         but if they are specified we should reset has_scalars to _FALSE_ before reading */
      ppt->has_scalars = _FALSE_;

      if ((strstr(string1, "s") != nullptr) || (strstr(string1, "S") != nullptr))
        ppt->has_scalars = _TRUE_;

      if ((strstr(string1, "v") != nullptr) || (strstr(string1, "V") != nullptr))
        ppt->has_vectors = _TRUE_;

      if ((strstr(string1, "t") != nullptr) || (strstr(string1, "T") != nullptr))
        ppt->has_tensors = _TRUE_;

      class_test(class_none_of_three(ppt->has_scalars, ppt->has_vectors, ppt->has_tensors),
                 errmsg,
                 "You wrote: modes='%s'. Could not identify any of the modes ('s', 'v', 't') in "
                 "such input",
                 string1);
    }

    if (ppt->has_scalars == _TRUE_) {
      class_call(parser_read_string(pfc, "ic", &string1, &flag1, errmsg), errmsg, errmsg);

      if (flag1 == _TRUE_) {
        /* if no initial conditions are specified, the default is has_ad=_TRUE_;
           but if they are specified we should reset has_ad to _FALSE_ before reading */
        ppt->has_ad = _FALSE_;

        if ((strstr(string1, "ad") != nullptr) || (strstr(string1, "AD") != nullptr))
          ppt->has_ad = _TRUE_;

        if ((strstr(string1, "bi") != nullptr) || (strstr(string1, "BI") != nullptr))
          ppt->has_bi = _TRUE_;

        if ((strstr(string1, "cdi") != nullptr) || (strstr(string1, "CDI") != nullptr))
          ppt->has_cdi = _TRUE_;

        if ((strstr(string1, "nid") != nullptr) || (strstr(string1, "NID") != nullptr))
          ppt->has_nid = _TRUE_;

        if ((strstr(string1, "niv") != nullptr) || (strstr(string1, "NIV") != nullptr))
          ppt->has_niv = _TRUE_;

        class_test(ppt->has_ad == _FALSE_ && ppt->has_bi == _FALSE_ && ppt->has_cdi == _FALSE_ &&
                       ppt->has_nid == _FALSE_ && ppt->has_niv == _FALSE_,
                   errmsg,
                   "You wrote: ic='%s'. Could not identify any of the initial conditions ('ad', "
                   "'bi', 'cdi', 'nid', 'niv') in such input",
                   string1);
      }
    }

    else {
      class_test(ppt->has_cl_cmb_lensing_potential == _TRUE_,
                 errmsg,
                 "Inconsistency: you want C_l's for cmb lensing potential, but no scalar modes\n");

      class_test(ppt->has_pk_matter == _TRUE_,
                 errmsg,
                 "Inconsistency: you want P(k) of matter, but no scalar modes\n");
    }

    if (ppt->has_vectors == _TRUE_) {
      class_test(
          (ppt->has_cl_cmb_temperature == _FALSE_) && (ppt->has_cl_cmb_polarization == _FALSE_),
          errmsg,
          "inconsistent input: you asked for vectors, so you should have at least one non-zero "
          "tensor source type (temperature or polarization). Please adjust your input.");
    }

    if (ppt->has_tensors == _TRUE_) {
      class_test(
          (ppt->has_cl_cmb_temperature == _FALSE_) && (ppt->has_cl_cmb_polarization == _FALSE_),
          errmsg,
          "inconsistent input: you asked for tensors, so you should have at least one non-zero "
          "tensor source type (temperature or polarization). Please adjust your input.");
    }
  }

  /** (d) define the primordial spectrum */

  class_call(parser_read_string(pfc, "P_k_ini type", &string1, &flag1, errmsg), errmsg, errmsg);

  if (flag1 == _TRUE_) {
    flag2 = _FALSE_;
    if (strcmp(string1, "analytic_Pk") == 0) {
      ppm->primordial_spec_type = analytic_Pk;
      flag2                     = _TRUE_;
    }
    if (strcmp(string1, "two_scales") == 0) {
      ppm->primordial_spec_type = two_scales;
      flag2                     = _TRUE_;
    }
    if (strcmp(string1, "inflation_V") == 0) {
      ppm->primordial_spec_type = inflation_V;
      flag2                     = _TRUE_;
    }
    if (strcmp(string1, "inflation_H") == 0) {
      ppm->primordial_spec_type = inflation_H;
      flag2                     = _TRUE_;
    }
    if (strcmp(string1, "inflation_V_end") == 0) {
      ppm->primordial_spec_type = inflation_V_end;
      flag2                     = _TRUE_;
    }
    if (strcmp(string1, "external_Pk") == 0) {
      ppm->primordial_spec_type = external_Pk;
      flag2                     = _TRUE_;
    }
    class_test(flag2 == _FALSE_,
               errmsg,
               "could not identify primordial spectrum type, check that it is one of "
               "'analytic_pk', 'two_scales', 'inflation_V', 'inflation_H', 'external_Pk'...");
  }

  class_read_double("k_pivot", ppm->k_pivot);

  if (ppm->primordial_spec_type == two_scales) {
    class_read_double("k1", k1);
    class_read_double("k2", k2);
    class_test(k1 <= 0., errmsg, "enter strictly positive scale k1");
    class_test(k2 <= 0., errmsg, "enter strictly positive scale k2");

    if (ppt->has_scalars == _TRUE_) {
      class_read_double("P_{RR}^1", prr1);
      class_read_double("P_{RR}^2", prr2);
      class_test(prr1 <= 0., errmsg, "enter strictly positive scale P_{RR}^1");
      class_test(prr2 <= 0., errmsg, "enter strictly positive scale P_{RR}^2");

      ppm->n_s = log(prr2 / prr1) / log(k2 / k1) + 1.;
      ppm->A_s = prr1 * exp((ppm->n_s - 1.) * log(ppm->k_pivot / k1));

      if ((ppt->has_bi == _TRUE_) || (ppt->has_cdi == _TRUE_) || (ppt->has_nid == _TRUE_) ||
          (ppt->has_niv == _TRUE_)) {
        class_read_double("P_{II}^1", pii1);
        class_read_double("P_{II}^2", pii2);
        class_read_double("P_{RI}^1", pri1);
        class_read_double("|P_{RI}^2|", pri2);

        class_test(pii1 <= 0.,
                   errmsg,
                   "since you request iso modes, you should have P_{ii}^1 strictly positive");
        class_test(pii2 < 0.,
                   errmsg,
                   "since you request iso modes, you should have P_{ii}^2 positive or eventually "
                   "null");
        class_test(pri2 < 0.,
                   errmsg,
                   "by definition, you should have |P_{ri}^2| positive or eventually null");

        flag1 = _FALSE_;

        class_call(parser_read_string(pfc, "special iso", &string1, &flag1, errmsg),
                   errmsg,
                   errmsg);

        /* axion case, only one iso parameter: piir1  */
        if ((flag1 == _TRUE_) && (strstr(string1, "axion") != nullptr)) {
          n_iso = 1.;
          n_cor = 0.;
          c_cor = 0.;
        }
        /* curvaton case, only one iso parameter: piir1  */
        else if ((flag1 == _TRUE_) && (strstr(string1, "anticurvaton") != nullptr)) {
          n_iso = ppm->n_s;
          n_cor = 0.;
          c_cor = 1.;
        }
        /* inverted-correlation-curvaton case, only one iso parameter: piir1  */
        else if ((flag1 == _TRUE_) && (strstr(string1, "curvaton") != nullptr)) {
          n_iso = ppm->n_s;
          n_cor = 0.;
          c_cor = -1.;
        }
        /* general case, but if pii2 or pri2=0 the code interprets it
           as a request for n_iso=n_ad or n_cor=0 respectively */
        else {
          if (pii2 == 0.) {
            n_iso = ppm->n_s;
          }
          else {
            class_test((pii1 == 0.) || (pii2 == 0.) || (pii1 * pii2 < 0.),
                       errmsg,
                       "should NEVER happen");
            n_iso = log(pii2 / pii1) / log(k2 / k1) + 1.;
          }
          class_test(pri1 == 0,
                     errmsg,
                     "the general isocurvature case requires a non-zero P_{RI}^1");
          if (pri2 == 0.) {
            n_cor = 0.;
          }
          else {
            class_test((pri1 == 0.) || (pri2 <= 0.) || (pii1 * pii2 < 0),
                       errmsg,
                       "should NEVER happen");
            n_cor = log(pri2 / fabs(pri1)) / log(k2 / k1) - 0.5 * (ppm->n_s + n_iso - 2.);
          }
          class_test((pii1 * prr1 <= 0.), errmsg, "should NEVER happen");
          class_test(fabs(pri1) / sqrt(pii1 * prr1) > 1,
                     errmsg,
                     "too large ad-iso cross-correlation in k1");
          class_test(fabs(pri1) / sqrt(pii1 * prr1) * exp(n_cor * log(k2 / k1)) > 1,
                     errmsg,
                     "too large ad-iso cross-correlation in k2");
          c_cor = -pri1 / sqrt(pii1 * prr1) * exp(n_cor * log(ppm->k_pivot / k1));
        }
        /* formula for f_iso valid in all cases */
        class_test((pii1 == 0.) || (prr1 == 0.) || (pii1 * prr1 < 0.),
                   errmsg,
                   "should NEVER happen");
        f_iso = sqrt(pii1 / prr1) * exp(0.5 * (n_iso - ppm->n_s) * log(ppm->k_pivot / k1));
      }

      if (ppt->has_bi == _TRUE_) {
        ppm->f_bi    = f_iso;
        ppm->n_bi    = n_iso;
        ppm->c_ad_bi = c_cor;
        ppm->n_ad_bi = n_cor;
      }

      if (ppt->has_cdi == _TRUE_) {
        ppm->f_cdi    = f_iso;
        ppm->n_cdi    = n_iso;
        ppm->c_ad_cdi = c_cor;
        ppm->n_ad_cdi = n_cor;
      }

      if (ppt->has_nid == _TRUE_) {
        ppm->f_nid    = f_iso;
        ppm->n_nid    = n_iso;
        ppm->c_ad_nid = c_cor;
        ppm->n_ad_nid = n_cor;
      }

      if (ppt->has_niv == _TRUE_) {
        ppm->f_niv    = f_iso;
        ppm->n_niv    = n_iso;
        ppm->c_ad_niv = c_cor;
        ppm->n_ad_niv = n_cor;
      }
    }

    ppm->primordial_spec_type = analytic_Pk;
  }

  else if (ppm->primordial_spec_type == analytic_Pk) {
    if (ppt->has_scalars == _TRUE_) {
      int flag4;
      double param4;
      class_call(parser_read_double(pfc, "A_s", &param1, &flag1, errmsg), errmsg, errmsg);
      class_call(parser_read_double(pfc, "ln10^{10}A_s", &param2, &flag2, errmsg), errmsg, errmsg);
      class_call(parser_read_double(pfc, "sigma8", &param3, &flag3, errmsg), errmsg, errmsg);
      class_call(parser_read_double(pfc, "S8", &param4, &flag4, errmsg), errmsg, errmsg);
      class_test(class_at_least_two_of_four(flag1, flag2, flag3, flag4),
                 errmsg,
                 "In input file, you can only enter one of A_s, ln10^{10}A_s, sigma8 and S8, "
                 "choose one");
      if (flag1 == _TRUE_)
        ppm->A_s = param1;
      else if (flag2 == _TRUE_)
        ppm->A_s = exp(param2) * 1.e-10;
      else if (flag3 == _TRUE_) {
        ppm->sigma8 = param3;
        class_test(param3 < 0., errmsg, "sigma8 should be non-negative");
      }
      else if (flag4 == _TRUE_) {
        // CDM read from the resolved budget (pba->Omega0_cdm is no longer
        // written by the coupled-species parser block).
        const double Omega0_cdm_for_S8 = omega_budget_.cdm.value_or(0.);
        ppm->sigma8 = param4 / pow((pba->Omega0_b + Omega0_cdm_for_S8) / 0.3, 0.5);
        class_test(param4 < 0., errmsg, "S8 should be non-negative");
      }

      if (ppt->has_ad == _TRUE_) {
        class_read_double("n_s", ppm->n_s);
        class_read_double("alpha_s", ppm->alpha_s);
      }

      if (ppt->has_bi == _TRUE_) {
        class_read_double("f_bi", ppm->f_bi);
        class_read_double("n_bi", ppm->n_bi);
        class_read_double("alpha_bi", ppm->alpha_bi);
      }

      if (ppt->has_cdi == _TRUE_) {
        class_read_double("f_cdi", ppm->f_cdi);
        class_read_double("n_cdi", ppm->n_cdi);
        class_read_double("alpha_cdi", ppm->alpha_cdi);
      }

      if (ppt->has_nid == _TRUE_) {
        class_read_double("f_nid", ppm->f_nid);
        class_read_double("n_nid", ppm->n_nid);
        class_read_double("alpha_nid", ppm->alpha_nid);
      }

      if (ppt->has_niv == _TRUE_) {
        class_read_double("f_niv", ppm->f_niv);
        class_read_double("n_niv", ppm->n_niv);
        class_read_double("alpha_niv", ppm->alpha_niv);
      }

      if ((ppt->has_ad == _TRUE_) && (ppt->has_bi == _TRUE_)) {
        class_read_double_one_of_two("c_ad_bi", "c_bi_ad", ppm->c_ad_bi);
        class_read_double_one_of_two("n_ad_bi", "n_bi_ad", ppm->n_ad_bi);
        class_read_double_one_of_two("alpha_ad_bi", "alpha_bi_ad", ppm->alpha_ad_bi);
      }

      if ((ppt->has_ad == _TRUE_) && (ppt->has_cdi == _TRUE_)) {
        class_read_double_one_of_two("c_ad_cdi", "c_cdi_ad", ppm->c_ad_cdi);
        class_read_double_one_of_two("n_ad_cdi", "n_cdi_ad", ppm->n_ad_cdi);
        class_read_double_one_of_two("alpha_ad_cdi", "alpha_cdi_ad", ppm->alpha_ad_cdi);
      }

      if ((ppt->has_ad == _TRUE_) && (ppt->has_nid == _TRUE_)) {
        class_read_double_one_of_two("c_ad_nid", "c_nid_ad", ppm->c_ad_nid);
        class_read_double_one_of_two("n_ad_nid", "n_nid_ad", ppm->n_ad_nid);
        class_read_double_one_of_two("alpha_ad_nid", "alpha_nid_ad", ppm->alpha_ad_nid);
      }

      if ((ppt->has_ad == _TRUE_) && (ppt->has_niv == _TRUE_)) {
        class_read_double_one_of_two("c_ad_niv", "c_niv_ad", ppm->c_ad_niv);
        class_read_double_one_of_two("n_ad_niv", "n_niv_ad", ppm->n_ad_niv);
        class_read_double_one_of_two("alpha_ad_niv", "alpha_niv_ad", ppm->alpha_ad_niv);
      }

      if ((ppt->has_bi == _TRUE_) && (ppt->has_cdi == _TRUE_)) {
        class_read_double_one_of_two("c_bi_cdi", "c_cdi_bi", ppm->c_bi_cdi);
        class_read_double_one_of_two("n_bi_cdi", "n_cdi_bi", ppm->n_bi_cdi);
        class_read_double_one_of_two("alpha_bi_cdi", "alpha_cdi_bi", ppm->alpha_bi_cdi);
      }

      if ((ppt->has_bi == _TRUE_) && (ppt->has_nid == _TRUE_)) {
        class_read_double_one_of_two("c_bi_nid", "c_nid_bi", ppm->c_bi_nid);
        class_read_double_one_of_two("n_bi_nid", "n_nid_bi", ppm->n_bi_nid);
        class_read_double_one_of_two("alpha_bi_nid", "alpha_nid_bi", ppm->alpha_bi_nid);
      }

      if ((ppt->has_bi == _TRUE_) && (ppt->has_niv == _TRUE_)) {
        class_read_double_one_of_two("c_bi_niv", "c_niv_bi", ppm->c_bi_niv);
        class_read_double_one_of_two("n_bi_niv", "n_niv_bi", ppm->n_bi_niv);
        class_read_double_one_of_two("alpha_bi_niv", "alpha_niv_bi", ppm->alpha_bi_niv);
      }

      if ((ppt->has_cdi == _TRUE_) && (ppt->has_nid == _TRUE_)) {
        class_read_double_one_of_two("c_cdi_nid", "c_nid_cdi", ppm->c_cdi_nid);
        class_read_double_one_of_two("n_cdi_nid", "n_nid_cdi", ppm->n_cdi_nid);
        class_read_double_one_of_two("alpha_cdi_nid", "alpha_nid_cdi", ppm->alpha_cdi_nid);
      }

      if ((ppt->has_cdi == _TRUE_) && (ppt->has_niv == _TRUE_)) {
        class_read_double_one_of_two("c_cdi_niv", "c_niv_cdi", ppm->c_cdi_niv);
        class_read_double_one_of_two("n_cdi_niv", "n_niv_cdi", ppm->n_cdi_niv);
        class_read_double_one_of_two("alpha_cdi_niv", "alpha_niv_cdi", ppm->alpha_cdi_niv);
      }

      if ((ppt->has_nid == _TRUE_) && (ppt->has_niv == _TRUE_)) {
        class_read_double_one_of_two("c_nid_niv", "c_niv_nid", ppm->c_nid_niv);
        class_read_double_one_of_two("n_nid_niv", "n_niv_nid", ppm->n_nid_niv);
        class_read_double_one_of_two("alpha_nid_niv", "alpha_niv_nid", ppm->alpha_nid_niv);
      }
    }

    if (ppt->has_tensors == _TRUE_) {
      class_read_double("r", ppm->r);

      if (ppt->has_scalars == _FALSE_) {
        class_read_double("A_s", ppm->A_s);
      }

      if (ppm->r <= 0) {
        ppt->has_tensors = _FALSE_;
      }
      else {
        class_call(parser_read_string(pfc, "n_t", &string1, &flag1, errmsg), errmsg, errmsg);

        if ((flag1 == _TRUE_) &&
            !((strstr(string1, "SCC") != nullptr) || (strstr(string1, "scc") != nullptr))) {
          class_read_double("n_t", ppm->n_t);
        }
        else {
          /* enforce single slow-roll self-consistency condition (order 2 in slow-roll) */
          ppm->n_t = -ppm->r / 8. * (2. - ppm->r / 8. - ppm->n_s);
        }

        class_call(parser_read_string(pfc, "alpha_t", &string1, &flag1, errmsg), errmsg, errmsg);

        if ((flag1 == _TRUE_) &&
            !((strstr(string1, "SCC") != nullptr) || (strstr(string1, "scc") != nullptr))) {
          class_read_double("alpha_t", ppm->alpha_t);
        }
        else {
          /* enforce single slow-roll self-consistency condition (order 2 in slow-roll) */
          ppm->alpha_t = ppm->r / 8. * (ppm->r / 8. + ppm->n_s - 1.);
        }
      }
    }
  }

  else if ((ppm->primordial_spec_type == inflation_V) ||
           (ppm->primordial_spec_type == inflation_H)) {
    double R0, R1, R2, R3, R4;
    double PSR0, PSR1, PSR2, PSR3, PSR4;
    double HSR0, HSR1, HSR2, HSR3, HSR4;

    if (ppm->primordial_spec_type == inflation_V) {
      class_call(parser_read_string(pfc, "potential", &string1, &flag1, errmsg), errmsg, errmsg);

      /* only polynomial coded so far: no need to interpret string1 **/

      class_call(parser_read_string(pfc, "PSR_0", &string1, &flag1, errmsg), errmsg, errmsg);

      if (flag1 == _TRUE_) {
        PSR0 = 0.;
        PSR1 = 0.;
        PSR2 = 0.;
        PSR3 = 0.;
        PSR4 = 0.;

        class_read_double("PSR_0", PSR0);
        class_read_double("PSR_1", PSR1);
        class_read_double("PSR_2", PSR2);
        class_read_double("PSR_3", PSR3);
        class_read_double("PSR_4", PSR4);

        class_test(PSR0 <= 0.,
                   errmsg,
                   "inconsistent parametrization of polynomial inflation potential");
        class_test(PSR1 <= 0.,
                   errmsg,
                   "inconsistent parametrization of polynomial inflation potential");

        R0 = PSR0;
        R1 = PSR1 * 16. * _PI_;
        R2 = PSR2 * 8. * _PI_;
        R3 = PSR3 * pow(8. * _PI_, 2);
        R4 = PSR4 * pow(8. * _PI_, 3);

        ppm->V0 = R0 * R1 * 3. / 128. / _PI_;
        ppm->V1 = -sqrt(R1) * ppm->V0;
        ppm->V2 = R2 * ppm->V0;
        ppm->V3 = R3 * ppm->V0 * ppm->V0 / ppm->V1;
        ppm->V4 = R4 * ppm->V0 / R1;
      }

      else {
        class_call(parser_read_string(pfc, "R_0", &string1, &flag1, errmsg), errmsg, errmsg);

        if (flag1 == _TRUE_) {
          R0 = 0.;
          R1 = 0.;
          R2 = 0.;
          R3 = 0.;
          R4 = 0.;

          class_read_double("R_0", R0);
          class_read_double("R_1", R1);
          class_read_double("R_2", R2);
          class_read_double("R_3", R3);
          class_read_double("R_4", R4);

          class_test(R0 <= 0.,
                     errmsg,
                     "inconsistent parametrization of polynomial inflation potential");
          class_test(R1 <= 0.,
                     errmsg,
                     "inconsistent parametrization of polynomial inflation potential");

          ppm->V0 = R0 * R1 * 3. / 128. / _PI_;
          ppm->V1 = -sqrt(R1) * ppm->V0;
          ppm->V2 = R2 * ppm->V0;
          ppm->V3 = R3 * ppm->V0 * ppm->V0 / ppm->V1;
          ppm->V4 = R4 * ppm->V0 / R1;
        }

        else {
          class_read_double("V_0", ppm->V0);
          class_read_double("V_1", ppm->V1);
          class_read_double("V_2", ppm->V2);
          class_read_double("V_3", ppm->V3);
          class_read_double("V_4", ppm->V4);
        }
      }
    }

    else {
      class_call(parser_read_string(pfc, "HSR_0", &string1, &flag1, errmsg), errmsg, errmsg);

      if (flag1 == _TRUE_) {
        HSR0 = 0.;
        HSR1 = 0.;
        HSR2 = 0.;
        HSR3 = 0.;
        HSR4 = 0.;

        class_read_double("HSR_0", HSR0);
        class_read_double("HSR_1", HSR1);
        class_read_double("HSR_2", HSR2);
        class_read_double("HSR_3", HSR3);
        class_read_double("HSR_4", HSR4);

        ppm->H0 = sqrt(HSR0 * HSR1 * _PI_);
        ppm->H1 = -sqrt(4. * _PI_ * HSR1) * ppm->H0;
        ppm->H2 = 4. * _PI_ * HSR2 * ppm->H0;
        ppm->H3 = 4. * _PI_ * HSR3 * ppm->H0 * ppm->H0 / ppm->H1;
        ppm->H4 = 4. * _PI_ * HSR4 * ppm->H0 * ppm->H0 * ppm->H0 / ppm->H1 / ppm->H1;
      }
      else {
        class_read_double("H_0", ppm->H0);
        class_read_double("H_1", ppm->H1);
        class_read_double("H_2", ppm->H2);
        class_read_double("H_3", ppm->H3);
        class_read_double("H_4", ppm->H4);
      }

      class_test(ppm->H0 <= 0.,
                 errmsg,
                 "inconsistent parametrization of polynomial inflation potential");
    }
  }

  else if (ppm->primordial_spec_type == inflation_V_end) {
    class_call(parser_read_string(pfc, "full_potential", &string1, &flag1, errmsg), errmsg, errmsg);

    if (flag1 == _TRUE_) {
      if (strcmp(string1, "polynomial") == 0) {
        ppm->potential = polynomial;
      }
      else if (strcmp(string1, "higgs_inflation") == 0) {
        ppm->potential = higgs_inflation;
      }
      else {
        class_stop(errmsg,
                   "did not recognize input parameter 'potential': should be one of 'polynomial' "
                   "or 'higgs_inflation'");
      }
    }

    class_read_double("phi_end", ppm->phi_end);
    class_read_double("Vparam0", ppm->V0);
    class_read_double("Vparam1", ppm->V1);
    class_read_double("Vparam2", ppm->V2);
    class_read_double("Vparam3", ppm->V3);
    class_read_double("Vparam4", ppm->V4);

    class_call(parser_read_string(pfc, "ln_aH_ratio", &string1, &flag1, errmsg), errmsg, errmsg);

    class_call(parser_read_string(pfc, "N_star", &string2, &flag2, errmsg), errmsg, errmsg);

    class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
               errmsg,
               "In input file, you can only enter one of ln_aH_ratio or N_star, the two are not "
               "compatible");

    if (flag1 == _TRUE_) {
      if ((strstr(string1, "auto") != nullptr) || (strstr(string1, "AUTO") != nullptr)) {
        ppm->phi_pivot_method = ln_aH_ratio_auto;
      }
      else {
        ppm->phi_pivot_method = ln_aH_ratio;
        class_read_double("ln_aH_ratio", ppm->phi_pivot_target);
      }
    }

    if (flag2 == _TRUE_) {
      ppm->phi_pivot_method = N_star;
      class_read_double("N_star", ppm->phi_pivot_target);
    }

    class_call(parser_read_string(pfc, "inflation_behavior", &string1, &flag1, errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {
      if (strstr(string1, "numerical") != nullptr) {
        ppm->behavior = numerical;
      }
      else if (strstr(string1, "analytical") != nullptr) {
        ppm->behavior = analytical;
      }
      else {
        class_stop(errmsg, "Your entry for 'inflation behavior' could not be understood");
      }
    }
  }

  else if (ppm->primordial_spec_type == external_Pk) {
    class_call(parser_read_string(pfc, "command", &(string1), &(flag1), errmsg), errmsg, errmsg);
    class_test(strlen(string1) == 0, errmsg, "You omitted to write a command for the external Pk");

    ppm->command = string1;
    class_read_double("custom1", ppm->custom1);
    class_read_double("custom2", ppm->custom2);
    class_read_double("custom3", ppm->custom3);
    class_read_double("custom4", ppm->custom4);
    class_read_double("custom5", ppm->custom5);
    class_read_double("custom6", ppm->custom6);
    class_read_double("custom7", ppm->custom7);
    class_read_double("custom8", ppm->custom8);
    class_read_double("custom9", ppm->custom9);
    class_read_double("custom10", ppm->custom10);
  }

  /* Tests moved from primordial module: */
  if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_H) ||
      (ppm->primordial_spec_type == inflation_V_end)) {
    class_test(ppt->has_scalars == _FALSE_,
               errmsg,
               "inflationary module cannot work if you do not ask for scalar modes");

    class_test(ppt->has_vectors == _TRUE_,
               errmsg,
               "inflationary module cannot work if you ask for vector modes");

    class_test(ppt->has_tensors == _FALSE_,
               errmsg,
               "inflationary module cannot work if you do not ask for tensor modes");

    class_test(ppt->has_bi == _TRUE_ || ppt->has_cdi == _TRUE_ || ppt->has_nid == _TRUE_ ||
                   ppt->has_niv == _TRUE_,
               errmsg,
               "inflationary module cannot work if you ask for isocurvature modes");
  }

  /** (e) parameters for final spectra */

  if (ppt->has_cls == _TRUE_) {
    if (ppt->has_scalars == _TRUE_) {
      if ((ppt->has_cl_cmb_temperature == _TRUE_) || (ppt->has_cl_cmb_polarization == _TRUE_) ||
          (ppt->has_cl_cmb_lensing_potential == _TRUE_))
        class_read_double("l_max_scalars", ppt->l_scalar_max);

      if ((ppt->has_cl_lensing_potential == _TRUE_) || (ppt->has_cl_number_count == _TRUE_))
        class_read_double("l_max_lss", ppt->l_lss_max);
    }

    if (ppt->has_vectors == _TRUE_) {
      class_read_double("l_max_vectors", ppt->l_vector_max);
    }

    if (ppt->has_tensors == _TRUE_) {
      class_read_double("l_max_tensors", ppt->l_tensor_max);
    }
  }

  class_call(parser_read_string(pfc, "lensing", &(string1), &(flag1), errmsg), errmsg, errmsg);

  if ((flag1 == _TRUE_) &&
      ((strstr(string1, "y") != nullptr) || (strstr(string1, "Y") != nullptr))) {
    if ((ppt->has_scalars == _TRUE_) &&
        ((ppt->has_cl_cmb_temperature == _TRUE_) || (ppt->has_cl_cmb_polarization == _TRUE_)) &&
        (ppt->has_cl_cmb_lensing_potential == _TRUE_)) {
      ple->has_lensed_cls = _TRUE_;
    }
    else {
      class_stop(errmsg,
                 "you asked for lensed CMB Cls, but this requires a minimal number of options: "
                 "'modes' should include 's', 'output' should include 'tCl' and/or 'pCL', and "
                 "also, importantly, 'lCl', the CMB lensing potential spectrum. You forgot one of "
                 "those in your input.");
    }
  }

  if ((ppt->has_scalars == _TRUE_) && (ppt->has_cl_cmb_lensing_potential == _TRUE_)) {
    class_read_double("lcmb_rescale", ptr->lcmb_rescale);
    class_read_double("lcmb_tilt", ptr->lcmb_tilt);
    class_read_double("lcmb_pivot", ptr->lcmb_pivot);
  }

  if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) ||
      (ppt->has_velocity_transfers == _TRUE_)) {
    class_call(parser_read_double(pfc, "P_k_max_h/Mpc", &param1, &flag1, errmsg), errmsg, errmsg);
    class_call(parser_read_double(pfc, "P_k_max_1/Mpc", &param2, &flag2, errmsg), errmsg, errmsg);
    class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
               errmsg,
               "In input file, you cannot enter both P_k_max_h/Mpc and P_k_max_1/Mpc, choose one");
    if (flag1 == _TRUE_) {
      ppt->k_max_for_pk = param1 * pba->h;
    }
    if (flag2 == _TRUE_) {
      ppt->k_max_for_pk = param2;
    }

    std::vector<double> zPkValues;
    class_call(readDoubleList(pfc, "z_pk", zPkValues, &flag1, errmsg), errmsg, errmsg);

    if (flag1 == _TRUE_) {
      int1 = static_cast<int>(zPkValues.size());
      class_test(int1 > _Z_PK_NUM_MAX_,
                 errmsg,
                 "you want to write some output for %d different values of z, hence you should "
                 "increase _Z_PK_NUM_MAX_ in include/output.h to at least this number",
                 int1);
      pop->z_pk_num = int1;
      for (i = 0; i < int1; i++) {
        pop->z_pk[i] = zPkValues[i];
      }
    }
  }

  /** Do we want density and velocity transfer functions in Nbody gauge? */
  if ((ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_)) {
    class_call(parser_read_string(pfc, "Nbody gauge transfer functions", &string1, &flag1, errmsg),
               errmsg,
               errmsg);

    if ((flag1 == _TRUE_) &&
        ((strstr(string1, "y") != nullptr) || (strstr(string1, "y") != nullptr))) {
      ppt->has_Nbody_gauge_transfers = _TRUE_;
    }
  }

  /* deal with selection functions */
  if ((ppt->has_cl_number_count == _TRUE_) || (ppt->has_cl_lensing_potential == _TRUE_)) {
    class_call(parser_read_string(pfc, "selection", &(string1), &(flag1), errmsg), errmsg, errmsg);

    if (flag1 == _TRUE_) {
      if (strstr(string1, "gaussian") != nullptr) {
        ppt->selection = gaussian;
      }
      else if (strstr(string1, "tophat") != nullptr) {
        ppt->selection = tophat;
      }
      else if (strstr(string1, "dirac") != nullptr) {
        ppt->selection = dirac;
      }
      else {
        class_stop(errmsg, "In selection function input: type '%s' is unclear", string1);
      }
    }

    std::vector<double> selectionValues;
    class_call(readDoubleList(pfc, "selection_mean", selectionValues, &flag1, errmsg),
               errmsg,
               errmsg);

    if ((flag1 == _TRUE_) && !selectionValues.empty()) {
      int1 = static_cast<int>(selectionValues.size());

      class_test(int1 > _SELECTION_NUM_MAX_,
                 errmsg,
                 "you want to compute density Cl's for %d different bins, hence you should "
                 "increase _SELECTION_NUM_MAX_ in include/perturbations.h to at least this number",
                 int1);

      ppt->selection_num = int1;
      for (i = 0; i < int1; i++) {
        class_test((selectionValues[i] < 0.) || (selectionValues[i] > 1000.),
                   errmsg,
                   "input of selection functions: you asked for a mean redshift equal to %e, "
                   "sounds odd",
                   selectionValues[i]);
        ppt->selection_mean[i] = selectionValues[i];
      }
      /* first set all widths to default; correct eventually later */
      for (i = 1; i < int1; i++) {
        class_test(ppt->selection_mean[i] <= ppt->selection_mean[i - 1],
                   errmsg,
                   "input of selection functions: the list of mean redshifts must be passed in "
                   "growing order; you entered %e before %e",
                   ppt->selection_mean[i - 1],
                   ppt->selection_mean[i]);
        ppt->selection_width[i]              = ppt->selection_width[0];
        ptr->selection_bias[i]               = ptr->selection_bias[0];
        ptr->selection_magnification_bias[i] = ptr->selection_magnification_bias[0];
      }

      selectionValues.clear();
      class_call(readDoubleList(pfc, "selection_width", selectionValues, &flag1, errmsg),
                 errmsg,
                 errmsg);

      if ((flag1 == _TRUE_) && !selectionValues.empty()) {
        int1 = static_cast<int>(selectionValues.size());

        if (int1 == 1) {
          for (i = 0; i < ppt->selection_num; i++) {
            ppt->selection_width[i] = selectionValues[0];
          }
        }
        else if (int1 == ppt->selection_num) {
          for (i = 0; i < int1; i++) {
            ppt->selection_width[i] = selectionValues[i];
          }
        }
        else {
          class_stop(errmsg,
                     "In input for selection function, you asked for %d bin centers and %d bin "
                     "widths; number of bins unclear; you should pass either one bin width (common "
                     "to all bins) or %d bin widths",
                     ppt->selection_num,
                     int1,
                     ppt->selection_num);
        }
      }

      selectionValues.clear();
      class_call(readDoubleList(pfc, "selection_bias", selectionValues, &flag1, errmsg),
                 errmsg,
                 errmsg);

      if ((flag1 == _TRUE_) && !selectionValues.empty()) {
        int1 = static_cast<int>(selectionValues.size());

        if (int1 == 1) {
          for (i = 0; i < ppt->selection_num; i++) {
            ptr->selection_bias[i] = selectionValues[0];
          }
        }
        else if (int1 == ppt->selection_num) {
          for (i = 0; i < int1; i++) {
            ptr->selection_bias[i] = selectionValues[i];
          }
        }
        else {
          class_stop(errmsg,
                     "In input for selection function, you asked for %d bin centers and %d bin "
                     "biases; number of bins unclear; you should pass either one bin bias (common "
                     "to all bins) or %d bin biases",
                     ppt->selection_num,
                     int1,
                     ppt->selection_num);
        }
      }

      selectionValues.clear();
      class_call(readDoubleList(pfc,
                                "selection_magnification_bias",
                                selectionValues,
                                &flag1,
                                errmsg),
                 errmsg,
                 errmsg);

      if ((flag1 == _TRUE_) && !selectionValues.empty()) {
        int1 = static_cast<int>(selectionValues.size());

        if (int1 == 1) {
          for (i = 0; i < ppt->selection_num; i++) {
            ptr->selection_magnification_bias[i] = selectionValues[0];
          }
        }
        else if (int1 == ppt->selection_num) {
          for (i = 0; i < int1; i++) {
            ptr->selection_magnification_bias[i] = selectionValues[i];
          }
        }
        else {
          class_stop(errmsg,
                     "In input for selection function, you asked for %d bin centers and %d bin "
                     "biases; number of bins unclear; you should pass either one bin bias (common "
                     "to all bins) or %d bin biases",
                     ppt->selection_num,
                     int1,
                     ppt->selection_num);
        }
      }
    }

    if (ppt->selection_num > 1) {
      class_read_int("non_diagonal", psp->non_diag);
      if ((psp->non_diag < 0) || (psp->non_diag >= ppt->selection_num))
        class_stop(errmsg,
                   "Input for non_diagonal is %d, while it is expected to be between 0 and %d\n",
                   psp->non_diag,
                   ppt->selection_num - 1);
    }

    class_call(parser_read_string(pfc, "dNdz_selection", &(string1), &(flag1), errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {
      if ((strstr(string1, "analytic") != nullptr)) {
        ptr->has_nz_analytic = _TRUE_;
      }
      else {
        ptr->has_nz_file = _TRUE_;
        class_read_string("dNdz_selection", ptr->nz_file_name);
      }
    }

    class_call(parser_read_string(pfc, "dNdz_evolution", &(string1), &(flag1), errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {
      if ((strstr(string1, "analytic") != nullptr)) {
        ptr->has_nz_evo_analytic = _TRUE_;
      }
      else {
        ptr->has_nz_evo_file = _TRUE_;
        class_read_string("dNdz_evolution", ptr->nz_evo_file_name);
      }
    }

    flag1 = _FALSE_;
    class_call(parser_read_double(pfc, "bias", &param1, &flag1, errmsg), errmsg, errmsg);
    class_test(flag1 == _TRUE_,
               errmsg,
               "the input parameter 'bias' is obsolete, because you can now pass an independent "
               "light-to-mass bias for each bin/selection function. The new input name is "
               "'selection_bias'. It can be set to a single number (common bias for all bins) or "
               "as many numbers as the number of bins");

    flag1 = _FALSE_;
    class_call(parser_read_double(pfc, "s_bias", &param1, &flag1, errmsg), errmsg, errmsg);
    class_test(flag1 == _TRUE_,
               errmsg,
               "the input parameter 's_bias' is obsolete, because you can now pass an independent "
               "magnitude bias for each bin/selection function. The new input name is "
               "'selection_magnitude_bias'. It can be set to a single number (common magnitude "
               "bias for all bins) or as many numbers as the number of bins");
  }
  /* end of selection function section */

  /* deal with z_max issues */
  if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) ||
      (ppt->has_velocity_transfers == _TRUE_) || (ppt->has_cl_number_count == _TRUE_) ||
      (ppt->has_cl_lensing_potential == _TRUE_)) {
    class_call(parser_read_double(pfc, "z_max_pk", &param1, &flag1, errmsg), errmsg, errmsg);

    if (flag1 == _TRUE_) {
      ppt->z_max_pk = param1;
    }
    else {
      ppt->z_max_pk = 0.;

      if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) ||
          (ppt->has_velocity_transfers == _TRUE_)) {
        for (i = 0; i < pop->z_pk_num; i++) {
          ppt->z_max_pk = MAX(ppt->z_max_pk, pop->z_pk[i]);
        }
      }

      if ((ppt->has_cl_number_count == _TRUE_) || (ppt->has_cl_lensing_potential == _TRUE_)) {
        for (int bin = 0; bin < ppt->selection_num; bin++) {
          /* the few lines below should be consistent with their counterpart in transfer.c, in transfer_selection_times() */
          if (ppt->selection == gaussian) {
            z_max = ppt->selection_mean[bin] +
                    ppt->selection_width[bin] * ppr->selection_cut_at_sigma;
          }
          if (ppt->selection == tophat) {
            z_max = ppt->selection_mean[bin] +
                    (1. + ppr->selection_cut_at_sigma * ppr->selection_tophat_edge) *
                        ppt->selection_width[bin];
          }
          if (ppt->selection == dirac) {
            z_max = ppt->selection_mean[bin];
          }
          ppt->z_max_pk = MAX(ppt->z_max_pk, z_max);
        }
      }
    }
    psp->z_max_pk = ppt->z_max_pk;
  }
  /* end of z_max section */

  class_call(parser_read_string(pfc, "root", &string1, &flag1, errmsg), errmsg, errmsg);
  if (flag1 == _TRUE_) {
    class_test(strlen(string1) > _FILENAMESIZE_ - 32,
               errmsg,
               "Root directory name is too long. Please install in other directory, or increase "
               "_FILENAMESIZE_ in common.h");
    strcpy(pop->root, string1);
  }

  class_call(parser_read_string(pfc, "headers", &(string1), &(flag1), errmsg), errmsg, errmsg);

  if ((flag1 == _TRUE_) &&
      ((strstr(string1, "y") == nullptr) && (strstr(string1, "Y") == nullptr))) {
    pop->write_header = _FALSE_;
  }

  class_call(parser_read_string(pfc, "format", &string1, &flag1, errmsg), errmsg, errmsg);

  if (flag1 == _TRUE_) {
    if ((strstr(string1, "class") != nullptr) || (strstr(string1, "CLASS") != nullptr))
      pop->output_format = class_format;
    else {
      if ((strstr(string1, "camb") != nullptr) || (strstr(string1, "CAMB") != nullptr))
        pop->output_format = camb_format;
      else
        class_stop(errmsg,
                   "You wrote: format='%s'. Could not identify any of the possible formats "
                   "('class', 'CLASS', 'camb', 'CAMB')",
                   string1);
    }
  }

  /** (f) parameter related to the non-linear spectra computation */

  class_call(parser_read_string(pfc, "non linear", &(string1), &(flag1), errmsg), errmsg, errmsg);

  if (flag1 == _TRUE_) {
    class_test(ppt->has_perturbations == _FALSE_,
               errmsg,
               "You requested non linear computation but no linear computation. You must set "
               "output to tCl or similar.");

    if ((strstr(string1, "halofit") != nullptr) || (strstr(string1, "Halofit") != nullptr) ||
        (strstr(string1, "HALOFIT") != nullptr)) {
      pnl->method       = nl_halofit;
      ppt->k_max_for_pk = MAX(ppt->k_max_for_pk,
                              MAX(ppr->halofit_min_k_max, ppr->nonlinear_min_k_max));
      ppt->has_nl_corrections_based_on_delta_m = _TRUE_;
    }
    if ((strstr(string1, "hmcode") != nullptr) || (strstr(string1, "HMCODE") != nullptr) ||
        (strstr(string1, "HMcode") != nullptr) || (strstr(string1, "Hmcode") != nullptr)) {
      pnl->method       = nl_HMcode;
      ppt->k_max_for_pk = MAX(ppt->k_max_for_pk,
                              MAX(ppr->hmcode_min_k_max, ppr->nonlinear_min_k_max));
      ppt->has_nl_corrections_based_on_delta_m = _TRUE_;
      class_read_int("extrapolation_method", pnl->extrapolation_method);

      class_call(parser_read_string(pfc, "feedback model", &(string1), &(flag1), errmsg),
                 errmsg,
                 errmsg);

      if (flag1 == _TRUE_) {
        if (strstr(string1, "emu_dmonly") != nullptr) {
          pnl->feedback = nl_emu_dmonly;
        }
        if (strstr(string1, "owls_dmonly") != nullptr) {
          pnl->feedback = nl_owls_dmonly;
        }
        if (strstr(string1, "owls_ref") != nullptr) {
          pnl->feedback = nl_owls_ref;
        }
        if (strstr(string1, "owls_agn") != nullptr) {
          pnl->feedback = nl_owls_agn;
        }
        if (strstr(string1, "owls_dblim") != nullptr) {
          pnl->feedback = nl_owls_dblim;
        }
      }

      class_call(parser_read_double(pfc, "eta_0", &param2, &flag2, errmsg), errmsg, errmsg);
      class_call(parser_read_double(pfc, "c_min", &param3, &flag3, errmsg), errmsg, errmsg);

      class_test(((flag1 == _TRUE_) && ((flag2 == _TRUE_) || (flag3 == _TRUE_))),
                 errmsg,
                 "In input file, you cannot enter both a baryonic feedback model and a choice of "
                 "baryonic feedback parameters, choose one of both methods");

      if ((flag2 == _TRUE_) && (flag3 == _TRUE_)) {
        pnl->feedback = nl_user_defined;
        class_read_double("eta_0", pnl->eta_0);
        class_read_double("c_min", pnl->c_min);
      }
      else if ((flag2 == _TRUE_) && (flag3 == _FALSE_)) {
        pnl->feedback = nl_user_defined;
        class_read_double("eta_0", pnl->eta_0);
        pnl->c_min = (0.98 - pnl->eta_0) / 0.12;
      }
      else if ((flag2 == _FALSE_) && (flag3 == _TRUE_)) {
        pnl->feedback = nl_user_defined;
        class_read_double("c_min", pnl->c_min);
        pnl->eta_0 = 0.98 - 0.12 * pnl->c_min;
      }

      class_call(parser_read_double(pfc, "z_infinity", &param1, &flag1, errmsg), errmsg, errmsg);

      if (flag1 == _TRUE_) {
        class_read_double("z_infinity", pnl->z_infinity);
      }
    }
  }

  /** (g) amount of information sent to standard output (none if all set to zero) */

  class_read_int("background_verbose", pba->background_verbose);

  class_read_int("thermodynamics_verbose", pth->thermodynamics_verbose);

  class_read_int("perturbations_verbose", ppt->perturbations_verbose);

  class_read_int("transfer_verbose", ptr->transfer_verbose);

  class_read_int("primordial_verbose", ppm->primordial_verbose);

  class_read_int("spectra_verbose", psp->spectra_verbose);

  class_read_int("nonlinear_verbose", pnl->nonlinear_verbose);

  class_read_int("lensing_verbose", ple->lensing_verbose);

  class_read_int("output_verbose", pop->output_verbose);

  if (ppt->has_tensors == _TRUE_) {
    /** - ---> Include ur and ncdm shear in tensor computation? */
    class_call(parser_read_string(pfc, "tensor method", &string1, &flag1, errmsg), errmsg, errmsg);
    if (flag1 == _TRUE_) {
      if (strstr(string1, "photons") != nullptr)
        ppt->tensor_method = tm_photons_only;
      if (strstr(string1, "massless") != nullptr)
        ppt->tensor_method = tm_massless_approximation;
      if (strstr(string1, "exact") != nullptr)
        ppt->tensor_method = tm_exact;
    }
  }

  /** - ---> derivatives of baryon sound speed only computed if some non-minimal tight-coupling schemes is requested */
  if ((ppr->tight_coupling_approximation == (int) first_order_CLASS) ||
      (ppr->tight_coupling_approximation == (int) second_order_CLASS)) {
    pth->compute_cb2_derivatives = _TRUE_;
  }

  class_test(ppr->ur_fluid_trigger_tau_over_tau_k ==
                 ppr->radiation_streaming_trigger_tau_over_tau_k,
             errmsg,
             "please choose different values for precision parameters "
             "ur_fluid_trigger_tau_over_tau_k and radiation_streaming_trigger_tau_over_tau_k, in "
             "order to avoid switching two approximation schemes at the same time");

  if (omega_budget_.idr.value_or(0.) != 0.) {
    class_test(ppr->idr_streaming_trigger_tau_over_tau_k ==
                   ppr->radiation_streaming_trigger_tau_over_tau_k,
               errmsg,
               "please choose different values for precision parameters "
               "dark_radiation_trigger_tau_over_tau_k and "
               "radiation_streaming_trigger_tau_over_tau_k, in order to avoid switching two "
               "approximation schemes at the same time");

    class_test(ppr->idr_streaming_trigger_tau_over_tau_k == ppr->ur_fluid_trigger_tau_over_tau_k,
               errmsg,
               "please choose different values for precision parameters "
               "dark_radiation_trigger_tau_over_tau_k and ur_fluid_trigger_tau_over_tau_k, in "
               "order to avoid switching two approximation schemes at the same time");

    class_test(ppr->idr_streaming_trigger_tau_over_tau_k == ppr->ncdm_fluid_trigger_tau_over_tau_k,
               errmsg,
               "please choose different values for precision parameters "
               "dark_radiation_trigger_tau_over_tau_k and ncdm_fluid_trigger_tau_over_tau_k, in "
               "order to avoid switching two approximation schemes at the same time");
  }

  /**
   * Here we can place all obsolete (deprecated) names for the precision parameters,
   * so they will still get read.
   * The new parameter names should be used preferrably
   * */
  class_read_double(
      "k_scalar_min_tau0",
      ppr->k_min_tau0);  // obsolete precision parameter: read for compatibility with old precision files
  class_read_double(
      "k_scalar_max_tau0_over_l_max",
      ppr->k_max_tau0_over_l_max);  // obsolete precision parameter: read for compatibility with old precision files
  class_read_double(
      "k_scalar_step_sub",
      ppr->k_step_sub);  // obsolete precision parameter: read for compatibility with old precision files
  class_read_double(
      "k_scalar_step_super",
      ppr->k_step_super);  // obsolete precision parameter: read for compatibility with old precision files
  class_read_double(
      "k_scalar_step_transition",
      ppr->k_step_transition);  // obsolete precision parameter: read for compatibility with old precision files
  class_read_double(
      "k_scalar_k_per_decade_for_pk",
      ppr->k_per_decade_for_pk);  // obsolete precision parameter: read for compatibility with old precision files
  class_read_double(
      "k_scalar_k_per_decade_for_bao",
      ppr->k_per_decade_for_bao);  // obsolete precision parameter: read for compatibility with old precision files
  class_read_double(
      "k_scalar_bao_center",
      ppr->k_bao_center);  // obsolete precision parameter: read for compatibility with old precision files
  class_read_double(
      "k_scalar_bao_width",
      ppr->k_bao_width);  // obsolete precision parameter: read for compatibility with old precision files

  class_read_double(
      "k_step_trans_scalars",
      ppr->q_linstep);  // obsolete precision parameter: read for compatibility with old precision files
  class_read_double(
      "k_step_trans_tensors",
      ppr->q_linstep);  // obsolete precision parameter: read for compatibility with old precision files
  class_read_double(
      "k_step_trans",
      ppr->q_linstep);  // obsolete precision parameter: read for compatibility with old precision files
  class_read_double(
      "q_linstep_trans",
      ppr->q_linstep);  // obsolete precision parameter: read for compatibility with old precision files
  class_read_double(
      "q_logstep_trans",
      ppr->q_logstep_spline);  // obsolete precision parameter: read for compatibility with old precision files

  class_call(parser_read_string(pfc,
                                "l_switch_limber_for_cl_density_over_z",
                                &string1,
                                &flag1,
                                errmsg),
             errmsg,
             errmsg);

  class_test(flag1 == _TRUE_,
             errmsg,
             "You passed in input a precision parameter called "
             "l_switch_limber_for_cl_density_over_z. This syntax is deprecated since v2.5.0. "
             "Please use instead the two precision parameters l_switch_limber_for_nc_local_over_z, "
             "l_switch_limber_for_nc_los_over_z, defined in include/common.h, and allowing for "
             "better performance.");

  /** (i) Write values in file */
  if (ple->has_lensed_cls == _TRUE_)
    ppt->l_scalar_max += ppr->delta_l_max;

  /** - (i.1.) shall we write background quantities in a file? */

  class_call(parser_read_string(pfc, "write background", &string1, &flag1, errmsg), errmsg, errmsg);

  if ((flag1 == _TRUE_) &&
      ((strstr(string1, "y") != nullptr) || (strstr(string1, "Y") != nullptr))) {
    pop->write_background = _TRUE_;
  }

  /** - (i.2.) shall we write thermodynamics quantities in a file? */

  class_call(parser_read_string(pfc, "write thermodynamics", &string1, &flag1, errmsg),
             errmsg,
             errmsg);

  if ((flag1 == _TRUE_) &&
      ((strstr(string1, "y") != nullptr) || (strstr(string1, "Y") != nullptr))) {
    pop->write_thermodynamics = _TRUE_;
  }

  /** - (i.3.) shall we write perturbation quantities in files? */

  std::vector<double> kOutputValues;
  class_call(readDoubleList(pfc, "k_output_values", kOutputValues, &flag1, errmsg), errmsg, errmsg);

  if (flag1 == _TRUE_) {
    int1 = static_cast<int>(kOutputValues.size());
    class_test(int1 > _MAX_NUMBER_OF_K_FILES_,
               errmsg,
               "you want to write some output for %d different values of k, hence you should "
               "increase _MAX_NUMBER_OF_K_FILES_ in include/perturbations.h to at least this "
               "number",
               int1);
    ppt->k_output_values_num = int1;

    for (i = 0; i < int1; i++) {
      ppt->k_output_values[i] = kOutputValues[i];
    }

    /* Sort k_output_values ascending */
    std::sort(ppt->k_output_values, ppt->k_output_values + ppt->k_output_values_num);

    ppt->store_perturbations = _TRUE_;
    pop->write_perturbations = _TRUE_;
  }

  /** - (i.4.) shall we write primordial spectra in a file? */

  class_call(parser_read_string(pfc, "write primordial", &string1, &flag1, errmsg), errmsg, errmsg);

  if ((flag1 == _TRUE_) &&
      ((strstr(string1, "y") != nullptr) || (strstr(string1, "Y") != nullptr))) {
    pop->write_primordial = _TRUE_;
  }

  /* flags for calling the interpolation routine */
  pba->short_info  = 0;
  pba->normal_info = 1;
  pba->long_info   = 2;

  pba->inter_normal  = 0;
  pba->inter_closeby = 1;

  /** - (i.5) special steps if we want Halofit with wa_fld non-zero:
      so-called "Pk_equal method" of 0810.0190 and 1601.07230 */

  /* Species aren't built yet at this point — peek wa_fld from pfc directly.
     Note: this read is informational (no class_read_double-style write target),
     so we just probe the file content via parser_read_double.

     The original (pre-refactor) parser block only wrote pba->wa_fld in the CLP
     branch, so the Pk_equal gate must only engage for CLP. Probe
     fluid_equation_of_state first; default-when-absent is CLP (matching
     FluidSpecies::CreateAll). */
  bool fluid_eos_is_clp = true;
  class_call(parser_read_string(pfc, "fluid_equation_of_state", &string1, &flag1, errmsg),
             errmsg,
             errmsg);
  if (flag1 == _TRUE_) {
    if ((strstr(string1, "EDE") != nullptr) || (strstr(string1, "ede") != nullptr)) {
      fluid_eos_is_clp = false;
    }
    else if ((strstr(string1, "CLP") != nullptr) || (strstr(string1, "clp") != nullptr)) {
      fluid_eos_is_clp = true;
    }
    /* Other strings: leave as CLP default (FluidSpecies::CreateAll will error
       later on an unrecognized value). */
  }
  double wa_fld_peek = 0.;
  class_call(parser_read_double(pfc, "wa_fld", &param1, &flag1, errmsg), errmsg, errmsg);
  if (flag1 == _TRUE_)
    wa_fld_peek = param1;
  if ((pnl->method == nl_halofit) && fluid_present_pfc && fluid_eos_is_clp && (wa_fld_peek != 0.)) {
    class_call(parser_read_string(pfc, "pk_eq", &string1, &flag1, errmsg), errmsg, errmsg);

    if ((flag1 == _TRUE_) &&
        ((strstr(string1, "y") != nullptr) || (strstr(string1, "Y") != nullptr))) {
      pnl->has_pk_eq = _TRUE_;
    }
  }

  /* ── Precision-consistency tests for perturbation-hierarchy l_max values.
     Moved here from perturb_vector_init: these checks depend only on ppr
     (plus ppt->idr_nature for the IDR test), so they belong at input-parse
     time, not inside the per-(k, approximation) hot path.  Tests run
     unconditionally — a too-low l_max is a user-config error whether or not
     the species ends up active. */
  class_test(ppr->l_max_g < 4,
             errmsg,
             "ppr->l_max_g should be at least 4, i.e. we must integrate at least over photon "
             "density, velocity, shear, third and fourth momentum");
  class_test(ppr->l_max_pol_g < 4, errmsg, "ppr->l_max_pol_g should be at least 4");
  class_test(ppr->l_max_ur < 4,
             errmsg,
             "ppr->l_max_ur should be at least 4, i.e. we must integrate at least over "
             "neutrino/relic density, velocity, shear, third and fourth momentum");
  class_test(ppr->l_max_dr < 4,
             errmsg,
             "ppr->l_max_dr should be at least 4, i.e. we must integrate at least over "
             "neutrino/relic density, velocity, shear, third and fourth momentum");
  class_test((ppr->l_max_idr < 4) && (ppt->idr_nature == idr_free_streaming),
             errmsg,
             "ppr->l_max_idr should be at least 4, i.e. we must integrate at least over "
             "interacting dark radiation density, velocity, shear, third and fourth momentum");
  class_test(ppr->l_max_g_ten < 4,
             errmsg,
             "ppr->l_max_g_ten should be at least 4, i.e. we must integrate at least over photon "
             "density, velocity, shear, third momentum");
  class_test(ppr->l_max_pol_g_ten < 4, errmsg, "ppr->l_max_pol_g_ten should be at least 4");
  class_test(ppr->l_max_ncdm < 4,
             errmsg,
             "ppr->l_max_ncdm=%d should be at least 4, i.e. we must integrate at least "
             "over first four momenta of non-cold dark matter perturbed phase-space "
             "distribution",
             ppr->l_max_ncdm);

  return _SUCCESS_;
}

/**
 * All default parameter values (for input parameters)
 *
 * @return the error status
 */

int InputModule::input_default_params() {
  background* pba = &background_; /* for cosmological background */
  primordial* ppm = &primordial_; /* for primordial spectra */
  transfers* ptr  = &transfers_;  /* for transfer functions */
  spectra* psp    = &spectra_;    /* for output spectra */
  output* pop     = &output_;

  double sigma_B = 2. * pow(_PI_, 5) * pow(_k_B_, 4) / 15. / pow(_h_P_, 3) / pow(_c_, 2);

  /** Define computed default parameter values that depend on other defaults.
      Simple constant defaults are now set as in-struct default member initializers. */

  /** - background structure: computed defaults */

  /* 5.10.2014: default parameters matched to Planck 2013 + WP
     best-fitting model, with ones small difference: the published
     Planck 2013 + WP bestfit is with h=0.6704 and one massive
     neutrino species with m_ncdm=0.06eV; here we assume only massless
     neutrinos in the default model; for the CMB, taking m_ncdm = 0 or
     0.06 eV makes practically no difference, provided that we adapt
     the value of h in order ot get the same peak scale, i.e. the same
     100*theta_s. The Planck 2013 + WP best-fitting model with
     h=0.6704 gives 100*theta_s = 1.042143 (or equivalently
     100*theta_MC=1.04119). By taking only massless neutrinos, one
     gets the same 100*theta_s provided that h is increased to
     0.67556. Hence, we take h=0.67556, N_ur=3.046, N_ncdm=0, and all
     other parameters from the Planck2013 Cosmological Parameter
     paper. */

  pba->H0       = pba->h * 1.e5 / _c_;
  pba->Omega0_g = (4. * sigma_B / _c_ * pow(pba->T_cmb, 4.)) /
                  (3. * _c_ * _c_ * 1.e10 * pba->h * pba->h / _Mpc_over_m_ / _Mpc_over_m_ / 8. /
                   _PI_ / _G_);
  pba->Omega0_b = 0.022032 / pow(pba->h, 2);
  // pba->Omega0_{ur,cdm,lambda,...} are no longer stored: those defaults live
  // on the species (read from pfc or applied by ReadCoupledOmegaBudget).

  /** - primordial structure: computed defaults */

  ppm->n_t     = -ppm->r / 8. * (2. - ppm->r / 8. - ppm->n_s);
  ppm->alpha_t = ppm->r / 8. * (ppm->r / 8. + ppm->n_s - 1.);

  /** - transfer structure: array element defaults */

  ptr->selection_bias[0]               = 1.;
  ptr->selection_magnification_bias[0] = 0.;

  /** - spectra structure: computed default */

  psp->z_max_pk = pop->z_pk[0];

  /** - perturbation structure: array element defaults */

  perturbations_.selection_mean[0]  = 1.;
  perturbations_.selection_width[0] = 0.1;

  return _SUCCESS_;
}

/**
 * Initialize the precision parameter structure.
 *
 * All precision parameters used in the other modules are listed here
 * and assigned here a default value.
 *
 * @return the error status
 *
 */

int InputModule::input_default_precision() {
  precision* ppr = &precision_; /* for precision parameters */

  /** Numeric and type precision parameters now have in-struct defaults
      (via precisions.h macros). Only string parameters need runtime
      path concatenation with class_dir. */

  class_test(ppr->smallest_allowed_variation < 0,
             ppr->error_message,
             "smallest_allowed_variation = %e < 0",
             ppr->smallest_allowed_variation);

  /* String parameters require runtime path concatenation */
  strncpy(ppr->sBBN_file, ppr->class_dir, sizeof(ppr->sBBN_file));
  strcat(ppr->sBBN_file, "/bbn/sBBN_2017.dat");

  strncpy(ppr->hyrec_Alpha_inf_file, ppr->class_dir, sizeof(ppr->hyrec_Alpha_inf_file));
  strcat(ppr->hyrec_Alpha_inf_file, "/hyrec/Alpha_inf.dat");

  strncpy(ppr->hyrec_R_inf_file, ppr->class_dir, sizeof(ppr->hyrec_R_inf_file));
  strcat(ppr->hyrec_R_inf_file, "/hyrec/R_inf.dat");

  strncpy(ppr->hyrec_two_photon_tables_file,
          ppr->class_dir,
          sizeof(ppr->hyrec_two_photon_tables_file));
  strcat(ppr->hyrec_two_photon_tables_file, "/hyrec/two_photon_tables.dat");

  return _SUCCESS_;
}

/** Overloaded helpers for type-dispatched precision parameter reading. */
namespace {
void read(const FileContent& fc, const char* name, double& v) {
  fc.read_double(name, v);
}
void read(const FileContent& fc, const char* name, int& v) {
  fc.read_int(name, v);
}
void read(const FileContent& fc, const char* name, FileName& v) {
  std::string s;
  if (fc.read_string(name, s))
    strncpy(v, s.c_str(), sizeof(FileName) - 1);
}
template <typename E>
void read_enum(const FileContent& fc, const char* name, E& v) {
  int tmp;
  if (fc.read_int(name, tmp))
    v = static_cast<E>(tmp);
}
}  // anonymous namespace

void precision::parse(const FileContent& fc) {
  /* Background */
  read(fc, "a_ini_over_a_today_default", a_ini_over_a_today_default);
  read(fc, "back_integration_stepsize", back_integration_stepsize);
  read(fc, "tol_background_integration", tol_background_integration);
  read(fc, "tol_initial_Omega_r", tol_initial_Omega_r);
  read(fc, "tol_M_ncdm", tol_M_ncdm);
  read(fc, "tol_ncdm", tol_ncdm);
  read(fc, "tol_ncdm_synchronous", tol_ncdm_synchronous);
  read(fc, "tol_ncdm_newtonian", tol_ncdm_newtonian);
  read(fc, "tol_ncdm_bg", tol_ncdm_bg);
  read(fc, "tol_ncdm_initial_w", tol_ncdm_initial_w);
  read(fc, "tol_tau_eq", tol_tau_eq);
  read(fc, "Omega0_cdm_min_synchronous", Omega0_cdm_min_synchronous);
  read(fc, "sBBN file", sBBN_file);

  /* Thermodynamics */
  read(fc, "recfast_z_initial", recfast_z_initial);
  read(fc, "recfast_Nz0", recfast_Nz0);
  read(fc, "thermo_z_initial_idm_dr", thermo_z_initial_idm_dr);
  read(fc, "thermo_Nz1_idm_dr", thermo_Nz1_idm_dr);
  read(fc, "thermo_Nz2_idm_dr", thermo_Nz2_idm_dr);
  read(fc, "tol_thermo_integration", tol_thermo_integration);
  read(fc, "recfast_Heswitch", recfast_Heswitch);
  read(fc, "recfast_fudge_He", recfast_fudge_He);
  read(fc, "recfast_Hswitch", recfast_Hswitch);
  read(fc, "recfast_fudge_H", recfast_fudge_H);
  read(fc, "recfast_delta_fudge_H", recfast_delta_fudge_H);
  read(fc, "recfast_AGauss1", recfast_AGauss1);
  read(fc, "recfast_AGauss2", recfast_AGauss2);
  read(fc, "recfast_zGauss1", recfast_zGauss1);
  read(fc, "recfast_zGauss2", recfast_zGauss2);
  read(fc, "recfast_wGauss1", recfast_wGauss1);
  read(fc, "recfast_wGauss2", recfast_wGauss2);
  read(fc, "recfast_z_He_1", recfast_z_He_1);
  read(fc, "recfast_delta_z_He_1", recfast_delta_z_He_1);
  read(fc, "recfast_z_He_2", recfast_z_He_2);
  read(fc, "recfast_delta_z_He_2", recfast_delta_z_He_2);
  read(fc, "recfast_z_He_3", recfast_z_He_3);
  read(fc, "recfast_delta_z_He_3", recfast_delta_z_He_3);
  read(fc, "recfast_x_He0_trigger", recfast_x_He0_trigger);
  read(fc, "recfast_x_He0_trigger2", recfast_x_He0_trigger2);
  read(fc, "recfast_x_He0_trigger_delta", recfast_x_He0_trigger_delta);
  read(fc, "recfast_x_H0_trigger", recfast_x_H0_trigger);
  read(fc, "recfast_x_H0_trigger2", recfast_x_H0_trigger2);
  read(fc, "recfast_x_H0_trigger_delta", recfast_x_H0_trigger_delta);
  read(fc, "recfast_H_frac", recfast_H_frac);
  read(fc, "reionization_z_start_max", reionization_z_start_max);
  read(fc, "reionization_sampling", reionization_sampling);
  read(fc, "reionization_optical_depth_tol", reionization_optical_depth_tol);
  read(fc, "reionization_start_factor", reionization_start_factor);
  read(fc, "thermo_rate_smoothing_radius", thermo_rate_smoothing_radius);
  read(fc, "Alpha_inf hyrec file", hyrec_Alpha_inf_file);
  read(fc, "R_inf hyrec file", hyrec_R_inf_file);
  read(fc, "two_photon_tables hyrec file", hyrec_two_photon_tables_file);

  /* Perturbations */
  read(fc, "k_min_tau0", k_min_tau0);
  read(fc, "k_max_tau0_over_l_max", k_max_tau0_over_l_max);
  read(fc, "k_step_sub", k_step_sub);
  read(fc, "k_step_super", k_step_super);
  read(fc, "k_step_transition", k_step_transition);
  read(fc, "k_step_super_reduction", k_step_super_reduction);
  read(fc, "k_per_decade_for_pk", k_per_decade_for_pk);
  read(fc, "idmdr_boost_k_per_decade_for_pk", idmdr_boost_k_per_decade_for_pk);
  read(fc, "k_per_decade_for_bao", k_per_decade_for_bao);
  read(fc, "k_bao_center", k_bao_center);
  read(fc, "k_bao_width", k_bao_width);
  read(fc, "start_small_k_at_tau_c_over_tau_h", start_small_k_at_tau_c_over_tau_h);
  read(fc, "start_large_k_at_tau_h_over_tau_k", start_large_k_at_tau_h_over_tau_k);
  read(fc, "tight_coupling_trigger_tau_c_over_tau_h", tight_coupling_trigger_tau_c_over_tau_h);
  read(fc, "tight_coupling_trigger_tau_c_over_tau_k", tight_coupling_trigger_tau_c_over_tau_k);
  read(fc, "start_sources_at_tau_c_over_tau_h", start_sources_at_tau_c_over_tau_h);
  read(fc, "tight_coupling_approximation", tight_coupling_approximation);
  read(fc,
       "idm_dr_tight_coupling_trigger_tau_c_over_tau_k",
       idm_dr_tight_coupling_trigger_tau_c_over_tau_k);
  read(fc,
       "idm_dr_tight_coupling_trigger_tau_c_over_tau_h",
       idm_dr_tight_coupling_trigger_tau_c_over_tau_h);
  read(fc, "idm_drmd_tight_coupling_trigger_G_over_aH", idm_drmd_tight_coupling_trigger_G_over_aH);
  read(fc, "l_max_g", l_max_g);
  read(fc, "l_max_pol_g", l_max_pol_g);
  read(fc, "l_max_dr", l_max_dr);
  read(fc, "l_max_dr_col", l_max_dr_col);
  read(fc, "l_max_ur", l_max_ur);
  read(fc, "l_max_idr", l_max_idr);
  read(fc, "l_max_ncdm", l_max_ncdm);
  read(fc, "l_max_g_ten", l_max_g_ten);
  read(fc, "l_max_pol_g_ten", l_max_pol_g_ten);
  read(fc, "curvature_ini", curvature_ini);
  read(fc, "entropy_ini", entropy_ini);
  read(fc, "gw_ini", gw_ini);
  read(fc, "perturb_integration_stepsize", perturb_integration_stepsize);
  read(fc, "perturb_sampling_stepsize", perturb_sampling_stepsize);
  read(fc, "tol_perturb_integration", tol_perturb_integration);
  read(fc, "c_gamma_k_H_square_max", c_gamma_k_H_square_max);
  read(fc, "tol_tau_approx", tol_tau_approx);
  read(fc, "radiation_streaming_approximation", radiation_streaming_approximation);
  read(fc,
       "radiation_streaming_trigger_tau_over_tau_k",
       radiation_streaming_trigger_tau_over_tau_k);
  read(fc,
       "radiation_streaming_trigger_tau_c_over_tau",
       radiation_streaming_trigger_tau_c_over_tau);
  read(fc, "idr_streaming_approximation", idr_streaming_approximation);
  read(fc, "idr_streaming_trigger_tau_over_tau_k", idr_streaming_trigger_tau_over_tau_k);
  read(fc, "idr_streaming_trigger_tau_c_over_tau", idr_streaming_trigger_tau_c_over_tau);
  read(fc, "ur_fluid_approximation", ur_fluid_approximation);
  read(fc, "ur_fluid_trigger_tau_over_tau_k", ur_fluid_trigger_tau_over_tau_k);
  read(fc, "ncdm_fluid_approximation", ncdm_fluid_approximation);
  read(fc, "ncdm_fluid_trigger_tau_over_tau_k", ncdm_fluid_trigger_tau_over_tau_k);
  read(fc, "neglect_CMB_sources_below_visibility", neglect_CMB_sources_below_visibility);
  read_enum(fc, "evolver", evolver);

  /* Primordial */
  read(fc, "k_per_decade_primordial", k_per_decade_primordial);
  read(fc, "primordial_inflation_ratio_min", primordial_inflation_ratio_min);
  read(fc, "primordial_inflation_ratio_max", primordial_inflation_ratio_max);
  read(fc, "primordial_inflation_phi_ini_maxit", primordial_inflation_phi_ini_maxit);
  read(fc, "primordial_inflation_pt_stepsize", primordial_inflation_pt_stepsize);
  read(fc, "primordial_inflation_bg_stepsize", primordial_inflation_bg_stepsize);
  read(fc, "primordial_inflation_tol_integration", primordial_inflation_tol_integration);
  read(fc,
       "primordial_inflation_attractor_precision_pivot",
       primordial_inflation_attractor_precision_pivot);
  read(fc,
       "primordial_inflation_attractor_precision_initial",
       primordial_inflation_attractor_precision_initial);
  read(fc, "primordial_inflation_attractor_maxit", primordial_inflation_attractor_maxit);
  read(fc, "primordial_inflation_tol_curvature", primordial_inflation_tol_curvature);
  read(fc, "primordial_inflation_aH_ini_target", primordial_inflation_aH_ini_target);
  read(fc, "primordial_inflation_end_dphi", primordial_inflation_end_dphi);
  read(fc, "primordial_inflation_end_logstep", primordial_inflation_end_logstep);
  read(fc, "primordial_inflation_small_epsilon", primordial_inflation_small_epsilon);
  read(fc, "primordial_inflation_small_epsilon_tol", primordial_inflation_small_epsilon_tol);
  read(fc, "primordial_inflation_extra_efolds", primordial_inflation_extra_efolds);

  /* Transfer */
  read(fc, "l_linstep", l_linstep);
  read(fc, "l_logstep", l_logstep);
  read(fc, "hyper_x_min", hyper_x_min);
  read(fc, "hyper_sampling_flat", hyper_sampling_flat);
  read(fc, "hyper_sampling_curved_low_nu", hyper_sampling_curved_low_nu);
  read(fc, "hyper_sampling_curved_high_nu", hyper_sampling_curved_high_nu);
  read(fc, "hyper_nu_sampling_step", hyper_nu_sampling_step);
  read(fc, "hyper_phi_min_abs", hyper_phi_min_abs);
  read(fc, "hyper_x_tol", hyper_x_tol);
  read(fc, "hyper_flat_approximation_nu", hyper_flat_approximation_nu);
  read(fc, "q_linstep", q_linstep);
  read(fc, "q_logstep_spline", q_logstep_spline);
  read(fc, "q_logstep_open", q_logstep_open);
  read(fc, "q_logstep_trapzd", q_logstep_trapzd);
  read(fc, "q_numstep_transition", q_numstep_transition);
  read(fc, "transfer_neglect_delta_k_S_t0", transfer_neglect_delta_k_S_t0);
  read(fc, "transfer_neglect_delta_k_S_t1", transfer_neglect_delta_k_S_t1);
  read(fc, "transfer_neglect_delta_k_S_t2", transfer_neglect_delta_k_S_t2);
  read(fc, "transfer_neglect_delta_k_S_e", transfer_neglect_delta_k_S_e);
  read(fc, "transfer_neglect_delta_k_V_t1", transfer_neglect_delta_k_V_t1);
  read(fc, "transfer_neglect_delta_k_V_t2", transfer_neglect_delta_k_V_t2);
  read(fc, "transfer_neglect_delta_k_V_e", transfer_neglect_delta_k_V_e);
  read(fc, "transfer_neglect_delta_k_V_b", transfer_neglect_delta_k_V_b);
  read(fc, "transfer_neglect_delta_k_T_t2", transfer_neglect_delta_k_T_t2);
  read(fc, "transfer_neglect_delta_k_T_e", transfer_neglect_delta_k_T_e);
  read(fc, "transfer_neglect_delta_k_T_b", transfer_neglect_delta_k_T_b);
  read(fc, "transfer_neglect_late_source", transfer_neglect_late_source);
  read(fc, "l_switch_limber", l_switch_limber);
  read(fc, "l_switch_limber_for_nc_local_over_z", l_switch_limber_for_nc_local_over_z);
  read(fc, "l_switch_limber_for_nc_los_over_z", l_switch_limber_for_nc_los_over_z);
  read(fc, "selection_cut_at_sigma", selection_cut_at_sigma);
  read(fc, "selection_sampling", selection_sampling);
  read(fc, "selection_sampling_bessel", selection_sampling_bessel);
  read(fc, "selection_sampling_bessel_los", selection_sampling_bessel_los);
  read(fc, "selection_tophat_edge", selection_tophat_edge);

  /* Nonlinear */
  read(fc, "sigma_k_per_decade", sigma_k_per_decade);
  read(fc, "nonlinear_min_k_max", nonlinear_min_k_max);
  read(fc, "halofit_min_k_nonlinear", halofit_min_k_nonlinear);
  read(fc, "halofit_min_k_max", halofit_min_k_max);
  read(fc, "halofit_k_per_decade", halofit_k_per_decade);
  read(fc, "halofit_sigma_precision", halofit_sigma_precision);
  read(fc, "halofit_tol_sigma", halofit_tol_sigma);
  read(fc, "pk_eq_z_max", pk_eq_z_max);
  read(fc, "pk_eq_tol", pk_eq_tol);
  read(fc, "hmcode_max_k_extra", hmcode_max_k_extra);
  read(fc, "hmcode_min_k_max", hmcode_min_k_max);
  read(fc, "hmcode_tol_sigma", hmcode_tol_sigma);
  read(fc, "n_hmcode_tables", n_hmcode_tables);
  read(fc, "rmin_for_sigtab", rmin_for_sigtab);
  read(fc, "rmax_for_sigtab", rmax_for_sigtab);
  read(fc, "ainit_for_growtab", ainit_for_growtab);
  read(fc, "amax_for_growtab", amax_for_growtab);
  read(fc, "nsteps_for_p1h_integral", nsteps_for_p1h_integral);
  read(fc, "mmin_for_p1h_integral", mmin_for_p1h_integral);
  read(fc, "mmax_for_p1h_integral", mmax_for_p1h_integral);

  /* Lensing */
  read(fc, "accurate_lensing", accurate_lensing);
  read(fc, "num_mu_minus_lmax", num_mu_minus_lmax);
  read(fc, "delta_l_max", delta_l_max);
  read(fc, "tol_gauss_legendre", tol_gauss_legendre);
}

int class_version(char* version) {
  snprintf(version, 12, "%s", _VERSION_);
  return _SUCCESS_;
}

// ── Hook-based shooting (the species-owned replacement for the enum dispatch) ──

int InputModule::ShootingResidual(
    double* x, int x_size, void* pworkspace, double* output, ErrorMsg /*error_message*/) {
  auto* w = static_cast<ShootingWorkspace*>(pworkspace);

  // Write each trial unknown into the file content.
  for (int i = 0; i < x_size; ++i) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%.17g", x[i]);
    w->fc.set(w->targets[i].unknown_param, buf);
  }
  // Mark as a shooting build so the lazily-evaluated module does not re-enter DoShooting.
  w->fc.is_shooting = true;

  Cosmology cosmology{std::make_unique<InputModule>(w->fc)};
  BackgroundModulePtr bgm  = cosmology.GetBackgroundModule();
  const double* bg_today   = bgm->background_table_.data() + (bgm->bt_size_ - 1) * bgm->bg_size_;
  const InputModulePtr& im = cosmology.GetInputModule();

  bool need_thermo = false;
  for (const auto& t : w->targets)
    if (t.target_name == "100*theta_s")
      need_thermo = true;
  ThermodynamicsModulePtr thm;
  if (need_thermo)
    thm = cosmology.GetThermodynamicsModule();

  // Assemble residuals in workspace order. theta_s is module-level (from thermo rs/ra);
  // each species target is routed to its owning species via the recorded key, passing the
  // authoritative target collected at discovery (the species must not re-derive it — in
  // this iteration build DoShooting has set the unknown, so for species with overlapping
  // target/unknown keys the user's target is not recoverable from the file content).
  ShootingResidualContext ctx{&im->background_, bg_today};
  for (size_t i = 0; i < w->targets.size(); ++i) {
    const ShootingTarget& t = w->targets[i];
    if (t.target_name == "100*theta_s") {
      output[i] = 100. * thm->rs_rec_ / thm->ra_rec_ - t.target_value;
    }
    else {
      const auto& sp = im->all_species_.at(w->target_species_keys[i]);
      output[i]      = sp->ComputeShootingResidual(ctx, t);
    }
  }
  return _SUCCESS_;
}

InputModulePtr InputModule::DoShooting(InputModulePtr input_module) {
  // Guard: never shoot from within a shooting build (the residual marks its fc this way).
  if (input_module->file_content_.is_shooting)
    return input_module;

  FileContent& fc = input_module->file_content_;
  ErrorMsg errmsg = "";  // initialized so a failure that skips writing it can't surface garbage

  ShootingWorkspace w(fc);
  std::vector<double> xguess, dxdF;

  // Cosmological target: 100*theta_s varies h (residual from thermo rs/ra). Module-level.
  {
    double tv;
    int flag = _FALSE_;
    parser_read_double(&fc, "100*theta_s", &tv, &flag, errmsg);
    if (flag == _TRUE_) {
      w.targets.push_back({"100*theta_s", "h", tv});
      w.target_species_keys.emplace_back();  // module-level: no owning species
      xguess.push_back(3.54 * tv * tv - 5.455 * tv + 2.548);
      dxdF.push_back(7.08 * tv - 5.455);
    }
  }

  // Per-species targets (all_species_ lex order). The discovery module's species already
  // guessed their unknowns at construction; query them for {target, guess} here.
  // SpeciesBuildContext::ncdm_settings is non-null by contract (some guesses fall back to
  // it when ppr is unavailable); build it from the module as ConstructSpecies does.
  NcdmSettings ncdm_settings;
  ncdm_settings.h           = input_module->background_.h;
  ncdm_settings.T_cmb       = input_module->background_.T_cmb;
  ncdm_settings.tol_ncdm    = input_module->precision_.tol_ncdm;
  ncdm_settings.tol_ncdm_bg = input_module->precision_.tol_ncdm_bg;
  ncdm_settings.tol_M_ncdm  = input_module->precision_.tol_M_ncdm;
  const SpeciesBuildContext gctx{&fc,
                                 &input_module->background_,
                                 &input_module->precision_,
                                 &ncdm_settings,
                                 /*bgm=*/nullptr,
                                 /*all_species=*/&input_module->all_species_,
                                 /*omega_budget=*/&input_module->omega_budget_};
  for (const auto& [key, sp] : input_module->all_species_) {
    std::vector<ShootingTarget> tgts = sp->GetShootingTargets();
    if (tgts.empty())
      continue;
    std::vector<double> g, d;
    sp->ComputeShootingGuess(gctx, g, d);
    // A species must report one guess + one Jacobian seed per target (same order), or the
    // flattened unknown vector below desyncs. Guard against a mis-implemented hook.
    if (g.size() != tgts.size() || d.size() != tgts.size()) {
      throw std::runtime_error("species '" + key + "' reported " + std::to_string(tgts.size()) +
                               " shooting target(s) but " + std::to_string(g.size()) +
                               " guess(es) / " + std::to_string(d.size()) + " Jacobian seed(s)");
    }
    for (size_t j = 0; j < tgts.size(); ++j) {
      w.targets.push_back(tgts[j]);
      w.target_species_keys.push_back(key);
      xguess.push_back(g[j]);
      dxdF.push_back(d[j]);
    }
  }

  if (w.targets.empty())
    return input_module;  // nothing to shoot — the common path

  // Solve (fzero_Newton handles n>=1 and writes the solution back into xguess).
  fc.is_shooting = true;
  int fevals     = 0;
  if (fzero_Newton(ShootingResidual,
                   xguess.data(),
                   dxdF.data(),
                   static_cast<int>(w.targets.size()),
                   1e-3,
                   1e-3,
                   &w,
                   &fevals,
                   errmsg) != _SUCCESS_) {
    throw std::runtime_error(std::string("Shooting (DoShooting) failed: ") + errmsg);
  }

  // Write the resolved unknowns; build and return a fresh, fully-resolved module.
  for (size_t i = 0; i < w.targets.size(); ++i) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%.17g", xguess[i]);
    fc.set(w.targets[i].unknown_param, buf);
  }
  return std::make_shared<InputModule>(fc);
}

int input_prepare_pk_eq(const struct precision* ppr_input,
                        const struct background* pba_input,
                        const struct thermo* pth_input,
                        struct nonlinear* pnl,
                        int input_verbose,
                        ErrorMsg errmsg) {
  return _SUCCESS_;
}
