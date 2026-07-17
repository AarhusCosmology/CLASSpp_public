/** @file input_module.cpp
 * Parses CLASS input and constructs the shared input data and species set.
 */

#include "input_module.h"

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <sstream>

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
#include "../species/ncdm_family.h"
#include "../species/species_input.h"

namespace {

void readDoubleList(FileContent* pfc, const char* name, std::vector<double>& values, int* found) {
  try {
    if (auto l = pfc->get<std::vector<double>>(name)) {
      values = *l;
      *found = true;
    }
    else
      *found = false;
  }
  catch (const std::exception& e) {
    class_stop("%s", e.what());
  }
}

// Number of the given optionals that hold a value.
template <class... Opt>
int n_present(const Opt&... opts) {
  return (0 + ... + (opts.has_value() ? 1 : 0));
}

// Value of whichever of two names is present; error if BOTH present; nullopt if neither.
template <class T>
std::optional<T> read_one_of(const FileContent& fc, const char* n1, const char* n2) {
  auto a = fc.get<T>(n1);
  auto b = fc.get<T>(n2);
  class_test(a && b, "In input file, you can only enter one of %s, %s, choose one", n1, n2);
  return a ? a : b;
}

}  // namespace
/**
 * Use this routine to extract initial parameters from files 'xxx.ini'
 * and/or 'xxx.pre'. They can be the arguments of the main() routine.
 *
 * If class is embedded into another code, you will probably prefer to
 * populate a 'file_content' structure directly and pass it to the
 * InputModule constructor, which then runs input_read_precisions,
 * ReadContext, ConstructSpecies, ReadDerived,
 * and WriteParameterFiles in sequence.
 */

void InputModule::file_content_from_arguments(int argc, char** argv, FileContent& fc) {
  /** Summary: */

  /** - define local variables */
  FileContent fc_input;     /** - --> a temporary structure with all input parameters */
  FileContent fc_precision; /** - --> a temporary structure with all precision parameters */
  FileContent fc_root;      /** - --> a temporary structure with only the root name */
  FileContent fc_inputroot; /** - --> sum of fc_inoput and fc_root */
  FileContent* pfc_input;   /** - --> a pointer to either fc_root or fc_inputroot */

  std::string input_file;
  std::string precision_file;
  std::string tmp_file;
  std::string extension;
  std::string inifilename;
  int filenum;
  bool root_set = false;

  pfc_input = &fc_input;

  /** - Initialize the two file_content structures (for input
      parameters and precision parameters) to some null content. If no
      arguments are passed, they will remain null and inform
      init_params() that all parameters take default values. */

  fc           = FileContent();
  fc_input     = FileContent();
  fc_precision = FileContent();
  input_file.clear();
  precision_file.clear();

  /** - If some arguments are passed, identify eventually some 'xxx.ini'
      and 'xxx.pre' files, and store their name. */

  if (argc > 1) {
    for (int i = 1; i < argc; i++) {
      {
        size_t arglen = strlen(argv[i]);
        extension     = (arglen >= 4) ? std::string(argv[i] + arglen - 4) : std::string(argv[i]);
      }
      if (extension == ".ini") {
        class_test(!input_file.empty(),
                   "You have passed more than one input file with extension '.ini', choose one.");
        input_file = argv[i];
      }
      else if (extension == ".pre") {
        class_test(!precision_file.empty(),
                   "You have passed more than one precision with extension '.pre', choose one.");
        precision_file = argv[i];
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

  if (!input_file.empty()) {
    parser_read_file(input_file.c_str(), &fc_input);

    /** - check whether a root name has been set */

    // Only presence matters here (was "root" set?); the value is read again in ReadDerived.
    root_set = fc_input.get<std::string>("root").has_value();

    /** - if root has not been set, use root=output/inputfilennameN_ */

    if (!root_set) {
      inifilename = input_file.substr(0, input_file.size() - 4);
      for (filenum = 0; filenum < 100; filenum++) {
        std::error_code ec;
        {
          std::ostringstream oss;
          oss << "output/" << inifilename << std::setw(2) << std::setfill('0') << filenum
              << "_cl.dat";
          tmp_file = oss.str();
        }
        if (std::filesystem::exists(tmp_file, ec))
          continue;
        {
          std::ostringstream oss;
          oss << "output/" << inifilename << std::setw(2) << std::setfill('0') << filenum
              << "_pk.dat";
          tmp_file = oss.str();
        }
        if (std::filesystem::exists(tmp_file, ec))
          continue;
        {
          std::ostringstream oss;
          oss << "output/" << inifilename << std::setw(2) << std::setfill('0') << filenum
              << "_tk.dat";
          tmp_file = oss.str();
        }
        if (std::filesystem::exists(tmp_file, ec))
          continue;
        {
          std::ostringstream oss;
          oss << "output/" << inifilename << std::setw(2) << std::setfill('0') << filenum
              << "_parameters.ini";
          tmp_file = oss.str();
        }
        if (std::filesystem::exists(tmp_file, ec))
          continue;
        break;
      }
      {
        std::ostringstream oss;
        oss << "output/" << inifilename << std::setw(2) << std::setfill('0') << filenum << "_";
        tmp_file = oss.str();
      }
      fc_root.set("root", tmp_file);
      parser_cat(&fc_input, &fc_root, &fc_inputroot);
      pfc_input = &fc_inputroot;
    }
  }

  /** - if there is an 'xxx.pre' file, read it and store its content. */

  if (!precision_file.empty())
    parser_read_file(precision_file.c_str(), &fc_precision);

  /** - if one or two files were read, merge their contents in a
      single 'file_content' structure. */

  if (!input_file.empty() || !precision_file.empty())
    parser_cat(pfc_input, &fc_precision, &fc);
}

InputModule::InputModule(FileContent& fc) : file_content_(fc) {
  file_content_.mark_all_unread();
  try {
    // Translate dot-syntax for single-instance legacy species (bary.Omega ->
    // Omega_b) before any consumer reads the file.
    TranslateSingleInstanceDotSyntax(&file_content_);
    input_read_precisions();
    ReadContext();
    ConstructSpecies();
    ReadDerived();
    WriteParameterFiles();
  }
  catch (const std::runtime_error& e) {
    throw std::invalid_argument(e.what());
  }
}

void InputModule::ConstructSpecies() {
  background* pba = &background_;  // need non-const to write closure Omega

  NcdmSettings ncdm_settings;
  ncdm_settings.h     = pba->h;
  ncdm_settings.T_cmb = pba->T_cmb;
  // The ncdm perturbation momentum grid is gauge-dependent: synchronous gauge
  // converges with far fewer q-bins than Newtonian (synchronous reaches sub-0.1%
  // P(k) at ~5 bins, Newtonian needs ~11). Use the gauge-appropriate tolerance.
  ncdm_settings.tol_ncdm    = (perturbations_.gauge == possible_gauges::newtonian)
                                  ? precision_.tol_ncdm_newtonian
                                  : precision_.tol_ncdm_synchronous;
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
  ctx.coupled_inputs = &coupled_inputs_;

  // Read input_verbose for the closure verbose message.
  int input_verbose = file_content_.get_or("input_verbose", 0);

  const std::string_view closure_name = ClosureSpeciesName(pba->closure_species);

  // Pass 1: build every non-closure species, summing Omega0 contributions.
  double omega0_sum = 0.0;
  for (const auto& entry : kAllSpeciesFactories) {
    if (entry.name == closure_name)
      continue;
    auto produced = entry.create_all(ctx);
    for (auto& e : produced) {
      omega0_sum += e.species->GetOmega0();
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
  // species is present (any NCDMBaseSpecies, or the DNCDM_DR composite).
  const precision* ppr = &precision_;
  if (all_species_.has_ncdm()) {
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

void InputModule::ReadCoupledCluster() {
  FileContent* pfc    = &file_content_;
  precision* ppr      = &precision_;
  background* pba     = &background_;
  thermo* /*pth*/ pth = &thermodynamics_;
  perturbs* ppt       = &perturbations_;

  (void) pth;  // future-proof; the budget computation itself doesn't touch pth.

  int input_verbose = pfc->get_or("input_verbose", 0);

  // ── IDR: stat_f_idr * (T_idr / T_cmb)^4 * Omega0_g ───────────────────────────
  // T_idr is a physics param owned by IDRSpecies; compute a local copy here
  // solely to derive the omega budget entry for IDR.
  double stat_f_idr = 7. / 8.;
  stat_f_idr        = pfc->get_or("stat_f_idr", stat_f_idr);

  auto N_idr  = pfc->get<double>("N_idr");
  auto N_dg   = pfc->get<double>("N_dg");
  auto xi_idr = pfc->get<double>("xi_idr");
  class_test(n_present(N_idr, N_dg, xi_idr) > 1,
             "In input file, you can only enter one of N_idr, N_dg or xi_idr, choose one");

  double T_idr_local = 0.;
  if (N_idr) {
    T_idr_local = pow(*N_idr / stat_f_idr * (7. / 8.) / pow(11. / 4., (4. / 3.)), (1. / 4.)) *
                  pba->T_cmb;
    if (input_verbose > 1)
      printf(
          "You passed N_idr = N_dg = %e, this is equivalent to xi_idr = %e in the ETHOS notation. "
          "\n",
          *N_idr,
          T_idr_local / pba->T_cmb);
  }
  else if (N_dg) {
    T_idr_local = pow(*N_dg / stat_f_idr * (7. / 8.) / pow(11. / 4., (4. / 3.)), (1. / 4.)) *
                  pba->T_cmb;
    if (input_verbose > 2)
      printf(
          "You passed N_dg = N_idr = %e, this is equivalent to xi_idr = %e in the ETHOS notation. "
          "\n",
          *N_dg,
          T_idr_local / pba->T_cmb);
  }
  else if (xi_idr) {
    T_idr_local = *xi_idr * pba->T_cmb;
    if (input_verbose > 1)
      printf(
          "You passed xi_idr = %e, this is equivalent to N_idr = N_dg = %e in the NADM notation. "
          "\n",
          *xi_idr,
          stat_f_idr * pow(*xi_idr, 4.) / (7. / 8.) * pow(11. / 4., (4. / 3.)));
  }

  // Store IDR intermediates into coupled_inputs_ for factory reuse.
  coupled_inputs_.stat_f_idr = stat_f_idr;
  if (N_idr || N_dg || xi_idr)
    coupled_inputs_.T_idr = T_idr_local;

  // Mark the budget slot present iff any of the three IDR-temperature inputs was given.
  // (Matches the legacy semantics: pba->Omega0_idr was always written, but only became
  // nonzero when T_idr was set by one of these inputs.)
  if (N_idr || N_dg || xi_idr) {
    omega_budget_.idr = stat_f_idr * pow(T_idr_local / pba->T_cmb, 4.) * pba->Omega0_g;
  }

  // ── CDM: parser value (or default), then synchronous-gauge minimum ───────────
  // Default CDM fallback, frozen at the default h=0.67556 (matches the historical
  // classyref value 0.12038/h^2). Must NOT use pba->h here: ReadCoupledCluster
  // runs after the user's h/H0 is read, so dividing by the live h would let the
  // default drift when h is set or 100*theta_s is shot
  // (the omega_b/Omega0_g defaults are likewise frozen at default h in background.h).
  double omega0_cdm = 0.12038 / (0.67556 * 0.67556);
  bool cdm_user_set = false;

  auto Omega_cdm = pfc->get<double>("Omega_cdm");
  auto omega_cdm = pfc->get<double>("omega_cdm");
  class_test((Omega_cdm && omega_cdm),
             "In input file, you can only enter one of Omega_cdm or omega_cdm, choose one");
  if (Omega_cdm) {
    omega0_cdm   = *Omega_cdm;
    cdm_user_set = true;
  }
  if (omega_cdm) {
    omega0_cdm   = *omega_cdm / pba->h / pba->h;
    cdm_user_set = true;
  }

  // Synchronous-gauge minimum: if CDM ends up zero in synchronous gauge, bump it
  // to Omega0_cdm_min_synchronous. Track presence: the gauge minimum kicks in
  // even if the user supplied nothing, so the budget slot becomes present.
  bool cdm_present = cdm_user_set || omega0_cdm != 0.;
  if ((ppt->gauge == possible_gauges::synchronous) && (omega0_cdm == 0.)) {
    omega0_cdm  = ppr->Omega0_cdm_min_synchronous;
    cdm_present = true;
  }
  if (cdm_present)
    omega_budget_.cdm = omega0_cdm;

  // ── IDM_DR: Omega_idm_dr / omega_idm_dr / f_idm_dr ───────────────────────────
  auto Omega_idm_dr = pfc->get<double>("Omega_idm_dr");
  auto omega_idm_dr = pfc->get<double>("omega_idm_dr");
  auto f_idm_dr     = pfc->get<double>("f_idm_dr");
  class_test(n_present(Omega_idm_dr, omega_idm_dr, f_idm_dr) > 1,
             "In input file, you can only enter one of Omega_idm_dr, omega_idm_dr or f_idm_dr, "
             "choose one");

  if (Omega_idm_dr)
    omega_budget_.idm_dr = *Omega_idm_dr;
  if (omega_idm_dr)
    omega_budget_.idm_dr = *omega_idm_dr / pba->h / pba->h;

  if (f_idm_dr) {
    class_test((*f_idm_dr < 0.) || (*f_idm_dr > 1.),
               "The fraction of interacting DM with DR must be between 0 and 1, you asked for "
               "f_idm_dr=%e",
               *f_idm_dr);
    const double cdm_for_frac = omega_budget_.cdm.value_or(0.);
    class_test((*f_idm_dr > 0.) && (cdm_for_frac == 0.),
               "If you want a fraction of interacting DM with DR, to be consistent, you should not "
               "set the fraction of CDM to zero");

    double new_idm_dr = *f_idm_dr * cdm_for_frac;
    double new_cdm    = cdm_for_frac - new_idm_dr;
    // avoid CDM=0 in synchronous gauge after the subtraction
    if ((ppt->gauge == possible_gauges::synchronous) && (new_cdm == 0.)) {
      new_cdm    += ppr->Omega0_cdm_min_synchronous;
      new_idm_dr -= ppr->Omega0_cdm_min_synchronous;
    }
    omega_budget_.idm_dr = new_idm_dr;
    omega_budget_.cdm    = new_cdm;
  }

  // ── DCDM_DR: Omega_dcdmdr / omega_dcdmdr ────────────────────────────────────
  // Gamma_dcdm and Omega_ini_dcdm are physics params owned by DCDMSpecies;
  // they are parsed in DCDM_DR_Species::CreateAll, not here.
  auto Omega_dcdmdr = pfc->get<double>("Omega_dcdmdr");
  auto omega_dcdmdr = pfc->get<double>("omega_dcdmdr");
  class_test((Omega_dcdmdr && omega_dcdmdr),
             "In input file, you can only enter one of Omega_dcdmdr or omega_dcdmdr, choose one");
  if (Omega_dcdmdr)
    omega_budget_.dcdmdr = *Omega_dcdmdr;
  if (omega_dcdmdr)
    omega_budget_.dcdmdr = *omega_dcdmdr / pba->h / pba->h;

  // ── DRMD: z_stop, G_over_aH_drmd_ini, f_idm_drmd, delta_Neff_drmd ────────────
  // These fields now live on IDM_DRMD_IDR_DRMD_Species; we parse them locally here
  // only for the budget math (delta_Neff_drmd → omega_budget_.idr_drmd;
  // f_idm_drmd → omega_budget_.idm_drmd with the CDM subtraction).
  auto z_stop          = pfc->get<double>("z_stop");
  auto G_over_aH_drmd  = pfc->get<double>("G_over_aH_drmd_ini");
  auto f_idm_drmd      = pfc->get<double>("f_idm_drmd");
  auto delta_Neff_drmd = pfc->get<double>("delta_Neff_drmd");

  const double z_stop_drmd           = z_stop ? *z_stop : 0.;
  const double f_idm_drmd_local      = f_idm_drmd ? *f_idm_drmd : 0.;
  const double delta_Neff_drmd_local = delta_Neff_drmd ? *delta_Neff_drmd : 0.;

  const int any_drmd = z_stop || G_over_aH_drmd || f_idm_drmd || delta_Neff_drmd;
  const int all_drmd = z_stop && G_over_aH_drmd && f_idm_drmd && delta_Neff_drmd;
  class_test(any_drmd && !all_drmd,
             "If any DRMD parameter is set, all of them must be non-zero.\nDRMD parameters are "
             "'z_stop', 'G_over_aH_drmd_ini', 'f_idm_drmd' and 'delta_Neff_drmd'.");

  // Store DRMD intermediates into coupled_inputs_ for factory reuse.
  coupled_inputs_.z_stop             = z_stop_drmd;
  coupled_inputs_.f_idm_drmd         = f_idm_drmd_local;
  coupled_inputs_.delta_Neff_drmd    = delta_Neff_drmd_local;
  coupled_inputs_.G_over_aH_drmd_ini = G_over_aH_drmd ? *G_over_aH_drmd : 0.;

  if (delta_Neff_drmd_local > 0.) {
    omega_budget_.idr_drmd = delta_Neff_drmd_local * 7. / 8. * pow(4. / 11., 4. / 3.) *
                             pba->Omega0_g;
    if (f_idm_drmd_local > 0) {
      class_test((z_stop_drmd > 200000.),
                 "z_stop is chosen too large. If you want to probe z_stop > 1000000 you need to "
                 "start evolving perturbations earlier in CLASS by changing the precision "
                 "settings. Also you should check that the exponential suppression factor does not "
                 "lead to numerical problems.");
    }
  }

  if (f_idm_drmd_local > 0) {
    class_test((f_idm_drmd_local > 1.),
               "The fraction of interacting DM with DR must be between 0 and 1, you asked for "
               "f_idm_drmd=%e",
               f_idm_drmd_local);
    const double cdm_for_drmd = omega_budget_.cdm.value_or(0.);
    class_test((cdm_for_drmd == 0.),
               "If you want a fraction of interacting DM with DRMD, to be consistent, you should "
               "not set the fraction of CDM to zero");

    double new_idm_drmd = f_idm_drmd_local * cdm_for_drmd;
    double new_cdm      = cdm_for_drmd - new_idm_drmd;
    if ((ppt->gauge == possible_gauges::synchronous) && (new_cdm == 0.)) {
      new_cdm      += ppr->Omega0_cdm_min_synchronous;
      new_idm_drmd -= ppr->Omega0_cdm_min_synchronous;
    }
    omega_budget_.idm_drmd = new_idm_drmd;
    omega_budget_.cdm      = new_cdm;
  }
}

/**
 * Write the read parameters to a file, the unread parameters to another
 * file, and warnings about unread parameters. Runs after ConstructSpecies
 * and ReadDerived, so species inputs read inside species factories are
 * correctly classified as read.
 *
 */

void InputModule::WriteParameterFiles() {
  FileContent* pfc = &file_content_;

  /** - eventually write all the read parameters in a file, unread parameters in another file, and warnings about unread parameters */

  if (auto write_parameters = pfc->get<std::string>("write parameters");
      write_parameters && ((write_parameters->find("y") != std::string::npos) ||
                           (write_parameters->find("Y") != std::string::npos))) {
    output* pop                   = &output_;
    std::string param_output_name = pop->root + "parameters.ini";
    std::string param_unused_name = pop->root + "unused_parameters";

    FILE* param_output;
    FILE* param_unused;
    class_open(param_output, param_output_name.c_str(), "w");
    class_open(param_unused, param_unused_name.c_str(), "w");

    fprintf(param_output, "# List of input/precision parameters actually read\n");
    fprintf(param_output, "# (all other parameters set to default values)\n");
    fprintf(param_output, "# Obtained with CLASS %s\n", _VERSION_);
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

  if (auto write_warnings = pfc->get<std::string>("write warnings");
      write_warnings && ((write_warnings->find("y") != std::string::npos) ||
                         (write_warnings->find("Y") != std::string::npos))) {
    pfc->for_each([](const std::string& name, const std::string& value, bool read) {
      if (!read)
        fprintf(stdout,
                "[WARNING: input line not recognized and not taken into account: '%s=%s']\n",
                name.c_str(),
                value.c_str());
    });
  }
}
void InputModule::input_read_precisions() {
  precision* ppr = &precision_;

  if (auto class_dir = file_content_.get<std::string>("class_dir"))
    ppr->class_dir = *class_dir;
  else
    ppr->class_dir = __CLASSDIR__;

  /** - resolve runtime data-file paths against class_dir (defaults are relative) */
  ppr->ResolveDataPaths();

  SynthesiseNcdmFluidApproximation(&file_content_);

  /** - parse all precision parameters from config file */
  ppr->parse(file_content_);

  class_test(ppr->smallest_allowed_variation < 0,
             "smallest_allowed_variation = %e < 0",
             ppr->smallest_allowed_variation);
}
/**
 * Phase i: read the non-species inputs that building the species needs
 * (gauge, h, photon/baryon densities, the coupled-Omega budget, closure
 * selection). Runs before ConstructSpecies.
 */
void InputModule::ReadContext() {
  /** Summary: */

  /** - define local variables */
  FileContent* pfc = &file_content_;
  precision* ppr   = &precision_;      /* for precision parameters */
  background* pba  = &background_;     /* for cosmological background */
  thermo* pth      = &thermodynamics_; /* for thermodynamics */
  perturbs* ppt    = &perturbations_;  /* for source functions */

  int input_verbose = 0;

  /** - if entries passed in file_content structure, carefully read
      and interpret each of them, and tune the relevant input
      parameters accordingly*/

  input_verbose = pfc->get_or("input_verbose", input_verbose);
  if (input_verbose > 0)
    printf("Reading input parameters\n");

  if (auto threads = pfc->get<int>("threads")) {
    pba->number_of_threads = *threads;
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

  if (auto gauge = pfc->get<std::string>("gauge")) {
    if ((gauge->find("newtonian") != std::string::npos) ||
        (gauge->find("Newtonian") != std::string::npos) ||
        (gauge->find("new") != std::string::npos)) {
      ppt->gauge = possible_gauges::newtonian;
    }

    if ((gauge->find("synchronous") != std::string::npos) ||
        (gauge->find("sync") != std::string::npos) ||
        (gauge->find("Synchronous") != std::string::npos)) {
      ppt->gauge = possible_gauges::synchronous;
    }
  }

  /** (a) background parameters */

  /** - h (dimensionless) and [\f$ H_0/c\f$] in \f$ Mpc^{-1} = h / 2997.9... = h * 10^5 / c \f$ */
  auto H0 = pfc->get<double>("H0");
  auto h  = pfc->get<double>("h");
  class_test(H0 && h, "In input file, you cannot enter both h and H0, choose one");
  if (H0) {
    pba->H0 = *H0 * 1.e3 / _c_;
    pba->h  = *H0 / 100.;
  }
  if (h) {
    pba->H0 = *h * 1.e5 / _c_;
    pba->h  = *h;
  }

  /** - 100*theta_s is a module-level shooting target (varies h; resolved later by
   *  DoShooting). Consume it here so the unread-parameter check passes, and reject
   *  combining it with a direct h/H0 — except inside a shooting build, where DoShooting
   *  has itself set the trial h. h keeps its default; DoShooting seeds and solves it. */
  auto theta_s = pfc->get<double>("100*theta_s");
  class_test(theta_s && (H0 || h) && !pfc->is_shooting,
             "In input file, you cannot enter both 100*theta_s and h (or H0), choose one");

  /** - Omega_0_g (photons) and T_cmb */
  auto T_cmb   = pfc->get<double>("T_cmb");
  auto Omega_g = pfc->get<double>("Omega_g");
  auto omega_g = pfc->get<double>("omega_g");
  class_test(n_present(T_cmb, Omega_g, omega_g) > 1,
             "In input file, you can only enter one of T_cmb, Omega_g or omega_g, choose one");

  if (n_present(T_cmb, Omega_g, omega_g) == 0) {
    pba->Omega0_g = Omega0gFromTcmb(pba->T_cmb, pba->h);
  }
  else {
    if (T_cmb) {
      /** - Omega0_g = rho_g / rho_c0, each of them expressed in \f$ Kg/m/s^2 \f$*/
      /** - rho_g = (4 sigma_B / c) \f$ T^4 \f$*/
      /** - rho_c0 \f$ = 3 c^2 H_0^2 / (8 \pi G) \f$*/
      pba->Omega0_g = Omega0gFromTcmb(*T_cmb, pba->h);
      pba->T_cmb    = *T_cmb;
    }

    if (Omega_g) {
      pba->Omega0_g = *Omega_g;
      pba->T_cmb    = TcmbFromOmega0g(pba->Omega0_g, pba->h);
    }

    if (omega_g) {
      pba->Omega0_g = *omega_g / pba->h / pba->h;
      pba->T_cmb    = TcmbFromOmega0g(pba->Omega0_g, pba->h);
    }
  }

  /** - Omega_0_b (baryons) */
  auto Omega_b = pfc->get<double>("Omega_b");
  auto omega_b = pfc->get<double>("omega_b");
  class_test((Omega_b && omega_b),
             "In input file, you can only enter one of Omega_b or omega_b, choose one");
  if (Omega_b)
    pba->Omega0_b = *Omega_b;
  if (omega_b)
    pba->Omega0_b = *omega_b / pba->h / pba->h;

  // NOTE: N_ur / N_eff / Omega_ur / omega_ur parsing has been moved to
  // UltraRelativisticSpecies::CreateAll, which reads pfc directly.

  if (auto ceff2_ur = pfc->get<double>("ceff2_ur"))
    ppt->three_ceff2_ur = 3. * *ceff2_ur;

  if (auto cvis2_ur = pfc->get<double>("cvis2_ur"))
    ppt->three_cvis2_ur = 3. * *cvis2_ur;

  auto G_eff_ur       = pfc->get<double>("G_eff_ur");
  auto log10_G_eff_ur = pfc->get<double>("log10_G_eff_ur");
  if (G_eff_ur)
    ppt->G_eff_ur = *G_eff_ur;
  if (log10_G_eff_ur) {
    ppt->G_eff_ur = *log10_G_eff_ur;
    class_test(G_eff_ur.has_value(),
               "In input file, you cannot enter both log10_G_eff_ur and G_eff_ur, choose one");
    ppt->G_eff_ur = pow(10.0, ppt->G_eff_ur);
  }

  // Coupled-species cluster (CDM, IDR, IDM_DR, DCDM_DR, IDM_DRMD, IDR_DRMD)
  // ReadCoupledCluster computes the coupled-Omega budget and fills
  // CoupledClusterInputs; IDM/IDR interaction physics parameters
  // (a_idm_dr/nindex/idr_nature/m_idm/b_idr/alpha/beta) are parsed
  // exclusively in IDM_DR_IDR_Species::CreateAll.
  ReadCoupledCluster();

  // The idm/idr interaction parameters (a_idm_dr/nindex/idr_nature/m_idm/b_idr/
  // alpha/beta) and the "non-zero IDR density" validation are now owned and
  // parsed entirely by IDM_DR_IDR_Species::CreateAll. The transitional pth/ppt
  // copies have been removed.

  // Omega_dcdmdr parsing has been moved to ReadCoupledCluster
  // (Omega_dcdmdr → omega_budget_.dcdmdr); Gamma_dcdm and Omega_ini_dcdm are
  // owned by DCDMSpecies and parsed in DCDM_DR_Species::CreateAll.
  // T_idr and l_max_idr are owned by IDRSpecies and parsed in
  // IDM_DR_IDR_Species::CreateAll.

  /** - Omega_0_k (effective fractional density of curvature) */
  pba->Omega0_k = pfc->get_or("Omega_k", pba->Omega0_k);
  /** - Set curvature parameter K */
  pba->K = -pba->Omega0_k * pow(pba->H0, 2);
  /** - Set curvature sign */
  if (pba->K > 0.)
    pba->sgnK = 1;
  else if (pba->K < 0.)
    pba->sgnK = -1;

  // DRMD parameter block (z_stop, G_over_aH_drmd_ini, f_idm_drmd, delta_Neff_drmd,
  // and the resulting Omega0_idr_drmd / Omega0_idm_drmd contributions) has been
  // moved to ReadCoupledCluster.

  /** - Omega_0_lambda (cosmological constant), Omega0_fld (dark energy fluid), Omega0_scf (scalar field) */

  auto Omega_Lambda = pfc->get<double>("Omega_Lambda");
  auto Omega_fld    = pfc->get<double>("Omega_fld");
  auto Omega_scf    = pfc->get<double>("Omega_scf");

  class_test(Omega_Lambda && Omega_fld && (!Omega_scf || (*Omega_scf >= 0.)),
             "In input file, either Omega_Lambda or Omega_fld must be left unspecified, except if "
             "Omega_scf is set and <0.0, in which case the contribution from the scalar field will "
             "be the free parameter.");

  /** - --> (!Omega_scf || (*Omega_scf >= 0.)) explained:
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
  if (!Omega_Lambda) {
    pba->closure_species = ClosureSpecies::Lambda;
  }
  else if (!Omega_fld) {
    pba->closure_species = ClosureSpecies::Fluid;
  }
  else if (Omega_scf && (*Omega_scf < 0.)) {
    pba->closure_species = ClosureSpecies::ScalarField;
  }
  else {
    pba->closure_species = ClosureSpecies::None;
  }

  /** - Test that the user have not specified Omega_scf = -1 but left either
      Omega_lambda or Omega_fld unspecified:*/
  class_test((!Omega_Lambda || !Omega_fld) && (Omega_scf && (*Omega_scf < 0.)),
             "It looks like you want to fulfil the closure relation sum Omega = 1 using the scalar "
             "field, so you have to specify both Omega_lambda and Omega_fld in the .ini file");

  /* Fluid physics params (fluid_equation_of_state, w0_fld, wa_fld, cs2_fld,
     Omega_EDE) are parsed inside FluidSpecies::CreateAll; PPF-only params
     (use_ppf, c_gamma_over_c_fld) are parsed there too but live on PpfFluid.
     No per-key writes to pba here.

     Scalar field physics params (scf_parameters, scf_tuning_index,
     scf_shooting_parameter, attractor_ic_scf, phi_ini_scf, phi_prime_ini_scf)
     are parsed inside ScalarFieldSpecies::CreateAll directly from pfc; no
     per-key writes to pba here. */
}

/**
 * Phase iii: read everything that does not feed species construction —
 * thermodynamics, reionization, energy injection, output / perturbation /
 * primordial / transfer / spectra / lensing configuration, and the
 * species-dependent S8 and halofit reads. Runs after ConstructSpecies.
 */
void InputModule::ReadDerived() {
  /** - define local variables */
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

  int n_list;
  double k1    = 0.;
  double k2    = 0.;
  double prr1  = 0.;
  double prr2  = 0.;
  double pii1  = 0.;
  double pii2  = 0.;
  double pri1  = 0.;
  double pri2  = 0.;
  double n_iso = 0.;
  double f_iso = 0.;
  double n_cor = 0.;
  double c_cor = 0.;

  int i;

  double z_max = 0.;

  /** (b) assign values to thermodynamics cosmological parameters */

  /** - primordial helium fraction */
  if (auto YHe = pfc->get<std::string>("YHe")) {
    if ((YHe->find("BBN") != std::string::npos) || (YHe->find("bbn") != std::string::npos)) {
      pth->YHe = _BBN_;
    }
    else {
      pth->YHe = pfc->get_or("YHe", pth->YHe);
    }
  }

  /** - recombination parameters */
  if (auto recombination = pfc->get<std::string>("recombination")) {
    if ((recombination->find("HYREC") != std::string::npos) ||
        (recombination->find("hyrec") != std::string::npos) ||
        (recombination->find("HyRec") != std::string::npos)) {
      pth->recombination = hyrec;
    }
  }

  /** - reionization parametrization */
  if (auto reio_parametrization = pfc->get<std::string>("reio_parametrization")) {
    bool recognized = false;
    if (*reio_parametrization == "reio_none") {
      pth->reio_parametrization = reio_none;
      recognized                = true;
    }
    if (*reio_parametrization == "reio_camb") {
      pth->reio_parametrization = reio_camb;
      recognized                = true;
    }
    if (*reio_parametrization == "reio_bins_tanh") {
      pth->reio_parametrization = reio_bins_tanh;
      recognized                = true;
    }
    if (*reio_parametrization == "reio_half_tanh") {
      pth->reio_parametrization = reio_half_tanh;
      recognized                = true;
    }
    if (*reio_parametrization == "reio_many_tanh") {
      pth->reio_parametrization = reio_many_tanh;
      recognized                = true;
    }
    if (*reio_parametrization == "reio_inter") {
      pth->reio_parametrization = reio_inter;
      recognized                = true;
    }

    class_test(!recognized,
               "could not identify reionization_parametrization value, check that it is one of "
               "'reio_none', 'reio_camb', 'reio_bins_tanh', 'reio_half_tanh', 'reio_many_tanh', "
               "'reio_inter'...");
  }

  /** - reionization parameters if reio_parametrization=reio_camb */
  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {
    auto z_reio   = pfc->get<double>("z_reio");
    auto tau_reio = pfc->get<double>("tau_reio");
    class_test((z_reio && tau_reio),
               "In input file, you can only enter one of z_reio or tau_reio, choose one");
    if (z_reio) {
      pth->z_reio        = *z_reio;
      pth->reio_z_or_tau = reio_z;
    }
    if (tau_reio) {
      pth->tau_reio      = *tau_reio;
      pth->reio_z_or_tau = reio_tau;
    }

    pth->reionization_exponent = pfc->get_or("reionization_exponent", pth->reionization_exponent);
    pth->reionization_width    = pfc->get_or("reionization_width", pth->reionization_width);
    pth->helium_fullreio_redshift = pfc->get_or("helium_fullreio_redshift",
                                                pth->helium_fullreio_redshift);
    pth->helium_fullreio_width = pfc->get_or("helium_fullreio_width", pth->helium_fullreio_width);
  }

  /** - reionization parameters if reio_parametrization=reio_bins_tanh */
  if (pth->reio_parametrization == reio_bins_tanh) {
    int found;
    pth->binned_reio_num = pfc->get_or("binned_reio_num", pth->binned_reio_num);
    readDoubleList(pfc, "binned_reio_z", pth->binned_reio_z, &found);
    class_test(!found || static_cast<int>(pth->binned_reio_z.size()) != pth->binned_reio_num,
               "Number of entries in binned_reio_z does not match expected number, %d.",
               pth->binned_reio_num);
    readDoubleList(pfc, "binned_reio_xe", pth->binned_reio_xe, &found);
    class_test(!found || static_cast<int>(pth->binned_reio_xe.size()) != pth->binned_reio_num,
               "Number of entries in binned_reio_xe does not match expected number, %d.",
               pth->binned_reio_num);
    pth->binned_reio_step_sharpness = pfc->get_or("binned_reio_step_sharpness",
                                                  pth->binned_reio_step_sharpness);
  }

  /** - reionization parameters if reio_parametrization=reio_many_tanh */
  if (pth->reio_parametrization == reio_many_tanh) {
    int found;
    pth->many_tanh_num = pfc->get_or("many_tanh_num", pth->many_tanh_num);
    readDoubleList(pfc, "many_tanh_z", pth->many_tanh_z, &found);
    class_test(!found || static_cast<int>(pth->many_tanh_z.size()) != pth->many_tanh_num,
               "Number of entries in many_tanh_z does not match expected number, %d.",
               pth->many_tanh_num);
    readDoubleList(pfc, "many_tanh_xe", pth->many_tanh_xe, &found);
    class_test(!found || static_cast<int>(pth->many_tanh_xe.size()) != pth->many_tanh_num,
               "Number of entries in many_tanh_xe does not match expected number, %d.",
               pth->many_tanh_num);
    pth->many_tanh_width = pfc->get_or("many_tanh_width", pth->many_tanh_width);
  }

  /** - reionization parameters if reio_parametrization=reio_many_tanh */
  if (pth->reio_parametrization == reio_inter) {
    int found;
    pth->reio_inter_num = pfc->get_or("reio_inter_num", pth->reio_inter_num);
    readDoubleList(pfc, "reio_inter_z", pth->reio_inter_z, &found);
    class_test(!found || static_cast<int>(pth->reio_inter_z.size()) != pth->reio_inter_num,
               "Number of entries in reio_inter_z does not match expected number, %d.",
               pth->reio_inter_num);
    readDoubleList(pfc, "reio_inter_xe", pth->reio_inter_xe, &found);
    class_test(!found || static_cast<int>(pth->reio_inter_xe.size()) != pth->reio_inter_num,
               "Number of entries in reio_inter_xe does not match expected number, %d.",
               pth->reio_inter_num);
  }

  /** - energy injection parameters from CDM annihilation/decay */

  pth->annihilation = pfc->get_or("annihilation", pth->annihilation);

  if (pth->annihilation > 0.) {
    pth->annihilation_variation = pfc->get_or("annihilation_variation",
                                              pth->annihilation_variation);
    pth->annihilation_z         = pfc->get_or("annihilation_z", pth->annihilation_z);
    pth->annihilation_zmax      = pfc->get_or("annihilation_zmax", pth->annihilation_zmax);
    pth->annihilation_zmin      = pfc->get_or("annihilation_zmin", pth->annihilation_zmin);
    pth->annihilation_f_halo    = pfc->get_or("annihilation_f_halo", pth->annihilation_f_halo);
    pth->annihilation_z_halo    = pfc->get_or("annihilation_z_halo", pth->annihilation_z_halo);

    if (auto on_the_spot = pfc->get<std::string>("on the spot")) {
      if ((on_the_spot->find("y") != std::string::npos) ||
          (on_the_spot->find("Y") != std::string::npos)) {
        pth->has_on_the_spot = true;
      }
      else {
        if ((on_the_spot->find("n") != std::string::npos) ||
            (on_the_spot->find("N") != std::string::npos)) {
          pth->has_on_the_spot = false;
        }
        else {
          class_stop("incomprehensible input '%s' for the field 'on the spot'",
                     on_the_spot->c_str());
        }
      }
    }
  }

  pth->decay = pfc->get_or("decay", pth->decay);

  if (auto compute_damping_scale = pfc->get<std::string>("compute damping scale")) {
    if ((compute_damping_scale->find("y") != std::string::npos) ||
        (compute_damping_scale->find("Y") != std::string::npos)) {
      pth->compute_damping_scale = true;
    }
    else {
      if ((compute_damping_scale->find("n") != std::string::npos) ||
          (compute_damping_scale->find("N") != std::string::npos)) {
        pth->compute_damping_scale = false;
      }
      else {
        class_stop("incomprehensible input '%s' for the field 'compute damping scale'",
                   compute_damping_scale->c_str());
      }
    }
  }

  /** (c) define which perturbations and sources should be computed, and down to which scale */

  ppt->has_perturbations = false;
  ppt->has_cls           = false;

  if (auto output = pfc->get<std::string>("output")) {
    if ((output->find("tCl") != std::string::npos) || (output->find("TCl") != std::string::npos) ||
        (output->find("TCL") != std::string::npos)) {
      ppt->has_cl_cmb_temperature = true;
      ppt->has_perturbations      = true;
      ppt->has_cls                = true;
    }

    if ((output->find("pCl") != std::string::npos) || (output->find("PCl") != std::string::npos) ||
        (output->find("PCL") != std::string::npos)) {
      ppt->has_cl_cmb_polarization = true;
      ppt->has_perturbations       = true;
      ppt->has_cls                 = true;
    }

    if ((output->find("lCl") != std::string::npos) || (output->find("LCl") != std::string::npos) ||
        (output->find("LCL") != std::string::npos)) {
      ppt->has_cl_cmb_lensing_potential = true;
      ppt->has_perturbations            = true;
      ppt->has_cls                      = true;
    }

    if ((output->find("nCl") != std::string::npos) || (output->find("NCl") != std::string::npos) ||
        (output->find("NCL") != std::string::npos) || (output->find("dCl") != std::string::npos) ||
        (output->find("DCl") != std::string::npos) || (output->find("DCL") != std::string::npos)) {
      ppt->has_cl_number_count = true;
      ppt->has_perturbations   = true;
      ppt->has_cls             = true;
    }

    if ((output->find("sCl") != std::string::npos) || (output->find("SCl") != std::string::npos) ||
        (output->find("SCL") != std::string::npos)) {
      ppt->has_cl_lensing_potential = true;
      ppt->has_perturbations        = true;
      ppt->has_cls                  = true;
    }

    if ((output->find("mPk") != std::string::npos) || (output->find("MPk") != std::string::npos) ||
        (output->find("MPK") != std::string::npos)) {
      ppt->has_pk_matter     = true;
      ppt->has_perturbations = true;
    }

    if ((output->find("mTk") != std::string::npos) || (output->find("MTk") != std::string::npos) ||
        (output->find("MTK") != std::string::npos) || (output->find("dTk") != std::string::npos) ||
        (output->find("DTk") != std::string::npos) || (output->find("DTK") != std::string::npos)) {
      ppt->has_density_transfers = true;
      ppt->has_perturbations     = true;
    }

    if ((output->find("vTk") != std::string::npos) || (output->find("VTk") != std::string::npos) ||
        (output->find("VTK") != std::string::npos)) {
      ppt->has_velocity_transfers = true;
      ppt->has_perturbations      = true;
    }
  }

  if (ppt->has_density_transfers) {
    if (auto extra = pfc->get<std::string>("extra metric transfer functions");
        extra &&
        ((extra->find("y") != std::string::npos) || (extra->find("y") != std::string::npos))) {
      ppt->has_metricpotential_transfers = true;
    }
  }

  if (ppt->has_cl_cmb_temperature) {
    if (auto temperature_contributions = pfc->get<std::string>("temperature contributions")) {
      const std::string& tc = *temperature_contributions;
      ppt->switch_sw        = 0;
      ppt->switch_eisw      = 0;
      ppt->switch_lisw      = 0;
      ppt->switch_dop       = 0;
      ppt->switch_pol       = 0;

      if ((tc.find("tsw") != std::string::npos) || (tc.find("TSW") != std::string::npos))
        ppt->switch_sw = 1;
      if ((tc.find("eisw") != std::string::npos) || (tc.find("EISW") != std::string::npos))
        ppt->switch_eisw = 1;
      if ((tc.find("lisw") != std::string::npos) || (tc.find("LISW") != std::string::npos))
        ppt->switch_lisw = 1;
      if ((tc.find("dop") != std::string::npos) || (tc.find("Dop") != std::string::npos))
        ppt->switch_dop = 1;
      if ((tc.find("pol") != std::string::npos) || (tc.find("Pol") != std::string::npos))
        ppt->switch_pol = 1;

      class_test((ppt->switch_sw == 0) && (ppt->switch_eisw == 0) && (ppt->switch_lisw == 0) &&
                     (ppt->switch_dop == 0) && (ppt->switch_pol == 0),
                 "In the field 'output', you selected CMB temperature, but in the field "
                 "'temperature contributions', you removed all contributions");

      ppt->eisw_lisw_split_z = pfc->get_or("early/late isw redshift", ppt->eisw_lisw_split_z);
    }
  }

  if (ppt->has_cl_number_count) {
    if (auto number_count_contributions = pfc->get<std::string>("number count contributions")) {
      const std::string& ncc = *number_count_contributions;
      if (ncc.find("density") != std::string::npos)
        ppt->has_nc_density = true;
      if (ncc.find("rsd") != std::string::npos)
        ppt->has_nc_rsd = true;
      if (ncc.find("lensing") != std::string::npos)
        ppt->has_nc_lens = true;
      if (ncc.find("gr") != std::string::npos)
        ppt->has_nc_gr = true;

      class_test((!ppt->has_nc_density) && (!ppt->has_nc_rsd) && (!ppt->has_nc_lens) &&
                     (!ppt->has_nc_gr),
                 "In the field 'output', you selected number count Cl's, but in the field 'number "
                 "count contributions', you removed all contributions");
    }

    else {
      /* default: only the density contribution */
      ppt->has_nc_density = true;
    }
  }

  if (ppt->has_perturbations) {
    /* perturbed recombination */
    if (auto perturbed_recombination = pfc->get<std::string>("perturbed recombination");
        perturbed_recombination && ((perturbed_recombination->find("y") != std::string::npos) ||
                                    (perturbed_recombination->find("Y") != std::string::npos))) {
      ppt->has_perturbed_recombination = true;
    }

    /* modes */
    if (auto modes = pfc->get<std::string>("modes")) {
      /* if no modes are specified, the default is has_scalars=true;
         but if they are specified we should reset has_scalars to false before reading */
      ppt->has_scalars = false;

      if ((modes->find("s") != std::string::npos) || (modes->find("S") != std::string::npos))
        ppt->has_scalars = true;

      if ((modes->find("v") != std::string::npos) || (modes->find("V") != std::string::npos))
        ppt->has_vectors = true;

      if ((modes->find("t") != std::string::npos) || (modes->find("T") != std::string::npos))
        ppt->has_tensors = true;

      class_test(!ppt->has_scalars && !ppt->has_vectors && !ppt->has_tensors,
                 "You wrote: modes='%s'. Could not identify any of the modes ('s', 'v', 't') in "
                 "such input",
                 modes->c_str());
    }

    if (ppt->has_scalars) {
      if (auto ic = pfc->get<std::string>("ic")) {
        /* if no initial conditions are specified, the default is has_ad=true;
           but if they are specified we should reset has_ad to false before reading */
        ppt->has_ad = false;

        if ((ic->find("ad") != std::string::npos) || (ic->find("AD") != std::string::npos))
          ppt->has_ad = true;

        if ((ic->find("bi") != std::string::npos) || (ic->find("BI") != std::string::npos))
          ppt->has_bi = true;

        if ((ic->find("cdi") != std::string::npos) || (ic->find("CDI") != std::string::npos))
          ppt->has_cdi = true;

        if ((ic->find("nid") != std::string::npos) || (ic->find("NID") != std::string::npos))
          ppt->has_nid = true;

        if ((ic->find("niv") != std::string::npos) || (ic->find("NIV") != std::string::npos))
          ppt->has_niv = true;

        class_test(!ppt->has_ad && !ppt->has_bi && !ppt->has_cdi && !ppt->has_nid && !ppt->has_niv,
                   "You wrote: ic='%s'. Could not identify any of the initial conditions ('ad', "
                   "'bi', 'cdi', 'nid', 'niv') in such input",
                   ic->c_str());
      }
    }

    else {
      class_test(ppt->has_cl_cmb_lensing_potential,
                 "Inconsistency: you want C_l's for cmb lensing potential, but no scalar modes\n");

      class_test(ppt->has_pk_matter,
                 "Inconsistency: you want P(k) of matter, but no scalar modes\n");
    }

    if (ppt->has_vectors) {
      class_test((!ppt->has_cl_cmb_temperature) && (!ppt->has_cl_cmb_polarization),
                 "inconsistent input: you asked for vectors, so you should have at least one "
                 "non-zero "
                 "tensor source type (temperature or polarization). Please adjust your input.");
    }

    if (ppt->has_tensors) {
      class_test((!ppt->has_cl_cmb_temperature) && (!ppt->has_cl_cmb_polarization),
                 "inconsistent input: you asked for tensors, so you should have at least one "
                 "non-zero "
                 "tensor source type (temperature or polarization). Please adjust your input.");
    }
  }

  /** (d) define the primordial spectrum */

  if (auto P_k_ini_type = pfc->get<std::string>("P_k_ini type")) {
    bool recognized = false;
    if (*P_k_ini_type == "analytic_Pk") {
      ppm->primordial_spec_type = analytic_Pk;
      recognized                = true;
    }
    if (*P_k_ini_type == "two_scales") {
      ppm->primordial_spec_type = two_scales;
      recognized                = true;
    }
    if (*P_k_ini_type == "inflation_V") {
      ppm->primordial_spec_type = inflation_V;
      recognized                = true;
    }
    if (*P_k_ini_type == "inflation_H") {
      ppm->primordial_spec_type = inflation_H;
      recognized                = true;
    }
    if (*P_k_ini_type == "inflation_V_end") {
      ppm->primordial_spec_type = inflation_V_end;
      recognized                = true;
    }
    if (*P_k_ini_type == "external_Pk") {
      ppm->primordial_spec_type = external_Pk;
      recognized                = true;
    }
    class_test(!recognized,
               "could not identify primordial spectrum type, check that it is one of "
               "'analytic_pk', 'two_scales', 'inflation_V', 'inflation_H', 'external_Pk'...");
  }

  ppm->k_pivot = pfc->get_or("k_pivot", ppm->k_pivot);

  if (ppm->primordial_spec_type == two_scales) {
    k1 = pfc->get_or("k1", k1);
    k2 = pfc->get_or("k2", k2);
    class_test(k1 <= 0., "enter strictly positive scale k1");
    class_test(k2 <= 0., "enter strictly positive scale k2");

    if (ppt->has_scalars) {
      prr1 = pfc->get_or("P_{RR}^1", prr1);
      prr2 = pfc->get_or("P_{RR}^2", prr2);
      class_test(prr1 <= 0., "enter strictly positive scale P_{RR}^1");
      class_test(prr2 <= 0., "enter strictly positive scale P_{RR}^2");

      ppm->n_s = log(prr2 / prr1) / log(k2 / k1) + 1.;
      ppm->A_s = prr1 * exp((ppm->n_s - 1.) * log(ppm->k_pivot / k1));

      if ((ppt->has_bi) || (ppt->has_cdi) || (ppt->has_nid) || (ppt->has_niv)) {
        pii1 = pfc->get_or("P_{II}^1", pii1);
        pii2 = pfc->get_or("P_{II}^2", pii2);
        pri1 = pfc->get_or("P_{RI}^1", pri1);
        pri2 = pfc->get_or("|P_{RI}^2|", pri2);

        class_test(pii1 <= 0.,
                   "since you request iso modes, you should have P_{ii}^1 strictly positive");
        class_test(pii2 < 0.,
                   "since you request iso modes, you should have P_{ii}^2 positive or eventually "
                   "null");
        class_test(pri2 < 0.,
                   "by definition, you should have |P_{ri}^2| positive or eventually null");

        auto special_iso = pfc->get<std::string>("special iso");

        /* axion case, only one iso parameter: piir1  */
        if (special_iso && (special_iso->find("axion") != std::string::npos)) {
          n_iso = 1.;
          n_cor = 0.;
          c_cor = 0.;
        }
        /* curvaton case, only one iso parameter: piir1  */
        else if (special_iso && (special_iso->find("anticurvaton") != std::string::npos)) {
          n_iso = ppm->n_s;
          n_cor = 0.;
          c_cor = 1.;
        }
        /* inverted-correlation-curvaton case, only one iso parameter: piir1  */
        else if (special_iso && (special_iso->find("curvaton") != std::string::npos)) {
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
            class_test((pii1 == 0.) || (pii2 == 0.) || (pii1 * pii2 < 0.), "should NEVER happen");
            n_iso = log(pii2 / pii1) / log(k2 / k1) + 1.;
          }
          class_test(pri1 == 0, "the general isocurvature case requires a non-zero P_{RI}^1");
          if (pri2 == 0.) {
            n_cor = 0.;
          }
          else {
            class_test((pri1 == 0.) || (pri2 <= 0.) || (pii1 * pii2 < 0), "should NEVER happen");
            n_cor = log(pri2 / fabs(pri1)) / log(k2 / k1) - 0.5 * (ppm->n_s + n_iso - 2.);
          }
          class_test((pii1 * prr1 <= 0.), "should NEVER happen");
          class_test(fabs(pri1) / sqrt(pii1 * prr1) > 1,
                     "too large ad-iso cross-correlation in k1");
          class_test(fabs(pri1) / sqrt(pii1 * prr1) * exp(n_cor * log(k2 / k1)) > 1,
                     "too large ad-iso cross-correlation in k2");
          c_cor = -pri1 / sqrt(pii1 * prr1) * exp(n_cor * log(ppm->k_pivot / k1));
        }
        /* formula for f_iso valid in all cases */
        class_test((pii1 == 0.) || (prr1 == 0.) || (pii1 * prr1 < 0.), "should NEVER happen");
        f_iso = sqrt(pii1 / prr1) * exp(0.5 * (n_iso - ppm->n_s) * log(ppm->k_pivot / k1));
      }

      if (ppt->has_bi) {
        ppm->f_bi    = f_iso;
        ppm->n_bi    = n_iso;
        ppm->c_ad_bi = c_cor;
        ppm->n_ad_bi = n_cor;
      }

      if (ppt->has_cdi) {
        ppm->f_cdi    = f_iso;
        ppm->n_cdi    = n_iso;
        ppm->c_ad_cdi = c_cor;
        ppm->n_ad_cdi = n_cor;
      }

      if (ppt->has_nid) {
        ppm->f_nid    = f_iso;
        ppm->n_nid    = n_iso;
        ppm->c_ad_nid = c_cor;
        ppm->n_ad_nid = n_cor;
      }

      if (ppt->has_niv) {
        ppm->f_niv    = f_iso;
        ppm->n_niv    = n_iso;
        ppm->c_ad_niv = c_cor;
        ppm->n_ad_niv = n_cor;
      }
    }

    ppm->primordial_spec_type = analytic_Pk;
  }

  else if (ppm->primordial_spec_type == analytic_Pk) {
    if (ppt->has_scalars) {
      auto A_s    = pfc->get<double>("A_s");
      auto ln_A_s = pfc->get<double>("ln10^{10}A_s");
      auto sigma8 = pfc->get<double>("sigma8");
      auto S8     = pfc->get<double>("S8");
      class_test(n_present(A_s, ln_A_s, sigma8, S8) > 1,
                 "In input file, you can only enter one of A_s, ln10^{10}A_s, sigma8 and S8, "
                 "choose one");
      if (A_s)
        ppm->A_s = *A_s;
      else if (ln_A_s)
        ppm->A_s = exp(*ln_A_s) * 1.e-10;
      else if (sigma8) {
        ppm->sigma8 = *sigma8;
        class_test(*sigma8 < 0., "sigma8 should be non-negative");
      }
      else if (S8) {
        // CDM read from the built species (ReadDerived runs after ConstructSpecies).
        const double Omega0_cdm_for_S8 = all_species_.count("CDM")
                                             ? all_species_.at("CDM")->GetOmega0()
                                             : 0.0;
        ppm->sigma8                    = *S8 / pow((pba->Omega0_b + Omega0_cdm_for_S8) / 0.3, 0.5);
        class_test(*S8 < 0., "S8 should be non-negative");
      }

      if (ppt->has_ad) {
        ppm->n_s     = pfc->get_or("n_s", ppm->n_s);
        ppm->alpha_s = pfc->get_or("alpha_s", ppm->alpha_s);
      }

      if (ppt->has_bi) {
        ppm->f_bi     = pfc->get_or("f_bi", ppm->f_bi);
        ppm->n_bi     = pfc->get_or("n_bi", ppm->n_bi);
        ppm->alpha_bi = pfc->get_or("alpha_bi", ppm->alpha_bi);
      }

      if (ppt->has_cdi) {
        ppm->f_cdi     = pfc->get_or("f_cdi", ppm->f_cdi);
        ppm->n_cdi     = pfc->get_or("n_cdi", ppm->n_cdi);
        ppm->alpha_cdi = pfc->get_or("alpha_cdi", ppm->alpha_cdi);
      }

      if (ppt->has_nid) {
        ppm->f_nid     = pfc->get_or("f_nid", ppm->f_nid);
        ppm->n_nid     = pfc->get_or("n_nid", ppm->n_nid);
        ppm->alpha_nid = pfc->get_or("alpha_nid", ppm->alpha_nid);
      }

      if (ppt->has_niv) {
        ppm->f_niv     = pfc->get_or("f_niv", ppm->f_niv);
        ppm->n_niv     = pfc->get_or("n_niv", ppm->n_niv);
        ppm->alpha_niv = pfc->get_or("alpha_niv", ppm->alpha_niv);
      }

      if ((ppt->has_ad) && (ppt->has_bi)) {
        if (auto v = read_one_of<double>(*pfc, "c_ad_bi", "c_bi_ad"))
          ppm->c_ad_bi = *v;
        if (auto v = read_one_of<double>(*pfc, "n_ad_bi", "n_bi_ad"))
          ppm->n_ad_bi = *v;
        if (auto v = read_one_of<double>(*pfc, "alpha_ad_bi", "alpha_bi_ad"))
          ppm->alpha_ad_bi = *v;
      }

      if ((ppt->has_ad) && (ppt->has_cdi)) {
        if (auto v = read_one_of<double>(*pfc, "c_ad_cdi", "c_cdi_ad"))
          ppm->c_ad_cdi = *v;
        if (auto v = read_one_of<double>(*pfc, "n_ad_cdi", "n_cdi_ad"))
          ppm->n_ad_cdi = *v;
        if (auto v = read_one_of<double>(*pfc, "alpha_ad_cdi", "alpha_cdi_ad"))
          ppm->alpha_ad_cdi = *v;
      }

      if ((ppt->has_ad) && (ppt->has_nid)) {
        if (auto v = read_one_of<double>(*pfc, "c_ad_nid", "c_nid_ad"))
          ppm->c_ad_nid = *v;
        if (auto v = read_one_of<double>(*pfc, "n_ad_nid", "n_nid_ad"))
          ppm->n_ad_nid = *v;
        if (auto v = read_one_of<double>(*pfc, "alpha_ad_nid", "alpha_nid_ad"))
          ppm->alpha_ad_nid = *v;
      }

      if ((ppt->has_ad) && (ppt->has_niv)) {
        if (auto v = read_one_of<double>(*pfc, "c_ad_niv", "c_niv_ad"))
          ppm->c_ad_niv = *v;
        if (auto v = read_one_of<double>(*pfc, "n_ad_niv", "n_niv_ad"))
          ppm->n_ad_niv = *v;
        if (auto v = read_one_of<double>(*pfc, "alpha_ad_niv", "alpha_niv_ad"))
          ppm->alpha_ad_niv = *v;
      }

      if ((ppt->has_bi) && (ppt->has_cdi)) {
        if (auto v = read_one_of<double>(*pfc, "c_bi_cdi", "c_cdi_bi"))
          ppm->c_bi_cdi = *v;
        if (auto v = read_one_of<double>(*pfc, "n_bi_cdi", "n_cdi_bi"))
          ppm->n_bi_cdi = *v;
        if (auto v = read_one_of<double>(*pfc, "alpha_bi_cdi", "alpha_cdi_bi"))
          ppm->alpha_bi_cdi = *v;
      }

      if ((ppt->has_bi) && (ppt->has_nid)) {
        if (auto v = read_one_of<double>(*pfc, "c_bi_nid", "c_nid_bi"))
          ppm->c_bi_nid = *v;
        if (auto v = read_one_of<double>(*pfc, "n_bi_nid", "n_nid_bi"))
          ppm->n_bi_nid = *v;
        if (auto v = read_one_of<double>(*pfc, "alpha_bi_nid", "alpha_nid_bi"))
          ppm->alpha_bi_nid = *v;
      }

      if ((ppt->has_bi) && (ppt->has_niv)) {
        if (auto v = read_one_of<double>(*pfc, "c_bi_niv", "c_niv_bi"))
          ppm->c_bi_niv = *v;
        if (auto v = read_one_of<double>(*pfc, "n_bi_niv", "n_niv_bi"))
          ppm->n_bi_niv = *v;
        if (auto v = read_one_of<double>(*pfc, "alpha_bi_niv", "alpha_niv_bi"))
          ppm->alpha_bi_niv = *v;
      }

      if ((ppt->has_cdi) && (ppt->has_nid)) {
        if (auto v = read_one_of<double>(*pfc, "c_cdi_nid", "c_nid_cdi"))
          ppm->c_cdi_nid = *v;
        if (auto v = read_one_of<double>(*pfc, "n_cdi_nid", "n_nid_cdi"))
          ppm->n_cdi_nid = *v;
        if (auto v = read_one_of<double>(*pfc, "alpha_cdi_nid", "alpha_nid_cdi"))
          ppm->alpha_cdi_nid = *v;
      }

      if ((ppt->has_cdi) && (ppt->has_niv)) {
        if (auto v = read_one_of<double>(*pfc, "c_cdi_niv", "c_niv_cdi"))
          ppm->c_cdi_niv = *v;
        if (auto v = read_one_of<double>(*pfc, "n_cdi_niv", "n_niv_cdi"))
          ppm->n_cdi_niv = *v;
        if (auto v = read_one_of<double>(*pfc, "alpha_cdi_niv", "alpha_niv_cdi"))
          ppm->alpha_cdi_niv = *v;
      }

      if ((ppt->has_nid) && (ppt->has_niv)) {
        if (auto v = read_one_of<double>(*pfc, "c_nid_niv", "c_niv_nid"))
          ppm->c_nid_niv = *v;
        if (auto v = read_one_of<double>(*pfc, "n_nid_niv", "n_niv_nid"))
          ppm->n_nid_niv = *v;
        if (auto v = read_one_of<double>(*pfc, "alpha_nid_niv", "alpha_niv_nid"))
          ppm->alpha_nid_niv = *v;
      }
    }

    if (ppt->has_tensors) {
      ppm->r = pfc->get_or("r", ppm->r);

      if (!ppt->has_scalars) {
        ppm->A_s = pfc->get_or("A_s", ppm->A_s);
      }

      if (ppm->r <= 0) {
        ppt->has_tensors = false;
      }
      else {
        auto n_t = pfc->get<std::string>("n_t");

        if (n_t &&
            !((n_t->find("SCC") != std::string::npos) || (n_t->find("scc") != std::string::npos))) {
          ppm->n_t = pfc->get_or("n_t", ppm->n_t);
        }
        else {
          /* enforce single slow-roll self-consistency condition (order 2 in slow-roll) */
          ppm->n_t = -ppm->r / 8. * (2. - ppm->r / 8. - ppm->n_s);
        }

        auto alpha_t = pfc->get<std::string>("alpha_t");

        if (alpha_t && !((alpha_t->find("SCC") != std::string::npos) ||
                         (alpha_t->find("scc") != std::string::npos))) {
          ppm->alpha_t = pfc->get_or("alpha_t", ppm->alpha_t);
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
      pfc->get<std::string>("potential");

      /* only polynomial coded so far: no need to interpret the value **/

      if (pfc->get<double>("PSR_0")) {
        PSR0 = 0.;
        PSR1 = 0.;
        PSR2 = 0.;
        PSR3 = 0.;
        PSR4 = 0.;

        PSR0 = pfc->get_or("PSR_0", PSR0);
        PSR1 = pfc->get_or("PSR_1", PSR1);
        PSR2 = pfc->get_or("PSR_2", PSR2);
        PSR3 = pfc->get_or("PSR_3", PSR3);
        PSR4 = pfc->get_or("PSR_4", PSR4);

        class_test(PSR0 <= 0., "inconsistent parametrization of polynomial inflation potential");
        class_test(PSR1 <= 0., "inconsistent parametrization of polynomial inflation potential");

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
        if (pfc->get<double>("R_0")) {
          R0 = 0.;
          R1 = 0.;
          R2 = 0.;
          R3 = 0.;
          R4 = 0.;

          R0 = pfc->get_or("R_0", R0);
          R1 = pfc->get_or("R_1", R1);
          R2 = pfc->get_or("R_2", R2);
          R3 = pfc->get_or("R_3", R3);
          R4 = pfc->get_or("R_4", R4);

          class_test(R0 <= 0., "inconsistent parametrization of polynomial inflation potential");
          class_test(R1 <= 0., "inconsistent parametrization of polynomial inflation potential");

          ppm->V0 = R0 * R1 * 3. / 128. / _PI_;
          ppm->V1 = -sqrt(R1) * ppm->V0;
          ppm->V2 = R2 * ppm->V0;
          ppm->V3 = R3 * ppm->V0 * ppm->V0 / ppm->V1;
          ppm->V4 = R4 * ppm->V0 / R1;
        }

        else {
          ppm->V0 = pfc->get_or("V_0", ppm->V0);
          ppm->V1 = pfc->get_or("V_1", ppm->V1);
          ppm->V2 = pfc->get_or("V_2", ppm->V2);
          ppm->V3 = pfc->get_or("V_3", ppm->V3);
          ppm->V4 = pfc->get_or("V_4", ppm->V4);
        }
      }
    }

    else {
      if (pfc->get<double>("HSR_0")) {
        HSR0 = 0.;
        HSR1 = 0.;
        HSR2 = 0.;
        HSR3 = 0.;
        HSR4 = 0.;

        HSR0 = pfc->get_or("HSR_0", HSR0);
        HSR1 = pfc->get_or("HSR_1", HSR1);
        HSR2 = pfc->get_or("HSR_2", HSR2);
        HSR3 = pfc->get_or("HSR_3", HSR3);
        HSR4 = pfc->get_or("HSR_4", HSR4);

        ppm->H0 = sqrt(HSR0 * HSR1 * _PI_);
        ppm->H1 = -sqrt(4. * _PI_ * HSR1) * ppm->H0;
        ppm->H2 = 4. * _PI_ * HSR2 * ppm->H0;
        ppm->H3 = 4. * _PI_ * HSR3 * ppm->H0 * ppm->H0 / ppm->H1;
        ppm->H4 = 4. * _PI_ * HSR4 * ppm->H0 * ppm->H0 * ppm->H0 / ppm->H1 / ppm->H1;
      }
      else {
        ppm->H0 = pfc->get_or("H_0", ppm->H0);
        ppm->H1 = pfc->get_or("H_1", ppm->H1);
        ppm->H2 = pfc->get_or("H_2", ppm->H2);
        ppm->H3 = pfc->get_or("H_3", ppm->H3);
        ppm->H4 = pfc->get_or("H_4", ppm->H4);
      }

      class_test(ppm->H0 <= 0., "inconsistent parametrization of polynomial inflation potential");
    }
  }

  else if (ppm->primordial_spec_type == inflation_V_end) {
    if (auto full_potential = pfc->get<std::string>("full_potential")) {
      if (*full_potential == "polynomial") {
        ppm->potential = polynomial;
      }
      else if (*full_potential == "higgs_inflation") {
        ppm->potential = higgs_inflation;
      }
      else {
        class_stop(
            "did not recognize input parameter 'potential': should be one of 'polynomial' "
            "or 'higgs_inflation'");
      }
    }

    ppm->phi_end = pfc->get_or("phi_end", ppm->phi_end);
    ppm->V0      = pfc->get_or("Vparam0", ppm->V0);
    ppm->V1      = pfc->get_or("Vparam1", ppm->V1);
    ppm->V2      = pfc->get_or("Vparam2", ppm->V2);
    ppm->V3      = pfc->get_or("Vparam3", ppm->V3);
    ppm->V4      = pfc->get_or("Vparam4", ppm->V4);

    auto ln_aH_ratio_str = pfc->get<std::string>("ln_aH_ratio");
    auto N_star_str      = pfc->get<std::string>("N_star");

    class_test(ln_aH_ratio_str && N_star_str,
               "In input file, you can only enter one of ln_aH_ratio or N_star, the two are not "
               "compatible");

    if (ln_aH_ratio_str) {
      if ((ln_aH_ratio_str->find("auto") != std::string::npos) ||
          (ln_aH_ratio_str->find("AUTO") != std::string::npos)) {
        ppm->phi_pivot_method = ln_aH_ratio_auto;
      }
      else {
        ppm->phi_pivot_method = ln_aH_ratio;
        ppm->phi_pivot_target = pfc->get_or("ln_aH_ratio", ppm->phi_pivot_target);
      }
    }

    if (N_star_str) {
      ppm->phi_pivot_method = N_star;
      ppm->phi_pivot_target = pfc->get_or("N_star", ppm->phi_pivot_target);
    }

    if (auto inflation_behavior = pfc->get<std::string>("inflation_behavior")) {
      if (inflation_behavior->find("numerical") != std::string::npos) {
        ppm->behavior = numerical;
      }
      else if (inflation_behavior->find("analytical") != std::string::npos) {
        ppm->behavior = analytical;
      }
      else {
        class_stop("Your entry for 'inflation behavior' could not be understood");
      }
    }
  }

  else if (ppm->primordial_spec_type == external_Pk) {
    auto command = pfc->get<std::string>("command");
    class_test(!command || command->empty(), "You omitted to write a command for the external Pk");

    ppm->command  = *command;
    ppm->custom1  = pfc->get_or("custom1", ppm->custom1);
    ppm->custom2  = pfc->get_or("custom2", ppm->custom2);
    ppm->custom3  = pfc->get_or("custom3", ppm->custom3);
    ppm->custom4  = pfc->get_or("custom4", ppm->custom4);
    ppm->custom5  = pfc->get_or("custom5", ppm->custom5);
    ppm->custom6  = pfc->get_or("custom6", ppm->custom6);
    ppm->custom7  = pfc->get_or("custom7", ppm->custom7);
    ppm->custom8  = pfc->get_or("custom8", ppm->custom8);
    ppm->custom9  = pfc->get_or("custom9", ppm->custom9);
    ppm->custom10 = pfc->get_or("custom10", ppm->custom10);
  }

  /* Tests moved from primordial module: */
  if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_H) ||
      (ppm->primordial_spec_type == inflation_V_end)) {
    class_test(!ppt->has_scalars,
               "inflationary module cannot work if you do not ask for scalar modes");

    class_test(ppt->has_vectors, "inflationary module cannot work if you ask for vector modes");

    class_test(!ppt->has_tensors,
               "inflationary module cannot work if you do not ask for tensor modes");

    class_test(ppt->has_bi || ppt->has_cdi || ppt->has_nid || ppt->has_niv,
               "inflationary module cannot work if you ask for isocurvature modes");
  }

  /** (e) parameters for final spectra */

  if (ppt->has_cls) {
    if (ppt->has_scalars) {
      if ((ppt->has_cl_cmb_temperature) || (ppt->has_cl_cmb_polarization) ||
          (ppt->has_cl_cmb_lensing_potential))
        ppt->l_scalar_max = pfc->get_or("l_max_scalars", ppt->l_scalar_max);

      if ((ppt->has_cl_lensing_potential) || (ppt->has_cl_number_count))
        ppt->l_lss_max = pfc->get_or("l_max_lss", ppt->l_lss_max);
    }

    if (ppt->has_vectors) {
      ppt->l_vector_max = pfc->get_or("l_max_vectors", ppt->l_vector_max);
    }

    if (ppt->has_tensors) {
      ppt->l_tensor_max = pfc->get_or("l_max_tensors", ppt->l_tensor_max);
    }
  }

  if (auto lensing = pfc->get<std::string>("lensing");
      lensing &&
      ((lensing->find("y") != std::string::npos) || (lensing->find("Y") != std::string::npos))) {
    if ((ppt->has_scalars) && ((ppt->has_cl_cmb_temperature) || (ppt->has_cl_cmb_polarization)) &&
        (ppt->has_cl_cmb_lensing_potential)) {
      ple->has_lensed_cls = true;
    }
    else {
      class_stop(
          "you asked for lensed CMB Cls, but this requires a minimal number of options: "
          "'modes' should include 's', 'output' should include 'tCl' and/or 'pCL', and "
          "also, importantly, 'lCl', the CMB lensing potential spectrum. You forgot one of "
          "those in your input.");
    }
  }

  if ((ppt->has_scalars) && (ppt->has_cl_cmb_lensing_potential)) {
    ptr->lcmb_rescale = pfc->get_or("lcmb_rescale", ptr->lcmb_rescale);
    ptr->lcmb_tilt    = pfc->get_or("lcmb_tilt", ptr->lcmb_tilt);
    ptr->lcmb_pivot   = pfc->get_or("lcmb_pivot", ptr->lcmb_pivot);
  }

  if (auto full_limber = pfc->get<std::string>("want_lcmb_full_limber")) {
    class_test(full_limber->empty() || (((*full_limber)[0] != 'y') && ((*full_limber)[0] != 'Y') &&
                                        ((*full_limber)[0] != 'n') && ((*full_limber)[0] != 'N')),
               "In input file, want_lcmb_full_limber must begin with 'y' or 'n', not '%s'",
               full_limber->c_str());
    ppt->want_lcmb_full_limber = ((*full_limber)[0] == 'y') || ((*full_limber)[0] == 'Y');
  }

  if ((ppt->has_pk_matter) || (ppt->has_density_transfers) || (ppt->has_velocity_transfers)) {
    auto P_k_max_h_Mpc = pfc->get<double>("P_k_max_h/Mpc");
    auto P_k_max_1_Mpc = pfc->get<double>("P_k_max_1/Mpc");
    class_test(P_k_max_h_Mpc && P_k_max_1_Mpc,
               "In input file, you cannot enter both P_k_max_h/Mpc and P_k_max_1/Mpc, choose one");
    if (P_k_max_h_Mpc) {
      ppt->k_max_for_pk = *P_k_max_h_Mpc * pba->h;
    }
    if (P_k_max_1_Mpc) {
      ppt->k_max_for_pk = *P_k_max_1_Mpc;
    }

    std::vector<double> zPkValues;
    int found;
    readDoubleList(pfc, "z_pk", zPkValues, &found);

    if (found) {
      n_list = static_cast<int>(zPkValues.size());
      class_test(n_list > _Z_PK_NUM_MAX_,
                 "you want to write some output for %d different values of z, hence you should "
                 "increase _Z_PK_NUM_MAX_ in include/output.h to at least this number",
                 n_list);
      pop->z_pk_num = n_list;
      for (i = 0; i < n_list; i++) {
        pop->z_pk[i] = zPkValues[i];
      }
    }
  }

  /** Do we want density and velocity transfer functions in Nbody gauge? */
  if ((ppt->has_density_transfers) || (ppt->has_velocity_transfers)) {
    if (auto nbody = pfc->get<std::string>("Nbody gauge transfer functions");
        nbody &&
        ((nbody->find("y") != std::string::npos) || (nbody->find("y") != std::string::npos))) {
      ppt->has_Nbody_gauge_transfers = true;
    }
  }

  /* deal with selection functions */
  if ((ppt->has_cl_number_count) || (ppt->has_cl_lensing_potential)) {
    int found;
    if (auto selection = pfc->get<std::string>("selection")) {
      if (selection->find("gaussian") != std::string::npos) {
        ppt->selection = gaussian;
      }
      else if (selection->find("tophat") != std::string::npos) {
        ppt->selection = tophat;
      }
      else if (selection->find("dirac") != std::string::npos) {
        ppt->selection = dirac;
      }
      else {
        class_stop("In selection function input: type '%s' is unclear", selection->c_str());
      }
    }

    std::vector<double> selectionValues;
    readDoubleList(pfc, "selection_mean", selectionValues, &found);

    if ((found) && !selectionValues.empty()) {
      n_list = static_cast<int>(selectionValues.size());

      class_test(n_list > _SELECTION_NUM_MAX_,
                 "you want to compute density Cl's for %d different bins, hence you should "
                 "increase _SELECTION_NUM_MAX_ in include/perturbations.h to at least this number",
                 n_list);

      ppt->selection_num = n_list;
      for (i = 0; i < n_list; i++) {
        class_test((selectionValues[i] < 0.) || (selectionValues[i] > 1000.),
                   "input of selection functions: you asked for a mean redshift equal to %e, "
                   "sounds odd",
                   selectionValues[i]);
        ppt->selection_mean[i] = selectionValues[i];
      }
      /* first set all widths to default; correct eventually later */
      for (i = 1; i < n_list; i++) {
        class_test(ppt->selection_mean[i] <= ppt->selection_mean[i - 1],
                   "input of selection functions: the list of mean redshifts must be passed in "
                   "growing order; you entered %e before %e",
                   ppt->selection_mean[i - 1],
                   ppt->selection_mean[i]);
        ppt->selection_width[i]              = ppt->selection_width[0];
        ptr->selection_bias[i]               = ptr->selection_bias[0];
        ptr->selection_magnification_bias[i] = ptr->selection_magnification_bias[0];
      }

      selectionValues.clear();
      readDoubleList(pfc, "selection_width", selectionValues, &found);

      if ((found) && !selectionValues.empty()) {
        n_list = static_cast<int>(selectionValues.size());

        if (n_list == 1) {
          for (i = 0; i < ppt->selection_num; i++) {
            ppt->selection_width[i] = selectionValues[0];
          }
        }
        else if (n_list == ppt->selection_num) {
          for (i = 0; i < n_list; i++) {
            ppt->selection_width[i] = selectionValues[i];
          }
        }
        else {
          class_stop(
              "In input for selection function, you asked for %d bin centers and %d bin "
              "widths; number of bins unclear; you should pass either one bin width (common "
              "to all bins) or %d bin widths",
              ppt->selection_num,
              n_list,
              ppt->selection_num);
        }
      }

      selectionValues.clear();
      readDoubleList(pfc, "selection_bias", selectionValues, &found);

      if ((found) && !selectionValues.empty()) {
        n_list = static_cast<int>(selectionValues.size());

        if (n_list == 1) {
          for (i = 0; i < ppt->selection_num; i++) {
            ptr->selection_bias[i] = selectionValues[0];
          }
        }
        else if (n_list == ppt->selection_num) {
          for (i = 0; i < n_list; i++) {
            ptr->selection_bias[i] = selectionValues[i];
          }
        }
        else {
          class_stop(
              "In input for selection function, you asked for %d bin centers and %d bin "
              "biases; number of bins unclear; you should pass either one bin bias (common "
              "to all bins) or %d bin biases",
              ppt->selection_num,
              n_list,
              ppt->selection_num);
        }
      }

      selectionValues.clear();
      readDoubleList(pfc, "selection_magnification_bias", selectionValues, &found);

      if ((found) && !selectionValues.empty()) {
        n_list = static_cast<int>(selectionValues.size());

        if (n_list == 1) {
          for (i = 0; i < ppt->selection_num; i++) {
            ptr->selection_magnification_bias[i] = selectionValues[0];
          }
        }
        else if (n_list == ppt->selection_num) {
          for (i = 0; i < n_list; i++) {
            ptr->selection_magnification_bias[i] = selectionValues[i];
          }
        }
        else {
          class_stop(
              "In input for selection function, you asked for %d bin centers and %d bin "
              "biases; number of bins unclear; you should pass either one bin bias (common "
              "to all bins) or %d bin biases",
              ppt->selection_num,
              n_list,
              ppt->selection_num);
        }
      }
    }

    if (ppt->selection_num > 1) {
      psp->non_diag = pfc->get_or("non_diagonal", psp->non_diag);
      if ((psp->non_diag < 0) || (psp->non_diag >= ppt->selection_num))
        class_stop("Input for non_diagonal is %d, while it is expected to be between 0 and %d\n",
                   psp->non_diag,
                   ppt->selection_num - 1);
    }

    if (auto dNdz_selection = pfc->get<std::string>("dNdz_selection")) {
      if ((dNdz_selection->find("analytic") != std::string::npos)) {
        ptr->has_nz_analytic = true;
      }
      else {
        ptr->has_nz_file  = true;
        ptr->nz_file_name = pfc->get_or("dNdz_selection", ptr->nz_file_name);
      }
    }

    if (auto dNdz_evolution = pfc->get<std::string>("dNdz_evolution")) {
      if ((dNdz_evolution->find("analytic") != std::string::npos)) {
        ptr->has_nz_evo_analytic = true;
      }
      else {
        ptr->has_nz_evo_file  = true;
        ptr->nz_evo_file_name = pfc->get_or("dNdz_evolution", ptr->nz_evo_file_name);
      }
    }

    class_test(pfc->get<double>("bias").has_value(),
               "the input parameter 'bias' is obsolete, because you can now pass an independent "
               "light-to-mass bias for each bin/selection function. The new input name is "
               "'selection_bias'. It can be set to a single number (common bias for all bins) or "
               "as many numbers as the number of bins");

    class_test(pfc->get<double>("s_bias").has_value(),
               "the input parameter 's_bias' is obsolete, because you can now pass an independent "
               "magnitude bias for each bin/selection function. The new input name is "
               "'selection_magnitude_bias'. It can be set to a single number (common magnitude "
               "bias for all bins) or as many numbers as the number of bins");
  }
  /* end of selection function section */

  /* deal with z_max issues */
  if ((ppt->has_pk_matter) || (ppt->has_density_transfers) || (ppt->has_velocity_transfers) ||
      (ppt->has_cl_number_count) || (ppt->has_cl_lensing_potential)) {
    if (auto z_max_pk = pfc->get<double>("z_max_pk")) {
      ppt->z_max_pk = *z_max_pk;
    }
    else {
      ppt->z_max_pk = 0.;

      if ((ppt->has_pk_matter) || (ppt->has_density_transfers) || (ppt->has_velocity_transfers)) {
        for (i = 0; i < pop->z_pk_num; i++) {
          ppt->z_max_pk = std::max(ppt->z_max_pk, pop->z_pk[i]);
        }
      }

      if ((ppt->has_cl_number_count) || (ppt->has_cl_lensing_potential)) {
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
          ppt->z_max_pk = std::max(ppt->z_max_pk, z_max);
        }
      }
    }
    psp->z_max_pk = ppt->z_max_pk;
  }
  /* end of z_max section */

  if (auto root = pfc->get<std::string>("root")) {
    pop->root = *root;
  }

  if (auto headers = pfc->get<std::string>("headers");
      headers &&
      ((headers->find("y") == std::string::npos) && (headers->find("Y") == std::string::npos))) {
    pop->write_header = false;
  }

  if (auto format = pfc->get<std::string>("format")) {
    if ((format->find("class") != std::string::npos) ||
        (format->find("CLASS") != std::string::npos))
      pop->output_format = file_format::class_format;
    else {
      if ((format->find("camb") != std::string::npos) ||
          (format->find("CAMB") != std::string::npos))
        pop->output_format = file_format::camb_format;
      else
        class_stop(
            "You wrote: format='%s'. Could not identify any of the possible formats "
            "('class', 'CLASS', 'camb', 'CAMB')",
            format->c_str());
    }
  }

  /** (f) parameter related to the non-linear spectra computation */

  if (auto non_linear = pfc->get<std::string>("non linear")) {
    class_test(!ppt->has_perturbations,
               "You requested non linear computation but no linear computation. You must set "
               "output to tCl or similar.");

    if ((non_linear->find("halofit") != std::string::npos) ||
        (non_linear->find("Halofit") != std::string::npos) ||
        (non_linear->find("HALOFIT") != std::string::npos)) {
      pnl->method       = nl_halofit;
      ppt->k_max_for_pk = std::max(ppt->k_max_for_pk,
                                   std::max(ppr->halofit_min_k_max, ppr->nonlinear_min_k_max));
      ppt->has_nl_corrections_based_on_delta_m = true;
    }
    if ((non_linear->find("hmcode") != std::string::npos) ||
        (non_linear->find("HMCODE") != std::string::npos) ||
        (non_linear->find("HMcode") != std::string::npos) ||
        (non_linear->find("Hmcode") != std::string::npos)) {
      pnl->method       = nl_HMcode;
      ppt->k_max_for_pk = std::max(ppt->k_max_for_pk,
                                   std::max(ppr->hmcode_min_k_max, ppr->nonlinear_min_k_max));
      ppt->has_nl_corrections_based_on_delta_m = true;
      pnl->extrapolation_method =
          (enum source_extrapolation) pfc->get_or<int>("extrapolation_method",
                                                       pnl->extrapolation_method);

      auto feedback_model = pfc->get<std::string>("feedback model");

      if (feedback_model) {
        if (feedback_model->find("emu_dmonly") != std::string::npos) {
          pnl->feedback = nl_emu_dmonly;
        }
        if (feedback_model->find("owls_dmonly") != std::string::npos) {
          pnl->feedback = nl_owls_dmonly;
        }
        if (feedback_model->find("owls_ref") != std::string::npos) {
          pnl->feedback = nl_owls_ref;
        }
        if (feedback_model->find("owls_agn") != std::string::npos) {
          pnl->feedback = nl_owls_agn;
        }
        if (feedback_model->find("owls_dblim") != std::string::npos) {
          pnl->feedback = nl_owls_dblim;
        }
      }

      auto eta_0 = pfc->get<double>("eta_0");
      auto c_min = pfc->get<double>("c_min");

      class_test((feedback_model && (eta_0 || c_min)),
                 "In input file, you cannot enter both a baryonic feedback model and a choice of "
                 "baryonic feedback parameters, choose one of both methods");

      if (eta_0 && c_min) {
        pnl->feedback = nl_user_defined;
        pnl->eta_0    = pfc->get_or("eta_0", pnl->eta_0);
        pnl->c_min    = pfc->get_or("c_min", pnl->c_min);
      }
      else if (eta_0 && !c_min) {
        pnl->feedback = nl_user_defined;
        pnl->eta_0    = pfc->get_or("eta_0", pnl->eta_0);
        pnl->c_min    = (0.98 - pnl->eta_0) / 0.12;
      }
      else if (!eta_0 && c_min) {
        pnl->feedback = nl_user_defined;
        pnl->c_min    = pfc->get_or("c_min", pnl->c_min);
        pnl->eta_0    = 0.98 - 0.12 * pnl->c_min;
      }

      if (pfc->get<double>("z_infinity")) {
        pnl->z_infinity = pfc->get_or("z_infinity", pnl->z_infinity);
      }
    }
  }

  /* The nonlinear module deals only with scalar perturbations, so it is a
     no-op when no scalars are requested. This is a deterministic function of
     the input, so we resolve it here rather than mutating pnl->method from
     inside the module (see issue #20). With no scalars, has_pk_matter is
     already false (enforced above), so clearing method makes the nonlinear
     module skip via its generic "nothing requested" gate. */
  if (!ppt->has_scalars) {
    pnl->method = nl_none;
  }

  /** (g) amount of information sent to standard output (none if all set to zero) */

  // get_or<int> then narrow to short; matches the old class_read_int narrowing cast.
  pba->background_verbose = pfc->get_or<int>("background_verbose", pba->background_verbose);

  pth->thermodynamics_verbose = pfc->get_or<int>("thermodynamics_verbose",
                                                 pth->thermodynamics_verbose);

  ppt->perturbations_verbose = pfc->get_or<int>("perturbations_verbose",
                                                ppt->perturbations_verbose);

  ptr->transfer_verbose = pfc->get_or<int>("transfer_verbose", ptr->transfer_verbose);

  ppm->primordial_verbose = pfc->get_or<int>("primordial_verbose", ppm->primordial_verbose);

  psp->spectra_verbose = pfc->get_or<int>("spectra_verbose", psp->spectra_verbose);

  pnl->nonlinear_verbose = pfc->get_or<int>("nonlinear_verbose", pnl->nonlinear_verbose);

  ple->lensing_verbose = pfc->get_or<int>("lensing_verbose", ple->lensing_verbose);

  pop->output_verbose = pfc->get_or<int>("output_verbose", pop->output_verbose);

  if (ppt->has_tensors) {
    /** - ---> Include ur and ncdm shear in tensor computation? */
    if (auto tensor_method = pfc->get<std::string>("tensor method")) {
      if (tensor_method->find("photons") != std::string::npos)
        ppt->tensor_method = tm_photons_only;
      if (tensor_method->find("massless") != std::string::npos)
        ppt->tensor_method = tm_massless_approximation;
      if (tensor_method->find("exact") != std::string::npos)
        ppt->tensor_method = tm_exact;
    }
  }

  /** - ---> derivatives of baryon sound speed only computed if some non-minimal tight-coupling schemes is requested */
  if ((ppr->tight_coupling_approximation == static_cast<int>(tca_method::first_order_CLASS)) ||
      (ppr->tight_coupling_approximation == static_cast<int>(tca_method::second_order_CLASS))) {
    pth->compute_cb2_derivatives = true;
  }

  class_test(ppr->ur_fluid_trigger_tau_over_tau_k ==
                 ppr->radiation_streaming_trigger_tau_over_tau_k,
             "please choose different values for precision parameters "
             "ur_fluid_trigger_tau_over_tau_k and radiation_streaming_trigger_tau_over_tau_k, in "
             "order to avoid switching two approximation schemes at the same time");

  if (omega_budget_.idr.value_or(0.) != 0.) {
    class_test(ppr->idr_streaming_trigger_tau_over_tau_k ==
                   ppr->radiation_streaming_trigger_tau_over_tau_k,
               "please choose different values for precision parameters "
               "dark_radiation_trigger_tau_over_tau_k and "
               "radiation_streaming_trigger_tau_over_tau_k, in order to avoid switching two "
               "approximation schemes at the same time");

    class_test(ppr->idr_streaming_trigger_tau_over_tau_k == ppr->ur_fluid_trigger_tau_over_tau_k,
               "please choose different values for precision parameters "
               "dark_radiation_trigger_tau_over_tau_k and ur_fluid_trigger_tau_over_tau_k, in "
               "order to avoid switching two approximation schemes at the same time");

    class_test(ppr->idr_streaming_trigger_tau_over_tau_k == ppr->ncdm_fluid_trigger_tau_over_tau_k,
               "please choose different values for precision parameters "
               "dark_radiation_trigger_tau_over_tau_k and ncdm_fluid_trigger_tau_over_tau_k, in "
               "order to avoid switching two approximation schemes at the same time");
  }

  /**
   * Here we can place all obsolete (deprecated) names for the precision parameters,
   * so they will still get read.
   * The new parameter names should be used preferrably
   * */
  // obsolete precision parameter: read for compatibility with old precision files
  ppr->k_min_tau0 = pfc->get_or("k_scalar_min_tau0", ppr->k_min_tau0);
  // obsolete precision parameter: read for compatibility with old precision files
  ppr->k_max_tau0_over_l_max = pfc->get_or("k_scalar_max_tau0_over_l_max",
                                           ppr->k_max_tau0_over_l_max);
  // obsolete precision parameter: read for compatibility with old precision files
  ppr->k_step_sub = pfc->get_or("k_scalar_step_sub", ppr->k_step_sub);
  // obsolete precision parameter: read for compatibility with old precision files
  ppr->k_step_super = pfc->get_or("k_scalar_step_super", ppr->k_step_super);
  // obsolete precision parameter: read for compatibility with old precision files
  ppr->k_step_transition = pfc->get_or("k_scalar_step_transition", ppr->k_step_transition);
  // obsolete precision parameter: read for compatibility with old precision files
  ppr->k_per_decade_for_pk = pfc->get_or("k_scalar_k_per_decade_for_pk", ppr->k_per_decade_for_pk);
  // obsolete precision parameter: read for compatibility with old precision files
  ppr->k_per_decade_for_bao = pfc->get_or("k_scalar_k_per_decade_for_bao",
                                          ppr->k_per_decade_for_bao);
  // obsolete precision parameter: read for compatibility with old precision files
  ppr->k_bao_center = pfc->get_or("k_scalar_bao_center", ppr->k_bao_center);
  // obsolete precision parameter: read for compatibility with old precision files
  ppr->k_bao_width = pfc->get_or("k_scalar_bao_width", ppr->k_bao_width);

  // obsolete precision parameter: read for compatibility with old precision files
  ppr->q_linstep = pfc->get_or("k_step_trans_scalars", ppr->q_linstep);
  // obsolete precision parameter: read for compatibility with old precision files
  ppr->q_linstep = pfc->get_or("k_step_trans_tensors", ppr->q_linstep);
  // obsolete precision parameter: read for compatibility with old precision files
  ppr->q_linstep = pfc->get_or("k_step_trans", ppr->q_linstep);
  // obsolete precision parameter: read for compatibility with old precision files
  ppr->q_linstep = pfc->get_or("q_linstep_trans", ppr->q_linstep);
  // obsolete precision parameter: read for compatibility with old precision files
  ppr->q_logstep_spline = pfc->get_or("q_logstep_trans", ppr->q_logstep_spline);

  class_test(pfc->get<std::string>("l_switch_limber_for_cl_density_over_z").has_value(),
             "You passed in input a precision parameter called "
             "l_switch_limber_for_cl_density_over_z. This syntax is deprecated since v2.5.0. "
             "Please use instead the two precision parameters l_switch_limber_for_nc_local_over_z, "
             "l_switch_limber_for_nc_los_over_z, defined in include/common.h, and allowing for "
             "better performance.");

  /** (i) Write values in file */
  if (ple->has_lensed_cls)
    ppt->l_scalar_max += ppr->delta_l_max;

  /** - (i.1.) shall we write background quantities in a file? */

  if (auto write_background = pfc->get<std::string>("write background");
      write_background && ((write_background->find("y") != std::string::npos) ||
                           (write_background->find("Y") != std::string::npos))) {
    pop->write_background = true;
  }

  /** - (i.2.) shall we write thermodynamics quantities in a file? */

  if (auto write_thermodynamics = pfc->get<std::string>("write thermodynamics");
      write_thermodynamics && ((write_thermodynamics->find("y") != std::string::npos) ||
                               (write_thermodynamics->find("Y") != std::string::npos))) {
    pop->write_thermodynamics = true;
  }

  /** - (i.3.) shall we write perturbation quantities in files? */

  std::vector<double> kOutputValues;
  int k_output_found;
  readDoubleList(pfc, "k_output_values", kOutputValues, &k_output_found);

  if (k_output_found) {
    n_list = static_cast<int>(kOutputValues.size());
    class_test(n_list > _MAX_NUMBER_OF_K_FILES_,
               "you want to write some output for %d different values of k, hence you should "
               "increase _MAX_NUMBER_OF_K_FILES_ in include/perturbations.h to at least this "
               "number",
               n_list);
    ppt->k_output_values_num = n_list;

    for (i = 0; i < n_list; i++) {
      ppt->k_output_values[i] = kOutputValues[i];
    }

    /* Sort k_output_values ascending */
    std::sort(ppt->k_output_values, ppt->k_output_values + ppt->k_output_values_num);

    ppt->store_perturbations = true;
    pop->write_perturbations = true;
  }

  /** - (i.4.) shall we write primordial spectra in a file? */

  if (auto write_primordial = pfc->get<std::string>("write primordial");
      write_primordial && ((write_primordial->find("y") != std::string::npos) ||
                           (write_primordial->find("Y") != std::string::npos))) {
    pop->write_primordial = true;
  }

  /* flags for calling the interpolation routine */
  pba->short_info  = 0;
  pba->normal_info = 1;
  pba->long_info   = 2;

  pba->inter_normal  = 0;
  pba->inter_closeby = 1;

  /** - (i.5) special steps if we want Halofit with wa_fld non-zero:
      so-called "Pk_equal method" of 0810.0190 and 1601.07230 */

  /* ReadDerived runs after ConstructSpecies, so query the built FluidSpecies
     directly instead of peeking the raw parameter file. The AND-chain matches
     the original gate: halofit enabled + fluid present + CLP EoS + wa_fld≠0. */
  if ((pnl->method == nl_halofit) && all_species_.count("Fluid")) {
    const auto& fluid = static_cast<const FluidSpecies&>(*all_species_.at("Fluid"));
    if ((fluid.fluid_eos() == CLP) && (fluid.wa_fld() != 0.)) {
      if (auto pk_eq = pfc->get<std::string>("pk_eq");
          pk_eq &&
          ((pk_eq->find("y") != std::string::npos) || (pk_eq->find("Y") != std::string::npos))) {
        pnl->has_pk_eq = true;
      }
    }
  }

  /* ── Precision-consistency tests for perturbation-hierarchy l_max values.
     Moved here from perturb_vector_init: these checks depend only on ppr
     (plus ppt->idr_nature for the IDR test), so they belong at input-parse
     time, not inside the per-(k, approximation) hot path.  Tests run
     unconditionally — a too-low l_max is a user-config error whether or not
     the species ends up active. */
  class_test(ppr->l_max_g < 4,
             "ppr->l_max_g should be at least 4, i.e. we must integrate at least over photon "
             "density, velocity, shear, third and fourth momentum");
  class_test(ppr->l_max_pol_g < 4, "ppr->l_max_pol_g should be at least 4");
  class_test(ppr->l_max_ur < 4,
             "ppr->l_max_ur should be at least 4, i.e. we must integrate at least over "
             "neutrino/relic density, velocity, shear, third and fourth momentum");
  class_test(ppr->l_max_dr < 4,
             "ppr->l_max_dr should be at least 4, i.e. we must integrate at least over "
             "neutrino/relic density, velocity, shear, third and fourth momentum");
  // idr_nature now lives on the IDM_DR_IDR species. Preserve the original
  // semantics exactly: with no IDM_DR_IDR species the legacy ppt->idr_nature
  // defaulted to idr_free_streaming, so the test fired iff l_max_idr < 4.
  const int eff_idr_nature = all_species_.count("IDM_DR_IDR")
                                 ? static_cast<const IDM_DR_IDR_Species&>(
                                       *all_species_.at("IDM_DR_IDR"))
                                       .idr()
                                       .idr_nature()
                                 : idr_free_streaming;
  class_test((ppr->l_max_idr < 4) && (eff_idr_nature == idr_free_streaming),
             "ppr->l_max_idr should be at least 4, i.e. we must integrate at least over "
             "interacting dark radiation density, velocity, shear, third and fourth momentum");
  class_test(ppr->l_max_g_ten < 4,
             "ppr->l_max_g_ten should be at least 4, i.e. we must integrate at least over photon "
             "density, velocity, shear, third momentum");
  class_test(ppr->l_max_pol_g_ten < 4, "ppr->l_max_pol_g_ten should be at least 4");
  class_test(ppr->l_max_ncdm < 4,
             "ppr->l_max_ncdm=%d should be at least 4, i.e. we must integrate at least "
             "over first four momenta of non-cold dark matter perturbed phase-space "
             "distribution",
             ppr->l_max_ncdm);
}

/** Overloaded helpers for type-dispatched precision parameter reading. */
namespace {
void read(const FileContent& fc, const char* name, double& v) {
  v = fc.get_or(name, v);
}
void read(const FileContent& fc, const char* name, int& v) {
  v = fc.get_or(name, v);
}
void read(const FileContent& fc, const char* name, std::string& v) {
  v = fc.get_or(name, v);
}
template <typename E>
void read_enum(const FileContent& fc, const char* name, E& v) {
  v = static_cast<E>(fc.get_or<int>(name, static_cast<int>(v)));
}
}  // anonymous namespace

void precision::ResolveDataPaths() {
  // Prepend the runtime class_dir to each field's relative-path default.
  sBBN_file                    = class_dir + sBBN_file;
  hyrec_Alpha_inf_file         = class_dir + hyrec_Alpha_inf_file;
  hyrec_R_inf_file             = class_dir + hyrec_R_inf_file;
  hyrec_two_photon_tables_file = class_dir + hyrec_two_photon_tables_file;
}

void precision::parse(const FileContent& fc) {
  /* Background */
  read(fc, "a_ini_over_a_today_default", a_ini_over_a_today_default);
  read(fc, "back_integration_stepsize", back_integration_stepsize);
  read(fc, "tol_background_integration", tol_background_integration);
  read(fc, "tol_initial_Omega_r", tol_initial_Omega_r);
  read(fc, "tol_M_ncdm", tol_M_ncdm);
  // "tol_ncdm" is a convenience input that sets BOTH gauge-specific tolerances;
  // the gauge-specific keys then override per gauge if also given.
  read(fc, "tol_ncdm", tol_ncdm_synchronous);
  read(fc, "tol_ncdm", tol_ncdm_newtonian);
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
  read(fc,
       "perturbations_sampling_boost_above_age_fraction",
       perturbations_sampling_boost_above_age_fraction);
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
  read(fc, "q_logstep_limber", q_logstep_limber);
  read(fc, "k_max_limber_over_l_max_scalars", k_max_limber_over_l_max_scalars);
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

// ── Hook-based shooting (the species-owned replacement for the enum dispatch) ──

void InputModule::ShootingResidual(double* x, int x_size, void* pworkspace, double* output) {
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
}

InputModulePtr InputModule::DoShooting(InputModulePtr input_module) {
  // Guard: never shoot from within a shooting build (the residual marks its fc this way).
  if (input_module->file_content_.is_shooting)
    return input_module;

  FileContent& fc = input_module->file_content_;

  ShootingWorkspace w(fc);
  std::vector<double> xguess, dxdF;

  // Cosmological target: 100*theta_s varies h (residual from thermo rs/ra). Module-level.
  if (auto theta_s = fc.get<double>("100*theta_s")) {
    const double tv = *theta_s;
    w.targets.push_back({"100*theta_s", "h", tv});
    w.target_species_keys.emplace_back();  // module-level: no owning species
    xguess.push_back(3.54 * tv * tv - 5.455 * tv + 2.548);
    dxdF.push_back(7.08 * tv - 5.455);
  }

  // Per-species targets (all_species_ lex order). The discovery module's species already
  // guessed their unknowns at construction; query them for {target, guess} here.
  // SpeciesBuildContext::ncdm_settings is non-null by contract (some guesses fall back to
  // it when ppr is unavailable); build it from the module as ConstructSpecies does.
  NcdmSettings ncdm_settings;
  ncdm_settings.h     = input_module->background_.h;
  ncdm_settings.T_cmb = input_module->background_.T_cmb;
  // Gauge-dependent ncdm momentum tolerance (see ConstructSpecies for rationale).
  ncdm_settings.tol_ncdm    = (input_module->perturbations_.gauge == possible_gauges::newtonian)
                                  ? input_module->precision_.tol_ncdm_newtonian
                                  : input_module->precision_.tol_ncdm_synchronous;
  ncdm_settings.tol_ncdm_bg = input_module->precision_.tol_ncdm_bg;
  ncdm_settings.tol_M_ncdm  = input_module->precision_.tol_M_ncdm;
  SpeciesBuildContext gctx{&fc,
                           &input_module->background_,
                           &input_module->precision_,
                           &ncdm_settings,
                           /*bgm=*/nullptr,
                           /*all_species=*/&input_module->all_species_,
                           /*omega_budget=*/&input_module->omega_budget_};
  gctx.coupled_inputs = &input_module->coupled_inputs_;
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
  fzero_Newton(ShootingResidual,
               xguess.data(),
               dxdF.data(),
               static_cast<int>(w.targets.size()),
               1e-3,
               1e-3,
               &w,
               &fevals);

  // Write the resolved unknowns; build and return a fresh, fully-resolved module.
  for (size_t i = 0; i < w.targets.size(); ++i) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%.17g", xguess[i]);
    fc.set(w.targets[i].unknown_param, buf);
  }
  return std::make_shared<InputModule>(fc);
}
