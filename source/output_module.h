#ifndef OUTPUT_MODULE_H
#define OUTPUT_MODULE_H

#include <string>

#include "base_module.h"
#include "input_module.h"

class OutputModule : public BaseModule {
 public:
  OutputModule(InputModulePtr input_module,
               BackgroundModulePtr background_module,
               ThermodynamicsModulePtr thermodynamics_module,
               PerturbationsModulePtr perturbations_module,
               PrimordialModulePtr primordial_module,
               NonlinearModulePtr nonlinear_module,
               SpectraModulePtr spectra_module,
               LensingModulePtr lensing_module);

 private:
  void output_total_cl_at_l(int l, double* cl);
  void output_init();
  void output_cl();
  void output_pk(enum pk_outputs pk_output);
  void output_tk();
  void output_background();
  void output_thermodynamics();
  void output_perturbations();
  void output_primordial();
  void output_print_data(FILE* out, const std::string& titles, const double* dataptr, int tau_size);
  void output_open_cl_file(FILE** clfile,
                           const std::string& filename,
                           const char* first_line,
                           int lmax);
  void output_one_line_of_cl(FILE* clfile, double l, double* cl, int ct_size);
  void output_open_pk_file(FILE** pkfile,
                           const std::string& filename,
                           const char* first_line,
                           double z);
  void output_one_line_of_pk(FILE* tkfile, double one_k, double one_pk);

  BackgroundModulePtr background_module_;
  ThermodynamicsModulePtr thermodynamics_module_;
  PerturbationsModulePtr perturbations_module_;
  PrimordialModulePtr primordial_module_;
  NonlinearModulePtr nonlinear_module_;
  SpectraModulePtr spectra_module_;
  LensingModulePtr lensing_module_;
};

#endif  //OUTPUT_MODULE_H
