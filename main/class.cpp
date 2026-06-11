/** @file class.c
 * Julien Lesgourgues, 17.04.2011
 */

#include "class.h"

#include "cosmology.h"
#include "output_module.h"

int main(int argc, char** argv) {
  try {
    FileContent fc;
    InputModule::file_content_from_arguments(argc, argv, fc);

    Cosmology cosmology{fc};

    OutputModule output_module(cosmology.GetInputModule(),
                               cosmology.GetBackgroundModule(),
                               cosmology.GetThermodynamicsModule(),
                               cosmology.GetPerturbationsModule(),
                               cosmology.GetPrimordialModule(),
                               cosmology.GetNonlinearModule(),
                               cosmology.GetSpectraModule(),
                               cosmology.GetLensingModule());
  }
  catch (const std::exception& e) {
    printf("\n\nError in Class:\n%s\n", e.what());
    return _FAILURE_;
  }

  return _SUCCESS_;
}
