// classy._core, the compiled half of classy.Class, in nanobind. It holds the C++ Cosmology's
// lifetime and parameter gate, one method per C++ call, exception translation and loops over
// input arrays, and nothing else: physics and checks are in C++, composition in the plain-Python
// package (_class.py subclasses PyCosmology).
//
// It replaced a Cython core with the same interface, call for call
// (docs/superpowers/specs/2026-11-21-classy-nanobind-results.md); the method and argument names,
// return types and exceptions are classy's public API, which Cobaya, MontePython and the tests
// use.
#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <new>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#include "background_module.h"
#include "bessel_table_cache.h"
#include "cosmology.h"
#include "hyperspherical.h"
#include "input_module.h"
#include "lensing_module.h"
#include "nonlinear_module.h"
#include "output_table.h"
#include "parser.h"
#include "perturbations_module.h"
#include "primordial_module.h"
#include "spectra_module.h"
#include "thermodynamics_module.h"
#include "version.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace {

using Vector1d  = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using Vector3d  = nb::ndarray<const double, nb::ndim<3>, nb::c_contig, nb::device::cpu>;
using IntVector = nb::ndarray<const int, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

// Set when the module is imported, and kept for the life of the interpreter.
PyObject* severe_error      = nullptr;  // <pkg>._errors.CosmoSevereError
PyObject* computation_error = nullptr;  // <pkg>._errors.CosmoComputationError
std::string data_dir;                   // <pkg>/data: bbn/, hyrec/, the default class_dir

// The C++ errors as classy's exceptions (STYLE.md, error conventions): by the dynamic type's
// name, std::invalid_argument (class_test_severe) is CosmoSevereError, std::runtime_error
// (class_test) is CosmoComputationError, and anything else is CosmoSevereError, because a
// sampler must not take an unexpected failure for a bad point. The exact type, not a base: a
// std::overflow_error is severe too. The message goes as bytes; CosmoError decodes it.
void TranslateException(const std::exception_ptr& exception, void* /*payload*/) {
  try {
    std::rethrow_exception(exception);
  }
  catch (const std::exception& e) {
    const std::string type = typeid(e).name();
    PyObject* error        = severe_error;
    if (type.find("invalid_argument") == std::string::npos &&
        type.find("runtime_error") != std::string::npos) {
      error = computation_error;
    }
    PyErr_SetObject(error, nb::bytes(e.what()).ptr());
  }
}

// The Bessel functions are not a cosmology: their C++ errors become Python's built-in
// exceptions, invalid_argument a ValueError and so on, as they always did.
[[noreturn]] void RaiseAsBuiltin(const std::exception& e) {
  PyObject* type = PyExc_RuntimeError;
  if (dynamic_cast<const std::bad_alloc*>(&e) != nullptr) {
    PyErr_NoMemory();
    throw nb::python_error();
  }
  if (dynamic_cast<const std::bad_cast*>(&e) != nullptr ||
      dynamic_cast<const std::bad_typeid*>(&e) != nullptr) {
    type = PyExc_TypeError;
  }
  else if (dynamic_cast<const std::domain_error*>(&e) != nullptr ||
           dynamic_cast<const std::invalid_argument*>(&e) != nullptr) {
    type = PyExc_ValueError;
  }
  else if (dynamic_cast<const std::out_of_range*>(&e) != nullptr) {
    type = PyExc_IndexError;
  }
  else if (dynamic_cast<const std::overflow_error*>(&e) != nullptr) {
    type = PyExc_OverflowError;
  }
  else if (dynamic_cast<const std::range_error*>(&e) != nullptr ||
           dynamic_cast<const std::underflow_error*>(&e) != nullptr) {
    type = PyExc_ArithmeticError;
  }
  PyErr_SetString(type, e.what());
  throw nb::python_error();
}

bool Truthy(nb::handle value) {
  const int truth = PyObject_IsTrue(value.ptr());
  if (truth < 0) {
    throw nb::python_error();
  }
  return truth != 0;
}

bool Equal(nb::handle a, nb::handle b) {
  const int equal = PyObject_RichCompareBool(a.ptr(), b.ptr(), Py_EQ);
  if (equal < 0) {
    throw nb::python_error();
  }
  return equal != 0;
}

// An integer argument as classy has always taken one: an int, or anything with __int__, a
// float truncated (raw_cl(2.5e3) is raw_cl(2500)). Out of the type's range is an OverflowError.
template <class Int>
Int IntArg(nb::handle value) {
  nb::object number;
  if (PyLong_Check(value.ptr())) {
    number = nb::borrow(value);
  }
  else {
    PyNumberMethods* methods = Py_TYPE(value.ptr())->tp_as_number;
    if (methods == nullptr || methods->nb_int == nullptr) {
      PyErr_Format(PyExc_TypeError,
                   "an integer is required, not '%s'",
                   Py_TYPE(value.ptr())->tp_name);
      throw nb::python_error();
    }
    number = nb::steal(methods->nb_int(value.ptr()));
    if (!number.is_valid()) {
      throw nb::python_error();
    }
  }
  const long long result = PyLong_AsLongLong(number.ptr());
  if (result == -1 && PyErr_Occurred() != nullptr) {
    throw nb::python_error();
  }
  if (result < std::numeric_limits<Int>::min() || result > std::numeric_limits<Int>::max()) {
    PyErr_SetString(PyExc_OverflowError, "value too large to convert to int");
    throw nb::python_error();
  }
  return static_cast<Int>(result);
}

// A str's or bytes' content up to its first NUL, as a char* hands it to std::string.
std::string Bytes(const nb::object& bytes) {
  const char* data = PyBytes_AsString(bytes.ptr());
  if (data == nullptr) {
    throw nb::python_error();
  }
  return data;
}

nb::module_ Numpy() {
  return nb::module_::import_("numpy");
}

nb::object Builtin(const char* name) {
  return nb::module_::import_("builtins").attr(name);
}

// A numpy array that owns `values`, moved into a capsule: nothing is copied.
nb::object ToArray(std::vector<double>&& values, std::vector<size_t> shape = {}) {
  if (shape.empty()) {
    shape.push_back(values.size());
  }
  auto* owned = new std::vector<double>(std::move(values));
  owned->reserve(1);  // a valid pointer for an empty array
  nb::capsule owner(owned, [](void* p) noexcept { delete static_cast<std::vector<double>*>(p); });
  return nb::ndarray<nb::numpy, double>(owned->data(), shape.size(), shape.data(), owner).cast();
}

// An output table as a dict of arrays, in column order. The table moves into a capsule and
// every column is a view of it: nothing is copied.
nb::dict TableDict(OutputTable&& table) {
  auto* owned = new OutputTable(std::move(table));
  nb::capsule owner(owned, [](void* p) noexcept { delete static_cast<OutputTable*>(p); });
  nb::dict out;
  for (size_t index = 0; index < owned->names().size(); ++index) {
    const OutputColumn column = owned->column(index);
    const nb::str name(owned->names()[index].c_str());
    if (column.size() == 0) {
      out[name] = ToArray({});
      continue;
    }
    const size_t shape[1] = {column.size()};
    // The table is ours, so its arrays may be written to, as any numpy array; the column is a
    // read-only view only to the table's other users.
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
    auto* data = const_cast<double*>(column.data());
    out[name]  = nb::ndarray<nb::numpy, double>(data, 1, shape, owner).cast();
  }
  return out;
}

nb::dict ClDict(std::map<std::string, std::vector<double>> cl) {
  nb::dict out;
  for (auto& [name, values] : cl) {
    out[nb::str(name.c_str())] = ToArray(std::move(values));
  }
  return out;
}

void CheckSizes(const char* what, std::initializer_list<std::pair<Py_ssize_t, size_t>> sizes) {
  int axis = 0;
  for (const auto& [asked, size] : sizes) {
    if (asked < 0) {
      throw nb::value_error("negative dimensions are not allowed");
    }
    if (static_cast<size_t>(asked) > size) {
      throw nb::index_error((std::string("Out of bounds on buffer access (") + what + ", axis " +
                             std::to_string(axis) + ")")
                                .c_str());
    }
    ++axis;
  }
}

class PyCosmology {
 public:
  explicit PyCosmology(nb::handle input_parameters) {
    // Copy: never alias the caller's dict, so external mutation cannot bypass the
    // parameters_changed gate (and the class_dir default does not leak into the caller's dict).
    pars_ = nb::borrow<nb::dict>(
        Builtin("dict")(input_parameters.is_none() ? nb::dict() : nb::borrow(input_parameters)));
    if (!pars_.contains("class_dir")) {
      pars_["class_dir"] = data_dir;
    }
    Reset();
  }

  nb::dict pars() const {
    // Snapshot, not the live dict: a caller mutating the returned dict would otherwise bypass
    // the parameters_changed gate and stale-serve. set() is the mutation API.
    return nb::steal<nb::dict>(PyDict_Copy(pars_.ptr()));
  }

  void Empty() {
    pars_               = nb::dict();
    parameters_changed_ = true;
  }

  // Keywords as well as a dict: Cobaya calls set(**args). Keys that are not Python
  // identifiers ('100*theta_s', 'nu1.m') can only travel in the dict. False when nothing
  // changed.
  bool Set(nb::handle input_parameters, const nb::kwargs& kwargs) {
    nb::object merged = Builtin("dict")(Truthy(input_parameters) ? nb::borrow(input_parameters)
                                                                 : nb::dict());
    if (PyDict_Update(merged.ptr(), kwargs.ptr()) < 0) {
      throw nb::python_error();
    }
    const int subset =
        PyObject_RichCompareBool(merged.attr("items")().ptr(), pars_.attr("items")().ptr(), Py_LE);
    if (subset < 0) {
      throw nb::python_error();
    }
    if (subset != 0) {
      return false;
    }
    if (PyDict_Update(pars_.ptr(), merged.ptr()) < 0) {
      throw nb::python_error();
    }
    parameters_changed_ = true;
    return true;
  }

  // Single gate to the C++ object: rebuilds the Cosmology if the parameter dict changed since
  // it was last built. Everything reaches the C++ object through here, so results always
  // reflect the current parameters.
  Cosmology& cosmo() {
    if (parameters_changed_ || !cosmology_) {
      Reset();
    }
    return *cosmology_;
  }

  // The input structs live inside the InputModule, and DoShooting replaces the InputModule
  // instance on first use; fetching it on every access is the only way a cached pointer
  // cannot go stale.
  const perturbs& pt() {
    return cosmo().GetInputModule()->perturbations_;
  }
  const nonlinear& nl() {
    return cosmo().GetInputModule()->nonlinear_;
  }
  // The computed modules, each built on first use.
  const BackgroundModule& background() {
    return *cosmo().GetBackgroundModule();
  }
  const ThermodynamicsModule& thermo() {
    return *cosmo().GetThermodynamicsModule();
  }
  const PrimordialModule& primordial() {
    return *cosmo().GetPrimordialModule();
  }
  const NonlinearModule& nonlinear_module() {
    return *cosmo().GetNonlinearModule();
  }
  const PerturbationsModule& perturbations() {
    return *cosmo().GetPerturbationsModule();
  }

  pk_outputs DefaultPkOutput() {
    // Which P(k) a bare pk()/get_pk() call means: non-linear if one was asked for.
    return nl().method == nl_none ? pk_linear : pk_nonlinear;
  }

  double Pk(double k, double z, int index_pk, pk_outputs output) {
    double pk = 0.;
    cosmo().GetNonlinearModule()->nonlinear_pk_at_k_and_z(output, k, z, index_pk, &pk, nullptr);
    return pk;
  }

  nb::object GetPk(const Vector3d& k,
                   const Vector1d& z,
                   nb::handle k_size_arg,
                   nb::handle z_size_arg,
                   nb::handle mu_size_arg,
                   int index_pk,
                   pk_outputs output) {
    const auto k_size  = IntArg<Py_ssize_t>(k_size_arg);
    const auto z_size  = IntArg<Py_ssize_t>(z_size_arg);
    const auto mu_size = IntArg<Py_ssize_t>(mu_size_arg);
    CheckSizes("k", {{k_size, k.shape(0)}, {z_size, k.shape(1)}, {mu_size, k.shape(2)}});
    CheckSizes("z", {{z_size, z.shape(0)}});
    std::vector<double> pk(static_cast<size_t>(k_size * z_size * mu_size));
    for (Py_ssize_t index_k = 0; index_k < k_size; index_k++) {
      for (Py_ssize_t index_z = 0; index_z < z_size; index_z++) {
        for (Py_ssize_t index_mu = 0; index_mu < mu_size; index_mu++) {
          pk[(index_k * z_size + index_z) * mu_size + index_mu] =
              Pk(k(index_k, index_z, index_mu), z(index_z), index_pk, output);
        }
      }
    }
    return ToArray(std::move(pk),
                   {static_cast<size_t>(k_size),
                    static_cast<size_t>(z_size),
                    static_cast<size_t>(mu_size)});
  }

  nb::dict GeneralCl(int lmax, bool is_lensed) {
    std::map<std::string, std::vector<double>> cl;
    if (is_lensed) {
      const auto lensing = cosmo().GetLensingModule();
      if (lmax == -1) {
        lmax = lensing->l_lensed_max_;
      }
      cl = lensing->cl_output(lmax);
    }
    else {
      const auto spectra = cosmo().GetSpectraModule();
      if (lmax == -1) {
        lmax = spectra->l_max_tot_;
      }
      cl = spectra->cl_output(lmax);
    }
    nb::dict out = ClDict(std::move(cl));
    out["ell"]   = Numpy().attr("arange")(lmax + 1);
    return out;
  }

  nb::dict RawClNoCopy(int lmax) {
    const auto spectra = cosmo().GetSpectraModule();
    if (lmax == -1) {
      lmax = spectra->l_max_tot_;
    }
    CheckSizes("lmax + 1", {{static_cast<Py_ssize_t>(lmax) + 1, SIZE_MAX}});
    const std::map<std::string, int> index_map = spectra->cl_output_index_map();
    std::vector<double*> output_pointers(spectra->ct_size_, nullptr);
    std::map<std::string, std::vector<double>> cl;
    for (const auto& [name, index] : index_map) {
      std::vector<double>& values = cl[name];
      values.resize(lmax + 1);
      output_pointers[index] = values.data();
    }
    spectra->cl_output_no_copy(lmax, output_pointers);
    nb::dict out = ClDict(std::move(cl));
    out["ell"]   = Numpy().attr("arange")(lmax + 1);
    return out;
  }

  nb::object TransferGeneral(double z, int index_tau, nb::handle output_format) {
    if (!pt().has_density_transfers && !pt().has_velocity_transfers) {
      return nb::dict();
    }
    const file_format format = Equal(output_format, nb::str("camb")) ? file_format::camb_format
                                                                     : file_format::class_format;
    const PerturbationsModule& module = perturbations();
    std::vector<OutputTable> tables   = index_tau < 0
                                            ? module.GetTransferTables(format, z)
                                            : module.GetTransferTablesAtIndexTau(format, index_tau);
    if (tables.size() == 1) {
      return TableDict(std::move(tables[0]));
    }
    nb::dict transfers;
    for (size_t index_ic = 0; index_ic < tables.size(); index_ic++) {
      std::string ic_info;
      std::string ic_suffix;
      module.perturb_output_firstline_and_ic_suffix(static_cast<int>(index_ic), ic_info, ic_suffix);
      transfers[nb::str(ic_suffix.c_str())] = TableDict(std::move(tables[index_ic]));
    }
    return transfers;
  }

  nb::dict GetPerturbations() {
    nb::dict perturbations_out;
    if (pt().k_output_values_num < 1) {
      return perturbations_out;
    }
    const PerturbationsModule& module = perturbations();
    const struct {
      const char* name;
      bool has_mode;
      int index_md;
    } modes[] = {{"scalar", pt().has_scalars, module.index_md_scalars_},
                 {"vector", pt().has_vectors, module.index_md_vectors_},
                 {"tensor", pt().has_tensors, module.index_md_tensors_}};
    for (const auto& mode : modes) {
      if (!mode.has_mode) {
        continue;
      }
      std::vector<OutputTable> tables = module.GetPerturbationTables(mode.index_md);
      nb::list per_k;
      for (OutputTable& table : tables) {
        per_k.append(TableDict(std::move(table)));
      }
      perturbations_out[mode.name] = per_k;
    }
    return perturbations_out;
  }

  nb::object NonlinearScales(const Vector1d& z, nb::handle z_size_arg, bool cb) {
    const auto z_size = IntArg<Py_ssize_t>(z_size_arg);
    CheckSizes("z", {{z_size, z.shape(0)}});
    const NonlinearModule& module = nonlinear_module();
    std::vector<double> k_nl(static_cast<size_t>(z_size));
    for (Py_ssize_t index_z = 0; index_z < z_size; index_z++) {
      double k_nl_value    = 0.;
      double k_nl_cb_value = 0.;
      module.nonlinear_k_nl_at_z(z(index_z), &k_nl_value, &k_nl_cb_value);
      k_nl[index_z] = cb ? k_nl_cb_value : k_nl_value;
    }
    return ToArray(std::move(k_nl));
  }

  // P(k, z) and P_cb(k, z) on a k and a z array, flat; P_cb is P when the module stores no
  // separate cdm+baryon spectrum.
  std::pair<nb::object, nb::object> GetPkArrays(const Vector1d& k,
                                                const Vector1d& z,
                                                nb::handle k_size_arg,
                                                nb::handle z_size_arg,
                                                nb::handle nonlinear) {
    const int k_size        = IntArg<int>(k_size_arg);
    const int z_size        = IntArg<int>(z_size_arg);
    const pk_outputs output = Equal(nonlinear, nb::int_(0)) ? pk_linear : pk_nonlinear;
    // Empty k and z arrays are refused, as classy always did, and so are sizes beyond them.
    CheckSizes("k", {{std::max<Py_ssize_t>(k_size, 1), k.shape(0)}});
    CheckSizes("z", {{std::max<Py_ssize_t>(z_size, 1), z.shape(0)}});
    CheckSizes("pk", {{static_cast<Py_ssize_t>(k_size) * z_size, SIZE_MAX}});
    std::vector<double> pk(static_cast<size_t>(k_size) * z_size);
    std::vector<double> pk_cb(pk.size());
    const auto module = cosmo().GetNonlinearModule();
    module->nonlinear_pks_at_kvec_and_zvec(output,
                                           k.data(),
                                           k_size,
                                           z.data(),
                                           z_size,
                                           pk.data(),
                                           pk_cb.data());
    nb::object pk_array = ToArray(std::move(pk));
    // The module stores P_cb apart only when it differs from P_m.
    nb::object pk_cb_array = module->has_pk_cb_ ? ToArray(std::move(pk_cb)) : pk_array;
    return {pk_array, pk_cb_array};
  }

 private:
  void Reset() {
    UpdateFileContent();
    if (!bessel_cache_) {
      bessel_cache_ = std::make_shared<BesselTableCache>();
    }
    // Cosmology(FileContent&) rejects parameters that nothing read
    // (InputModule::RejectUnreadParameters), unless the input sets allow_unread_parameters.
    cosmology_ = std::make_unique<Cosmology>(fc_, bessel_cache_);
    // Only clear the dirty flag after all throwing operations have succeeded, so that a failed
    // construction is retried and re-raised on the next call instead of serving stale results.
    parameters_changed_ = false;
  }

  void UpdateFileContent() {
    fc_ = FileContent();
    // empty() and set() can leave the parameters without class_dir; the C++ default would then
    // be the source tree this module was built from, which a wheel need not have.
    if (!pars_.contains("class_dir")) {
      fc_.set("class_dir", data_dir);
    }
    // The key encoded (a key that is not a str fails there, with an AttributeError) and the
    // value's str() encoded, as classy has always passed them.
    for (const auto& [name, value] : pars_) {
      fc_.set(Bytes(name.attr("encode")()), Bytes(nb::str(value).attr("encode")()));
    }
  }

  // Before cosmology_, so that it outlives it: the InputModule keeps a reference to it.
  FileContent fc_;
  std::unique_ptr<Cosmology> cosmology_;
  // The transfer module's Bessel table, kept across the Cosmology objects this wrapper builds
  // one after another (include/bessel_table_cache.h). A cached table never changes a result.
  std::shared_ptr<BesselTableCache> bessel_cache_;
  nb::dict pars_;
  bool parameters_changed_ = false;
};

// One background quantity at z, by its index in the background vector.
double BackgroundAtZ(PyCosmology& self, double z, int BackgroundModule::* index) {
  const BackgroundModule& module = self.background();
  return module.GetBackgroundAtZ(z)[module.*index];
}

}  // namespace

NB_MODULE(_core, m) {  // NOLINT(cppcoreguidelines-pro-type-cstyle-cast): nanobind's macro
  m.doc() =
      "The compiled core of classy: the Cosmology's lifetime, one method per C++ call, "
      "exception translation and loops over input arrays (nanobind).";

  // The exceptions are plain Python (_errors.py), defined once for the whole package, so that
  // `except classy.CosmoSevereError` catches what this module raises. Imported by the
  // package's runtime name, which is not "classy" in a build renamed with
  // CLASS_PYTHON_MODULE_NAME.
  const std::string name     = nb::cast<std::string>(m.attr("__name__"));
  const std::string package  = name.substr(0, name.rfind('.'));
  nb::module_ errors         = nb::module_::import_((package + "._errors").c_str());
  severe_error               = nb::object(errors.attr("CosmoSevereError")).release().ptr();
  computation_error          = nb::object(errors.attr("CosmoComputationError")).release().ptr();
  m.attr("CosmoError")       = errors.attr("CosmoError");
  m.attr("CosmoSevereError") = errors.attr("CosmoSevereError");
  m.attr("CosmoComputationError") = errors.attr("CosmoComputationError");
  // The package's data (bbn/, hyrec/), next to its __init__.py and to this module.
  nb::object path = nb::module_::import_("os").attr("path");
  nb::object init = nb::module_::import_("sys").attr("modules")[package.c_str()].attr("__file__");
  data_dir        = nb::cast<std::string>(
      path.attr("join")(path.attr("dirname")(path.attr("abspath")(init)), "data"));
  m.attr("_DATA_DIR") = data_dir;
  // The version CMake stamps into version.h from pyproject.toml. Cobaya refuses a classy
  // without it.
  m.attr("__version__") = _VERSION_;

  nb::register_exception_translator(TranslateException);
  // Cosmologies alive at interpreter exit are normal (tests keep them in module globals).
  nb::set_leak_warnings(false);

  m.def(
      "_bessel_direct",
      [](int K, double beta, const IntVector& l, const Vector1d& x) {
        CheckSizes("l", {{1, l.shape(0)}});
        CheckSizes("x", {{1, x.shape(0)}});
        try {
          std::vector<int> lvec(l.data(), l.data() + l.shape(0));
          std::vector<double> phi(l.shape(0) * x.shape(0));
          hyperspherical_bessel_direct_vector(K,
                                              beta,
                                              lvec.data(),
                                              static_cast<int>(l.shape(0)),
                                              x.data(),
                                              static_cast<int>(x.shape(0)),
                                              phi.data());
          return ToArray(std::move(phi), {l.shape(0), x.shape(0)});
        }
        catch (const nb::python_error&) {
          throw;
        }
        catch (const std::exception& e) {
          RaiseAsBuiltin(e);
        }
      },
      "K"_a,
      "beta"_a,
      "l"_a,
      "x"_a,
      R"doc(Phi_l^beta(x) by direct recurrence, as an (n_l, n_x) array (see classy._bessel).)doc");

  m.def(
      "_bessel_interpolate",
      [](int K,
         double beta,
         const IntVector& l_unique,
         const Vector1d& x,
         double xmin,
         double xmax,
         double sampling,
         int l_WKB,
         double phiminabs,
         nb::handle derivatives) {
        // The interpolation table needs an l, and the evaluation an x.
        CheckSizes("l_unique", {{1, l_unique.shape(0)}});
        CheckSizes("x", {{1, x.shape(0)}});
        try {
          const bool with_derivatives = Truthy(derivatives);
          const size_t nl             = l_unique.shape(0);
          const size_t nx             = x.shape(0);
          std::vector<int> lvec(l_unique.data(), l_unique.data() + nl);
          std::vector<double> phi(nl * nx);
          std::vector<double> dphi(with_derivatives ? nl * nx : 0);
          std::vector<double> d2phi(with_derivatives ? nl * nx : 0);
          const HyperInterpStruct his(K,
                                      beta,
                                      static_cast<int>(nl),
                                      lvec.data(),
                                      xmin,
                                      xmax,
                                      sampling,
                                      l_WKB,
                                      phiminabs);
          for (size_t i = 0; i < nl; i++) {
            if (with_derivatives) {
              hyperspherical_Hermite_interpolation<6, true, true, true>(&his,
                                                                        static_cast<int>(nx),
                                                                        static_cast<int>(i),
                                                                        x.data(),
                                                                        phi.data() + i * nx,
                                                                        dphi.data() + i * nx,
                                                                        d2phi.data() + i * nx);
            }
            else {
              hyperspherical_Hermite_interpolation<6, true, false, false>(&his,
                                                                          static_cast<int>(nx),
                                                                          static_cast<int>(i),
                                                                          x.data(),
                                                                          phi.data() + i * nx,
                                                                          nullptr,
                                                                          nullptr);
            }
          }
          nb::object Phi = ToArray(std::move(phi), {nl, nx});
          if (!with_derivatives) {
            return nb::make_tuple(Phi, nb::none(), nb::none());
          }
          return nb::make_tuple(Phi,
                                ToArray(std::move(dphi), {nl, nx}),
                                ToArray(std::move(d2phi), {nl, nx}));
        }
        catch (const nb::python_error&) {
          throw;
        }
        catch (const std::exception& e) {
          RaiseAsBuiltin(e);
        }
      },
      "K"_a,
      "beta"_a,
      "l_unique"_a,
      "x"_a,
      "xmin"_a,
      "xmax"_a,
      "sampling"_a,
      "l_WKB"_a,
      "phiminabs"_a,
      "derivatives"_a.none(),
      "Phi_l^beta(x) on a grid over [xmin, xmax] with Hermite interpolation, one row per l in\n"
      "l_unique (sorted, unique); with derivatives, also dPhi and d2Phi, else None for both.");

  m.def(
      "_bessel_wkb",
      [](int K, double beta, const IntVector& l, const Vector1d& x) {
        std::vector<double> phi(l.shape(0) * x.shape(0));
        for (size_t il = 0; il < l.shape(0); il++) {
          for (size_t ix = 0; ix < x.shape(0); ix++) {
            bool ok = false;
            try {
              ok = hyperspherical_WKB(K, l(il), beta, x(ix), &phi[il * x.shape(0) + ix]);
            }
            catch (const std::exception& e) {
              RaiseAsBuiltin(e);
            }
            if (!ok) {
              PyErr_SetString(severe_error, "hyperspherical_WKB failed");
              throw nb::python_error();
            }
          }
        }
        return ToArray(std::move(phi), {l.shape(0), x.shape(0)});
      },
      "K"_a,
      "beta"_a,
      "l"_a,
      "x"_a,
      R"doc(Phi_l^beta(x) by the WKB (Airy) approximation, as an (n_l, n_x) array.)doc");

  nb::class_<PyCosmology> cls(m,
                              "PyCosmology",
                              "The C++ Cosmology: its parameters, its lifetime and one method per "
                              "C++ call. classy.Class "
                              "subclasses it.");

  cls.def(
         "__init__",
         [](PyCosmology* self, nb::handle input_parameters) {
           new (self) PyCosmology(input_parameters);
         },
         "input_parameters"_a.none() = nb::none())
      .def_prop_ro("pars", &PyCosmology::pars)
      .def_prop_ro("state", [](PyCosmology&) { return true; })
      .def_prop_ro("nonlinear_method",
                   [](PyCosmology& self) { return static_cast<int>(self.nl().method); })
      // struct_cleanup(), empty(), set() and compute() are not necessary, but they are here to
      // support legacy code and MontePython.
      .def("struct_cleanup", [](PyCosmology&) {})
      .def("empty",
           [](PyCosmology& self) {
             self.Empty();
             return nb::find(self);
           })
      .def(
          "set",
          [](PyCosmology& self, nb::handle input_parameters, const nb::kwargs& kwargs) {
            return self.Set(input_parameters, kwargs) ? nb::find(self) : nb::none();
          },
          "input_parameters"_a.none() = nb::none(),
          "kwargs"_a)
      .def(
          "compute",
          [](PyCosmology& self, nb::handle level) {
            // Legacy shim (MontePython): modules are built lazily on first access, so
            // compute() only pre-triggers the requested level. Its remaining value is that
            // computation errors surface here, where callers expect to catch them.
            // The level is compared as Python compares it; an unknown one builds nothing.
            const nb::object final_level = level.is_none() ? nb::object(nb::str("lensing"))
                                                           : level[0].attr("lower")();
            using Build                  = void (*)(Cosmology&);
            static const std::pair<const char*, Build> levels[] =
                {{"background", [](Cosmology& c) { c.GetBackgroundModule(); }},
                 {"thermodynamics", [](Cosmology& c) { c.GetThermodynamicsModule(); }},
                 {"perturb", [](Cosmology& c) { c.GetPerturbationsModule(); }},
                 {"primordial", [](Cosmology& c) { c.GetPrimordialModule(); }},
                 {"nonlinear", [](Cosmology& c) { c.GetNonlinearModule(); }},
                 {"transfer", [](Cosmology& c) { c.GetTransferModule(); }},
                 {"spectra", [](Cosmology& c) { c.GetSpectraModule(); }},
                 {"lensing", [](Cosmology& c) { c.GetLensingModule(); }}};
            for (const auto& [name, build] : levels) {
              if (Equal(final_level, nb::str(name))) {
                build(self.cosmo());
                break;
              }
            }
            return nb::find(self);
          },
          "level"_a.none() = nb::none())
      .def(
          "raw_cl",
          [](PyCosmology& self, nb::handle lmax) {
            return self.GeneralCl(IntArg<int>(lmax), false);
          },
          "lmax"_a = -1)
      .def(
          "lensed_cl",
          [](PyCosmology& self, nb::handle lmax) {
            return self.GeneralCl(IntArg<int>(lmax), true);
          },
          "lmax"_a = -1)
      .def(
          "raw_cl_no_copy",
          [](PyCosmology& self, nb::handle lmax) { return self.RawClNoCopy(IntArg<int>(lmax)); },
          "lmax"_a = -1)
      .def("lensed_cl_computed",
           [](PyCosmology& self) {
             nb::dict out = ClDict(self.cosmo().GetLensingModule()->cl_output_computed());
             out["ell"]   = out["ell"].attr("astype")(Builtin("int"));
             return out;
           })
      .def(
          "z_of_r",
          [](PyCosmology& self, const Vector1d& z_array) {
            const BackgroundModule& module = self.background();
            std::vector<double> r(z_array.shape(0));
            std::vector<double> dzdr(z_array.shape(0));
            for (size_t i = 0; i < z_array.shape(0); i++) {
              const std::vector<double> pvecback = module.GetBackgroundAtZ(z_array(i));
              r[i]                               = pvecback[module.index_bg_conf_distance_];
              // dz/dr = H
              dzdr[i] = pvecback[module.index_bg_H_];
            }
            return nb::make_tuple(ToArray(std::move(r)), ToArray(std::move(dzdr)));
          },
          "z_array"_a)
      .def(
          "luminosity_distance",
          [](PyCosmology& self, double z) {
            return BackgroundAtZ(self, z, &BackgroundModule::index_bg_lum_distance_);
          },
          "z"_a,
          R"doc(luminosity_distance(z))doc")
      .def(
          "pk",
          [](PyCosmology& self, double k, double z) {
            return self.Pk(k, z, self.nonlinear_module().index_pk_m_, self.DefaultPkOutput());
          },
          "k"_a,
          "z"_a,
          R"doc(Gives the total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z (will be non linear if requested to Class, linear otherwise)

.. note::

    there is an additional check that output contains `mPk`,
    because otherwise a segfault will occur)doc")
      .def(
          "pk_cb",
          [](PyCosmology& self, double k, double z) {
            return self.Pk(k, z, self.nonlinear_module().index_pk_cluster_, self.DefaultPkOutput());
          },
          "k"_a,
          "z"_a,
          R"doc(Gives the cdm+b pk (in Mpc**3) for a given k (in 1/Mpc) and z (will be non linear if requested to Class, linear otherwise)

.. note::

    there is an additional check that output contains `mPk`,
    because otherwise a segfault will occur)doc")
      .def(
          "pk_lin",
          [](PyCosmology& self, double k, double z) {
            return self.Pk(k, z, self.nonlinear_module().index_pk_m_, pk_linear);
          },
          "k"_a,
          "z"_a,
          R"doc(Gives the linear total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z

.. note::

    there is an additional check that output contains `mPk`,
    because otherwise a segfault will occur)doc")
      .def(
          "pk_cb_lin",
          [](PyCosmology& self, double k, double z) {
            return self.Pk(k, z, self.nonlinear_module().index_pk_cluster_, pk_linear);
          },
          "k"_a,
          "z"_a,
          R"doc(Gives the linear cdm+b pk (in Mpc**3) for a given k (in 1/Mpc) and z

.. note::

    there is an additional check that output contains `mPk`,
    because otherwise a segfault will occur)doc")
      .def(
          "get_pk",
          [](PyCosmology& self,
             const Vector3d& k,
             const Vector1d& z,
             nb::handle k_size,
             nb::handle z_size,
             nb::handle mu_size) {
            return self.GetPk(k,
                              z,
                              k_size,
                              z_size,
                              mu_size,
                              self.nonlinear_module().index_pk_m_,
                              self.DefaultPkOutput());
          },
          "k"_a,
          "z"_a,
          "k_size"_a,
          "z_size"_a,
          "mu_size"_a,
          R"doc(Fast function to get the power spectrum on a k and z array )doc")
      .def(
          "get_pk_cb",
          [](PyCosmology& self,
             const Vector3d& k,
             const Vector1d& z,
             nb::handle k_size,
             nb::handle z_size,
             nb::handle mu_size) {
            return self.GetPk(k,
                              z,
                              k_size,
                              z_size,
                              mu_size,
                              self.nonlinear_module().index_pk_cluster_,
                              self.DefaultPkOutput());
          },
          "k"_a,
          "z"_a,
          "k_size"_a,
          "z_size"_a,
          "mu_size"_a,
          R"doc(Fast function to get the power spectrum on a k and z array )doc")
      .def(
          "get_pk_lin",
          [](PyCosmology& self,
             const Vector3d& k,
             const Vector1d& z,
             nb::handle k_size,
             nb::handle z_size,
             nb::handle mu_size) {
            return self.GetPk(k,
                              z,
                              k_size,
                              z_size,
                              mu_size,
                              self.nonlinear_module().index_pk_m_,
                              pk_linear);
          },
          "k"_a,
          "z"_a,
          "k_size"_a,
          "z_size"_a,
          "mu_size"_a,
          R"doc(Fast function to get the linear power spectrum on a k and z array )doc")
      .def(
          "get_pk_cb_lin",
          [](PyCosmology& self,
             const Vector3d& k,
             const Vector1d& z,
             nb::handle k_size,
             nb::handle z_size,
             nb::handle mu_size) {
            return self.GetPk(k,
                              z,
                              k_size,
                              z_size,
                              mu_size,
                              self.nonlinear_module().index_pk_cluster_,
                              pk_linear);
          },
          "k"_a,
          "z"_a,
          "k_size"_a,
          "z_size"_a,
          "mu_size"_a,
          R"doc(Fast function to get the linear power spectrum on a k and z array )doc")
      .def(
          "_z_grid_for_output",
          [](PyCosmology& self) { return ToArray(self.perturbations().GetLateTimeRedshifts()); },
          R"doc(The redshifts P(k,z) and T_i(k,z) are tabulated at, decreasing to 0.

One owner for this axis. Both tables live on the perturbation module's
ln_tau grid -- the Fourier module holds a copy -- so anything that
combines them, as get_Weyl_pk_and_k_and_z does, is aligned by
construction rather than by coincidence.)doc")
      .def(
          "_pk_grid",
          [](PyCosmology& self, nb::handle nonlinear, nb::handle only_clustering_species) {
            const NonlinearModule& module = self.nonlinear_module();
            std::vector<double> k;
            std::vector<double> pk;
            module.GetPkGrid(Truthy(nonlinear) ? pk_nonlinear : pk_linear,
                             Truthy(only_clustering_species) ? module.index_pk_cluster_
                                                             : module.index_pk_m_,
                             k,
                             pk);
            return nb::make_tuple(ToArray(std::move(pk)), ToArray(std::move(k)));
          },
          "nonlinear"_a.none(),
          "only_clustering_species"_a.none(),
          R"doc(The stored P(k, z) grid, flat and z-major, and its k in 1/Mpc.)doc")
      .def(
          "_scalar_transfer_k",
          [](PyCosmology& self) { return ToArray(self.perturbations().GetScalarTransferK()); },
          R"doc(The k axis of the scalar transfer grids, in 1/Mpc; refuses runs it does not fit.)doc")
      .def(
          "_transfer_at_node",
          [](PyCosmology& self, int index_tau, nb::handle output_format) {
            return self.TransferGeneral(0., index_tau, output_format);
          },
          "index_tau"_a,
          "output_format"_a.none(),
          R"doc(The transfer functions at a node of the stored late-time sampling.)doc")
      .def(
          "_sigma",
          [](PyCosmology& self, double R, double z, nb::handle cb) {
            const NonlinearModule& module = self.nonlinear_module();
            return module.GetSigma(R,
                                   z,
                                   Truthy(cb) ? module.index_pk_cluster_ : module.index_pk_m_);
          },
          "R"_a,
          "z"_a,
          "cb"_a.none(),
          R"doc(sigma(R, z) of all matter, or of cdm + baryons with cb; R in Mpc.)doc")
      .def(
          "pk_tilt",
          [](PyCosmology& self, double k, double z) {
            const auto module = self.cosmo().GetNonlinearModule();
            double pk_tilt    = 0.;
            module->nonlinear_pk_tilt_at_k_and_z(pk_linear, k, z, module->index_pk_m_, &pk_tilt);
            return pk_tilt;
          },
          "k"_a,
          "z"_a,
          R"doc(Gives effective logarithmic slope of P_L(k,z) (total matter) for a given k and z
(k is the wavenumber in units of 1/Mpc, z is the redshift, the output is dimensionless)

.. note::

    there is an additional check to verify whether output contains `mPk` and whether k is in the right range)doc")
      .def(
          "angular_distance",
          [](PyCosmology& self, double z) {
            return BackgroundAtZ(self, z, &BackgroundModule::index_bg_ang_distance_);
          },
          "z"_a,
          R"doc(angular_distance(z)

Return the angular diameter distance (exactly, the quantity defined by Class
as index_bg_ang_distance in the background module)

Parameters
----------
z : float
        Desired redshift)doc")
      .def(
          "angular_distance_from_to",
          [](PyCosmology& self, double z1, double z2) {
            return self.background().GetAngularDistanceFromTo(z1, z2);
          },
          "z1"_a,
          "z2"_a,
          R"doc(Angular diameter distance of an object at z2 as seen by an observer at z1,
that is, sin_K((chi2 - chi1)*sqrt(|K|))/sqrt(|K|)/(1 + z2). Zero if z1 >= z2.

Parameters
----------
z1 : float
        Observer redshift
z2 : float
        Source redshift

Returns
-------
d_A(z1, z2) in Mpc)doc")
      .def(
          "scale_independent_growth_factor",
          [](PyCosmology& self, double z) {
            return BackgroundAtZ(self, z, &BackgroundModule::index_bg_D_);
          },
          "z"_a,
          R"doc(scale_independent_growth_factor(z)

Return the scale invariant growth factor D(a) for CDM perturbations
(exactly, the quantity defined by Class as index_bg_D in the background module)

Parameters
----------
z : float
        Desired redshift)doc")
      .def(
          "scale_independent_growth_factor_f",
          [](PyCosmology& self, double z) {
            return BackgroundAtZ(self, z, &BackgroundModule::index_bg_f_);
          },
          "z"_a,
          R"doc(scale_independent_growth_factor_f(z)

Return the scale invariant growth factor f(z)=d ln D / d ln a for CDM perturbations
(exactly, the quantity defined by Class as index_bg_f in the background module)

Parameters
----------
z : float
        Desired redshift)doc")
      .def(
          "z_of_tau",
          [](PyCosmology& self, double tau) { return self.background().GetZOfTau(tau); },
          "tau"_a,
          R"doc(Redshift corresponding to a given conformal time.

Parameters
----------
tau : float
        Conformal time)doc")
      .def(
          "Hubble",
          [](PyCosmology& self, double z) {
            return BackgroundAtZ(self, z, &BackgroundModule::index_bg_H_);
          },
          "z"_a,
          R"doc(Hubble(z)

Return the Hubble rate (exactly, the quantity defined by Class as index_bg_H
in the background module)

Parameters
----------
z : float
        Desired redshift)doc")
      .def(
          "Om_m",
          [](PyCosmology& self, double z) {
            return BackgroundAtZ(self, z, &BackgroundModule::index_bg_Omega_m_);
          },
          "z"_a,
          R"doc(Omega_m(z)

Return the matter density fraction (exactly, the quantity defined by Class as index_bg_Omega_m
in the background module)

Parameters
----------
z : float
        Desired redshift)doc")
      .def(
          "ionization_fraction",
          [](PyCosmology& self, double z) {
            const ThermodynamicsModule& module = self.thermo();
            return module.GetThermodynamicsAtZ(z)[module.index_th_xe_];
          },
          "z"_a,
          R"doc(ionization_fraction(z)

Return the ionization fraction for a given redshift z

Parameters
----------
z : float
        Desired redshift)doc")
      .def(
          "baryon_temperature",
          [](PyCosmology& self, double z) {
            const ThermodynamicsModule& module = self.thermo();
            return module.GetThermodynamicsAtZ(z)[module.index_th_Tb_];
          },
          "z"_a,
          R"doc(baryon_temperature(z)

Give the baryon temperature for a given redshift z

Parameters
----------
z : float
        Desired redshift)doc")
      .def(
          "get_background",
          [](PyCosmology& self) { return TableDict(self.background().GetOutputTable()); },
          R"doc(Return an array of the background quantities at all times.

Parameters
----------

Returns
-------
background : dictionary containing background.)doc")
      .def(
          "get_thermodynamics",
          [](PyCosmology& self) { return TableDict(self.thermo().GetOutputTable()); },
          R"doc(Return the thermodynamics quantities.

Returns
-------
thermodynamics : dictionary containing thermodynamics.)doc")
      .def(
          "get_primordial",
          [](PyCosmology& self) { return TableDict(self.primordial().GetOutputTable()); },
          R"doc(Return the primordial scalar and/or tensor spectrum depending on 'modes'.
'output' must be set to something, e.g. 'tCl'.

Returns
-------
primordial : dictionary containing k-vector and primordial scalar and tensor P(k).)doc")
      .def(
          "get_perturbations",
          [](PyCosmology& self) { return self.GetPerturbations(); },
          R"doc(Return scalar, vector and/or tensor perturbations as arrays for requested
k-values.

.. note::

    you need to specify both 'k_output_values', and have some
    perturbations computed, for instance by setting 'output' to 'tCl'.

Returns
-------
perturbations : dict of array of dicts
        perturbations['scalar'] is an array of length 'k_output_values' of
        dictionary containing scalar perturbations.
        Similar for perturbations['vector'] and perturbations['tensor'].)doc")
      .def(
          "get_transfer",
          [](PyCosmology& self, double z, nb::handle output_format) {
            return self.TransferGeneral(z, -1, output_format);
          },
          "z"_a                    = 0.,
          "output_format"_a.none() = "class",
          R"doc(Return the density and/or velocity transfer functions for all initial
conditions today. You must include 'mTk' and/or 'vCTk' in the list of
'output'. The transfer functions can also be computed at higher redshift z
provided that 'z_pk' has been set and that 0<z<z_pk.

Parameters
----------
z  : redshift (default = 0)
output_format  : ('class' or 'camb') Format transfer functions according to
                 CLASS convention (default) or CAMB convention.

Returns
-------
tk : dictionary containing transfer functions.)doc")
      .def(
          "_derived_parameter",
          [](PyCosmology& self, const std::string& name) {
            return self.cosmo().GetDerivedParameter(name);
          },
          "name"_a,
          R"doc(One derived parameter, by name (Cosmology::GetDerivedParameter).)doc")
      .def(
          "nonlinear_scale",
          [](PyCosmology& self, const Vector1d& z, nb::handle z_size) {
            return self.NonlinearScales(z, z_size, false);
          },
          "z"_a,
          "z_size"_a,
          R"doc(nonlinear_scale(z, z_size)

Return the nonlinear scale for all the redshift specified in z, of size
z_size

Parameters
----------
z : numpy array
        Array of requested redshifts
z_size : int
        Size of the redshift array)doc")
      .def(
          "nonlinear_scale_cb",
          [](PyCosmology& self, const Vector1d& z, nb::handle z_size) {
            return self.NonlinearScales(z, z_size, true);
          },
          "z"_a,
          "z_size"_a,
          R"doc(nonlinear_scale_cb(z, z_size)

Return the nonlinear scale of the cdm+baryon spectrum for all the
redshift specified in z, of size z_size

Parameters
----------
z : numpy array
        Array of requested redshifts
z_size : int
        Size of the redshift array)doc")
      .def(
          "get_pk_array",
          [](PyCosmology& self,
             const Vector1d& k,
             const Vector1d& z,
             nb::handle k_size,
             nb::handle z_size,
             nb::handle nonlinear) {
            return self.GetPkArrays(k, z, k_size, z_size, nonlinear).first;
          },
          "k"_a,
          "z"_a,
          "k_size"_a,
          "z_size"_a,
          "nonlinear"_a.none())
      .def(
          "get_pk_cb_array",
          [](PyCosmology& self,
             const Vector1d& k,
             const Vector1d& z,
             nb::handle k_size,
             nb::handle z_size,
             nb::handle nonlinear) {
            return self.GetPkArrays(k, z, k_size, z_size, nonlinear).second;
          },
          "k"_a,
          "z"_a,
          "k_size"_a,
          "z_size"_a,
          "nonlinear"_a.none());
}
