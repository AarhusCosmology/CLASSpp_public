# CLASS: Cosmic Linear Anisotropy Solving System

Main authors: Julien Lesgourgues, Thomas Tram, Nils Schoeneberg

with several major inputs from other people, especially Benjamin
Audren, Simon Prunet, Deanna Hooper, Maria Archidiacono, Karim
Benabed, Thejs Brinckmann, Marius Millea, Simeon
Bird, Samuel Brieden, Jesus Torrado, Miguel Zumalacarregui, Francesco
Montanari, Jeppe Dakin, Emil Brinch Holm, etc.

For download and information, see [class-code.net](http://class-code.net),
[original Github](https://github.com/lesgourg/class_public) and [this fork](https://github.com/AarhusCosmology/CLASSpp).

## About this fork

In March 2020 a branch of the development version of CLASS v2.9 was started in order to start converting CLASS into C++. The reasons for doing this were many and are summarised below:

- C++ is as fast as C, but it has many features that enables one to do more with fewer lines of code.
- CLASS is often used through `classy`, a Python wrapper written in Cython. C++ allows a more seamless wrapping where e.g. exceptions can be thrown in C++ code and automatically handled in Cython/Python.
- C++ allows automatic (and safe) allocation and deallocation of resources (memory) using constructors and destructors which can eliminate the possibility of memory leaks.
- With C, there was always the possibility that a function was called while the corresponding module is in an incomplete state, often leading to segmentation fault. In CLASS++, if construction of a module fails, it is not possible to call its member methods.

## Why not Rust, Julia, Go, or something else?

The primary reason for C++ is that it allowed/allows the conversion to be gradual. At the beginning, a GitHub Actions workflow was established that tested the observables computed in many models against the main branch before allowing a PR to merge. If a new code had been written, it would have required a more thorough validation in the end.

That said, the benefits of Rust compared to C++ are primarily in multithreaded code where race-conditions can be discovered at compile-time. An Einstein--Boltzmann code in Julia would be interesting, but it might be dificult to make it sufficiently fast using just-in-time compilation.

Regarding other languages, anyone are encouraged to write a new Einstein--Boltzmann solver!

## Installation

To install CLASS for use exclusively through the Python wrapper `classy`, use one of the following options:

1. **Installation from PyPI:**  
```
pip install classy-community
```
2. **Installation from GitHub:**  
```
pip install 'classy-community @ git+https://github.com/AarhusCosmology/CLASSpp_public.git'
```  
or  
```
pip install 'classy-community @ git+https://github.com/AarhusCosmology/CLASSpp_public.git@5d8b27678fae945080080cdaeb3a40857fa49deb'
```  
for a specific version.
3. **Installation from sources:**  
```
git clone https://github.com/AarhusCosmology/CLASSpp_public
pip install .
```

On MacOS and Anaconda, the path to the archiver `ar` can be wrong. *If you get an error*, do

```bash
export AR=/usr/bin/ar
```
and try again.

## Installing the command line application:
If you want to run CLASS from the command line, you need to clone the repository and build the executable using one of the following two options:

1. **Build using GNU Make:**  
```
make -j
```  
Optionally modify the Makefile to specify compiler flags etc.
2. **Build using Xcode command-line tools on MacOS:**  
```
xcodebuild -configuration Release SYMROOT=$PWD
```  
or  
```
xcodebuild -configuration Debug SYMROOT=$PWD
```  
for the debug build. The executable `class` will be stored in `./Release/` and `./Debug/` respectively. Consider copying/moving it to the project root.
3. **Build using Xcode GUI:**  
The code can be build, profiled and debugged using Xcode by opening the project file.


To check that the code runs, type:

    ./class explanatory.ini

The explanatory.ini file is THE reference input file, containing and
explaining the use of all possible input parameters. We recommend to
read it, to keep it unchanged (for future reference), and to create
for your own purposes some shorter input files, containing only the
input lines which are useful for you. Input files must have a *.ini
extension.

If you want to play with the precision/speed of the code, you can use
one of the provided precision files (e.g. cl_permille.pre) or modify
one of them, and run with two input files, for instance:

    ./class test.ini cl_permille.pre

The files *.pre are suppposed to specify the precision parameters for
which you don't want to keep default values. If you find it more
convenient, you can pass these precision parameter values in your *.ini
file instead of an additional *.pre file.

The automatically-generated documentation is located in

    doc/manual/html/index.html
    doc/manual/CLASS_manual.pdf

On top of that, if you wish to modify the code, you will find lots of
comments directly in the files.

### Plotting utility

There are two outdated plotting scripts that may however still be useful from time to time.One os for Python, CPU.py, (Class Plotting Utility), written by Benjamin Audren and
Jesus Torrado. It can plot the Cl's, the P(k) or any other CLASS
output, for one or several models, as well as their ratio or percentage
difference. The syntax and list of available options is obtained by
typing 'pyhton CPU.py -h'. The other is for MATLAB. To use it, once in MATLAB, type 'help
plot_CLASS_output.m'

## Bug tracking and contributing

If you discover a bug or have a nice idea for a feature that you would like to see, open an issue at 
[https://github.com/AarhusCosmology/CLASSpp_public/issues](https://github.com/AarhusCosmology/CLASSpp_public/issues). *A bug report should provide clear instructions on how to reproduce the bug*.

If you want to contribute to CLASS, make a fork of (https://github.com/AarhusCosmology/CLASSpp_public)[https://github.com/AarhusCosmology/CLASSpp_public]. Before submitting the PR, make sure that an **issue** exists which clearly explains the bug or feature that the PR aims to solve.

## Using and citing the code

You can use CLASS freely, but it is expected that in your publications, you cite
at least the paper `CLASS II: Approximation schemes
<http://arxiv.org/abs/1104.2933>`. (See COPYING, LICENSE and AUTHORS
for formal details.) Feel free to cite more CLASS papers!

```
@article{Blas:2011rf,
    author = "Blas, Diego and Lesgourgues, Julien and Tram, Thomas",
    title = "{The Cosmic Linear Anisotropy Solving System (CLASS) II: Approximation schemes}",
    eprint = "1104.2933",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    reportNumber = "CERN-PH-TH-2011-082, LAPTH-010-11",
    doi = "10.1088/1475-7516/2011/07/034",
    journal = "JCAP",
    volume = "07",
    pages = "034",
    year = "2011"
}
```

```
@article{Lesgourgues:2011rh,
    author = "Lesgourgues, Julien and Tram, Thomas",
    title = "{The Cosmic Linear Anisotropy Solving System (CLASS) IV: efficient implementation of non-cold relics}",
    eprint = "1104.2935",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    reportNumber = "CERN-PH-TH-2011-084, LAPTH-012-11",
    doi = "10.1088/1475-7516/2011/09/032",
    journal = "JCAP",
    volume = "09",
    pages = "032",
    year = "2011"
}
```

## Naming conventions and code-style

We never explicitly defined a code-style for CLASS, but there is a set of implicit rules that we have mostly adhered to:

- Indentation uses two spaces, no tabs.
- All variable names are typeset using snake_case.
- Defined constants `_EXAMPLE_` use capital letters and has a starting and a trailing underscore.
- There is empty space around comparison operators and binary operators `+` and `-`.
- There is *no* empty space around binary operators `*` and `/`.
- There is empty space after every comma `,`.

The last three rules are not followed in all of the code-base, but they should be changed whenever a modification is neccessary. If-else statements are indented like this:

```cpp
if (a == 0) {
  t = 1 + 1.0/2.0;
}
else {
  t = pow(1, a); 
}
```

We will need to add a few more rules to support classes and objects:

- Names of classes, structs and enums are written in PascalCase.
- Methods are also written in PascalCase.
- Struct data members are written in snake_case.
- Class data members are written in snake_case without a trailing underscore.
- Constant data members (declared with the const specifier) are written as `kDaysInAWeek`. 
- Static, constant data members (declared static const ) are written as `sDaysInAWeek`.

These choices are inspired by Google's [C++ style guide](https://google.github.io/styleguide/cppguide.html). The trailing `_` on class data members is very important due to shadowing. Consider the following code:

```cpp
class ClModule() {
  int lmax = 2500; 
  int GetEllMax()
}

int ClModule::GetEllMax() {
  int lmax = 2;
  return lmax;
}
```
In this case `GetEllMax()` returns `2` since the local `lmax` variable *shadows* the class member variable `lmax` . This can lead to unexpected behaviour, especially in large classes. A trailing underscore effectively protects against shadowing, and together with the syntax highlighting of the editor makes class variables easy to recognise in the code.

For structs, this is not a problem, since they should atmost contain a simple constructor and destructor. (Otherwise, the struct in question should be converted into a class.)

## Scoped variables

One programming style which should be completely abandoned throughout the code-base is declaring a variable before it is needed and without initialising it. Since C99, this was not even neccessary to do in pure C, so that should have been avoided from the start. Here are a few examples trying to accomplish the same thing:

```cpp
// Worst possible way!
double a_prime_over_a = 0.;
...
if (a > 0.5) {
  printf("a_prime_over_a = %.3e\n", a_prime_over_a);
}
```
In this case we have a variable `a_prime_over_a` which is initialised to the wrong value! We **hope** that somewhere during ... , we correct the value, but this is not very robust: perhaps we only overwrite it under certain conditions.

```cpp
// Still not good!
double a_prime_over_a;
...
if (a > 0.5) {
  printf("a_prime_over_a = %.3e\n", a_prime_over_a);
}
```
In this example we are at least not initialising the variable to a **wrong** value. This means that we get a compiler warning (which we are grateful for, because it tells us that we are doing something stupid!) if we use it before we initialise it.

```cpp
// Getting better!
double a_prime_over_a = a*H;
...
if (a > 0.5) {
  printf("a_prime_over_a = %.3e\n", a_prime_over_a);
}
```
Now `a_prime_over_a` has the correct value for its complete life-time. I cannot compile the code if I am using it before it is initialised, which is good.

```cpp
// Best!
...
if (a > 0.5) {
  double a_prime_over_a = a*H;
  printf("a_prime_over_a = %.3e\n", a_prime_over_a);
}
```

This is the best solution: the life-time of the `a_prime_over_a` object is tied to the scope of the if-statement, so it cannot be changed from outside the if- statement. Furthermore, the code is much more readable, since the reader immediately sees the definition of `a_prime_over_a`. Thus, we shall adhere to the following rules:

- Variables should always be initialised when they are declared, except if it is not possible.
- Variables should be declared in as narrow a scope as possible.

A typical example of this is the counter inside if-statements:

```cpp
int i = 123; // Not changed by loop
for (int i = 0; i < 10; ++i) {
  ... 
}
```

This way, the variable `i` is local to the for-loop, so it does not change the value of the variable `i` defined outside the loop.

## The BaseModule class

There are some things that all modules (not including `InputModule`) have in common, and this is collected in the `BaseModule` class which all other modules inherit from. A `BaseModule` is initialised with an `InputModulePtr`, and upon construction it sets the usual pointers `ppr`, `pba`, etc. to point towards the "input structs" that are stored inside the `InputModule`.

```cpp
class BaseModule {
public:
  BaseModule(InputModulePtr input_module)
  : ncdm_(input_module->ncdm_)
  , dr_(input_module->dr_)
  , ppr(&input_module->precision_)
  , pba(&input_module->background_)
  , pth(&input_module->thermodynamics_)
  , ppt(&input_module->perturbations_)
  , ppm(&input_module->primordial_)
  , pnl(&input_module->nonlinear_)
  , ptr(&input_module->transfers_)
  , psp(&input_module->spectra_)
  , ple(&input_module->lensing_)
  , pop(&input_module->output_) {
    input_module_ = std::move(input_module);
    error_message_[0] = '\n';
  }
  BaseModule(const BaseModule&) = delete;
  
  mutable ErrorMsg error_message_;
public:
  const std::shared_ptr<NonColdDarkMatter> ncdm_;
  const std::shared_ptr<DarkRadiation> dr_;
protected:
  InputModulePtr input_module_;

  const precision* const ppr;
  const background* const pba;
  const thermo* const pth;
  const perturbs* const ppt;
  const primordial* const ppm;
  const nonlinear* const pnl;
  const transfers* const ptr;
  const spectra* const psp;
  const lensing* const ple;
  const output* const pop;
};
```

There are a few things to note here: First, the `BaseModule` contains a shared pointer to the `InputModule` object, so all pointers are guaranteed to be valid for the life-time of the `BaseModule`. Second, the input pointers are declared with two `const` keywords. This means that the pointer address cannot be changed, and the content that the pointer points toward can also not be changed.

Once a module inherits from `BaseModule`, it inherits all these pointers to input structs, (`pba`, `pth` etc.). Thus, whenever we are inside a module method, we can simply write `pba->` and `error_mesage_`, and we never need to pass those around as input arguments. This is a massive advantage compared to vanilla CLASS, and it simplifies writing new methods a lot. Note that we did not adhere to the `_` notation for these members. The reason is that we found the names sufficiently cannonical to warrant the omission.

## The Cosmology class

The full CLASS has become a `Cosmology` class. The main function reads:

```cpp
int main(int argc, char **argv) {

  FileContent fc;
  ErrorMsg error_message;
  if (InputModule::file_content_from_arguments(argc, argv, fc, error_message) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n", error_message);
    return _FAILURE_;
  }

  Cosmology cosmology{fc};

  OutputModule output_module(cosmology.GetInputModule(),
                             cosmology.GetBackgroundModule(),
                             cosmology.GetThermodynamicsModule(),
                             cosmology.GetPerturbationsModule(),
                             cosmology.GetPrimordialModule(),
                             cosmology.GetNonlinearModule(),
                             cosmology.GetSpectraModule(),
                             cosmology.GetLensingModule());

  return _SUCCESS_;
```

The `Cosmology` class itself is defined as this:

```cpp
class Cosmology {
public:
  Cosmology(FileContent& fc)
  : input_module_ptr_(InputModulePtr(new InputModule(fc))) {}
  Cosmology(std::unique_ptr<InputModule> input_module)
  : input_module_ptr_(std::move(input_module)) {}

  InputModulePtr& GetInputModule();
  BackgroundModulePtr& GetBackgroundModule();
  ThermodynamicsModulePtr& GetThermodynamicsModule();
  PerturbationsModulePtr& GetPerturbationsModule();
  PrimordialModulePtr& GetPrimordialModule();
  NonlinearModulePtr& GetNonlinearModule();
  TransferModulePtr& GetTransferModule();
  SpectraModulePtr& GetSpectraModule();
  LensingModulePtr& GetLensingModule();

private:
  InputModulePtr input_module_ptr_;
  BackgroundModulePtr background_module_ptr_;
  ThermodynamicsModulePtr thermodynamics_module_ptr_;
  PerturbationsModulePtr perturbations_module_ptr_;
  PrimordialModulePtr primordial_module_ptr_;
  NonlinearModulePtr nonlinear_module_ptr_;
  TransferModulePtr transfer_module_ptr_;
  SpectraModulePtr spectra_module_ptr_;
  LensingModulePtr lensing_module_ptr_;
};
```

There are two ways to construct a `Cosmology` class:

1. By passing a reference to a `FileContent` struct.
2. By passing a unique pointer to an input module.

The second way is useful when the input structure needs to be modified manually, which happens when we want to use non-default precision parameters during shooting.