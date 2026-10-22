#ifndef CLASS_FUNCTION_REF_H
#define CLASS_FUNCTION_REF_H

/** @file function_ref.h A non-owning reference to a callable: C++26's std::function_ref,
 *  until the code base moves to C++26.
 *
 *  Design, measurements and the alternatives rejected (std::function, templated evolvers):
 *  docs/superpowers/specs/2026-10-22-function-ref-callbacks-design.md
 *
 *  A FunctionRef is the (function pointer, void* context) pair that the C interfaces passed
 *  by hand, with the compiler writing the forwarding function. A call costs one indirect
 *  call, as before, and the callable's body is inlined into the forwarder:
 *
 *      auto derivs = [&](double tau, double* y, double* dy) { perturb_derivs(tau, y, dy, &ppaw); };
 *      evolver_ndf15(derivs, tau_ini, tau_end, y, y_size, options);
 *
 *  It refers to the callable and does not copy it, so the callable must outlive every call
 *  through the FunctionRef. It is never null: a callback that may be absent is a
 *  std::optional<FunctionRef<...>>, as it would be a std::optional<std::function_ref<...>>.
 *
 *  The interface is the part of std::function_ref the code base uses, spelled the same way,
 *  so moving to C++26 is a rename. Every program that compiles against FunctionRef compiles
 *  against std::function_ref and does the same. The one difference is deliberate:
 *  FunctionRef binds only lvalues. std::function_ref also binds a temporary, which is safe
 *  as a function argument but dangles once stored, and options structs store callbacks:
 *
 *      options.output = [&](double x, double* y, double* dy, int i) { ... };  // dangles
 *
 *  Here that line does not compile. Name the lambda first. The guard goes with the rename.
 */

#include <memory>
#include <type_traits>
#include <utility>

template <class Signature>
class FunctionRef;

template <class R, class... Args>
class FunctionRef<R(Args...)> {
 public:
  /** A free function, which is its own context. Must not be null. */
  template <
      class F,
      std::enable_if_t<std::is_function_v<F> && std::is_invocable_r_v<R, F&, Args...>, int> = 0>
  FunctionRef(F* function) noexcept : bound_(function), forward_(&ForwardFunction<F>) {}

  /** Any other callable, by reference. Lvalues only: see the file comment. */
  template <
      class F,
      std::enable_if_t<!std::is_same_v<std::remove_cv_t<F>, FunctionRef> &&
                           !std::is_function_v<F> && !std::is_pointer_v<F> &&
                           !std::is_member_pointer_v<F> && std::is_invocable_r_v<R, F&, Args...>,
                       int> = 0>
  FunctionRef(F& callable) noexcept
      : bound_(std::addressof(callable)), forward_(&ForwardObject<F>) {}

  FunctionRef(const FunctionRef&) noexcept = default;
  // The defaulted copy assigns the union. NOLINTNEXTLINE(cppcoreguidelines-pro-type-union-access)
  FunctionRef& operator=(const FunctionRef&) noexcept = default;

  /** As in std::function_ref, assigning a callable is an error, since nothing would keep
   *  it alive; assign a FunctionRef or a function pointer instead. */
  template <class T,
            std::enable_if_t<!std::is_same_v<T, FunctionRef> && !std::is_pointer_v<T>, int> = 0>
  FunctionRef& operator=(T) = delete;

  R operator()(Args... args) const {
    return forward_(bound_, std::forward<Args>(args)...);
  }

 private:
  /* The type erasure itself: a union holding what was erased, and the casts that restore
     it. Converting a function pointer to void* is only conditionally supported, hence the
     union; casting a function pointer to another function pointer type and back is well
     defined. */
  // NOLINTBEGIN(cppcoreguidelines-pro-type-union-access,cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast)
  union Bound {
    template <class F>
    explicit Bound(F* callable) noexcept : object(const_cast<std::remove_cv_t<F>*>(callable)) {}
    template <class F, class... A>
    explicit Bound(F (*fn)(A...)) noexcept : function(reinterpret_cast<void (*)()>(fn)) {}

    void* object;
    void (*function)();
  };

  template <class F>
  static R ForwardObject(Bound bound, Args... args) {
    F& callable = *static_cast<F*>(bound.object);
    if constexpr (std::is_void_v<R>) {
      callable(std::forward<Args>(args)...);
    }
    else {
      return callable(std::forward<Args>(args)...);
    }
  }

  template <class F>
  static R ForwardFunction(Bound bound, Args... args) {
    F* function = reinterpret_cast<F*>(bound.function);
    if constexpr (std::is_void_v<R>) {
      function(std::forward<Args>(args)...);
    }
    else {
      return function(std::forward<Args>(args)...);
    }
  }
  // NOLINTEND(cppcoreguidelines-pro-type-union-access,cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast)

  Bound bound_;
  R (*forward_)(Bound, Args...);
};

#endif  // CLASS_FUNCTION_REF_H
