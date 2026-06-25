#pragma once

#include <algorithm>
#include <cassert>
#include <memory>
#include <string>
#include <vector>

#include "base_species.h"

/**
 * Ordered, named collection of cosmological species.
 *
 * Mutation is restricted to the construction phase (InputModule::ConstructSpecies).
 * Each species is inserted with an explicit lookup key — today this is the
 * hard-coded class identifier (e.g. "Photons", "CDM", "IDM_DR_IDR"); under the
 * object-based input syntax planned in issue #246 it will be the user-supplied
 * instance name (e.g. "nu1", "ncdm__1"). The key is NOT derived from
 * BaseSpecies::name() because post-#246 several instances may share the same
 * type name (e.g. two "massive_ncdm" instances).
 *
 * Call freeze() once all species have been inserted; after that the container
 * is sorted by key and the hot-path accessors photons()/baryons() are valid.
 *
 * Iteration yields `Entry&`, where Entry exposes .key and .species. The common
 * form
 *     for (auto& sp : all_species_) { sp->Foo(); }
 * is supported via an Entry that acts like a unique_ptr<BaseSpecies>& (operator->
 * and operator* forward to the owned species). Use `sp.key` when the lookup key
 * is needed; use `sp->name()` only when the BaseSpecies type name is wanted
 * (these differ once #246 lands).
 */
class SpeciesCollection {
 public:
  using Ptr = std::unique_ptr<BaseSpecies>;

  struct Entry {
    std::string key;
    Ptr species;

    // Forwarders mimic std::unique_ptr: const-on-Entry does NOT propagate to
    // the pointee, matching how `const std::unique_ptr<T>&` iteration behaved
    // under the previous std::map-backed storage.
    BaseSpecies* operator->() const {
      return species.get();
    }
    BaseSpecies& operator*() const {
      return *species;
    }
    BaseSpecies* get() const {
      return species.get();
    }
    explicit operator bool() const {
      return static_cast<bool>(species);
    }
  };

  using Container = std::vector<Entry>;

  SpeciesCollection()                                    = default;
  SpeciesCollection(const SpeciesCollection&)            = delete;
  SpeciesCollection& operator=(const SpeciesCollection&) = delete;
  SpeciesCollection(SpeciesCollection&&)                 = default;
  SpeciesCollection& operator=(SpeciesCollection&&)      = default;

  // ── Construction-phase mutation ─────────────────────────────────────────
  /** Insert a species under an explicit lookup key. Must be called before
   *  freeze(). Duplicate keys abort in debug builds. */
  void insert(std::string key, Ptr species);

  /** Sort by key and resolve cached Photons/Baryons pointers.
   *  Must be called exactly once, after all insert() calls. */
  void freeze();

  // ── Keyed lookup (replaces std::map API) ────────────────────────────────
  std::size_t count(const std::string& key) const;

  /** Returns the owning unique_ptr, matching std::map::at semantics.
   *  Throws std::out_of_range if absent. Call sites typically write
   *  `*all_species_.at("X")` (BaseSpecies&) or
   *  `all_species_.at("X")->Foo()` (forwards through unique_ptr). */
  std::unique_ptr<BaseSpecies>& at(const std::string& key);
  const std::unique_ptr<BaseSpecies>& at(const std::string& key) const;

  /** Pointer-or-nullptr lookup. Returns the owning unique_ptr (not its payload)
   *  so call sites can write `if (auto* p = all_species_.find("X")) (*p)->Foo();`
   *  (pairing with the at() shape). */
  std::unique_ptr<BaseSpecies>* find(const std::string& key);
  const std::unique_ptr<BaseSpecies>* find(const std::string& key) const;

  /** Lex-sorted index of an inserted species, or size() if absent.
   *  Valid only after freeze(); the lower_bound search relies on the
   *  container being sorted. Asserted in debug builds. */
  std::size_t index_of(const std::string& key) const {
    assert(frozen_);
    auto it = std::lower_bound(species_.begin(),
                               species_.end(),
                               key,
                               [](const Entry& e, const std::string& k) { return e.key < k; });
    if (it == species_.end() || it->key != key)
      return species_.size();
    return static_cast<std::size_t>(it - species_.begin());
  }

  // ── Hot-path cached accessors ───────────────────────────────────────────
  // Valid only after freeze(); asserted in debug builds.
  BaseSpecies& photons() {
    assert(frozen_ && photons_);
    return *photons_;
  }
  BaseSpecies& photons() const {
    assert(frozen_ && photons_);
    return *photons_;
  }
  BaseSpecies& baryons() {
    assert(frozen_ && baryons_);
    return *baryons_;
  }
  BaseSpecies& baryons() const {
    assert(frozen_ && baryons_);
    return *baryons_;
  }
  std::size_t photons_index() const {
    assert(frozen_);
    return photons_index_;
  }
  std::size_t baryons_index() const {
    assert(frozen_);
    return baryons_index_;
  }
  /** True iff the collection holds any NCDM-family species (an NCDMBaseSpecies,
   *  or the DNCDM_DR_Species composite). Valid only after freeze(). */
  bool has_ncdm() const {
    assert(frozen_);
    return has_ncdm_;
  }

  // ── Iteration / size ────────────────────────────────────────────────────
  Container::iterator begin() {
    return species_.begin();
  }
  Container::iterator end() {
    return species_.end();
  }
  Container::const_iterator begin() const {
    return species_.begin();
  }
  Container::const_iterator end() const {
    return species_.end();
  }

  bool empty() const {
    return species_.empty();
  }
  std::size_t size() const {
    return species_.size();
  }

  /** Index-based access for parallel-vector dispatch (e.g. with
   *  perturb_vector::species_layouts). i must be < size(). */
  BaseSpecies* operator[](std::size_t i) {
    assert(i < species_.size());
    return species_[i].species.get();
  }
  const BaseSpecies* operator[](std::size_t i) const {
    assert(i < species_.size());
    return species_[i].species.get();
  }

 private:
  Container species_;
  BaseSpecies* photons_      = nullptr;
  BaseSpecies* baryons_      = nullptr;
  std::size_t photons_index_ = 0;
  std::size_t baryons_index_ = 0;
  bool has_ncdm_             = false;
  bool frozen_               = false;
};
