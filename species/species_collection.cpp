#include "species_collection.h"

#include <algorithm>
#include <stdexcept>

#include "dncdm_dr_species.h"
#include "ncdm_base_species.h"

void SpeciesCollection::insert(std::string key, Ptr species) {
  if (frozen_) {
    throw std::logic_error("SpeciesCollection::insert: called after freeze()");
  }
  if (!species) {
    throw std::invalid_argument("SpeciesCollection::insert: null species");
  }
  for (const auto& e : species_) {
    if (e.key == key) {
      throw std::invalid_argument("SpeciesCollection::insert: duplicate key '" + key + "'");
    }
  }
  species_.push_back(Entry{std::move(key), std::move(species)});
}

void SpeciesCollection::freeze() {
  if (frozen_) {
    throw std::logic_error("SpeciesCollection::freeze: called twice");
  }
  std::sort(species_.begin(), species_.end(), [](const Entry& a, const Entry& b) {
    return a.key < b.key;
  });
  /* Stamp each species with its sorted index so PrintVariables can look up
     its own layout in ppw->pv->species_layouts[collection_index_]. */
  for (std::size_t i = 0; i < species_.size(); ++i)
    species_[i].species->collection_index_ = i;
  auto* photons_slot = find("Photons");
  auto* baryons_slot = find("Baryons");
  if (!photons_slot) {
    throw std::logic_error("SpeciesCollection::freeze: Photons species missing");
  }
  if (!baryons_slot) {
    throw std::logic_error("SpeciesCollection::freeze: Baryons species missing");
  }
  photons_ = photons_slot->get();
  baryons_ = baryons_slot->get();
  for (std::size_t i = 0; i < species_.size(); ++i) {
    if (species_[i].key == "Photons")
      photons_index_ = i;
    if (species_[i].key == "Baryons")
      baryons_index_ = i;
    if (dynamic_cast<NCDMBaseSpecies*>(species_[i].get()) ||
        dynamic_cast<DNCDM_DR_Species*>(species_[i].get()))
      has_ncdm_ = true;
  }
  frozen_ = true;
}

namespace {
template <class Vec>
auto find_impl(Vec& v, const std::string& key) -> decltype(&v[0].species) {
  for (auto& e : v)
    if (e.key == key)
      return &e.species;
  return nullptr;
}
}  // namespace

std::size_t SpeciesCollection::count(const std::string& key) const {
  return find(key) != nullptr ? 1U : 0U;
}

std::unique_ptr<BaseSpecies>* SpeciesCollection::find(const std::string& key) {
  return find_impl(species_, key);
}
const std::unique_ptr<BaseSpecies>* SpeciesCollection::find(const std::string& key) const {
  return find_impl(species_, key);
}

std::unique_ptr<BaseSpecies>& SpeciesCollection::at(const std::string& key) {
  if (auto* p = find(key))
    return *p;
  throw std::out_of_range("SpeciesCollection::at: no species with key '" + key + "'");
}
const std::unique_ptr<BaseSpecies>& SpeciesCollection::at(const std::string& key) const {
  if (auto* p = find(key))
    return *p;
  throw std::out_of_range("SpeciesCollection::at: no species with key '" + key + "'");
}
