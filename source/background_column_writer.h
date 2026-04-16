#pragma once
#include "common.h"  // class_store_columntitle, class_store_double, _TRUE_, _FALSE_

/**
 * Thin helper for background output — analogous to PerturbColumnWriter.
 * Construct in title mode (titles != nullptr) or data mode (dataptr != nullptr).
 * Call Add() once per column in both modes; the writer handles the branching.
 */
class BackgroundColumnWriter {
 public:
  // Title mode
  explicit BackgroundColumnWriter(char* titles) : titles_(titles) {}

  // Data mode
  BackgroundColumnWriter(double* dataptr, int& storeidx)
      : dataptr_(dataptr), storeidx_(&storeidx) {}

  bool IsTitleMode() const {
    return titles_ != nullptr;
  }

  /** Write one column: title in title mode, value in data mode.
   *  condition=false: column is omitted entirely — no title appended, no slot
   *  written and the store index is NOT advanced. This matches class_store_double
   *  semantics. Both WriteBackgroundColumnTitles and WriteBackgroundData must
   *  pass the same condition so title and data counts stay aligned. */
  void Add(const char* title, double value, bool condition = true);

 private:
  char* titles_    = nullptr;
  double* dataptr_ = nullptr;
  int* storeidx_   = nullptr;
};
