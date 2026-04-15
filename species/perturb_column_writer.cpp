#include "perturb_source_context.h"

void PerturbColumnWriter::Add(const char* title, int tp_index, bool active) {
  if (titles_) {
    class_store_columntitle(titles_, title, active ? _TRUE_ : _FALSE_);
  }
  else if (dataptr_ && tk_) {
    // Guard tp_index access: only read tk_[tp_index] when active, because
    // inactive columns have an uninitialized tp_index (the source was never
    // registered via class_define_index).
    // NOTE: (*storeidx_) — explicit parentheses required so the macro's
    // "dataindex++" expands to "(*storeidx_)++" (increment the int value)
    // rather than "*storeidx_++" (advance the pointer).
    if (active) {
      class_store_double(dataptr_, tk_[tp_index], _TRUE_, (*storeidx_));
    }
  }
}

void PerturbColumnWriter::Add(const char* title, double value, bool active) {
  if (titles_) {
    class_store_columntitle(titles_, title, active ? _TRUE_ : _FALSE_);
  }
  else if (dataptr_) {
    // Same parenthesis rule: (*storeidx_) not *storeidx_ to avoid pointer advance.
    class_store_double(dataptr_, value, active ? _TRUE_ : _FALSE_, (*storeidx_));
  }
}
