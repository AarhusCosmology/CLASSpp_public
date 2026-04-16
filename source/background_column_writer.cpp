#include "background_column_writer.h"

void BackgroundColumnWriter::Add(const char* title, double value, bool condition) {
  if (titles_) {
    class_store_columntitle(titles_, title, condition ? _TRUE_ : _FALSE_);
  }
  else if (dataptr_) {
    class_store_double(dataptr_, value, condition ? _TRUE_ : _FALSE_, (*storeidx_));
  }
}
