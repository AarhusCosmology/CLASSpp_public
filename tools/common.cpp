#include "common.h"

#include <string>

void class_protect_sprintf(char* dest, const char* tpl, ...) {
  va_list args;
  va_start(args, tpl);
  vsnprintf(dest, 2048, tpl, args);
  va_end(args);
}

void class_protect_fprintf(FILE* stream, char* tpl, ...) {
  va_list args;
  char dest[6000];
  va_start(args, tpl);
  vsnprintf(dest, 2048, tpl, args);
  va_end(args);
  fprintf(stream, "%s", dest);
}

void* class_protect_memcpy(void* dest, void* from, size_t sz) {
  return memcpy(dest, from, sz);
}

int get_number_of_titles(const std::string& titlestring) {
  int number_of_titles = 0;
  for (char c : titlestring)
    if (c == '\t')
      number_of_titles++;
  return number_of_titles;
}
