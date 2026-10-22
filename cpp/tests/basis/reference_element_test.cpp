#include <cassert>

#include "wasl/basis/reference_element.hpp"

int main() {
  static_assert(wasl::line_mode_count(0) == 1);
  static_assert(wasl::line_mode_count(1) == 2);
  static_assert(wasl::line_mode_count(2) == 3);
  static_assert(wasl::line_mode_count(3) == 4);
  static_assert(wasl::triangle_mode_count(0) == 1);
  static_assert(wasl::triangle_mode_count(1) == 3);
  static_assert(wasl::triangle_mode_count(2) == 6);
  static_assert(wasl::triangle_mode_count(3) == 10);
  assert(wasl::triangle_mode_count(4) == 15);
}
