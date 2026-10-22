#include <iostream>

#include "wasl/basis/reference_element.hpp"

int main() {
  std::cout << "Wasl CPU toolchain ready; P2 line modes: " << wasl::line_mode_count(2) << '\n';
}
