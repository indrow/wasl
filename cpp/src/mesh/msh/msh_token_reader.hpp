// Whitespace-token reader with descriptive MSH errors.

#pragma once

#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace wasl {

class MshTokenReader {
public:
  explicit MshTokenReader(const std::string &path) : input_{path} {
    if (!input_)
      throw std::runtime_error("MSH reader: cannot open " + path);
  }

  template <class T> T read(const char *context) {
    T value{};
    if (!(input_ >> value))
      throw std::runtime_error("MSH reader: invalid or missing " + std::string{context});
    return value;
  }

  std::string read_quoted(const char *context) {
    std::string value;
    if (!(input_ >> std::quoted(value)))
      throw std::runtime_error("MSH reader: invalid or missing " + std::string{context});
    return value;
  }

  // Read a nonnegative int-sized count.
  int count(const char *context) {
    const auto value = read<long long>(context);
    if (value < 0 || value > std::numeric_limits<int>::max())
      throw std::runtime_error("MSH reader: invalid " + std::string{context});
    return static_cast<int>(value);
  }

  std::vector<int> read_ints(int count, const char *context) {
    std::vector<int> values;
    values.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i)
      values.push_back(read<int>(context));
    return values;
  }

  void expect(const char *expected) {
    const auto token = read<std::string>(expected);
    if (token != expected)
      throw std::runtime_error("MSH reader: expected " + std::string{expected} + ", got " + token);
  }

  bool next(std::string &token) { return static_cast<bool>(input_ >> token); }

private:
  std::ifstream input_;
};

} // namespace wasl
