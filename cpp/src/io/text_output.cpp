#include "io/text_output.hpp"

#include <iomanip>
#include <stdexcept>

namespace wasl {

void create_parent_directories(const std::filesystem::path &path) {
  if (path.has_parent_path())
    std::filesystem::create_directories(path.parent_path());
}

std::ofstream open_text_output(const std::filesystem::path &path) {
  create_parent_directories(path);
  std::ofstream output(path);
  if (!output)
    throw std::runtime_error("file output: cannot open " + path.string());
  output << std::setprecision(17);
  return output;
}

std::string xml_escape(const std::string &value) {
  std::string result;
  for (const char character : value) {
    switch (character) {
    case '&':
      result += "&amp;";
      break;
    case '<':
      result += "&lt;";
      break;
    case '"':
      result += "&quot;";
      break;
    default:
      result += character;
    }
  }
  return result;
}

} // namespace wasl
