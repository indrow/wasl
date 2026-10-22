// File locations of one case invocation.

#pragma once

#include <filesystem>

namespace wasl {

struct CaseFiles {
  std::filesystem::path prepared_mesh;
  std::filesystem::path output_directory;
  std::filesystem::path restart_checkpoint; // empty for a fresh start
};

} // namespace wasl
