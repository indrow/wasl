// Order-independent key of a triangle edge: its sorted pair of global vertex
// IDs, so an edge shared by two elements maps to one slot.

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>

namespace wasl {

struct FaceKey {
  std::int64_t first;  // smaller vertex ID
  std::int64_t second; // larger vertex ID

  bool operator==(const FaceKey &other) const noexcept {
    return first == other.first && second == other.second;
  }
};

inline FaceKey make_face_key(std::int64_t a, std::int64_t b) noexcept {
  return a < b ? FaceKey{a, b} : FaceKey{b, a};
}

struct FaceKeyHash {
  std::size_t operator()(const FaceKey &key) const noexcept {
    const auto h1 = std::hash<std::int64_t>{}(key.first);
    const auto h2 = std::hash<std::int64_t>{}(key.second);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};

} // namespace wasl
