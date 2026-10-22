// Typed scalar attributes and one-dimensional datasets over the HDF5 C++ API.
// Private to the io subsystem: HDF5 types never reach public headers.

#pragma once

#include <H5Cpp.h>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace wasl {

template <typename T> const H5::PredType &h5_type();
template <> inline const H5::PredType &h5_type<double>() { return H5::PredType::NATIVE_DOUBLE; }
template <> inline const H5::PredType &h5_type<std::int64_t>() {
  return H5::PredType::NATIVE_INT64;
}
template <> inline const H5::PredType &h5_type<std::int32_t>() {
  return H5::PredType::NATIVE_INT32;
}
template <> inline const H5::PredType &h5_type<std::uint64_t>() {
  return H5::PredType::NATIVE_UINT64;
}
template <> inline const H5::PredType &h5_type<unsigned char>() {
  return H5::PredType::NATIVE_UCHAR;
}

template <typename T> void write_attribute(H5::H5Object &object, const char *name, const T &value) {
  H5::DataSpace scalar(H5S_SCALAR);
  object.createAttribute(name, h5_type<T>(), scalar).write(h5_type<T>(), &value);
}

template <typename T> T read_attribute(const H5::H5Object &object, const char *name) {
  T value{};
  object.openAttribute(name).read(h5_type<T>(), &value);
  return value;
}

template <typename T>
void write_vector(H5::H5File &file, const char *name, const std::vector<T> &values) {
  const hsize_t size = values.size();
  H5::DataSpace space(1, &size);
  auto dataset = file.createDataSet(name, h5_type<T>(), space);
  if (!values.empty())
    dataset.write(values.data(), h5_type<T>());
}

template <typename T> std::vector<T> read_vector(const H5::H5File &file, const char *name) {
  const auto dataset = file.openDataSet(name);
  const auto space = dataset.getSpace();
  if (space.getSimpleExtentNdims() != 1)
    throw std::runtime_error(std::string("HDF5 I/O: dataset ") + name + " must be one-dimensional");
  hsize_t size = 0;
  space.getSimpleExtentDims(&size);
  std::vector<T> values(static_cast<std::size_t>(size));
  if (!values.empty())
    dataset.read(values.data(), h5_type<T>());
  return values;
}

// Oldest file format version any reader accepts.
constexpr int earliest_supported_format_version = 1;

} // namespace wasl
