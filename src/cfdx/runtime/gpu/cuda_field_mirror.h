#pragma once

#include "cuda_device_buffer.h"
#include "cfdx/core/field/field.h"
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::runtime::gpu {

// Explicit field residency bridge. The core Field remains backend-neutral;
// this adapter owns device buffers and makes host/device transitions visible
// to the runtime instead of hiding transfers in numerical kernels.
template<cfdx::core::Location L>
class CudaFieldMirror {
public:
    CudaFieldMirror() = default;

    void allocate(const cfdx::core::Field<double,L>& field) {
        components_.clear();
        components_.resize(field.dimension());
        for (auto& b : components_) b.resize(field.size() * sizeof(double));
        size_ = field.size();
        dimension_ = field.dimension();
    }

    void upload(const cfdx::core::Field<double,L>& field) {
        if (field.size() != size_ || field.dimension() != dimension_)
            allocate(field);
        for (std::size_t c = 0; c < dimension_; ++c)
            components_[c].copy_from_host(field.component_data(c),
                                          field.size() * sizeof(double));
    }

    void download(cfdx::core::Field<double,L>& field) const {
        if (field.size() != size_ || field.dimension() != dimension_)
            throw std::invalid_argument("CudaFieldMirror: field shape mismatch");
        for (std::size_t c = 0; c < dimension_; ++c)
            components_[c].copy_to_host(field.component_data(c),
                                        field.size() * sizeof(double));
    }

    std::size_t size() const noexcept { return size_; }
    std::size_t dimension() const noexcept { return dimension_; }
    const CudaDeviceBuffer& component(std::size_t c) const { return components_.at(c); }
    CudaDeviceBuffer& component(std::size_t c) { return components_.at(c); }

private:
    std::size_t size_ = 0;
    std::size_t dimension_ = 0;
    std::vector<CudaDeviceBuffer> components_;
};

} // namespace cfdx::runtime::gpu
