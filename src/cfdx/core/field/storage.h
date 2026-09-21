#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include "cfdx/core/field/field.h"

namespace cfdx::core {

enum class StorageState : std::uint8_t { HOST_ONLY=0, DEVICE_ONLY, SYNCHRONIZED, HOST_DIRTY, DEVICE_DIRTY };
inline const char* to_string(StorageState s){
    switch(s){case StorageState::HOST_ONLY:return "host_only";case StorageState::DEVICE_ONLY:return "device_only";case StorageState::SYNCHRONIZED:return "synchronized";case StorageState::HOST_DIRTY:return "host_dirty";case StorageState::DEVICE_DIRTY:return "device_dirty";}
    return "unknown";
}

struct HostStorage {
    std::vector<double,AlignedAllocator<double,64>> data;
    std::size_t size()const noexcept{return data.size();}
    bool empty()const noexcept{return data.empty();}
    void resize(std::size_t n){data.resize(n,0.0);}
    void clear(){data.clear();}
    double* data_ptr()noexcept{return data.data();}
    const double* data_ptr()const noexcept{return data.data();}
};

struct DeviceStorage {
    std::size_t size_bytes=0;
    std::size_t size()const noexcept{return size_bytes/sizeof(double);}
    bool available()const noexcept{return false;}
};

struct WorkingSetStorage {
    std::size_t capacity=0,used=0;
    bool available()const noexcept{return capacity>0;}
};

class StorageHandle {
public:
    StorageHandle()=default;
    HostStorage& host()noexcept{return host_;}
    const HostStorage& host()const noexcept{return host_;}
    DeviceStorage& device()noexcept{return device_;}
    const DeviceStorage& device()const noexcept{return device_;}
    StorageState state()const noexcept{return state_;}
    void set_state(StorageState s)noexcept{state_=s;}
    void mark_host_dirty(){if(state_==StorageState::SYNCHRONIZED)state_=StorageState::HOST_DIRTY;}
    void mark_device_dirty(){if(state_==StorageState::SYNCHRONIZED)state_=StorageState::DEVICE_DIRTY;}
    void sync_host_to_device(){throw std::runtime_error("StorageHandle: CUDA device synchronization is unavailable in the core storage layer; use the GPU runtime backend");}
    void sync_device_to_host(){throw std::runtime_error("StorageHandle: CUDA device synchronization is unavailable in the core storage layer; use the GPU runtime backend");}
    std::size_t size()const noexcept{return host_.size();}
    bool empty()const noexcept{return host_.empty();}
    void resize(std::size_t n){host_.resize(n);device_.size_bytes=n*sizeof(double);state_=StorageState::HOST_ONLY;}
    void clear(){host_.clear();device_.size_bytes=0;state_=StorageState::HOST_ONLY;}
private:
    HostStorage host_;
    DeviceStorage device_;
    StorageState state_=StorageState::HOST_ONLY;
};

} // namespace cfdx::core
