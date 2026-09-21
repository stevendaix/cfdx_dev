#pragma once
#include <cstddef>
#include <string>
#include <vector>
namespace cfdx::core::memory {
enum class Residency { Persistent, Cached, Ephemeral };
struct FieldLifetime { std::string name;std::size_t bytes=0;int birth=0;int death=0;Residency residency=Residency::Ephemeral; };
inline bool reusable(const FieldLifetime&a,const FieldLifetime&b){return a.residency!=Residency::Persistent&&b.residency!=Residency::Persistent&&(a.death<=b.birth||b.death<=a.birth);}
} // namespace cfdx::core::memory
