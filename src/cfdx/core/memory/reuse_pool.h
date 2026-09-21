#pragma once
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>
namespace cfdx::core::memory {

// Lifetime-independent scratch allocator. Buffers are retained and reused by
// size/location, avoiding repeated heap allocation in iterative CFD kernels.
class ReusePool {
public:
 struct Block { std::size_t offset=0,size=0; bool free=true; };
 explicit ReusePool(std::size_t initial=0):storage_(initial){}
 std::size_t acquire(std::size_t bytes) {
  if(bytes==0)return 0;
  for(auto& b:blocks_) if(b.free&&b.size>=bytes){b.free=false;return b.offset;}
  const std::size_t offset=storage_.size();storage_.resize(offset+bytes);blocks_.push_back({offset,bytes,false});return offset;
 }
 void release(std::size_t offset) {
  for(auto& b:blocks_)if(b.offset==offset){b.free=true;return;}
  throw std::out_of_range("unknown pool block");
 }
 void reset(){for(auto& b:blocks_)b.free=true;}
 std::size_t bytes()const noexcept{return storage_.size();}
 std::size_t allocated_bytes()const noexcept{std::size_t n=0;for(const auto&b:blocks_)if(!b.free)n+=b.size;return n;}
private:
 std::vector<std::byte> storage_;
 std::vector<Block> blocks_;
};
} // namespace cfdx::core::memory
