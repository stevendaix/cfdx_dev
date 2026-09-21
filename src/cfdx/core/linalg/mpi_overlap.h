#pragma once
#include <functional>
#include <utility>
namespace cfdx::core {

// Backend-neutral communication/computation overlap schedule.
// The MPI implementation can map begin_exchange/wait_exchange to
// MPI_Irecv/MPI_Isend/MPI_Waitall without making the core solver MPI-aware.
class HaloOverlap {
public:
 using Callback=std::function<void()>;
 HaloOverlap(Callback begin_exchange,Callback interior,Callback wait_exchange,Callback boundary)
 :begin_(std::move(begin_exchange)),interior_(std::move(interior)),wait_(std::move(wait_exchange)),boundary_(std::move(boundary)){}
 void execute()const{begin_();interior_();wait_();boundary_();}
private:Callback begin_,interior_,wait_,boundary_;
};
} // namespace cfdx::core
