#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace cfdx::core {

enum class Location : std::uint8_t { CELL=0, FACE, POINT, BOUNDARY, UNKNOWN };

inline const char* to_string(Location loc) {
    switch (loc) {
        case Location::CELL: return "cell";
        case Location::FACE: return "face";
        case Location::POINT: return "point";
        case Location::BOUNDARY: return "boundary";
        default: return "unknown";
    }
}
inline Location location_from_string(const std::string& s) {
    if (s=="cell") return Location::CELL;
    if (s=="face") return Location::FACE;
    if (s=="point") return Location::POINT;
    if (s=="boundary") return Location::BOUNDARY;
    return Location::UNKNOWN;
}

enum class Precision : std::uint8_t { FLOAT32=0, FLOAT64 };
inline const char* to_string(Precision p) { return p==Precision::FLOAT32 ? "float32" : "float64"; }
inline Precision precision_from_string(const std::string& s) {
    if (s=="float32" || s=="float") return Precision::FLOAT32;
    if (s=="float64" || s=="double") return Precision::FLOAT64;
    return Precision::FLOAT64;
}

struct FieldMetadata {
    std::string name;
    Location location = Location::UNKNOWN;
    std::size_t dimension = 1;
    std::string unit;
    Precision precision = Precision::FLOAT64;
    std::string storage = "host";
};

struct Vec3 {
    double x=0.0, y=0.0, z=0.0;
    Vec3()=default;
    Vec3(double v):x(v),y(v),z(v){}
    Vec3(double x_,double y_,double z_):x(x_),y(y_),z(z_){}
    Vec3 operator+(const Vec3&o)const{return{x+o.x,y+o.y,z+o.z};}
    Vec3& operator+=(const Vec3&o)noexcept{x+=o.x;y+=o.y;z+=o.z;return *this;}
    Vec3 operator-(const Vec3&o)const{return{x-o.x,y-o.y,z-o.z};}
    Vec3 operator*(double s)const{return{x*s,y*s,z*s};}
    Vec3& operator*=(double s)noexcept{x*=s;y*=s;z*=s;return *this;}
    double dot(const Vec3&o)const{return x*o.x+y*o.y+z*o.z;}
    Vec3 cross(const Vec3&o)const{return{y*o.z-z*o.y,z*o.x-x*o.z,x*o.y-y*o.x};}
    double mag()const{return std::sqrt(x*x+y*y+z*z);}
    double mag2()const{return x*x+y*y+z*z;}
    Vec3 normalized()const{const double m=mag();return m==0.0?Vec3{}:Vec3{x/m,y/m,z/m};}
};

template<class T,std::size_t Alignment=64>
class AlignedAllocator {
public:
    using value_type=T;
    using size_type=std::size_t;
    using difference_type=std::ptrdiff_t;
    template<class U> struct rebind { using other=AlignedAllocator<U,Alignment>; };
    AlignedAllocator() noexcept = default;
    template<class U> AlignedAllocator(const AlignedAllocator<U,Alignment>&) noexcept {}
    T* allocate(std::size_t n) {
        if(n>std::numeric_limits<std::size_t>::max()/sizeof(T)) throw std::bad_array_new_length();
        if(n==0) return nullptr;
        void* p=::operator new(n*sizeof(T),std::align_val_t(Alignment));
        return static_cast<T*>(p);
    }
    void deallocate(T* p,std::size_t) noexcept {
        ::operator delete(p,std::align_val_t(Alignment));
    }
};
template<class T,std::size_t A,class U,std::size_t B>
constexpr bool operator==(const AlignedAllocator<T,A>&,const AlignedAllocator<U,B>&) noexcept { return A==B; }
template<class T,std::size_t A,class U,std::size_t B>
constexpr bool operator!=(const AlignedAllocator<T,A>&a,const AlignedAllocator<U,B>&b) noexcept { return !(a==b); }

template<typename T=double, Location L=Location::CELL>
class Field {
    static_assert(std::is_arithmetic_v<T>, "Field<T> uses arithmetic component storage; use Field<Vec3,L> specialization for vectors.");
public:
    using ValueType=T;
    static constexpr Location location=L;
    Field()=default;
    explicit Field(std::size_t n,const std::string& name="",const std::string& unit="",std::size_t dim=1,Precision precision=Precision::FLOAT64)
      : meta_{name,L,dim?dim:1,unit,precision,"host"},n_(n),dim_(dim?dim:1) { allocate_storage(); }

    std::size_t size()const noexcept{return n_;}
    std::size_t dimension()const noexcept{return dim_;}
    void set_dimension(std::size_t dim){ if(dim==0) dim=1; if(dim!=dim_){dim_=dim;allocate_storage();} }
    bool empty()const noexcept{return n_==0;}
    Precision precision()const noexcept{return meta_.precision;}

    T operator()(std::size_t i)const{check_index(i);return data_[i];}
    T& operator()(std::size_t i){check_index(i);return data_[i];}
    void get(std::size_t i,T&x,T&y,T&z)const{check_index(i);require_vector();x=data_[i];y=data_[n_+i];z=data_[2*n_+i];}
    void set(std::size_t i,T x,T y,T z){check_index(i);require_vector();data_[i]=x;data_[n_+i]=y;data_[2*n_+i]=z;}
    T& operator()(std::size_t i,std::size_t comp){check_index(i);check_component(comp);return data_[comp*n_+i];}
    T operator()(std::size_t i,std::size_t comp)const{check_index(i);check_component(comp);return data_[comp*n_+i];}

    const FieldMetadata& metadata()const noexcept{return meta_;}
    FieldMetadata& metadata()noexcept{return meta_;}
    const std::string& name()const noexcept{return meta_.name;}
    Location loc()const noexcept{return meta_.location;}

    void fill(T value){std::fill(data_.begin(),data_.begin()+n_,value);if(dim_>1)for(std::size_t c=1;c<dim_;++c)std::fill(data_.begin()+c*n_,data_.begin()+(c+1)*n_,value);}
    void resize(std::size_t n){n_=n;allocate_storage();}
    void clear(){n_=0;data_.clear();component_ptrs_.clear();}
    const T* const* data()const noexcept{return component_ptrs_.data();}
    T* const* data()noexcept{return component_ptrs_.data();}
    const T* component_data(std::size_t comp)const noexcept{return comp<dim_?data_.data()+comp*n_:nullptr;}
    T* component_data(std::size_t comp)noexcept{return comp<dim_?data_.data()+comp*n_:nullptr;}
    bool is_valid()const{for(const T&v:data_)if(!std::isfinite(v))return false;return true;}

private:
    void check_index(std::size_t i)const{if(i>=n_)throw std::out_of_range("Field: index out of range");}
    void check_component(std::size_t c)const{if(c>=dim_)throw std::out_of_range("Field: component index out of range");}
    void require_vector()const{if(dim_<3)throw std::runtime_error("Field: dimension < 3");}
    void allocate_storage(){
        data_.assign(n_*dim_,T{});
        component_ptrs_.resize(dim_);
        for(std::size_t c=0;c<dim_;++c)component_ptrs_[c]=data_.data()+c*n_;
    }
    FieldMetadata meta_;
    std::size_t n_=0,dim_=1;
    std::vector<T,AlignedAllocator<T,64>> data_;
    std::vector<T*> component_ptrs_;
};

template<Location L>
class Field<Vec3,L> {
public:
    using ValueType=Vec3;
    static constexpr Location location=L;
    Field()=default;
    explicit Field(std::size_t n,const std::string& name="",const std::string& unit="",std::size_t dim=3,Precision precision=Precision::FLOAT64)
      : meta_{name,L,3,unit,precision,"host"},n_(n) { (void)dim; allocate_storage(); }
    std::size_t size()const noexcept{return n_;}
    std::size_t dimension()const noexcept{return 3;}
    void set_dimension(std::size_t dim){if(dim!=3)throw std::invalid_argument("Field<Vec3>: dimension must be 3");}
    bool empty()const noexcept{return n_==0;}
    Precision precision()const noexcept{return meta_.precision;}
    Vec3 operator()(std::size_t i)const{check(i);return{x_[i],y_[i],z_[i]};}
    Vec3& operator()(std::size_t)=delete;
    void get(std::size_t i,Vec3&v)const{v=(*this)(i);}
    void set(std::size_t i,double x,double y,double z){check(i);x_[i]=x;y_[i]=y;z_[i]=z;}
    double& operator()(std::size_t i,std::size_t c){check(i);if(c>2)throw std::out_of_range("Field<Vec3>: component index out of range");return c==0?x_[i]:(c==1?y_[i]:z_[i]);}
    double operator()(std::size_t i,std::size_t c)const{check(i);if(c>2)throw std::out_of_range("Field<Vec3>: component index out of range");return c==0?x_[i]:(c==1?y_[i]:z_[i]);}
    const FieldMetadata& metadata()const noexcept{return meta_;}
    FieldMetadata& metadata()noexcept{return meta_;}
    const std::string& name()const noexcept{return meta_.name;}
    Location loc()const noexcept{return meta_.location;}
    void fill(Vec3 v){std::fill(x_.begin(),x_.end(),v.x);std::fill(y_.begin(),y_.end(),v.y);std::fill(z_.begin(),z_.end(),v.z);}
    void resize(std::size_t n){n_=n;allocate_storage();}
    void clear(){n_=0;x_.clear();y_.clear();z_.clear();ptrs_.clear();}
    const double* const* data()const noexcept{return ptrs_.data();}
    double* const* data()noexcept{return ptrs_.data();}
    const double* component_data(std::size_t c)const noexcept{return c<3?ptrs_[c]:nullptr;}
    double* component_data(std::size_t c)noexcept{return c<3?ptrs_[c]:nullptr;}
    bool is_valid()const{for(std::size_t i=0;i<n_;++i)if(!std::isfinite(x_[i])||!std::isfinite(y_[i])||!std::isfinite(z_[i]))return false;return true;}
private:
    void check(std::size_t i)const{if(i>=n_)throw std::out_of_range("Field<Vec3>: index out of range");}
    void allocate_storage(){x_.assign(n_,0.0);y_.assign(n_,0.0);z_.assign(n_,0.0);ptrs_={x_.data(),y_.data(),z_.data()};}
    FieldMetadata meta_;
    std::size_t n_=0;
    std::vector<double,AlignedAllocator<double,64>> x_,y_,z_;
    std::vector<double*> ptrs_;
};

using ScalarCellField=Field<double,Location::CELL>;
using ScalarFaceField=Field<double,Location::FACE>;
using ScalarPointField=Field<double,Location::POINT>;
using Vec3CellField=Field<Vec3,Location::CELL>;
using Vec3FaceField=Field<Vec3,Location::FACE>;
using Float32CellField=Field<float,Location::CELL>;
using Float64CellField=Field<double,Location::CELL>;

} // namespace cfdx::core
