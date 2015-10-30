#include <Gripper/stl/stlGauntCoefficient.hpp>


////////////////////////////////////////
//                                    //
// Multipole::stl::Gaunt::Coefficient //
//                                    //
////////////////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::stl::Gaunt::Coefficient::Coefficient() {}


Multipole::stl::Gaunt::Coefficient::Coefficient(const Coefficient& in) : i1(in.i1), i2(in.i2), i3(in.i3), value(in.value) {}


Multipole::stl::Gaunt::Coefficient::Coefficient(Coefficient&& src) : i1(std::move(src.i1)), i2(std::move(src.i2)), i3(std::move(src.i3)), value(std::move(src.value)) {}


Multipole::stl::Gaunt::Coefficient::~Coefficient() {}


Multipole::stl::Gaunt::Coefficient& Multipole::stl::Gaunt::Coefficient::operator=(const Coefficient& rhs) { i1 = rhs.i1; i2 = rhs.i2; i3 = rhs.i3; value = rhs.value; return *this; }


Multipole::stl::Gaunt::Coefficient::Coefficient(const Multipole::Gaunt::Coefficient& in) : i1(in.l1m1), i2(in.l2m2), i3(in.l3m3), value(in.value) {}


Multipole::stl::Gaunt::Coefficient::Coefficient(const index_type& i1_in, const index_type& i2_in, const index_type& i3_in, const value_type& value_in) : i1(i1_in), i2(i2_in), i3(i3_in), value(value_in) {}


Multipole::stl::Gaunt::Coefficient::Coefficient(const Spherical::IndexPair& l1m1_in, const Spherical::IndexPair& l2m2_in, const Spherical::IndexPair& l3m3_in, const value_type& value_in) : i1(l1m1_in), i2(l2m2_in), i3(l3m3_in), value(value_in) {}