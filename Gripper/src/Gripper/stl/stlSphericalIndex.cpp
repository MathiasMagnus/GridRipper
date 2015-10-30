#include <Gripper/stl/stlSphericalIndex.hpp>


//////////////////////////////////////
//                                  //
// Multipole::stl::Spherical::Index //
//                                  //
//////////////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::stl::Spherical::Index::Index() {}


Multipole::stl::Spherical::Index::Index(const Index& in) : i(in.i) {}


Multipole::stl::Spherical::Index::Index(Index&& src) : i(std::move(src.i)) {}


Multipole::stl::Spherical::Index::~Index() {}


Multipole::stl::Spherical::Index& Multipole::stl::Spherical::Index::operator=(const Index& rhs) { i = rhs.i; return *this; }


Multipole::stl::Spherical::Index::Index(const value_type& in) : i(in) {}


Multipole::stl::Spherical::Index::operator value_type&() { return i; }

//////////////////////
// member operators //
//////////////////////

// Unary
Multipole::stl::Spherical::Index Multipole::stl::Spherical::Index::operator++(int rhs)   { Index result(i); i += 1; return result; }
Multipole::stl::Spherical::Index& Multipole::stl::Spherical::Index::operator++()         { i += 1; return *this; }
Multipole::stl::Spherical::Index Multipole::stl::Spherical::Index::operator--(int rhs)   { Index result(i); i -= 1; return result; }
Multipole::stl::Spherical::Index& Multipole::stl::Spherical::Index::operator--()         { i -= 1; return *this; }

// Binary
Multipole::stl::Spherical::Index& Multipole::stl::Spherical::Index::operator+=(const Index& rhs) { i += rhs.i; return *this; }
Multipole::stl::Spherical::Index& Multipole::stl::Spherical::Index::operator-=(const Index& rhs) { i -= rhs.i; return *this; }
Multipole::stl::Spherical::Index& Multipole::stl::Spherical::Index::operator*=(const Index& rhs) { i *= rhs.i; return *this; }
Multipole::stl::Spherical::Index& Multipole::stl::Spherical::Index::operator/=(const Index& rhs) { i /= rhs.i; return *this; }

//////////////////////////
// non-member operators //
//////////////////////////

// Unary 
Multipole::stl::Spherical::Index operator+(const Multipole::stl::Spherical::Index& rhs) { return Multipole::stl::Spherical::Index(rhs); }
Multipole::stl::Spherical::Index operator-(const Multipole::stl::Spherical::Index& rhs) { return Multipole::stl::Spherical::Index(-rhs.i); }

// Binary
Multipole::stl::Spherical::Index operator+(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs) { return Multipole::stl::Spherical::Index(lhs.i + rhs.i); }
Multipole::stl::Spherical::Index operator-(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs) { return Multipole::stl::Spherical::Index(lhs.i - rhs.i); }
Multipole::stl::Spherical::Index operator*(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs) { return Multipole::stl::Spherical::Index(lhs.i * rhs.i); }
Multipole::stl::Spherical::Index operator/(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs) { return Multipole::stl::Spherical::Index(lhs.i / rhs.i); }
Multipole::stl::Spherical::Index operator%(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs) { return Multipole::stl::Spherical::Index(lhs.i % rhs.i); }

bool operator< (const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs) { return lhs.i < rhs.i; }
bool operator>(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs) { return lhs.i > rhs.i; }
bool operator<=(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs) { return lhs.i <= rhs.i; }
bool operator>=(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs) { return lhs.i >= rhs.i; }
bool operator==(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs) { return lhs.i == rhs.i; }
bool operator!=(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs) { return lhs.i != rhs.i; }


//////////////////////////////////////////
//                                      //
// Multipole::stl::Spherical::IndexPair //
//                                      //
//////////////////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::stl::Spherical::IndexPair::IndexPair() {}


Multipole::stl::Spherical::IndexPair::IndexPair(const IndexPair& in) : l(in.l), m(in.m) {}


Multipole::stl::Spherical::IndexPair::IndexPair(const Index& l_in, const Index& m_in) : l(l_in), m(m_in) {}


Multipole::stl::Spherical::IndexPair::IndexPair(const Gaunt::Index& i_in)
{
    l = static_cast<Multipole::stl::Spherical::Index::value_type>(std::round((-1. + std::sqrt(1. + 4 * i_in.i)) / 2.));
    m = i_in.i - l*(l + 1);
}


Multipole::stl::Spherical::IndexPair::~IndexPair() {}

//////////////////////////
// non-member operators //
//////////////////////////

// Binary
bool operator==(const Multipole::stl::Spherical::IndexPair& lhs, const Multipole::stl::Spherical::IndexPair& rhs) { return (lhs.l == rhs.l) && (lhs.m == rhs.m); }
bool operator!=(const Multipole::stl::Spherical::IndexPair& lhs, const Multipole::stl::Spherical::IndexPair& rhs) { return (lhs.l != rhs.l) || (lhs.m != rhs.m); }