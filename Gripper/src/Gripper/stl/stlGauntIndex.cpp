#include <Gripper/stl/stlGauntIndex.hpp>


//////////////////////////////////
//                              //
// Multipole::stl::Gaunt::Index //
//                              //
//////////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::stl::Gaunt::Index::Index() {}


Multipole::stl::Gaunt::Index::Index(const Index& in) : i(in.i) {}


Multipole::stl::Gaunt::Index::Index(Index&& src) : i(std::move(src.i)) {}


Multipole::stl::Gaunt::Index::~Index() {}


Multipole::stl::Gaunt::Index& Multipole::stl::Gaunt::Index::operator=(const Index& rhs) { i = rhs.i; return *this; }


Multipole::stl::Gaunt::Index::Index(const value_type& i_in) : i(i_in) {}


Multipole::stl::Gaunt::Index::Index(const Spherical::IndexPair& lm_in) : i(static_cast<Spherical::Index::value_type>(lm_in.l * (lm_in.l + 1) + lm_in.m)) {}


Multipole::stl::Gaunt::Index::Index(const Multipole::Spherical::IndexPair& lm_in) : i(static_cast<value_type>(lm_in.l.i * (lm_in.l.i + 1u) + lm_in.m.i)) {}


Multipole::stl::Gaunt::Index::operator const Multipole::stl::Gaunt::Index::value_type&() { return i; }

//////////////////////
// member operators //
//////////////////////

// Unary
Multipole::stl::Gaunt::Index Multipole::stl::Gaunt::Index::operator++(int rhs)   { Index result(i); i += 1; return result; }
Multipole::stl::Gaunt::Index& Multipole::stl::Gaunt::Index::operator++()         { i += 1; return *this; }
Multipole::stl::Gaunt::Index Multipole::stl::Gaunt::Index::operator--(int rhs)   { Index result(i); i -= 1; return result; }
Multipole::stl::Gaunt::Index& Multipole::stl::Gaunt::Index::operator--()         { i -= 1; return *this; }

// Binary
Multipole::stl::Gaunt::Index& Multipole::stl::Gaunt::Index::operator+=(const Index& rhs) { i += rhs.i; return *this; }
Multipole::stl::Gaunt::Index& Multipole::stl::Gaunt::Index::operator-=(const Index& rhs) { i -= rhs.i; return *this; }
Multipole::stl::Gaunt::Index& Multipole::stl::Gaunt::Index::operator*=(const Index& rhs) { i *= rhs.i; return *this; }
Multipole::stl::Gaunt::Index& Multipole::stl::Gaunt::Index::operator/=(const Index& rhs) { i /= rhs.i; return *this; }

//////////////////////////
// non-member operators //
//////////////////////////

// Binary
Multipole::stl::Gaunt::Index operator+(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs) { return Multipole::stl::Gaunt::Index(lhs.i + rhs.i); }
Multipole::stl::Gaunt::Index operator-(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs) { return Multipole::stl::Gaunt::Index(lhs.i - rhs.i); }
Multipole::stl::Gaunt::Index operator*(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs) { return Multipole::stl::Gaunt::Index(lhs.i * rhs.i); }
Multipole::stl::Gaunt::Index operator/(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs) { return Multipole::stl::Gaunt::Index(lhs.i / rhs.i); }
Multipole::stl::Gaunt::Index operator%(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs) { return Multipole::stl::Gaunt::Index(lhs.i % rhs.i); }

bool operator< (const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs) { return lhs.i < rhs.i; }
bool operator> (const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs) { return lhs.i > rhs.i; }
bool operator<=(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs) { return lhs.i <= rhs.i; }
bool operator>=(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs) { return lhs.i >= rhs.i; }
bool operator==(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs) { return lhs.i == rhs.i; }
bool operator!=(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs) { return lhs.i != rhs.i; }