#include <Gripper/stl/stlRadialIndex.hpp>


///////////////////////////////////
//                               //
// Multipole::stl::Radial::Index //
//                               //
///////////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::stl::Radial::Index::Index() {}


Multipole::stl::Radial::Index::Index(const Index& in) : r(in.r) {}


Multipole::stl::Radial::Index::Index(Index&& src) : r(std::move(src.r)) {}


Multipole::stl::Radial::Index::~Index() {}


Multipole::stl::Radial::Index& Multipole::stl::Radial::Index::operator=(const Index& rhs) { r = rhs.r; return *this; }


Multipole::stl::Radial::Index::Index(const value_type& in) : r(in) {}


Multipole::stl::Radial::Index::operator value_type&() { return r; }

//////////////////////
// member operators //
//////////////////////

// Unary
Multipole::stl::Radial::Index Multipole::stl::Radial::Index::operator++(int rhs)   { Index result(r); r += 1; return result; }
Multipole::stl::Radial::Index& Multipole::stl::Radial::Index::operator++()         { r += 1; return *this; }
Multipole::stl::Radial::Index Multipole::stl::Radial::Index::operator--(int rhs)   { Index result(r); r -= 1; return result; }
Multipole::stl::Radial::Index& Multipole::stl::Radial::Index::operator--()         { r -= 1; return *this; }

// Binary
Multipole::stl::Radial::Index& Multipole::stl::Radial::Index::operator+=(const Index& rhs) { r += rhs.r; return *this; }
Multipole::stl::Radial::Index& Multipole::stl::Radial::Index::operator-=(const Index& rhs) { r *= rhs.r; return *this; }
Multipole::stl::Radial::Index& Multipole::stl::Radial::Index::operator*=(const Index& rhs) { r *= rhs.r; return *this; }
Multipole::stl::Radial::Index& Multipole::stl::Radial::Index::operator/=(const Index& rhs) { r /= rhs.r; return *this; }

//////////////////////////
// non-member operators //
//////////////////////////

// Unary 
Multipole::stl::Radial::Index operator+(const Multipole::stl::Radial::Index& rhs) { return Multipole::stl::Radial::Index(rhs); }
Multipole::stl::Radial::Index operator-(const Multipole::stl::Radial::Index& rhs) { return Multipole::stl::Radial::Index(-rhs.r); }
double Multipole::stl::Radial::pown(const Index& base, int power) { return std::pow(base.r, power); }

// Binary
Multipole::stl::Radial::Index operator+(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs) { return Multipole::stl::Radial::Index(lhs.r + rhs.r); }
Multipole::stl::Radial::Index operator-(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs) { return Multipole::stl::Radial::Index(lhs.r - rhs.r); }
Multipole::stl::Radial::Index operator*(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs) { return Multipole::stl::Radial::Index(lhs.r * rhs.r); }
Multipole::stl::Radial::Index operator/(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs) { return Multipole::stl::Radial::Index(lhs.r / rhs.r); }
Multipole::stl::Radial::Index operator%(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs) { return Multipole::stl::Radial::Index(lhs.r % rhs.r); }

bool operator< (const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs) { return lhs.r < rhs.r; }
bool operator>(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs) { return lhs.r > rhs.r; }
bool operator<=(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs) { return lhs.r <= rhs.r; }
bool operator>=(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs) { return lhs.r >= rhs.r; }
bool operator==(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs) { return lhs.r == rhs.r; }
bool operator!=(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs) { return lhs.r != rhs.r; }