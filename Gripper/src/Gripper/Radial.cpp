#include <Gripper/MultipoleTypes.hpp>


//////////////////////////////
//                          //
// Multipole::Radial::Index //
//                          //
//////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::Radial::Index::Index() {}


Multipole::Radial::Index::Index(const Multipole::Radial::Index& in) : r(in.r) {}


Multipole::Radial::Index::Index(const Multipole::Radial::IndexType& in) : r(in) {}


Multipole::Radial::Index::~Index() {}


Multipole::Radial::Index::operator const Multipole::Radial::IndexType&() { return r; }

//////////////////////
// member operators //
//////////////////////

// Unary
Multipole::Radial::Index Multipole::Radial::Index::operator++(int rhs)   { Index result(r); r += 1; return result; }
Multipole::Radial::Index& Multipole::Radial::Index::operator++()         { r += 1; return *this; }
Multipole::Radial::Index Multipole::Radial::Index::operator--(int rhs)   { Index result(r); r -= 1; return result; }
Multipole::Radial::Index& Multipole::Radial::Index::operator--()         { r -= 1; return *this; }

// Binary
Multipole::Radial::Index& Multipole::Radial::Index::operator+=(const Index& rhs) { r += rhs.r; return *this; }
Multipole::Radial::Index& Multipole::Radial::Index::operator-=(const Index& rhs) { r *= rhs.r; return *this; }
Multipole::Radial::Index& Multipole::Radial::Index::operator*=(const Index& rhs) { r *= rhs.r; return *this; }
Multipole::Radial::Index& Multipole::Radial::Index::operator/=(const Index& rhs) { r /= rhs.r; return *this; }

//////////////////////////
// non-member operators //
//////////////////////////

// Unary 
Multipole::Radial::Index operator+(const Multipole::Radial::Index& rhs) { return Multipole::Radial::Index(rhs); }
Multipole::Radial::Index operator-(const Multipole::Radial::Index& rhs) { return Multipole::Radial::Index(-rhs.r); }

// Binary
Multipole::Radial::Index operator+(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs) { return Multipole::Radial::Index(lhs.r + rhs.r); }
Multipole::Radial::Index operator-(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs) { return Multipole::Radial::Index(lhs.r - rhs.r); }
Multipole::Radial::Index operator*(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs) { return Multipole::Radial::Index(lhs.r * rhs.r); }
Multipole::Radial::Index operator/(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs) { return Multipole::Radial::Index(lhs.r / rhs.r); }
Multipole::Radial::Index operator%(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs) { return Multipole::Radial::Index(lhs.r % rhs.r); }

bool operator< (const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs) { return lhs.r < rhs.r; }
bool operator> (const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs) { return lhs.r > rhs.r; }
bool operator<=(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs) { return lhs.r <= rhs.r; }
bool operator>=(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs) { return lhs.r >= rhs.r; }
bool operator==(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs) { return lhs.r == rhs.r; }
bool operator!=(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs) { return lhs.r != rhs.r; }