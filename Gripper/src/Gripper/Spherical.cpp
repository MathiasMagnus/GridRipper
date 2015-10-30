#include <Gripper/MultipoleTypes.hpp>


/////////////////////////////////
//                             //
// Multipole::Spherical::Index //
//                             //
/////////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::Spherical::Index::Index() {}


Multipole::Spherical::Index::Index(const Multipole::Spherical::Index& in) : i(in.i) {}


Multipole::Spherical::Index::Index(const Multipole::Spherical::IndexType& in) : i(in) {}


Multipole::Spherical::Index::~Index() {}


Multipole::Spherical::Index::operator const Multipole::Spherical::IndexType&() { return i; }

//////////////////////
// member operators //
//////////////////////

// Unary
Multipole::Spherical::Index Multipole::Spherical::Index::operator++(int rhs)   { Index result(i); i += 1; return result; }
Multipole::Spherical::Index& Multipole::Spherical::Index::operator++()         { i += 1; return *this; }
Multipole::Spherical::Index Multipole::Spherical::Index::operator--(int rhs)   { Index result(i); i -= 1; return result; }
Multipole::Spherical::Index& Multipole::Spherical::Index::operator--()         { i -= 1; return *this; }

// Binary
Multipole::Spherical::Index& Multipole::Spherical::Index::operator+=(const Index& rhs) { i += rhs.i; return *this; }
Multipole::Spherical::Index& Multipole::Spherical::Index::operator-=(const Index& rhs) { i -= rhs.i; return *this; }
Multipole::Spherical::Index& Multipole::Spherical::Index::operator*=(const Index& rhs) { i *= rhs.i; return *this; }
Multipole::Spherical::Index& Multipole::Spherical::Index::operator/=(const Index& rhs) { i /= rhs.i; return *this; }

//////////////////////////
// non-member operators //
//////////////////////////

// Unary 
Multipole::Spherical::Index operator+(const Multipole::Spherical::Index& rhs) { return Multipole::Spherical::Index(rhs); }
Multipole::Spherical::Index operator-(const Multipole::Spherical::Index& rhs) { return Multipole::Spherical::Index(-rhs.i); }

// Binary
Multipole::Spherical::Index operator+(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs) { return Multipole::Spherical::Index(lhs.i + rhs.i); }
Multipole::Spherical::Index operator-(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs) { return Multipole::Spherical::Index(lhs.i - rhs.i); }
Multipole::Spherical::Index operator*(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs) { return Multipole::Spherical::Index(lhs.i * rhs.i); }
Multipole::Spherical::Index operator/(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs) { return Multipole::Spherical::Index(lhs.i / rhs.i); }
Multipole::Spherical::Index operator%(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs) { return Multipole::Spherical::Index(lhs.i % rhs.i); }

bool operator< (const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs) { return lhs.i < rhs.i; }
bool operator> (const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs) { return lhs.i > rhs.i; }
bool operator<=(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs) { return lhs.i <= rhs.i; }
bool operator>=(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs) { return lhs.i >= rhs.i; }
bool operator==(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs) { return lhs.i == rhs.i; }
bool operator!=(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs) { return lhs.i != rhs.i; }


/////////////////////////////////////
//                                 //
// Multipole::Spherical::IndexPair //
//                                 //
/////////////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::Spherical::IndexPair::IndexPair() {}


Multipole::Spherical::IndexPair::IndexPair(const IndexPair& in) : l(in.l), m(in.m) {}


Multipole::Spherical::IndexPair::IndexPair(const Index& l_in, const Index& m_in) : l(l_in), m(m_in) {}


Multipole::Spherical::IndexPair::IndexPair(const Gaunt::Index& i_in) { l = static_cast<Multipole::Spherical::IndexType>(std::round((-1. + std::sqrt(1. + 4 * i_in.i)) / 2.)); m = i_in.i - l*(l + 1); }


Multipole::Spherical::IndexPair::~IndexPair() {}

//////////////////////////
// non-member operators //
//////////////////////////

// Binary
bool operator==(const Multipole::Spherical::IndexPair& lhs, const Multipole::Spherical::IndexPair& rhs) { return (lhs.l == rhs.l) && (lhs.m == rhs.m); }
bool operator!=(const Multipole::Spherical::IndexPair& lhs, const Multipole::Spherical::IndexPair& rhs) { return (lhs.l != rhs.l) || (lhs.m != rhs.m); }