#include <Gripper/stl/stlSphericalExtent.hpp>


///////////////////////////////////////
//                                   //
// Multipole::stl::Spherical::Extent //
//                                   //
///////////////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::stl::Spherical::Extent::Extent() : m_initial(0), m_final(0) { /*CALLING*/ }


Multipole::stl::Spherical::Extent::Extent(const Extent& in) : m_initial(in.m_initial), m_final(in.m_final) { /*CALLING*/ }


Multipole::stl::Spherical::Extent::Extent(Extent&& src) : m_initial(std::move(src.m_initial)), m_final(std::move(src.m_final)) { /*CALLING*/ }


Multipole::stl::Spherical::Extent::Extent(const index_type& initial, const index_type& final) : m_initial(initial), m_final(final) { /*CALLING*/ }


Multipole::stl::Spherical::Extent::Extent(const Index& initial, const Index& final) : m_initial(IndexPair(initial, -initial)), m_final(index_type(IndexPair(final, final)) + 1) { /*CALLING*/ }


Multipole::stl::Spherical::Extent::~Extent() { /*CALLING*/ }


Multipole::stl::Spherical::Extent& Multipole::stl::Spherical::Extent::operator=(const Extent& rhs) { m_initial = rhs.m_initial; m_final = rhs.m_final; return *this; }

//////////////////////
// member functions //
//////////////////////

const Multipole::stl::Spherical::Extent::index_type& Multipole::stl::Spherical::Extent::initial() const { return m_initial; }


const Multipole::stl::Spherical::Extent::index_type& Multipole::stl::Spherical::Extent::final() const { return m_final; }


bool Multipole::stl::Spherical::Extent::contains(const index_type& index) const { return (index >= m_initial) && (index < m_final) ? true : false; }

//////////////////////////
// non-member operators //
//////////////////////////

// Binary
bool operator==(const Multipole::stl::Spherical::Extent& lhs, const Multipole::stl::Spherical::Extent& rhs) { return (lhs.final() == rhs.final()) && (lhs.initial() == rhs.initial()); }
bool operator!=(const Multipole::stl::Spherical::Extent& lhs, const Multipole::stl::Spherical::Extent& rhs) { return !((lhs.final() == rhs.final()) && (lhs.initial() == rhs.initial())); }