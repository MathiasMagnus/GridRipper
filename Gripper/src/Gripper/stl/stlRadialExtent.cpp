#include <Gripper/stl/stlRadialExtent.hpp>


////////////////////////////////////
//                                //
// Multipole::stl::Radial::Extent //
//                                //
////////////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::stl::Radial::Extent::Extent() : m_initial(0), m_final(0) { /*CALLING*/ }


Multipole::stl::Radial::Extent::Extent(const Extent& in) : m_initial(in.m_initial), m_final(in.m_final) { /*CALLING*/ }


Multipole::stl::Radial::Extent::Extent(Extent&& src) : m_initial(std::move(src.m_initial)), m_final(std::move(src.m_final)) { /*CALLING*/ }


Multipole::stl::Radial::Extent::Extent(const index_type& initial, const index_type& final) : m_initial(initial), m_final(final) { /*CALLING*/ }


Multipole::stl::Radial::Extent::~Extent() { /*CALLING*/ }


Multipole::stl::Radial::Extent& Multipole::stl::Radial::Extent::operator=(const Extent& rhs) { m_initial = rhs.m_initial; m_final = rhs.m_final; return *this; }

//////////////////////
// member functions //
//////////////////////

const Multipole::stl::Radial::Extent::index_type& Multipole::stl::Radial::Extent::initial() const { return m_initial; }


const Multipole::stl::Radial::Extent::index_type& Multipole::stl::Radial::Extent::final() const { return m_final; }


bool Multipole::stl::Radial::Extent::contains(const index_type& index) const { return (index >= m_initial) && (index < m_final); }

//////////////////////////
// non-member operators //
//////////////////////////

// Binary
bool operator==(const Multipole::stl::Radial::Extent& lhs, const Multipole::stl::Radial::Extent& rhs) { return (lhs.final() == rhs.final()) && (lhs.initial() == rhs.initial()); }
bool operator!=(const Multipole::stl::Radial::Extent& lhs, const Multipole::stl::Radial::Extent& rhs) { return !((lhs.final() == rhs.final()) && (lhs.initial() == rhs.initial())); }