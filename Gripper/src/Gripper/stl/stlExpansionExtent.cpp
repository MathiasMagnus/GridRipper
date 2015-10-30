#include <Gripper/stl/stlExpansionExtent.hpp>


////////////////////////////
//                        //
// Multipole::stl::Extent //
//                        //
////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::stl::Expansion::Extent::Extent() { /*CALLING*/ }


Multipole::stl::Expansion::Extent::Extent(const Extent& in) : m_radial(in.m_radial), m_spherical(in.m_spherical) { /*CALLING*/ }


Multipole::stl::Expansion::Extent::Extent(Extent&& src) : m_radial(std::move(src.m_radial)), m_spherical(std::move(src.m_spherical)) { /*CALLING*/ }


Multipole::stl::Expansion::Extent::Extent(const Radial::Extent& radial, const Spherical::Extent& spherical) : m_radial(radial), m_spherical(spherical) { /*CALLING*/ }


Multipole::stl::Expansion::Extent::~Extent() { /*CALLING*/ }


Multipole::stl::Expansion::Extent& Multipole::stl::Expansion::Extent::operator=(const Extent& rhs) { m_radial = rhs.m_radial; m_spherical = rhs.m_spherical; return *this; }

//////////////////////
// member functions //
//////////////////////

const Multipole::stl::Expansion::Extent::radial_extent_type& Multipole::stl::Expansion::Extent::radial() const { return m_radial; }


const Multipole::stl::Expansion::Extent::spherical_extent_type& Multipole::stl::Expansion::Extent::spherical() const { return m_spherical; }


bool Multipole::stl::Expansion::Extent::contains(const index_type& i) const { return m_radial.contains(i.radial()) && m_spherical.contains(i.spherical()); }

//////////////////////////
// non-member operators //
//////////////////////////

// Binary
bool operator==(const Multipole::stl::Expansion::Extent& lhs, const Multipole::stl::Expansion::Extent& rhs) { return (lhs.radial() == rhs.radial()) && (lhs.spherical() == rhs.spherical()); }
bool operator!=(const Multipole::stl::Expansion::Extent& lhs, const Multipole::stl::Expansion::Extent& rhs) { return !((lhs.radial() == rhs.radial()) && (lhs.spherical() == rhs.spherical())); }