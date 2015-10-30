#include <Gripper/stl/stlExpansionIndex.hpp>


//////////////////////////////////////
//                                  //
// Multipole::stl::Expansion::Index //
//                                  //
//////////////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::stl::Expansion::Index::Index() {}


Multipole::stl::Expansion::Index::Index(const Index& in) : m_radial(in.m_radial), m_spherical(in.m_spherical) {}


Multipole::stl::Expansion::Index::Index(Index&& src) : m_radial(std::move(src.m_radial)), m_spherical(std::move(src.m_spherical)) {}


Multipole::stl::Expansion::Index::Index(const radial_index_type& radial, const spherical_index_type& spherical) : m_radial(radial), m_spherical(spherical) {}


Multipole::stl::Expansion::Index::~Index() {}

//////////////////////
// member functions //
//////////////////////

const Multipole::stl::Expansion::Index::radial_index_type& Multipole::stl::Expansion::Index::radial() const { return m_radial; }


const Multipole::stl::Expansion::Index::spherical_index_type& Multipole::stl::Expansion::Index::spherical() const { return m_spherical; }