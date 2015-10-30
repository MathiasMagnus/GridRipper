#include <Gripper/stl/stlGauntMatrix.hpp>


///////////////////////////////////
//                               //
// Multipole::stl::Gaunt::Matrix //
//                               //
///////////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::stl::Gaunt::Matrix::Matrix() {}


Multipole::stl::Gaunt::Matrix::Matrix(const Matrix& in) : m_L_max(in.m_L_max), m_data(in.m_data), m_markers(in.m_markers) {}


Multipole::stl::Gaunt::Matrix::Matrix(Matrix&& in) : m_L_max(in.m_L_max), m_data(std::move(in.m_data)), m_markers(std::move(in.m_markers)) {}


Multipole::stl::Gaunt::Matrix::~Matrix() {}


Multipole::stl::Gaunt::Matrix::Matrix(const Multipole::Gaunt::Matrix& in) : m_L_max(in.getL().i)
{
    Gripper::clog << Gripper::LogLevel::DEBUG << "Entering " << __FUNCTION__;

    for (auto& coefficient : in) m_data.push_back(coefficient);

    for (Index i = 0; i < Index(Spherical::IndexPair(m_L_max, m_L_max))+1; ++i)
        m_markers.push_back(std::make_pair(std::find_if(m_data.begin(), m_data.end(), [=](Coefficient& coeff) {return coeff.i1 == i;}), std::find_if(m_data.begin(), m_data.end(), [=](Coefficient& coeff) {return coeff.i1 == i+1;})));
    /*
    for (Index i = 0; i <= Spherical::IndexPair(m_L_max, m_L_max); ++i)
    {
        auto first = std::find_if(m_data.begin(), m_data.end(), [&](value_type& coeff){ return coeff.i1 == i; });
        
        if (first != m_data.end())
        {
            auto last = std::find_if(first, m_data.end(), [&](value_type& coeff){ return coeff.i1 == i+1; });

            m_markers.push_back(std::make_pair(first, last));
        }
        else Gripper::clog << Gripper::DEBUG << "Combined i1 = " << i.i << " does not have coefficients";
    }
    */  
    Gripper::clog << Gripper::LogLevel::DEBUG << "Leaving " << __FUNCTION__;
}

//////////////////////
// member functions //
//////////////////////

std::pair<Multipole::stl::Gaunt::Matrix::iterator_type, Multipole::stl::Gaunt::Matrix::iterator_type>& Multipole::stl::Gaunt::Matrix::marker(const Index& index) { return m_markers.at(index.i); }

//////////////////////
// member operators //
//////////////////////

Multipole::stl::Gaunt::Matrix& Multipole::stl::Gaunt::Matrix::operator=(Matrix& in) { m_data = in.m_data; m_L_max = in.m_L_max; m_markers = in.m_markers; return *this; }
Multipole::stl::Gaunt::Matrix& Multipole::stl::Gaunt::Matrix::operator=(Matrix&& in) { m_data = std::move(in.m_data); m_L_max = in.m_L_max; m_markers = std::move(in.m_markers); return *this; }

////////////////////////
// iterator interface //
////////////////////////

Multipole::stl::Gaunt::Matrix::iterator_type        Multipole::stl::Gaunt::Matrix::begin() { return m_data.begin(); }
Multipole::stl::Gaunt::Matrix::const_iterator_type  Multipole::stl::Gaunt::Matrix::begin() const { return m_data.begin(); }
Multipole::stl::Gaunt::Matrix::const_iterator_type  Multipole::stl::Gaunt::Matrix::cbegin() const { return m_data.cbegin(); }

Multipole::stl::Gaunt::Matrix::iterator_type        Multipole::stl::Gaunt::Matrix::end() { return m_data.end(); }
Multipole::stl::Gaunt::Matrix::const_iterator_type  Multipole::stl::Gaunt::Matrix::end() const { return m_data.end(); }
Multipole::stl::Gaunt::Matrix::const_iterator_type  Multipole::stl::Gaunt::Matrix::cend() const { return m_data.cend(); }

Multipole::stl::Gaunt::Matrix::reverse_iterator_type        Multipole::stl::Gaunt::Matrix::rbegin() { return m_data.rbegin(); }
Multipole::stl::Gaunt::Matrix::const_reverse_iterator_type  Multipole::stl::Gaunt::Matrix::rbegin() const { return m_data.rbegin(); }
Multipole::stl::Gaunt::Matrix::const_reverse_iterator_type  Multipole::stl::Gaunt::Matrix::crbegin() const { return m_data.crbegin(); }

Multipole::stl::Gaunt::Matrix::reverse_iterator_type        Multipole::stl::Gaunt::Matrix::rend() { return m_data.rend(); }
Multipole::stl::Gaunt::Matrix::const_reverse_iterator_type  Multipole::stl::Gaunt::Matrix::rend() const { return m_data.rend(); }
Multipole::stl::Gaunt::Matrix::const_reverse_iterator_type  Multipole::stl::Gaunt::Matrix::crend() const { return m_data.crend(); }

Multipole::stl::Gaunt::Matrix::value_type           Multipole::stl::Gaunt::Matrix::at(size_type pos) { return m_data.at(pos); }
const Multipole::stl::Gaunt::Matrix::value_type&    Multipole::stl::Gaunt::Matrix::at(size_type pos) const { return m_data.at(pos); }
Multipole::stl::Gaunt::Matrix::value_type           Multipole::stl::Gaunt::Matrix::operator[](size_type pos) { return m_data[pos]; }
const Multipole::stl::Gaunt::Matrix::value_type&    Multipole::stl::Gaunt::Matrix::operator[](size_type pos) const { return m_data[pos]; }

Multipole::stl::Gaunt::Matrix::value_type*          Multipole::stl::Gaunt::Matrix::data() { return m_data.data(); }
const Multipole::stl::Gaunt::Matrix::value_type*    Multipole::stl::Gaunt::Matrix::data() const { return m_data.data(); }

Multipole::stl::Gaunt::Matrix::size_type     Multipole::stl::Gaunt::Matrix::size() const { return m_data.size(); }
void Multipole::stl::Gaunt::Matrix::clear() { m_data.clear(); }

namespace Gripper
{
    namespace stl
    {
        EXPORT Multipole::stl::Gaunt::Matrix gaunt;
    } // namespace stl

} // namespace Gripper