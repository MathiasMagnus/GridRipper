#include <Gripper/MultipoleTypes.hpp>
/*

/////////////////////////////
//                         //
// Multipole::Gaunt::Index //
//                         //
/////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::Gaunt::Index::Index() {}


Multipole::Gaunt::Index::Index(const Index& in) : i(in.i) {}


Multipole::Gaunt::Index::Index(const IndexType& i_in) : i(i_in) {}


Multipole::Gaunt::Index::Index(const Spherical::IndexPair& lm_in) : i(static_cast<Spherical::IndexType>(lm_in.l * (lm_in.l + 1) + lm_in.m)) {}


Multipole::Gaunt::Index::~Index() {}


Multipole::Gaunt::Index::operator const Multipole::Gaunt::IndexType() { return i; }

//////////////////////
// member operators //
//////////////////////

// Unary
Multipole::Gaunt::Index Multipole::Gaunt::Index::operator++(int rhs)   { Index result(i); i += 1; return result; }
Multipole::Gaunt::Index& Multipole::Gaunt::Index::operator++()         { i += 1; return *this; }
Multipole::Gaunt::Index Multipole::Gaunt::Index::operator--(int rhs)   { Index result(i); i -= 1; return result; }
Multipole::Gaunt::Index& Multipole::Gaunt::Index::operator--()         { i -= 1; return *this; }

// Binary
Multipole::Gaunt::Index& Multipole::Gaunt::Index::operator+=(const Index& rhs) { i += rhs.i; return *this; }
Multipole::Gaunt::Index& Multipole::Gaunt::Index::operator-=(const Index& rhs) { i -= rhs.i; return *this; }
Multipole::Gaunt::Index& Multipole::Gaunt::Index::operator*=(const Index& rhs) { i *= rhs.i; return *this; }
Multipole::Gaunt::Index& Multipole::Gaunt::Index::operator/=(const Index& rhs) { i /= rhs.i; return *this; }

//////////////////////////
// non-member operators //
//////////////////////////

// Binary
Multipole::Gaunt::Index operator+(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs) { return Multipole::Gaunt::Index(lhs.i + rhs.i); }
Multipole::Gaunt::Index operator-(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs) { return Multipole::Gaunt::Index(lhs.i - rhs.i); }
Multipole::Gaunt::Index operator*(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs) { return Multipole::Gaunt::Index(lhs.i * rhs.i); }
Multipole::Gaunt::Index operator/(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs) { return Multipole::Gaunt::Index(lhs.i / rhs.i); }
Multipole::Gaunt::Index operator%(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs) { return Multipole::Gaunt::Index(lhs.i % rhs.i); }

bool operator< (const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs) { return lhs.i < rhs.i; }
bool operator> (const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs) { return lhs.i > rhs.i; }
bool operator<=(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs) { return lhs.i <= rhs.i; }
bool operator>=(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs) { return lhs.i >= rhs.i; }
bool operator==(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs) { return lhs.i == rhs.i; }
bool operator!=(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs) { return lhs.i != rhs.i; }


///////////////////////////////////
//                               //
// Multipole::Gaunt::Coefficient //
//                               //
///////////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::Gaunt::Coefficient::Coefficient() {}


Multipole::Gaunt::Coefficient::Coefficient(const Coefficient& in) : l1m1(in.l1m1), l2m2(in.l2m2), l3m3(in.l3m3), value(in.value) {}


Multipole::Gaunt::Coefficient::Coefficient(const Index& i1_in, const Index& i2_in, const Index& i3_in, const ValueType& value_in) : l1m1(i1_in), l2m2(i2_in), l3m3(i3_in), value(value_in) {}


Multipole::Gaunt::Coefficient::Coefficient(const Spherical::IndexPair& l1m1_in, const Spherical::IndexPair& l2m2_in, const Spherical::IndexPair& l3m3_in, const ValueType& value_in) : l1m1(l1m1_in), l2m2(l2m2_in), l3m3(l3m3_in), value(value_in) {}


Multipole::Gaunt::Coefficient::~Coefficient() {}

//////////////////////////
// non-member operators //
//////////////////////////

bool operator< (const Multipole::Gaunt::Coefficient& lhs, const Multipole::Gaunt::Coefficient& rhs) { return lhs.value < rhs.value; }
bool operator> (const Multipole::Gaunt::Coefficient& lhs, const Multipole::Gaunt::Coefficient& rhs) { return lhs.value > rhs.value; }
bool operator<=(const Multipole::Gaunt::Coefficient& lhs, const Multipole::Gaunt::Coefficient& rhs) { return lhs.value <= rhs.value; }
bool operator>=(const Multipole::Gaunt::Coefficient& lhs, const Multipole::Gaunt::Coefficient& rhs) { return lhs.value >= rhs.value; }
bool operator==(const Multipole::Gaunt::Coefficient& lhs, const Multipole::Gaunt::Coefficient& rhs) { return lhs.value == rhs.value; }
bool operator!=(const Multipole::Gaunt::Coefficient& lhs, const Multipole::Gaunt::Coefficient& rhs) { return lhs.value != rhs.value; }


//////////////////////////////
//                          //
// Multipole::Gaunt::Matrix //
//                          //
//////////////////////////////

/////////////////////////////
// constructors/destructor //
/////////////////////////////

Multipole::Gaunt::Matrix::Matrix() {}


Multipole::Gaunt::Matrix::Matrix(const Matrix& in) : m_data(in.m_data), m_L_max(in.m_L_max) {}


Multipole::Gaunt::Matrix::Matrix(Matrix&& in) : m_data(std::move(in.m_data)), m_L_max(in.m_L_max) {}


Multipole::Gaunt::Matrix::Matrix(Spherical::Index& L_max_in) { m_L_max = L_max_in; compute(); }


Multipole::Gaunt::Matrix::~Matrix() {}

//////////////////////
// member operators //
//////////////////////

Multipole::Gaunt::Matrix& Multipole::Gaunt::Matrix::operator=(Matrix& in) { m_data = in.m_data; m_L_max = in.m_L_max; return *this; }
Multipole::Gaunt::Matrix& Multipole::Gaunt::Matrix::operator=(Matrix&& in) { m_data = std::move(in.m_data); m_L_max = in.m_L_max; return *this; }

////////////////////////
// iterator interface //
////////////////////////

Multipole::Gaunt::Matrix::iterator_type        Multipole::Gaunt::Matrix::begin() { return m_data.begin(); }
Multipole::Gaunt::Matrix::const_iterator_type  Multipole::Gaunt::Matrix::begin() const { return m_data.begin(); }
Multipole::Gaunt::Matrix::const_iterator_type  Multipole::Gaunt::Matrix::cbegin() const { return m_data.cbegin(); }

Multipole::Gaunt::Matrix::iterator_type        Multipole::Gaunt::Matrix::end() { return m_data.end(); }
Multipole::Gaunt::Matrix::const_iterator_type  Multipole::Gaunt::Matrix::end() const { return m_data.end(); }
Multipole::Gaunt::Matrix::const_iterator_type  Multipole::Gaunt::Matrix::cend() const { return m_data.cend(); }

Multipole::Gaunt::Matrix::reverse_iterator_type        Multipole::Gaunt::Matrix::rbegin() { return m_data.rbegin(); }
Multipole::Gaunt::Matrix::const_reverse_iterator_type  Multipole::Gaunt::Matrix::rbegin() const { return m_data.rbegin(); }
Multipole::Gaunt::Matrix::const_reverse_iterator_type  Multipole::Gaunt::Matrix::crbegin() const { return m_data.crbegin(); }

Multipole::Gaunt::Matrix::reverse_iterator_type        Multipole::Gaunt::Matrix::rend() { return m_data.rend(); }
Multipole::Gaunt::Matrix::const_reverse_iterator_type  Multipole::Gaunt::Matrix::rend() const { return m_data.rend(); }
Multipole::Gaunt::Matrix::const_reverse_iterator_type  Multipole::Gaunt::Matrix::crend() const { return m_data.crend(); }

Multipole::Gaunt::Coefficient           Multipole::Gaunt::Matrix::at(size_type pos) { return m_data.at(pos); }
const Multipole::Gaunt::Coefficient&    Multipole::Gaunt::Matrix::at(size_type pos) const { return m_data.at(pos); }
Multipole::Gaunt::Coefficient           Multipole::Gaunt::Matrix::operator[](size_type pos) { return m_data[pos]; }
const Multipole::Gaunt::Coefficient&    Multipole::Gaunt::Matrix::operator[](size_type pos) const { return m_data[pos]; }

Multipole::Gaunt::Coefficient*          Multipole::Gaunt::Matrix::data() { return m_data.data(); }
const Multipole::Gaunt::Coefficient*    Multipole::Gaunt::Matrix::data() const { return m_data.data(); }

Multipole::Gaunt::Matrix::size_type     Multipole::Gaunt::Matrix::size() const { return m_data.size(); }
void Multipole::Gaunt::Matrix::clear() { m_data.clear(); }

Multipole::Spherical::Index Multipole::Gaunt::Matrix::getL() const { return m_L_max; }

//////////////////////////////
// private member functions //
//////////////////////////////

void Multipole::Gaunt::Matrix::compute()
{
    using namespace Gripper;

    ENTERING

    clog << LogLevel::INFO << "Computing Gaunt coefficients up to L_max = " << m_L_max.i;

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<std::future<void>> jobs(m_L_max.i+1);
    std::vector<container_type> temp(m_L_max.i+1);

    for (Spherical::Index i = 0; i <= m_L_max; ++i) jobs.at(i.i) = std::async(std::launch::async, &Multipole::Gaunt::Matrix::computeL, this, i, &temp);
    for (auto& job : jobs) job.wait();

    clog << LogLevel::DEBUG << "Merging results";

    for (auto& tmp : temp) m_data.insert(m_data.end(), tmp.begin(), tmp.end());

    auto end = std::chrono::high_resolution_clock::now();

    clog << LogLevel::TIMING << "Computation took " << std::chrono::duration_cast<std::chrono::seconds>(end.time_since_epoch() - start.time_since_epoch()).count() << " seconds.";

    LEAVING
};


void Multipole::Gaunt::Matrix::computeL(Spherical::Index L1, std::vector<container_type>* result)
{
    using namespace Gripper;

    ENTERING

    clog << LogLevel::ALL << "Compute thread L1 = " << L1.i << " (id = " << std::this_thread::get_id() << ") started.";

    // Compute all coefficients with designated L1 value

    // NOTE 1: Running the index M1 to only 0 instead of L1 is due to the symmetry: Gl1m1l2m2l3m3=Gl1-m1l2-m2l3-m3
    //         This requires that every time an index structure is loaded, an extra multiplication need be made with negated m indices.
    //
    // NOTE 2: This optimization is a limited consequence of the symmetry: Gi1i2i3=G[i1i2i3] in other words Gaunt-coefficients for all
    //         permutations of indices. This means that when an index is found an extra multiplication needs to be made switched
    //         switched i2i3 indices.

    //	for(int M1 = -L1 ; M1 <= 0 ; ++M1) // NOTE 1! Symmetry optimization
    for (Spherical::Index M1 = -L1; M1 <= L1; ++M1)
    {
        for (Spherical::Index L2 = 0; L2 <= m_L_max; ++L2)
        {
            for (Spherical::Index M2 = -L2; M2 <= L2; ++M2)
            {
                for (Spherical::Index L3 = 0; L3 <= m_L_max; ++L3)
                {
                    for (Spherical::Index M3 = -L3; M3 <= L3; ++M3)
                    {
                        // Shortcuts
                        if ((L1>(L2 + L3)) || (L2>(L1 + L3)) || (L3>(L1 + L2))) continue;
                        if (static_cast<Spherical::IndexType>((L1 + L2 + L3) % 2)) continue;
                        if ((M1 + M2 + M3) != 0) continue;

                        // if(gauntIndexer(L3,M3)<gauntIndexer(L2,M2)) continue; // NOTE 2! Symmetry optimization

                        Spherical::IndexPair L1M1(L1 * 2, M1 * 2);
                        Spherical::IndexPair L2M2(L2 * 2, M2 * 2);
                        Spherical::IndexPair L3M3(L3 * 2, M3 * 2);
                        Spherical::IndexPair L1M0(L1 * 2, 0);
                        Spherical::IndexPair L2M0(L2 * 2, 0);
                        Spherical::IndexPair L3M0(L3 * 2, 0);
                        */
                        /*
                        ValueType temp1 = wigner3jSymbol(L1M1, L2M2, L3M3);
                        ValueType temp2 = wigner3jSymbol(L1M0, L2M0, L3M0);
                        */
/*
                        ValueType temp1 = gsl_sf_coupling_3j(L1M1.l.i, L2M2.l.i, L3M3.l.i, L1M1.m.i, L2M2.m.i, L3M3.m.i);
                        ValueType temp2 = gsl_sf_coupling_3j(L1M0.l.i, L2M0.l.i, L3M0.l.i, 0, 0, 0);
                        
                        ValueType res = std::sqrt(static_cast<Spherical::IndexType>((2 * L1 + 1) * (2 * L2 + 1) * (2 * L3 + 1))) * temp1 * temp2 / (2.0 * std::sqrt(M_PI));

                        if (std::fabs(res) > 1e-6)
                            result->at(L1.i).push_back(Coefficient(Spherical::IndexPair(L1, M1), Spherical::IndexPair(L2, M2), Spherical::IndexPair(L3, M3), res));
                    }
                }
            }
        }
    }

    clog << LogLevel::ALL << "Compute thread L = " << L1.i << " finished.";

    LEAVING
}
*/
/*
Multipole::Gaunt::ValueType Multipole::Gaunt::Matrix::wigner3jSymbol(Spherical::IndexPair& L1M1, Spherical::IndexPair& L2M2, Spherical::IndexPair& L3M3)
{
    using namespace Gripper;

    clog << LogLevel::INSANE << "Entering " << __FUNCTION__;

    if (L1M1.l < 0 || L2M2.l < 0 || L3M3.l < 0)
    {
        clog << LogLevel::ERR << __FUNCTION__ << " Gaunt matrix domain error with L1 = " << L1M1.l.i << " L2 = " << L2M2.l.i << " L3 = " << L3M3.l.i;
        return 0.0;
    }
    else if (triangleSelectionFail(L1M1, L2M2, L3M3) || mSelectionFail(L1M1, L2M2, L3M3))
    {
        return 0.0;
    }
    else
    {
        Spherical::Index jca = (-L1M1.l + L2M2.l + L3M3.l) / 2,
            jcb = (L1M1.l - L2M2.l + L3M3.l) / 2,
            jcc = (L1M1.l + L2M2.l - L3M3.l) / 2,
            jmma = (L1M1.l - L1M1.m) / 2,
            jmmb = (L2M2.l - L2M2.m) / 2,
            jmmc = (L3M3.l - L3M3.m) / 2,
            jpma = (L1M1.l + L1M1.m) / 2,
            jpmb = (L2M2.l + L2M2.m) / 2,
            jpmc = (L3M3.l + L3M3.m) / 2,
            jsum = (L1M1.l + L2M2.l + L3M3.l) / 2,
            //kmin = std::max( std::max( Spherical::Index(0), jpmb - jmmc), jmma - jpmc ),
            //kmax = std::min( std::min(jcc, jmma), jpmb ),
            kmin = std::max(std::max(0, jpmb.i - jmmc.i), jmma.i - jpmc.i),
            kmax = std::min(std::min(jcc.i, jmma.i), jpmb.i),
            sign = static_cast<Spherical::IndexType>(kmin - jpma + jmmb) % 2 ? -1 : 1;

        ValueType sum_pos = 0.0, sum_neg = 0.0, norm, term;
        value_with_error bc1, bc2, bc3, bcn1, bcn2, bcd1, bcd2, bcd3, bcd4, result;

        bool overflow = false;

        overflow |= choose(L1M1.l, jcc, bcn1);
        overflow |= choose(L2M2.l, jcc, bcn2);
        overflow |= choose(jsum + 1, jcc, bcd1);
        overflow |= choose(L1M1.l, jmma, bcd2);
        overflow |= choose(L2M2.l, jmmb, bcd3);
        overflow |= choose(L3M3.l, jpmc, bcd4);

        if (overflow != false) clog << LogLevel::ERR << __FUNCTION__ << " Gaunt matrix overflow error";

        norm = std::sqrt(bcn1.first * bcn2.first) / std::sqrt(bcd1.first * bcd2.first * bcd3.first * bcd4.first * (L3M3.l.i + 1.0));

        for (Spherical::Index k = kmin; k <= kmax; k++)
        {
            overflow |= choose(jcc, k, bc1);
            overflow |= choose(jcb, jmma - k, bc2);
            overflow |= choose(jca, jpmb - k, bc3);

            if (overflow != 0) clog << LogLevel::ERR << __FUNCTION__ << " Gaunt matrix overflow error";

            term = bc1.first * bc2.first * bc3.first;

            if (sign < 0) sum_neg += norm * term;
            else sum_pos += norm * term;

            sign = -sign;
        }

        result.first = sum_pos - sum_neg;
        result.second = 2.0 * std::numeric_limits<ValueType>::lowest() * (sum_pos + sum_neg);
        result.second += 2.0 * std::numeric_limits<ValueType>::lowest() * (kmax.i - kmin.i) * std::fabs(result.first);

        return result.first; // for the moment, we do not check the error
    }

    clog << LogLevel::INSANE << "Leaving " << __FUNCTION__;
}


bool Multipole::Gaunt::Matrix::triangleSelectionFail(Spherical::IndexPair& L1M1, Spherical::IndexPair& L2M2, Spherical::IndexPair& L3M3)
{
    return ((L2M2.l < std::abs(L1M1.l.i - L3M3.l.i)) || (L2M2.l > L1M1.l + L3M3.l));
}


bool Multipole::Gaunt::Matrix::mSelectionFail(Spherical::IndexPair& L1M1, Spherical::IndexPair& L2M2, Spherical::IndexPair& L3M3)
{
    return (std::abs(L1M1.m.i) > L1M1.l
        || std::abs(L2M2.m.i) > L2M2.l
        || std::abs(L3M3.m.i) > L3M3.l
        || (L1M1.l.i + L1M1.m.i) % 2
        || (L2M2.l.i + L2M2.m.i) % 2
        || (L3M3.l.i + L3M3.m.i) % 2
        || (L1M1.m + L2M2.m + L3M3.m) != 0 );
}


bool Multipole::Gaunt::Matrix::choose(Spherical::Index& n, Spherical::Index& m, value_with_error& result)
{
    using namespace Gripper;

    if (m > n)
    {
        clog << LogLevel::ERR << __FUNCTION__ << " Gaunt matrix domain error (m > n) with m = " << m.i << " n = " << n.i;
        result.first = 0.0;
        return EXIT_FAILURE;
    }
    else if (m == n || m == 0)
    {
        result.first = 1.0;
        result.second = 0.0;
        return EXIT_SUCCESS;
    }
    else if (n.i <= 297)
    {
        result.first = (factorial(n.i) / factorial(m.i)) / factorial(n.i - m.i);
        result.second = 6.0 * std::numeric_limits<ValueType>::lowest() * std::fabs(result.first);
        return EXIT_SUCCESS;
    }
    else
    {
        if (m * 2 < n) m = n - m;

        if (n - m < 64)  // compute product for a manageable number of terms
        {
            ValueType prod = 1.0;
            Spherical::Index k;

            for (k = n; k >= m + 1; k--)
            {
                ValueType tk = static_cast<ValueType>(k.i) / static_cast<Spherical::IndexType>(k - m);
                if (tk > std::numeric_limits<ValueType>::max() / prod)
                {
                    clog << LogLevel::ERR << __FUNCTION__ << " Gaunt matrix overflow error";
                    return EXIT_FAILURE;
                }
                prod *= tk;
            }
            result.first = prod;
            result.second = 2.0 * std::numeric_limits<ValueType>::lowest() * prod * std::abs(static_cast<Spherical::IndexType>(n - m));
            return EXIT_SUCCESS;
        }
        else
        {
            clog << LogLevel::CRIT << __FUNCTION__ << " We hoped this never happened";
            return EXIT_FAILURE;
        }
    }
}


Multipole::Gaunt::ValueType Multipole::Gaunt::Matrix::factorial(int i)
{
    switch (i)
    {
        case 0:   return static_cast<ValueType>(1.000000000000000000000000000); break;
        case 1:   return static_cast<ValueType>(1.000000000000000000000000000); break;
        case 2:   return static_cast<ValueType>(2.000000000000000000000000000); break;
        case 3:   return static_cast<ValueType>(3.000000000000000000000000000); break;
        case 4:   return static_cast<ValueType>(8.000000000000000000000000000); break;
        case 5:   return static_cast<ValueType>(15.00000000000000000000000000); break;
        case 6:   return static_cast<ValueType>(48.00000000000000000000000000); break;
        case 7:   return static_cast<ValueType>(105.0000000000000000000000000); break;
        case 8:   return static_cast<ValueType>(384.0000000000000000000000000); break;
        case 9:   return static_cast<ValueType>(945.0000000000000000000000000); break;
        case 10:  return static_cast<ValueType>(3840.000000000000000000000000); break;
        case 11:  return static_cast<ValueType>(10395.00000000000000000000000); break;
        case 12:  return static_cast<ValueType>(46080.00000000000000000000000); break;
        case 13:  return static_cast<ValueType>(135135.0000000000000000000000); break;
        case 14:  return static_cast<ValueType>(645120.00000000000000000000000); break;
        case 15:  return static_cast<ValueType>(2.02702500000000000000000000000e6); break;
        case 16:  return static_cast<ValueType>(1.03219200000000000000000000000e7); break;
        case 17:  return static_cast<ValueType>(3.4459425000000000000000000000e7); break;
        case 18:  return static_cast<ValueType>(1.85794560000000000000000000000e8); break;
        case 19:  return static_cast<ValueType>(6.5472907500000000000000000000e8); break;
        case 20:  return static_cast<ValueType>(3.7158912000000000000000000000e9); break;
        case 21:  return static_cast<ValueType>(1.37493105750000000000000000000e10); break;
        case 22:  return static_cast<ValueType>(8.1749606400000000000000000000e10); break;
        case 23:  return static_cast<ValueType>(3.1623414322500000000000000000e11); break;
        case 24:  return static_cast<ValueType>(1.96199055360000000000000000000e12); break;
        case 25:  return static_cast<ValueType>(7.9058535806250000000000000000e12); break;
        case 26:  return static_cast<ValueType>(5.1011754393600000000000000000e13); break;
        case 27:  return static_cast<ValueType>(2.13458046676875000000000000000e14); break;
        case 28:  return static_cast<ValueType>(1.42832912302080000000000000000e15); break;
        case 29:  return static_cast<ValueType>(6.1902833536293750000000000000e15); break;
        case 30:  return static_cast<ValueType>(4.2849873690624000000000000000e16); break;
        case 31:  return static_cast<ValueType>(1.91898783962510625000000000000e17); break;
        case 32:  return static_cast<ValueType>(1.37119595809996800000000000000e18); break;
        case 33:  return static_cast<ValueType>(6.3326598707628506250000000000e18); break;
        case 34:  return static_cast<ValueType>(4.6620662575398912000000000000e19); break;
        case 35:  return static_cast<ValueType>(2.21643095476699771875000000000e20); break;
        case 36:  return static_cast<ValueType>(1.67834385271436083200000000000e21); break;
        case 37:  return static_cast<ValueType>(8.2007945326378915593750000000e21); break;
        case 38:  return static_cast<ValueType>(6.3777066403145711616000000000e22); break;
        case 39:  return static_cast<ValueType>(3.1983098677287777081562500000e23); break;
        case 40:  return static_cast<ValueType>(2.55108265612582846464000000000e24); break;
        case 41:  return static_cast<ValueType>(1.31130704576879886034406250000e25); break;
        case 42:  return static_cast<ValueType>(1.07145471557284795514880000000e26); break;
        case 43:  return static_cast<ValueType>(5.6386202968058350994794687500e26); break;
        case 44:  return static_cast<ValueType>(4.7144007485205310026547200000e27); break;
        case 45:  return static_cast<ValueType>(2.53737913356262579476576093750e28); break;
        case 46:  return static_cast<ValueType>(2.16862434431944426122117120000e29); break;
        case 47:  return static_cast<ValueType>(1.19256819277443412353990764062e30); break;
        case 48:  return static_cast<ValueType>(1.04093968527333324538616217600e31); break;
        case 49:  return static_cast<ValueType>(5.8435841445947272053455474391e31); break;
        case 50:  return static_cast<ValueType>(5.2046984263666662269308108800e32); break;
        case 51:  return static_cast<ValueType>(2.98022791374331087472622919392e33); break;
        case 52:  return static_cast<ValueType>(2.70644318171066643800402165760e34); break;
        case 53:  return static_cast<ValueType>(1.57952079428395476360490147278e35); break;
        case 54:  return static_cast<ValueType>(1.46147931812375987652217169510e36); break;
        case 55:  return static_cast<ValueType>(8.6873643685617511998269581003e36); break;
        case 56:  return static_cast<ValueType>(8.1842841814930553085241614926e37); break;
        case 57:  return static_cast<ValueType>(4.9517976900801981839013661172e38); break;
        case 58:  return static_cast<ValueType>(4.7468848252659720789440136657e39); break;
        case 59:  return static_cast<ValueType>(2.92156063714731692850180600912e40); break;
        case 60:  return static_cast<ValueType>(2.84813089515958324736640819942e41); break;
        case 61:  return static_cast<ValueType>(1.78215198865986332638610166557e42); break;
        case 62:  return static_cast<ValueType>(1.76584115499894161336717308364e43); break;
        case 63:  return static_cast<ValueType>(1.12275575285571389562324404931e44); break;
        case 64:  return static_cast<ValueType>(1.13013833919932263255499077353e45); break;
        case 65:  return static_cast<ValueType>(7.2979123935621403215510863205e45); break;
        case 66:  return static_cast<ValueType>(7.4589130387155293748629391053e46); break;
        case 67:  return static_cast<ValueType>(4.8896013036866340154392278347e47); break;
        case 68:  return static_cast<ValueType>(5.0720608663265599749067985916e48); break;
        case 69:  return static_cast<ValueType>(3.3738248995437774706530672060e49); break;
        case 70:  return static_cast<ValueType>(3.5504426064285919824347590141e50); break;
        case 71:  return static_cast<ValueType>(2.39541567867608200416367771623e51); break;
        case 72:  return static_cast<ValueType>(2.55631867662858622735302649017e52); break;
        case 73:  return static_cast<ValueType>(1.74865344543353986303948473285e53); break;
        case 74:  return static_cast<ValueType>(1.89167582070515380824123960272e54); break;
        case 75:  return static_cast<ValueType>(1.31149008407515489727961354964e55); break;
        case 76:  return static_cast<ValueType>(1.43767362373591689426334209807e56); break;
        case 77:  return static_cast<ValueType>(1.00984736473786927090530243322e57); break;
        case 78:  return static_cast<ValueType>(1.12138542651401517752540683649e58); break;
        case 79:  return static_cast<ValueType>(7.9777941814291672401518892225e58); break;
        case 80:  return static_cast<ValueType>(8.9710834121121214202032546920e59); break;
        case 81:  return static_cast<ValueType>(6.4620132869576254645230302702e60); break;
        case 82:  return static_cast<ValueType>(7.3562883979319395645666688474e61); break;
        case 83:  return static_cast<ValueType>(5.3634710281748291355541151243e62); break;
        case 84:  return static_cast<ValueType>(6.1792822542628292342360018318e63); break;
        case 85:  return static_cast<ValueType>(4.5589503739486047652209978556e64); break;
        case 86:  return static_cast<ValueType>(5.3141827386660331414429615754e65); break;
        case 87:  return static_cast<ValueType>(3.9662868253352861457422681344e66); break;
        case 88:  return static_cast<ValueType>(4.6764808100261091644698061863e67); break;
        case 89:  return static_cast<ValueType>(3.5299952745484046697106186396e68); break;
        case 90:  return static_cast<ValueType>(4.2088327290234982480228255677e69); break;
        case 91:  return static_cast<ValueType>(3.2122956998390482494366629620e70); break;
        case 92:  return static_cast<ValueType>(3.8721261107016183881809995223e71); break;
        case 93:  return static_cast<ValueType>(2.98743500085031487197609655470e72); break;
        case 94:  return static_cast<ValueType>(3.6397985440595212848901395509e73); break;
        case 95:  return static_cast<ValueType>(2.83806325080779912837729172696e74); break;
        case 96:  return static_cast<ValueType>(3.4942066022971404334945339689e75); break;
        case 97:  return static_cast<ValueType>(2.75292135328356515452597297515e76); break;
        case 98:  return static_cast<ValueType>(3.4243224702511976248246432895e77); break;
        case 99:  return static_cast<ValueType>(2.72539213975072950298071324540e78); break;
        case 100: return static_cast<ValueType>(3.4243224702511976248246432895e79); break;
        case 101: return static_cast<ValueType>(2.75264606114823679801052037785e80); break;
        case 102: return static_cast<ValueType>(3.4928089196562215773211361553e81); break;
        case 103: return static_cast<ValueType>(2.83522544298268390195083598919e82); break;
        case 104: return static_cast<ValueType>(3.6325212764424704404139816015e83); break;
        case 105: return static_cast<ValueType>(2.97698671513181809704837778865e84); break;
        case 106: return static_cast<ValueType>(3.8504725530290186668388204976e85); break;
        case 107: return static_cast<ValueType>(3.1853757851910453638417642339e86); break;
        case 108: return static_cast<ValueType>(4.1585103572713401601859261374e87); break;
        case 109: return static_cast<ValueType>(3.4720596058582394465875230149e88); break;
        case 110: return static_cast<ValueType>(4.5743613929984741762045187512e89); break;
        case 111: return static_cast<ValueType>(3.8539861625026457857121505465e90); break;
        case 112: return static_cast<ValueType>(5.1232847601582910773490610013e91); break;
        case 113: return static_cast<ValueType>(4.3550043636279897378547301176e92); break;
        case 114: return static_cast<ValueType>(5.8405446265804518281779295415e93); break;
        case 115: return static_cast<ValueType>(5.0082550181721881985329396352e94); break;
        case 116: return static_cast<ValueType>(6.7750317668333241206863982681e95); break;
        case 117: return static_cast<ValueType>(5.8596583712614601922835393732e96); break;
        case 118: return static_cast<ValueType>(7.9945374848633224624099499564e97); break;
        case 119: return static_cast<ValueType>(6.9729934618011376288174118541e98); break;
        case 120: return static_cast<ValueType>(9.5934449818359869548919399477e99); break;
        case 121: return static_cast<ValueType>(8.4373220887793765308690683435e100); break;
        case 122: return static_cast<ValueType>(1.17040028778399040849681667362e102); break;
        case 123: return static_cast<ValueType>(1.03779061691986331329689540625e103); break;
        case 124: return static_cast<ValueType>(1.45129635685214810653605267528e104); break;
        case 125: return static_cast<ValueType>(1.29723827114982914162111925781e105); break;
        case 126: return static_cast<ValueType>(1.82863340963370661423542637086e106); break;
        case 127: return static_cast<ValueType>(1.64749260436028300985882145742e107); break;
        case 128: return static_cast<ValueType>(2.34065076433114446622134575470e108); break;
        case 129: return static_cast<ValueType>(2.12526545962476508271787968008e109); break;
        case 130: return static_cast<ValueType>(3.04284599363048780608774948111e110); break;
        case 131: return static_cast<ValueType>(2.78409775210844225836042238090e111); break;
        case 132: return static_cast<ValueType>(4.0165567115922439040358293151e112); break;
        case 133: return static_cast<ValueType>(3.7028500103042282036193617666e113); break;
        case 134: return static_cast<ValueType>(5.3821859935336068314080112822e114); break;
        case 135: return static_cast<ValueType>(4.9988475139107080748861383849e115); break;
        case 136: return static_cast<ValueType>(7.3197729512057052907148953438e116); break;
        case 137: return static_cast<ValueType>(6.8484210940576700625940095873e117); break;
        case 138: return static_cast<ValueType>(1.01012866726638733011865555744e119); break;
        case 139: return static_cast<ValueType>(9.5193053207401613870056733264e119); break;
        case 140: return static_cast<ValueType>(1.41418013417294226216611778042e121); break;
        case 141: return static_cast<ValueType>(1.34222205022436275556779993902e122); break;
        case 142: return static_cast<ValueType>(2.00813579052557801227588724819e123); break;
        case 143: return static_cast<ValueType>(1.91937753182083874046195391280e124); break;
        case 144: return static_cast<ValueType>(2.89171553835683233767727763739e125); break;
        case 145: return static_cast<ValueType>(2.78309742114021617366983317355e126); break;
        case 146: return static_cast<ValueType>(4.2219046860009752130088253506e127); break;
        case 147: return static_cast<ValueType>(4.0911532090761177752946547651e128); break;
        case 148: return static_cast<ValueType>(6.2484189352814433152530615189e129); break;
        case 149: return static_cast<ValueType>(6.0958182815234154851890356000e130); break;
        case 150: return static_cast<ValueType>(9.3726284029221649728795922783e131); break;
        case 151: return static_cast<ValueType>(9.2046856051003573826354437561e132); break;
        case 152: return static_cast<ValueType>(1.42463951724416907587769802630e134); break;
        case 153: return static_cast<ValueType>(1.40831689758035467954322289468e135); break;
        case 154: return static_cast<ValueType>(2.19394485655602037685165496051e136); break;
        case 155: return static_cast<ValueType>(2.18289119124954975329199548675e137); break;
        case 156: return static_cast<ValueType>(3.4225539762273917878885817384e138); break;
        case 157: return static_cast<ValueType>(3.4271391702617931126684329142e139); break;
        case 158: return static_cast<ValueType>(5.4076352824392790248639591467e140); break;
        case 159: return static_cast<ValueType>(5.4491512807162510491428083336e141); break;
        case 160: return static_cast<ValueType>(8.6522164519028464397823346347e142); break;
        case 161: return static_cast<ValueType>(8.7731335619531641891199214170e143); break;
        case 162: return static_cast<ValueType>(1.40165906520826112324473821082e145); break;
        case 163: return static_cast<ValueType>(1.43002077059836576282654719098e146); break;
        case 164: return static_cast<ValueType>(2.29872086694154824212137066574e147); break;
        case 165: return static_cast<ValueType>(2.35953427148730350866380286512e148); break;
        case 166: return static_cast<ValueType>(3.8158766391229700819214753051e149); break;
        case 167: return static_cast<ValueType>(3.9404222333837968594685507847e150); break;
        case 168: return static_cast<ValueType>(6.4106727537265897376280785126e151); break;
        case 169: return static_cast<ValueType>(6.6593135744186166925018508262e152); break;
        case 170: return static_cast<ValueType>(1.08981436813352025539677334714e154); break;
        case 171: return static_cast<ValueType>(1.13874262122558345441781649128e155); break;
        case 172: return static_cast<ValueType>(1.87448071318965483928245015709e156); break;
        case 173: return static_cast<ValueType>(1.97002473472025937614282252992e157); break;
        case 174: return static_cast<ValueType>(3.2615964409499994203514632733e158); break;
        case 175: return static_cast<ValueType>(3.4475432857604539082499394274e159); break;
        case 176: return static_cast<ValueType>(5.7404097360719989798185753611e160); break;
        case 177: return static_cast<ValueType>(6.1021516157960034176023927864e161); break;
        case 178: return static_cast<ValueType>(1.02179293302081581840770641427e163); break;
        case 179: return static_cast<ValueType>(1.09228513922748461175082830877e164); break;
        case 180: return static_cast<ValueType>(1.83922727943746847313387154568e165); break;
        case 181: return static_cast<ValueType>(1.97703610200174714726899923887e166); break;
        case 182: return static_cast<ValueType>(3.3473936485761926211036462131e167); break;
        case 183: return static_cast<ValueType>(3.6179760666631972795022686071e168); break;
        case 184: return static_cast<ValueType>(6.1592043133801944228307090322e169); break;
        case 185: return static_cast<ValueType>(6.6932557233269149670791969232e170); break;
        case 186: return static_cast<ValueType>(1.14561200228871616264651187999e172); break;
        case 187: return static_cast<ValueType>(1.25163882026213309884380982464e173); break;
        case 188: return static_cast<ValueType>(2.15375056430278638577544233437e174); break;
        case 189: return static_cast<ValueType>(2.36559737029543155681480056857e175); break;
        case 190: return static_cast<ValueType>(4.0921260721752941329733404353e176); break;
        case 191: return static_cast<ValueType>(4.5182909772642742735162690860e177); break;
        case 192: return static_cast<ValueType>(7.8568820585765647353088136358e178); break;
        case 193: return static_cast<ValueType>(8.7203015861200493478863993359e179); break;
        case 194: return static_cast<ValueType>(1.52423511936385355864990984535e181); break;
        case 195: return static_cast<ValueType>(1.70045880929340962283784787050e182); break;
        case 196: return static_cast<ValueType>(2.98750083395315297495382329688e183); break;
        case 197: return static_cast<ValueType>(3.3499038543080169569905603049e184); break;
        case 198: return static_cast<ValueType>(5.9152516512272428904085701278e185); break;
        case 199: return static_cast<ValueType>(6.6663086700729537444112150067e186); break;
        case 200: return static_cast<ValueType>(1.18305033024544857808171402556e188); break;
        case 201: return static_cast<ValueType>(1.33992804268466370262665421635e189); break;
        case 202: return static_cast<ValueType>(2.38976166709580612772506233164e190); break;
        case 203: return static_cast<ValueType>(2.72005392664986731633210805920e191); break;
        case 204: return static_cast<ValueType>(4.8751138008754445005591271565e192); break;
        case 205: return static_cast<ValueType>(5.5761105496322279984808215214e193); break;
        case 206: return static_cast<ValueType>(1.00427344298034156711518019425e195); break;
        case 207: return static_cast<ValueType>(1.15425488377387119568553005492e196); break;
        case 208: return static_cast<ValueType>(2.08888876139911045959957480403e197); break;
        case 209: return static_cast<ValueType>(2.41239270708739079898275781478e198); break;
        case 210: return static_cast<ValueType>(4.3866663989381319651591070885e199); break;
        case 211: return static_cast<ValueType>(5.0901486119543945858536189892e200); break;
        case 212: return static_cast<ValueType>(9.2997327657488397661373070276e201); break;
        case 213: return static_cast<ValueType>(1.08420165434628604678682084470e203); break;
        case 214: return static_cast<ValueType>(1.99014281187025170995338370390e204); break;
        case 215: return static_cast<ValueType>(2.33103355684451500059166481610e205); break;
        case 216: return static_cast<ValueType>(4.2987084736397436934993088004e206); break;
        case 217: return static_cast<ValueType>(5.0583428183525975512839126509e207); break;
        case 218: return static_cast<ValueType>(9.3711844725346412518284931849e208); break;
        case 219: return static_cast<ValueType>(1.10777707721921886373117687056e210); break;
        case 220: return static_cast<ValueType>(2.06166058395762107540226850068e211); break;
        case 221: return static_cast<ValueType>(2.44818734065447368884590088393e212); break;
        case 222: return static_cast<ValueType>(4.5768864963859187873930360715e213); break;
        case 223: return static_cast<ValueType>(5.4594577696594763261263589712e214); break;
        case 224: return static_cast<ValueType>(1.02522257519044580837604008002e216); break;
        case 225: return static_cast<ValueType>(1.22837799817338217337843076851e217); break;
        case 226: return static_cast<ValueType>(2.31700301993040752692985058084e218); break;
        case 227: return static_cast<ValueType>(2.78841805585357753356903784452e219); break;
        case 228: return static_cast<ValueType>(5.2827668854413291614000593243e220); break;
        case 229: return static_cast<ValueType>(6.3854773479046925518730966640e221); break;
        case 230: return static_cast<ValueType>(1.21503638365150570712201364459e223); break;
        case 231: return static_cast<ValueType>(1.47504526736598397948268532937e224); break;
        case 232: return static_cast<ValueType>(2.81888441007149324052307165546e225); break;
        case 233: return static_cast<ValueType>(3.4368554729627426721946568174e226); break;
        case 234: return static_cast<ValueType>(6.5961895195672941828239876738e227); break;
        case 235: return static_cast<ValueType>(8.0766103614624452796574435210e228); break;
        case 236: return static_cast<ValueType>(1.55670072661788142714646109101e230); break;
        case 237: return static_cast<ValueType>(1.91415665566659953127881411447e231); break;
        case 238: return static_cast<ValueType>(3.7049477293505577966085773966e232); break;
        case 239: return static_cast<ValueType>(4.5748344070431728797563657336e233); break;
        case 240: return static_cast<ValueType>(8.8918745504413387118605857518e234); break;
        case 241: return static_cast<ValueType>(1.10253509209740466402128414180e236); break;
        case 242: return static_cast<ValueType>(2.15183364120680396827026175195e237); break;
        case 243: return static_cast<ValueType>(2.67916027379669333357172046456e238); break;
        case 244: return static_cast<ValueType>(5.2504740845446016825794386748e239); break;
        case 245: return static_cast<ValueType>(6.5639426708018986672507151382e240); break;
        case 246: return static_cast<ValueType>(1.29161662479797201391454191399e242); break;
        case 247: return static_cast<ValueType>(1.62129383968806897081092663913e243); break;
        case 248: return static_cast<ValueType>(3.2032092294989705945080639467e244); break;
        case 249: return static_cast<ValueType>(4.0370216608232917373192073314e245); break;
        case 250: return static_cast<ValueType>(8.0080230737474264862701598667e246); break;
        case 251: return static_cast<ValueType>(1.01329243686664622606712104019e248); break;
        case 252: return static_cast<ValueType>(2.01802181458435147454008028642e249); break;
        case 253: return static_cast<ValueType>(2.56362986527261495194981623168e250); break;
        case 254: return static_cast<ValueType>(5.1257754090442527453318039275e251); break;
        case 255: return static_cast<ValueType>(6.5372561564451681274720313908e252); break;
        case 256: return static_cast<ValueType>(1.31219850471532870280494180544e254); break;
        case 257: return static_cast<ValueType>(1.68007483220640820876031206743e255); break;
        case 258: return static_cast<ValueType>(3.3854721421655480532367498580e256); break;
        case 259: return static_cast<ValueType>(4.3513938154145972606892082546e257); break;
        case 260: return static_cast<ValueType>(8.8022275696304249384155496309e258); break;
        case 261: return static_cast<ValueType>(1.13571378582320988503988335446e260); break;
        case 262: return static_cast<ValueType>(2.30618362324317133386487400329e261); break;
        case 263: return static_cast<ValueType>(2.98692725671504199765489322224e262); break;
        case 264: return static_cast<ValueType>(6.0883247653619723214032673687e263); break;
        case 265: return static_cast<ValueType>(7.9153572302948612937854670389e264); break;
        case 266: return static_cast<ValueType>(1.61949438758628463749326912007e266); break;
        case 267: return static_cast<ValueType>(2.11340038048872796544071969939e267); break;
        case 268: return static_cast<ValueType>(4.3402449587312428284819612418e268); break;
        case 269: return static_cast<ValueType>(5.6850470235146782270355359914e269); break;
        case 270: return static_cast<ValueType>(1.17186613885743556369012953528e271); break;
        case 271: return static_cast<ValueType>(1.54064774337247779952663025366e272); break;
        case 272: return static_cast<ValueType>(3.1874758976922247332371523360e273); break;
        case 273: return static_cast<ValueType>(4.2059683394068643927077005925e274); break;
        case 274: return static_cast<ValueType>(8.7336839596766957690697974006e275); break;
        case 275: return static_cast<ValueType>(1.15664129333688770799461766294e277); break;
        case 276: return static_cast<ValueType>(2.41049677287076803226326408256e278); break;
        case 277: return static_cast<ValueType>(3.2038963825431789511450909263e279); break;
        case 278: return static_cast<ValueType>(6.7011810285807351296918741495e280); break;
        case 279: return static_cast<ValueType>(8.9388709072954692736948036845e281); break;
        case 280: return static_cast<ValueType>(1.87633068800260583631372476186e283); break;
        case 281: return static_cast<ValueType>(2.51182272495002686590823983534e284); break;
        case 282: return static_cast<ValueType>(5.2912525401673484584047038284e285); break;
        case 283: return static_cast<ValueType>(7.1084583116085760305203187340e286); break;
        case 284: return static_cast<ValueType>(1.50271572140752696218693588728e288); break;
        case 285: return static_cast<ValueType>(2.02591061880844416869829083919e289); break;
        case 286: return static_cast<ValueType>(4.2977669632255271118546366376e290); break;
        case 287: return static_cast<ValueType>(5.8143634759802347641640947085e291); break;
        case 288: return static_cast<ValueType>(1.23775688540895180821413535163e293); break;
        case 289: return static_cast<ValueType>(1.68035104455828784684342337075e294); break;
        case 290: return static_cast<ValueType>(3.5894949676859602438209925197e295); break;
        case 291: return static_cast<ValueType>(4.8898215396646176343143620089e296); break;
        case 292: return static_cast<ValueType>(1.04813253056430039119572981576e298); break;
        case 293: return static_cast<ValueType>(1.43271771112173296685410806860e299); break;
        case 294: return static_cast<ValueType>(3.08150963985904315011544565835e300); break;
        case 295: return static_cast<ValueType>(4.2265172478091122522196188024e301); break;
        case 296: return static_cast<ValueType>(9.1212685339827677243417191487e302); break;
        case 297: return static_cast<ValueType>(1.25527562259930633890922678431e304); break;
        default: return static_cast<ValueType>(0.0); break;
    }
}
*/
/*
std::ostream& Multipole::Gaunt::operator<<(std::ostream& os, const Matrix& mat)
{
    using namespace Gripper;

    ENTERING

    clog << LogLevel::INFO << "Formatting Multipole::Gaunt::Matrix contents";

    std::ostringstream temp;

    for (auto& element : mat.m_data) temp << "{" << element.l1m1.l.i << "," << element.l1m1.m.i << "}\t{" << element.l2m2.l.i << "," << element.l2m2.m.i << "}\t{" << element.l3m3.l.i << "," << element.l3m3.m.i << "}\t" << element.value << "\n";

    LEAVING

    return os << temp.str();
}


namespace Gripper
{
    EXPORT Multipole::Gaunt::Matrix gaunt;
} // namespace Gripper

*/