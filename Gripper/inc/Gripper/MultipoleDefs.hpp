#ifndef MULTIPOLEDEFS_HPP
#define MULTIPOLEDEFS_HPP

// Standard C++ includes
#include <cstdint>


namespace Multipole
{
    namespace Radial
    {
        typedef int32_t IndexType;
        class Index;

    } // namespace Radial

    namespace Spherical
    {
        typedef int32_t IndexType;
        class Index;

        class IndexPair;

    } // namespace Spherical

	namespace Spin
	{
		typedef int32_t IndexType;
		class Index;

	} // namespace Spherical

    namespace Gaunt
    {
        typedef uint32_t IndexType;
        class Index;

        typedef double ValueType;

        class Coefficient;
        class Matrix;

    } // namespace Gaunt

} // namespace Multipole

#endif // MULTIPOLEDEFS_HPP