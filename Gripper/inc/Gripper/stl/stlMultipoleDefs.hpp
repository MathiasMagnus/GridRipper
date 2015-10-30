#ifndef STLMULTIPOLEDEFS_HPP
#define STLMULTIPOLEDEFS_HPP

// Gripper includes
#include <Gripper/MultipoleDefs.hpp> // Forward declaration of Multipole types

// Standard C++ includes
#include <cstdint>
#include <utility>


namespace Multipole
{
    namespace stl
    {
        namespace Gaunt
        {
            class Index;

            class Coefficient;
            class Matrix;

        } // namespace Gaunt

        namespace Radial
        {
            class Index;

            class Extent;

            template <typename VT> class Traits;
            template <typename ET, typename VT> class Expression;
            template <typename VT> class Vector;

        } // namespace Radial

        namespace Spherical
        {
            class Index;
            class IndexPair;

            class Extent;

            template <typename VT> class Traits;
            template <typename ET, typename VT> class Expression;
            template <typename VT> class Vector;

        } // namespace Spherical

        namespace Expansion
        {
            class Index;

            class Extent;

            namespace Constant
            {
                template <typename VT> class Traits;
                template <typename ET, typename VT> class Expression;
                template <typename VT> class Factor;

            } // namespace Constant

            template <typename VT> class Traits;
            template <typename ET, typename VT> class Expression;
            template <typename VT> class Field;
            /*
            enum DiffMethod
            {
                Exact,
                O421
            };

            enum DiffDirection
            {
                R,
                Theta,
                Phi
            };
            */
        } // namespace Expansion

    } // namespace stl

} // namespace Multipole

#endif // STLMULTIPOLEDEFS_HPP