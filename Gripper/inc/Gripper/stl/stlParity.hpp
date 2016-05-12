#pragma once

// ISO C++ includes
#include <cstdint>

namespace Multipole
{
    namespace stl
    {
        enum Parity : std::int16_t
        {
            Odd = -1,
            Even = 1
        };
    }
}
