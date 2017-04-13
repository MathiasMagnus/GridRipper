#pragma once

// ISO C++ includes
#include <cstdint>      // std::int16_t

namespace math
{
    namespace sws
    {
        enum parity : std::int16_t
        {
            odd = -1,
            even = 1
        };
    }
}
