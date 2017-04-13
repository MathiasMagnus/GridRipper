#pragma once

// Structured exception handling
#include <Windows.h>        // EXCEPTION_POINTERS
#include <eh.h>             // _set_se_translator

// Standard C++ includes
#include <exception>        // std::exception
#include <initializer_list> // std::initializer_list
#include <numeric>          // std::accumulate
#include <functional>       // std::bit_or
#include <float.h>          // _clearfp, _controlfp_s


namespace floating
{
    enum class control : unsigned int
    {
        denormal = _MCW_DN,
        interrupt = _MCW_EM,
        infinity = _MCW_IC,
        rounding = _MCW_RC,
        precision = _MCW_PC
    };

    enum class denormal_mode : unsigned int
    {
        save = _DN_SAVE,
        flush = _DN_FLUSH
    };

    enum exception_mask : unsigned int
    {
        inexact = _EM_INEXACT,
        underflow = _EM_UNDERFLOW,
        overflow = _EM_OVERFLOW,
        zero_divide = _EM_ZERODIVIDE,
        invalid = _EM_INVALID,
        denormal = _EM_DENORMAL
    };

    enum class infinity_mode : unsigned int
    {
        affine = _IC_AFFINE,
        projective = _IC_PROJECTIVE
    };

    enum class rounding_mode : unsigned int
    {
        chop = _RC_CHOP,
        up = _RC_UP,
        down = _RC_DOWN,
        nearest = _RC_NEAR
    };

    enum class precision_mode : unsigned int
    {
        _24 = _PC_24,
        _53 = _PC_53,
        _64 = _PC_64
    };

    class exception : public std::exception {};

    class denormal_operand : public exception {};
    class divide_by_zero : public exception {};
    class inexact_result : public exception {};
    class invalid_operation : public exception {};
    class overflow : public exception {};
    class stack_check : public exception {};
    class underflow : public exception {};

    namespace detail
    {
        void se_fe_translator_func(unsigned int u, EXCEPTION_POINTERS*)
        {
            switch (u)
            {
            case STATUS_FLOAT_DENORMAL_OPERAND:   throw denormal_operand{};
            case STATUS_FLOAT_DIVIDE_BY_ZERO:     throw divide_by_zero{};
            case STATUS_FLOAT_INEXACT_RESULT:     throw denormal_operand{};
            case STATUS_FLOAT_INVALID_OPERATION:  throw divide_by_zero{};
            case STATUS_FLOAT_OVERFLOW:           throw denormal_operand{};
            case STATUS_FLOAT_STACK_CHECK:        throw divide_by_zero{};
            case STATUS_FLOAT_UNDERFLOW:          throw denormal_operand{};
            };
        }

        struct fe_translator { fe_translator() { _set_se_translator(se_fe_translator_func); } } translator;

        //exception_mask bit_or(const exception_mask& lhs, const exception_mask& rhs) { return exception_mask{ static_cast<unsigned int>(lhs) | static_cast<unsigned int>(rhs) }; }
    }

    class scoped_exception_enabler
    {
    public:

        scoped_exception_enabler(std::initializer_list<exception_mask> list = { exception_mask::overflow,
                                                                                exception_mask::zero_divide,
                                                                                exception_mask::invalid })
        {
            // Save previous exception state
            _controlfp_s(&m_old_fp_config, 0, 0);

            // Clear exception state
            _clearfp();

            // Zero out the specified bits, leaving other bits alone.
            _controlfp_s(nullptr,
                         ~std::accumulate(list.begin(),
                                          list.end(),
                                          ~m_old_fp_config,
                                          std::bit_or<unsigned int>{}),
                         static_cast<unsigned int>(control::interrupt));
        }

        scoped_exception_enabler(const scoped_exception_enabler&) = delete;
        scoped_exception_enabler& operator=(const scoped_exception_enabler&) = delete;

        ~scoped_exception_enabler()
        {
            // Clear exception state
            _clearfp();

            // Restore exception previous exception state
            _controlfp_s(nullptr, m_old_fp_config, static_cast<unsigned int>(control::interrupt));
        }

    private:

        unsigned int m_old_fp_config;
    };

    class scoped_exception_disabler
    {
    public:

        scoped_exception_disabler(std::initializer_list<exception_mask> list)
        {
            // Save previous exception state
            _controlfp_s(&m_old_fp_config, 0, 0);

            // Clear exception state
            _clearfp();

            // Zero out the specified bits, leaving other bits alone.
            _controlfp_s(nullptr,
                            std::accumulate(list.begin(),
                                            list.end(),
                                            m_old_fp_config,
                                            std::bit_or<unsigned int>{}),
                            static_cast<unsigned int>(control::interrupt));
        }

        scoped_exception_disabler(const scoped_exception_disabler&) = delete;
        scoped_exception_disabler& operator=(const scoped_exception_disabler&) = delete;

        ~scoped_exception_disabler()
        {
            // Clear exception state
            _clearfp();

            // Restore exception previous exception state
            _controlfp_s(nullptr, m_old_fp_config, static_cast<unsigned int>(control::interrupt));
        }

    private:

        unsigned int m_old_fp_config;
    };
}
