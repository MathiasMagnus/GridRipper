#pragma once
/*
// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/Logger.hpp>
#include <Gripper/MultipoleTypes.hpp>

// Standard C++ includes
#include <memory>
#include <sstream>
#include <cstdlib>

// Command-line option parser include
#include <Gripper/AnyOption.hpp>


namespace Gripper
{
    enum FiniteDifferenceScheme
    {
        RK4
    };

    enum ArithmeticPrecision
    {
        SinglePrecision = 0,
        DoublePrecision = 1,
    };

    class EXPORT Runtime
    {
    public:

        Runtime();
        ~Runtime();

        virtual void initialize(int argc, char** argv);

        Multipole::Spherical::Index getL() const;
        ArithmeticPrecision getPrecision() const;

    protected:

        int m_argc;
        char** m_argv;

        std::unique_ptr<AnyOption> m_base_parser;

        std::string m_config_file;

        virtual void processCommandLineArgs();

        virtual void processConfigFile();

        virtual void logOptions();

        virtual void init();

    private:

        void applyOptions();

        Multipole::Spherical::Index m_L_max;
        ArithmeticPrecision m_precision;
    };

    namespace Private
    {
        // Global runtime instance
        EXPORT extern Runtime g_runtime;

    } // namespace Private


    // Global query functions

    EXPORT void initialize(int argc = 0, char** argv = nullptr);

    EXPORT Multipole::Spherical::Index getL();

    EXPORT ArithmeticPrecision getPrecision();

} // namespace Gripper
*/