#ifndef STLRUNTIME_HPP
#define STLRUNTIME_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/Runtime.hpp>
#include <Gripper/stl/stlLogger.hpp>
#include <Gripper/stl/stlGauntMatrix.hpp>


namespace Gripper
{
    namespace stl
    {
        class EXPORT Runtime : public Gripper::Runtime
        {
        public:

            Runtime();
            ~Runtime();

            virtual void initialize(int argc, char** argv) final;

            bool getCompilerParallelization() const;

        protected:

            std::unique_ptr<AnyOption> m_stl_parser;

            virtual void processCommandLineArgs() final;

            virtual void processConfigFile() final;

            virtual void logOptions() final;

            virtual void init() final;

        private:

            void applySTLOptions();

            bool m_compiler_parallelization;
        };

        namespace Private
        {
            // Global runtime instance

            EXPORT extern Runtime g_runtime;

        } // namespace Private

        EXPORT void initialize(int argc = 0, char** argv = nullptr);

        EXPORT Multipole::stl::Spherical::Index getL();

        EXPORT ArithmeticPrecision getPrecision();

        EXPORT bool getCompilerParallelization();

    } // namespace stl

} // namespace Gripper

#endif // STLRUNTIME_HPP