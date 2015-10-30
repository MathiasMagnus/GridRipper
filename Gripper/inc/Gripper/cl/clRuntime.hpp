#ifndef CLRUNTIME_HPP
#define CLRUNTIME_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/Runtime.hpp>
#include <Gripper/cl/clLogger.hpp>
#include <Gripper/cl/clMultipoleTypes.hpp>


namespace Gripper
{
    namespace cl
    {
        class EXPORT Runtime : Gripper::Runtime
        {
        public:

            Runtime();
            ~Runtime();

            virtual void initialize(int argc, char** argv) final;

        protected:

            std::unique_ptr<AnyOption> m_cl_parser;

            virtual void processCommandLineArgs() final;

            virtual void processConfigFile() final;

            virtual void logOptions() final;

            virtual void init() final;
        };

        namespace Private
        {
            // Global runtime instance

            EXPORT extern Runtime g_runtime;

        } // namespace Private

        EXPORT void initialize(int argc = 0, char** argv = nullptr);

    } // namespace cl

} // namespace Gripper

#endif // CLRUNTIME_HPP