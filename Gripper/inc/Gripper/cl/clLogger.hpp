///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Implementation of LoggerBase virtual class logging to std::cout,std::cerr //
//                                                                           //
// Author: Nagy-Egri Máté Ferenc                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef CLLOGGER_HPP
#define CLLOGGER_HPP


// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/Logger.hpp>

// Standard C++ includes
#include <array>
#include <string>
#include <sstream>

namespace Gripper
{
    namespace cl
    {
        class Logger;
        class LogBuffer;

        class EXPORT LogBuffer
        {
            friend EXPORT LogBuffer operator<< (Logger& lhs, LogLevel rhs);
            template <typename T> friend LogBuffer& operator<<(LogBuffer& lhs, const T& rhs);

            friend EXPORT Logger;

        public:

            ~LogBuffer();

        private:

            LogBuffer(const LogBuffer& in);
            LogBuffer(LogBuffer&& in);
            LogBuffer(Logger& logger_in, LogLevel log_lev_in);

            LogLevel m_log_level;       // Message level of verbosity

            std::stringstream m_data;   // Message contents
            Logger& m_output;           // Logger instance
        };

        class EXPORT Logger
        {
            friend EXPORT LogBuffer operator<<(Logger& lhs, LogLevel rhs);
            template <typename T> friend LogBuffer& operator<<(LogBuffer& lhs, const T& rhs);

            friend EXPORT LogBuffer;

        public:

            Logger(LogLevel log_lev_in, Gripper::Logger& dst);
            ~Logger();

            void setLogLevel(LogLevel target);

            void operator<<(LogBuffer& msg);

        private:

            LogLevel m_log_level;                       // Logger level of verbosity

            std::array<std::string, 9> m_log_strings;   // Loglevel strings
            Gripper::Logger& m_output;                  // Instance of base logger
        };

        EXPORT LogBuffer operator<< (Logger& lhs, LogLevel rhs);

        template <typename T> LogBuffer& operator<< (LogBuffer& lhs, const T& rhs) { if (lhs.m_log_level <= lhs.m_output.m_log_level) lhs.m_data << rhs; return lhs; }

        EXPORT extern Logger clog;

    } // namespace cl

} // namespace Gripper

#endif // CLLOGGER_HPP