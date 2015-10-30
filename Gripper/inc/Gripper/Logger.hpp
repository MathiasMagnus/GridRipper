//////////////////////////////////////////////////////////////////////
//                                                                  //
// Threadsafe logger class logging to std::cout,std::clog,std::cerr //
//                                                                  //
// Author: Nagy-Egri Máté Ferenc                                    //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#ifndef LOGGER_HPP
#define LOGGER_HPP


// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Standard C++ includes
#include <mutex>
#include <array>
#include <string>
#include <sstream>
#include <iostream>

namespace Gripper
{
    enum LogLevel
    {
        INSANE = 8, // everything that you don't shame
        ALL = 7, 	// the same as debug level but can be used for log inside loops (huge output)
        DEBUG = 6, 	// debug level for finding bugs, developing
        TIMING = 5,	// profiling level
        INFO = 4, 	// normal level, informing the user on what is happening
        WARN = 3, 	// logs only the warnings
        ERR = 2, 	// logs only error message
        CRIT = 1, 	// logs only critical error messages
        NONE = 0  	// no log output, silent run
    };

    class Logger;
    class LogBuffer;

    class EXPORT LogBuffer
    {
        friend EXPORT LogBuffer operator<< (Logger& lhs, LogLevel rhs);
        template <typename T> friend LogBuffer& operator<<(LogBuffer& lhs, const T& rhs);

        friend EXPORT Logger;

    public:

        LogBuffer(Logger& logger_in, LogLevel log_lev_in, std::string&& msg_in);    // Compute back-end logging interface
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

        Logger(LogLevel log_lev_in);
        ~Logger();

        void setLogLevel(LogLevel target);

        void operator<<(LogBuffer& msg);

    private:

        LogLevel m_log_level;                       // Logger level of verbosity

        std::mutex m_access;                        // Ownership mutex
        std::array<std::string, 9> m_log_strings;   // Loglevel strings
    };

    EXPORT LogBuffer operator<< (Logger& lhs, LogLevel rhs);

    template <typename T> LogBuffer& operator<< (LogBuffer& lhs, const T& rhs) { if (lhs.m_log_level <= lhs.m_output.m_log_level) lhs.m_data << rhs; return lhs; }

    EXPORT extern LogLevel log_level;
    EXPORT extern Logger clog;

#define ENTERING Gripper::clog << Gripper::ALL << "Entering " << __FUNCTION__;
#define LEAVING  Gripper::clog << Gripper::ALL << "Leaving "  << __FUNCTION__;
#define CALLING  Gripper::clog << Gripper::ALL << "Calling "  << __FUNCTION__;

} // namespace Gripper

#endif // LOGGER_HPP