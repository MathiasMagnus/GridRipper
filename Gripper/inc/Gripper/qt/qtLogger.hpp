///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Implementation of LoggerBase virtual class logging to std::cout,std::cerr //
//                                                                           //
// Author: Nagy-Egri Máté Ferenc                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOGGER_HPP
#define LOGGER_HPP

// Gripper includes
#include <Gripper/LoggerBase.hpp>

// Standard C++ includes
#include <string>
#include <iostream>

namespace Gripper
{
    namespace STL
    {
        class Logger : public LoggerBase
        {
        public:

            Logger();
            Logger(LogLevel loglevel);
            ~Logger();

            void log(const LogLevel logLev,		// the loglevel of the message
                const std::string trace,	// the location of the message
                const std::string msg,		// the message itself
                ...							// variable length of arguments, just as for printf
                );
        };

    } // namespace STL

} // namespace Gripper

#endif // LOGGER_HPP