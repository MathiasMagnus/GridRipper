#include <Gripper/stl/Logger.hpp>


Gripper::STL::Logger::Logger()
{
}


Gripper::STL::Logger::~Logger()
{
}


Gripper::STL::Logger::Logger(LogLevel logLev_in)
{
	logLevel = logLev_in;
	logStrings.at(0) = std::string("");
	logStrings.at(1) = std::string("CRIT  : ");
	logStrings.at(2) = std::string("ERROR : ");
	logStrings.at(3) = std::string("WARN  : ");
	logStrings.at(4) = std::string("INFO  : ");
	logStrings.at(5) = std::string("PROF  : ");
	logStrings.at(6) = std::string("DEBUG : ");
	logStrings.at(7) = std::string("ALL   : ");
	logStrings.at(8) = std::string("INSANE: ");
}


void Gripper::STL::Logger::log(const LogLevel logLev, const std::string trace, const std::string msg, ...)
{
	std::lock_guard<std::mutex> lock(access);
	if(logLev <= logLevel)
	{
		if(logLev <= LogLevel::WARN)
		{
			if(logLev != LogLevel::NONE) std::cout << logStrings.at(logLev) << msg << " from " << trace << std::endl;
		}
        else std::cerr << logStrings.at(logLev) << msg << " from " << trace << std::endl;
	}
	if(logLev == LogLevel::CRIT) exit(EXIT_FAILURE);
}
