#include <Gripper/stl/stlLogger.hpp>


Gripper::stl::Logger::Logger(LogLevel log_lev_in, Gripper::Logger& dst) : m_output(dst)
{
    m_log_level = log_lev_in;
    m_log_strings.at(0) = std::string("");
    m_log_strings.at(1) = std::string("CRIT  : ");
    m_log_strings.at(2) = std::string("ERROR : ");
    m_log_strings.at(3) = std::string("WARN  : ");
    m_log_strings.at(4) = std::string("INFO  : ");
    m_log_strings.at(5) = std::string("PROF  : ");
    m_log_strings.at(6) = std::string("DEBUG : ");
    m_log_strings.at(7) = std::string("ALL   : ");
    m_log_strings.at(8) = std::string("INSANE: ");
}


Gripper::stl::Logger::~Logger() {}


void Gripper::stl::Logger::setLogLevel(LogLevel target) { m_log_level = target; }


void Gripper::stl::Logger::operator<<(LogBuffer& msg) { Gripper::LogBuffer(m_output, msg.m_log_level, std::move(msg.m_data.str())); }


Gripper::stl::LogBuffer::~LogBuffer() { m_output << *this; }


Gripper::stl::LogBuffer::LogBuffer(const LogBuffer& in) : m_log_level(in.m_log_level), m_data(), m_output(in.m_output) {}


Gripper::stl::LogBuffer::LogBuffer(LogBuffer&& in) : m_log_level(in.m_log_level), m_data(std::move(in.m_data)), m_output(in.m_output) {}


Gripper::stl::LogBuffer::LogBuffer(Logger& logger_in, LogLevel log_lev_in) : m_log_level(log_lev_in), m_data(), m_output(logger_in) {}


Gripper::stl::LogBuffer Gripper::stl::operator<<(Logger& lhs, LogLevel rhs) { LogBuffer result(lhs, rhs); result << lhs.m_log_strings.at(rhs); return result; }


namespace Gripper
{
    namespace stl
    {
        EXPORT Logger clog(Gripper::log_level, Gripper::clog);

    } // namespace stl

} // namespace Gripper