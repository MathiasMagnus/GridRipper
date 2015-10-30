#include <Gripper/Logger.hpp>


Gripper::Logger::Logger(LogLevel log_lev_in)
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


Gripper::Logger::~Logger() {}


void Gripper::Logger::setLogLevel(LogLevel target) { m_log_level = target; }


void Gripper::Logger::operator<<(LogBuffer& msg)
{
    if (!msg.m_data.str().empty()) // Needed for the NRVO did not kick in at Gripper::operator<< (Logger& lhs, LogLevel rhs) and the moved local gets destroyed with empty brains (MSVC Debug mode only)
    {
        if (msg.m_log_level <= m_log_level)
        {
            std::lock_guard<std::mutex> lock(m_access);

            if (msg.m_log_level >= LogLevel::INFO)
            {
                if (msg.m_log_level != LogLevel::NONE) std::cout << msg.m_data.str() << std::endl;
            }
            else { std::cerr << msg.m_data.str() << std::endl; }
        }
        if (msg.m_log_level == LogLevel::CRIT)
            exit(EXIT_FAILURE);
    }
}


Gripper::LogBuffer::LogBuffer(Logger& logger_in, LogLevel log_lev_in, std::string&& msg_in) : m_log_level(log_lev_in), m_data(msg_in), m_output(logger_in) {}


Gripper::LogBuffer::LogBuffer(const LogBuffer& in) : m_log_level(in.m_log_level), m_data(), m_output(in.m_output) {}


Gripper::LogBuffer::LogBuffer(LogBuffer&& in) : m_log_level(in.m_log_level), m_data(std::move(in.m_data)), m_output(in.m_output) {}


Gripper::LogBuffer::LogBuffer(Logger& logger_in, LogLevel log_lev_in) : m_log_level(log_lev_in), m_data(), m_output(logger_in) {}


Gripper::LogBuffer::~LogBuffer() { m_output << *this; }


Gripper::LogBuffer Gripper::operator<< (Logger& lhs, LogLevel rhs)
{
    LogBuffer result(lhs, rhs); result << lhs.m_log_strings.at(rhs); return result;
}

namespace Gripper
{
    EXPORT LogLevel log_level = LogLevel::WARN;
    EXPORT Logger clog(log_level);

} // namespace Gripper