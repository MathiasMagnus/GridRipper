#include <Gripper/cl/clRuntime.hpp>

namespace Gripper
{
    namespace cl
    {
        namespace Private
        {
            // Global runtime instance

            EXPORT Runtime g_runtime;

        } // namespace Private

    } // namespace cl

} // namespace Gripper


Gripper::cl::Runtime::Runtime() :
    Gripper::Runtime(),
    m_cl_parser(new AnyOption())
{
    m_cl_parser->addUsage("");
    m_cl_parser->addUsage(" OpenCL compute back-end related switches");
    m_cl_parser->addUsage("");
    m_cl_parser->addUsage(" OpenCL options: <token> <value>");
    m_cl_parser->addUsage("");
    m_cl_parser->addUsage(" -d  --device      Possible values: cpu, gpu, all [default: gpu]");
    m_cl_parser->addUsage("");
}


Gripper::cl::Runtime::~Runtime()
{
}


void Gripper::cl::Runtime::initialize(int argc, char** argv)
{
    if (argc != 0)
    {
        m_argc = argc;
        m_argv = argv;
    }

    Gripper::Runtime::processCommandLineArgs();

    Gripper::Runtime::processConfigFile();

    processCommandLineArgs();

    processConfigFile();

    if (m_base_parser->getFlag('h'))
    {
        m_base_parser->printUsage();
        m_cl_parser->printUsage();

        exit(EXIT_SUCCESS);
    }

    Gripper::clog.setLogLevel(log_level);

    Gripper::clog << LogLevel::DEBUG << "Basic Logger instantiated";

    clog.setLogLevel(log_level);

    clog << LogLevel::DEBUG << "OpenCL Logger instantiated";

    Gripper::Runtime::logOptions();

    logOptions();

    clog << LogLevel::INFO << "Initializing Gripper base runtime";

    Gripper::Runtime::init();

    clog << LogLevel::INFO << "Initializing Gripper OpenCL runtime";

    init();
}


void Gripper::cl::Runtime::processCommandLineArgs()
{
    m_cl_parser->processCommandArgs(m_argc, m_argv);
}


void Gripper::cl::Runtime::processConfigFile()
{
    m_cl_parser->processFile(m_config_file.c_str());
}


void Gripper::cl::Runtime::logOptions()
{
    clog << LogLevel::DEBUG << "Entering " << __FUNCTION__;

    if (m_cl_parser->getFlag('t'))
        clog << LogLevel::INFO << "Using guided parallelization";
    else
        clog << LogLevel::INFO << "Using only compiler optimizations";

    clog << LogLevel::DEBUG << "Leaving " << __FUNCTION__;
}


void Gripper::stl::Runtime::init()
{
    clog << LogLevel::DEBUG << "Entering " << __FUNCTION__;

    clog << LogLevel::INFO << "Optimizing Gaunt coefficient matrix for STL back-end";

    gaunt = Multipole::Gaunt::Matrix(Gripper::gaunt);

    clog << LogLevel::DEBUG << "Leaving " << __FUNCTION__;
}


void Gripper::stl::initialize(int argc, char** argv)
{
    Private::g_runtime.initialize(argc, argv);
}