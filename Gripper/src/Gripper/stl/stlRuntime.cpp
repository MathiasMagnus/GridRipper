#include <Gripper/stl/stlRuntime.hpp>

namespace Gripper
{
    namespace stl
    {
        namespace Private
        {
            // Global runtime instance

            EXPORT Runtime g_runtime;

        } // namespace Private

    } // namespace stl

} // namespace Gripper


Gripper::stl::Runtime::Runtime() :
    Gripper::Runtime(),
    m_stl_parser(new AnyOption()),
    m_compiler_parallelization(false)
{
    m_stl_parser->addUsage("");
    m_stl_parser->addUsage(" STL compute back-end related switches");
    m_stl_parser->addUsage("");
    m_stl_parser->addUsage(" STL flags: <token>");
    m_stl_parser->addUsage("");
    m_stl_parser->addUsage(" -a  --auto        Auto-parallelized evaluation of arithmetic types");
    m_stl_parser->addUsage("");
}


Gripper::stl::Runtime::~Runtime()
{
}


void Gripper::stl::Runtime::initialize(int argc, char** argv)
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
        m_stl_parser->printUsage();

        exit(EXIT_SUCCESS);
    }

    Gripper::clog.setLogLevel(log_level);

    Gripper::clog << LogLevel::DEBUG << "Basic Logger instantiated";

    clog.setLogLevel(log_level);

    clog << LogLevel::DEBUG << "STL Logger instantiated";

    Gripper::Runtime::logOptions();

    logOptions();

    clog << LogLevel::INFO << "Initializing Gripper base runtime";

    Gripper::Runtime::init();

    clog << LogLevel::INFO << "Initializing Gripper STL runtime";

    init();
}


void Gripper::stl::Runtime::processCommandLineArgs()
{
    m_stl_parser->processCommandArgs(m_argc, m_argv);

    applySTLOptions();
}


void Gripper::stl::Runtime::processConfigFile()
{
    m_stl_parser->processFile(m_config_file.c_str());

    applySTLOptions();
}


void Gripper::stl::Runtime::logOptions()
{
    clog << LogLevel::DEBUG << "Entering " << __FUNCTION__;

    if (m_compiler_parallelization)
        clog << LogLevel::INFO << "Using hand-tuned parallelization";
    else
        clog << LogLevel::INFO << "Using compiler auto-parallelization";

    clog << LogLevel::DEBUG << "Leaving " << __FUNCTION__;
}


void Gripper::stl::Runtime::init()
{
    clog << LogLevel::DEBUG << "Entering " << __FUNCTION__;

    clog << LogLevel::INFO << "Optimizing Gaunt coefficient matrix for STL back-end";

    gaunt = Multipole::stl::Gaunt::Matrix(Gripper::gaunt);

    clog << LogLevel::DEBUG << "Leaving " << __FUNCTION__;
}


void Gripper::stl::Runtime::applySTLOptions()
{
    if (m_stl_parser->getFlag('a')) m_compiler_parallelization = true;
}


bool Gripper::stl::Runtime::getCompilerParallelization() const { return m_compiler_parallelization; }


////////////////////////////
//                        //
// Global query functions //
//                        //
////////////////////////////


void Gripper::stl::initialize(int argc, char** argv) { Private::g_runtime.initialize(argc, argv); }


Multipole::stl::Spherical::Index Gripper::stl::getL() { return Multipole::stl::Spherical::Index(Private::g_runtime.getL().i); }


Gripper::ArithmeticPrecision Gripper::stl::getPrecision() { return Gripper::getPrecision(); }


bool Gripper::stl::getCompilerParallelization() { return Private::g_runtime.getCompilerParallelization(); }