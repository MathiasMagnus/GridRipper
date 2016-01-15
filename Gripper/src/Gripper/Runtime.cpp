#include <Gripper/Runtime.hpp>
/*
namespace Gripper
{
    namespace Private
    {
        // Global runtime instance
        EXPORT Runtime g_runtime;

    } // namespace Private

} // namespace Gripper


Gripper::Runtime::Runtime() :
    m_base_parser(new AnyOption()),
    m_config_file("config.cfg"),
    m_L_max(3),
    m_precision(ArithmeticPrecision::DoublePrecision)
{
    m_base_parser->addUsage("");
    m_base_parser->addUsage("Usage: <executable> [SWITCHES...]");
    m_base_parser->addUsage("");
    m_base_parser->addUsage(" SWITCHES may be stand-alone switches (flags), or may have parameters (options).");
    m_base_parser->addUsage("");
    m_base_parser->addUsage(" Basic flags: <token>");
    m_base_parser->addUsage("");
    m_base_parser->addUsage(" -h  --help        Prints out this help");
    m_base_parser->addUsage(" -q  --quiet       Run without commmand-line output");
    m_base_parser->addUsage(" -v  --verbose     Run showing general information");
    m_base_parser->addUsage(" -p  --profiling   Run showing profiling information");
    m_base_parser->addUsage("     --debug       Run showing debug information");
    m_base_parser->addUsage("     --trace       Run showing trace information");
    m_base_parser->addUsage("     --insane      Run using insane log level");
    m_base_parser->addUsage("");
    m_base_parser->addUsage(" Basic options: <token> <value>");
    m_base_parser->addUsage("");
    m_base_parser->addUsage(" -c  --config      Filename to read configuration from [default: config.cfg]");
    m_base_parser->addUsage(" -o  --output      Filename to redirect log to");
    m_base_parser->addUsage(" -L  --L_max       L_max target in Gaunt-coefficient calculation");

    m_base_parser->setFlag("help", 'h');
    m_base_parser->setFlag("quiet", 'q');
    m_base_parser->setFlag("verbose", 'v');
    m_base_parser->setFlag("profiling", 'p');
    m_base_parser->setFlag("debug");
    m_base_parser->setFlag("trace");
    m_base_parser->setFlag("insane");

    m_base_parser->setOption("config", 'c');
    m_base_parser->setOption("output", 'o');
    m_base_parser->setOption("L_max", 'L');
}


Gripper::Runtime::~Runtime()
{
}


void Gripper::Runtime::initialize(int argc, char** argv)
{
    if (argc != 0)
    {
        m_argc = argc;
        m_argv = argv;
    }

    processCommandLineArgs();

    processConfigFile();

    if (m_base_parser->getFlag('h'))
    {
        m_base_parser->printUsage();

        exit(EXIT_SUCCESS);
    }

    clog.setLogLevel(log_level);

    clog << LogLevel::DEBUG << "Basic Logger initialized";

    logOptions();

    clog << LogLevel::INFO << "Initializing Gripper base runtime";

    init();
}



void Gripper::Runtime::processCommandLineArgs()
{
    m_base_parser->processCommandArgs(m_argc, m_argv);

    applyOptions();

}


void Gripper::Runtime::processConfigFile()
{
    m_base_parser->processFile(m_config_file.c_str());

    applyOptions();
}


void Gripper::Runtime::logOptions()
{
    ENTERING

    clog << LogLevel::INFO << "Using L_max = " << m_L_max.i;

    switch (m_precision)
    {
    case ArithmeticPrecision::SinglePrecision:
        clog << LogLevel::INFO << "Using precision = float";
        break;
    case ArithmeticPrecision::DoublePrecision:
        clog << LogLevel::INFO << "Using precision = double";
        break;
    default:
        clog << LogLevel::CRIT << "Using precision = UNKOWN";
        break;
    }

    LEAVING
}


void Gripper::Runtime::init()
{
    ENTERING

    clog << LogLevel::INFO << "Calculating Gaunt coefficient matrix";

    gaunt = Multipole::Gaunt::Matrix(m_L_max);

    LEAVING
}


void Gripper::Runtime::applyOptions()
{
    ENTERING

    if (m_base_parser->getFlag('q')) Gripper::log_level = Gripper::LogLevel::NONE;

    if (m_base_parser->getFlag('v')) Gripper::log_level = Gripper::LogLevel::INFO;

    if (m_base_parser->getFlag('p')) Gripper::log_level = Gripper::LogLevel::TIMING;

    if (m_base_parser->getFlag("debug")) Gripper::log_level = Gripper::LogLevel::DEBUG;

    if (m_base_parser->getFlag("trace")) Gripper::log_level = Gripper::LogLevel::ALL;

    if (m_base_parser->getFlag("insane")) Gripper::log_level = Gripper::LogLevel::INSANE;

    if (m_base_parser->getValue('c') != nullptr) m_config_file = std::string(m_base_parser->getValue('c'));

    if (m_base_parser->getValue('L') != nullptr) m_L_max = std::atoi(m_base_parser->getValue('L'));

    LEAVING
}


Multipole::Spherical::Index Gripper::Runtime::getL() const { return m_L_max; }


Gripper::ArithmeticPrecision Gripper::Runtime::getPrecision() const { return m_precision; }


////////////////////////////
//                        //
// Global query functions //
//                        //
////////////////////////////


void Gripper::initialize(int argc, char** argv) { Private::g_runtime.initialize(argc, argv); }


Multipole::Spherical::Index Gripper::getL() { return Private::g_runtime.getL(); }


Gripper::ArithmeticPrecision Gripper::getPrecision() { return Private::g_runtime.getPrecision(); }
*/