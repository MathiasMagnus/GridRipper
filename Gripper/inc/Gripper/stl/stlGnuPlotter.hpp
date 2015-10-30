///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Plotter class to create GnuPlot compatible datasets and plotting scripts  //
//                                                                           //
// Author: Nagy-Egri Máté Ferenc                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef STLGNUPLOTTER_HPP
#define STLGNUPLOTTER_HPP


// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/stl/stlLogger.hpp>
#include <Gripper/stl/stlRadialVector.hpp>
#include <Gripper/stl/stlSphericalVector.hpp>
#include <Gripper/stl/stlExpansionField.hpp>

// Standard C++ includes
#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include <sstream>

namespace Gripper
{
    namespace stl
    {
        namespace Graphics
        {
            enum Color
            {
                Undefined,
                Violet,
                Purple,
                Blue,
                Turquoise,
                Green,
                Yellow,
                Orange,
                Red,
            };

            enum Connectivity
            {
                None,
                Line,
                Dots,
                Segments,
                DottedLine
            };

            enum Marker
            {
                None,
                Point,
                X,
                Cross,
                Circle,
                Square,
                Triangle
            };


            class DataSet2D
            {
            public:

                DataSet2D() = default;
                DataSet2D(const DataSet2D& in) = default;
                DataSet2D(DataSet2D&& sry) = default;
                ~DataSet2D() = default;

                template <typename ET, typename VT> DataSet(const Multipole::stl::Radial::Expression<ET, VT>& in) { for (auto r = in.extent().initial(); in.extent().contains(r); ++r) m_data.push_back(std::make_pair(static_cast<double>(r.r), static_cast<double>(in.at(r)))); }
                template <typename ET, typename VT> DataSet(const Multipole::stl::Spherical::Expression<ET, VT>& in) { for (auto i = in.extent().initial(); in.extent().contains(i); ++i) m_data.push_back(std::make_pair(static_cast<double>(r.r), static_cast<double>(in.at(i)))); }

                template <typename ET, typename VT> void setData(const Multipole::stl::Radial::Expression<ET, VT>& in)
                {
                    m_data.clear();
                    for (auto r = in.extent().initial(); in.extent().contains(r); ++r)
                        m_data.push_back(std::make_pair(static_cast<double>(r.r), static_cast<double>(in.at(r))));
                }
                template <typename ET, typename VT> void setData(const Multipole::stl::Spherical::Expression<ET, VT>& in)
                {
                    m_data.clear();
                    for (auto i = in.extent().initial(); in.extent().contains(i); ++i)
                        m_data.push_back(std::make_pair(static_cast<double>(r.r), static_cast<double>(in.at(i))));
                }
                void setName(const char* name) { m_name = name; }
                void setColor(Color color) { m_color = color; }
                void setConnectivity(Connectivity conn) { m_conn = conn; }
                void setMarker(Marker marker) { m_marker = marker; }

                std::vector<std::pair<double, double>>& getData() { return m_data; }
                const std::vector<std::pair<double, double>>& getData() const { return m_data; }
                std::string getName() { return m_name; }
                Color getColor() { return m_color; }
                Connectivity getConnectivity() { return m_conn; }
                Marker getMarker() { return m_marker; }

            private:

                std::vector<std::pair<double, double>> m_data;
                std::string m_name;
                Color m_color;
                Connectivity m_conn;
                Marker m_marker;
            };


            class GnuPlot2D
            {
            public:

                GnuPlot2D() = default;
                GnuPlot2D(const GnuPlot2D& in) = default;
                GnuPlot2D(GnuPlot2D&& src) = default;
                ~GnuPlot2D() = default;

                GnuPlot2D(const char* prefix) : m_data(std::string(prefix) + ".dat"), m_script(std::string(prefix) + ".plt") {}

                bool is_open() { return m_data.is_open() && m_script.is_open(); }
                void open(const char* prefix)
                {
                    m_data.open(std::string(prefix) + ".dat");
                    m_script.open(std::string(prefix) + ".plt");
                }
                void close()
                {
                    // Sort DataSets according to their lengths
                    std::sort(m_data_sets.begin(), m_data_sets.end(), [](const DataSet2D& a, const DataSet2D& b){ return a.getData().size() > b.getData().size(); });

                    // Print them all to the datafile
                    for (auto line = 0; line < m_data_sets.at(0).getData().size(); ++line)
                    {
                        for (auto& set : m_data_sets)
                        {
                            if (line < set.getData().size())
                                m_data << set.getData().at(line).first << "\t" << set.getData().at(line).second << "\t";
                        }
                        m_data << std::endl;
                    }

                    // Find out plot ranges
                    double x_min, x_max, y_min, y_max;
                    for (auto& set : m_data_sets)
                        ;

                    m_data.close();
                    m_script.close();
                }

                GnuPlot2D& operator<<(const DataSet2D& in)
                {
                    m_data_sets.push_back(in);
                    return *this;
                }

            private:

                std::vector<DataSet2D> m_data_sets;

                std::ofstream m_data;
                std::ofstream m_script;
            };

        } // namespace Graphics

    } // namespace stl

} // namespace Gripper

#endif // STLGNUPLOTTER_HPP