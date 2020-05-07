/*
Tetrahedral Mesh Generator from SWC Morphology Data
Copyright (C) 2020 Okinawa Institute of Science and Technology, Japan.

Developer: Weiliang Chen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Timer.h>

#include <algorithm>
#include <fstream>
#include <functional>
#include <map>
#include <string>
#include <vector>
#include <tuple>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <CGAL/boost/graph/helpers.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <boost/geometry/index/rtree.hpp>


#include "utility.h"
#include "swc_morph.h"

namespace po = boost::program_options;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::FT FT;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;

typedef bg::model::point<double, 3, bg::cs::cartesian> sp_point;
typedef bg::model::box<sp_point> sp_box;
typedef std::pair<sp_box, uint> sp_value;

typedef std::function<FT(const Point&)> Function;

using namespace CGAL::parameters;

double distsqFromSurf(const Point & p, std::vector<PointR> & interpolation, bgi::rtree<sp_value, bgi::quadratic<16>> & rtree, double search_range)
{
  sp_box query(sp_point(p.x()-search_range, p.y()-search_range, p.z()-search_range), sp_point(p.x()+search_range, p.y()+search_range, p.z()+search_range));
  std::vector<sp_value> result_s;
  rtree.query(bgi::intersects(query), std::back_inserter(result_s));

  for(auto & r : result_s) {
    double distsq = interpolation[r.second].distsqFromSurf(p.x(), p.y(), p.z());
    if(distsq <= 0.0) {
      return distsq;
    }
  }
  return 1.0;
}

int main(int argc, char*argv[])
{
  #ifdef CGAL_CONCURRENT_MESH_3
    std::cout << "Running in concurrency mode.\n";
  #else
    std::cout << "Running in sequential mode.\n";
  #endif
  CGAL::Timer t_total;
  t_total.start();
  try {
    po::options_description desc("Options");
    desc.add_options()
        ("help,h", "Help message")
        ("input,i", po::value<std::string>(), "Input SWC morphology file")
        ("output,o", po::value<std::string>(), "Output tetrahedral mesh file")
        ("surfmesh", po::bool_switch()->default_value(false), "Also output surface mesh")
        ("interp-distance", po::value<double>()->default_value(0.05f, "0.05"), "Distance between two interpolation points")
        ("fc-angle", po::value<double>()->default_value(25.0f, "25.0"), "Facet criteria - Angle")
        ("fc-size", po::value<double>()->default_value(1.0f, "1.0"), "Facet criteria - Size")
        ("fc-distance", po::value<double>()->default_value(0.1f, "0.1"), "Facet criteria - Distance")
        ("cc-ratio", po::value<double>()->default_value(2.0f, "2.0"), "Cell criteria - Cell radius edge ratio")
        ("cc-size", po::value<double>()->default_value(1.0f, "1.0"), "Cell criteria - Size")
        ("odt", po::bool_switch()->default_value(false), "Enable ODT mesh optimization")
        ("odt-time", po::value<double>()->default_value(10.0f, "10.0"), "Time limit for ODT mesh optimization, in second")
        ("lloyd", po::bool_switch()->default_value(false), "Enable Lloyd mesh optimization")
        ("lloyd-time", po::value<double>()->default_value(10.0f, "10.0"), "Time limit for Lloyd mesh optimization, in second")
        ("perturb", po::bool_switch()->default_value(false), "Enable mesh sliver perturbation")
        ("perturb-time", po::value<double>()->default_value(10.0f, "10.0"), "Time limit for sliver perturbation, in second")
        ("perturb-bound", po::value<double>()->default_value(0.0f, "0.0"), "Targeted lower bound on dihedral angles of mesh cells for sliver perturbation, in degree")
        ("exude", po::bool_switch()->default_value(false), "Enable mesh sliver exudation")
        ("exude-time", po::value<double>()->default_value(10.0f, "10.0"), "Time limit for sliver exudation, in second")
        ("exude-bound", po::value<double>()->default_value(0.0f, "0.0"), "Targeted lower bound on dihedral angles of mesh cells for sliver exudation, in degree")
    ;

    po::positional_options_description p;
    p.add("input", 1);

    po::variables_map vm;
    
    po::store(po::command_line_parser(argc, argv).
              options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << "Usage: ./SWCMesher input [options]\n";
        std::cout << desc;
        return 0;
    }

    if (vm.count("input"))
    {
        std::cout << "Input SWC file is: "
              << vm["input"].as< std::string>() << "\n";
    }
    else {
      std::cerr << "Input SWC file is required.\n";
        std::cout << "Usage: ./SWCMesher input [options]\n";
      std::cerr << desc;
      return EXIT_FAILURE;
    }

    std::string output_file_name;
    std::string surfmesh_file_name;

    if (vm.count("output"))
    {
      output_file_name = vm["output"].as<std::string>();
      if(!ends_with(output_file_name, ".mesh")) {
        output_file_name += std::string(".mesh");
      }
    }
    else {
      output_file_name = vm["input"].as< std::string>();
      string_replace(output_file_name, ".swc", ".mesh");
    }
    if(vm["surfmesh"].as<bool>()) {
      surfmesh_file_name = output_file_name;
      string_replace(surfmesh_file_name, ".mesh", ".off");
      std::cout << "Surface mesh will be written to: " << surfmesh_file_name << "\n";
    }
    std::cout << "Tetrahedral mesh will be written to: " << output_file_name << "\n\n";

    std::cout << "\nReading SWC morphology...\n";
    MorphData morph(vm["input"].as< std::string>());

    std::cout << "Number of segments: " << morph.getNumSegments() << "\n";
    double length_min = morph.getMinSegmentLength();
    double length_max = morph.getMaxSegmentLength();
    std::cout << "Minimum segment length: " << length_min << "\n";
    std::cout << "Maximum segment length: " << length_max << "\n";

    double interp_distance = vm["interp-distance"].as<double>();
    std::cout << "Interpolation distance: " << interp_distance << "\n";
    if (interp_distance > length_max / 2)
    {
      std::cerr << "WARNING: Interpolation distance is larger than half of the maximum segment length.\n";
      std::cerr << "This may result in a poor interpolation.\n";
    }

    std::cout << "\nGenerating spherical interpolation...\n";
    std::vector<PointR> interpolation;
    morph.interpolate(interpolation, interp_distance);

    std::cout << "Number of spheres used for interpolation: " << interpolation.size() << "\n";

    std::cout << "\nCreating spatial index...\n";
    bgi::rtree<sp_value, bgi::quadratic<16>> rtree;

    size_t n_points = interpolation.size();
    for (auto i = 0; i < n_points; i++)
    {
      BoundingBox3 b = interpolation[i].getBBox();
      sp_box sp_b(sp_point(b.x_min, b.y_min, b.z_min), sp_point(b.x_max, b.y_max, b.z_max));
      rtree.insert(std::make_pair(sp_b, i));
    }

    std::function<FT(const Point &)> boundary_func = [&] (const Point & p) {
      return distsqFromSurf(p, interpolation, rtree, morph.getMaxRadius());
    };
    BoundingBox3 morph_bbox = morph.getBBox();
    Mesh_domain domain =
      Mesh_domain::create_implicit_mesh_domain(boundary_func,
                                               CGAL::Bbox_3(morph_bbox.x_min,
                                                            morph_bbox.y_min,
                                                            morph_bbox.z_min,
                                                            morph_bbox.x_max,
                                                            morph_bbox.y_max,
                                                            morph_bbox.z_max));

    std::cout << "\nMeshing...\n";
    Facet_criteria facet_criteria(vm["fc-angle"].as<double>(),vm["fc-size"].as<double>(), vm["fc-distance"].as<double>());
    Cell_criteria cell_criteria(vm["cc-ratio"].as<double>(),vm["cc-size"].as<double>());
    Mesh_criteria criteria(facet_criteria, cell_criteria);
    // Mesh generation

    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, 
                                        vm["odt"].as<bool>()?odt(time_limit=vm["odt-time"].as<double>()):no_odt(),
                                        vm["lloyd"].as<bool>()?lloyd(time_limit=vm["lloyd-time"].as<double>()):no_lloyd(),
                                        vm["perturb"].as<bool>()?perturb(time_limit=vm["perturb-time"].as<double>(), 
                                                                                      sliver_bound=vm["perturb-bound"].as<double>()):no_perturb(),
                                        vm["exude"].as<bool>()?exude(time_limit=vm["exude-time"].as<double>(), 
                                                                                      sliver_bound=vm["exude-bound"].as<double>()):no_exude(),
                                        CGAL::parameters::manifold());

    std::cout << "\nWriting mesh to output...\n";
  
    if(vm["surfmesh"].as<bool>()) {
      std::ofstream surf_file(surfmesh_file_name);
      c3t3.output_boundary_to_off(surf_file);
      std::cout << "Surface mesh has been written to " << surfmesh_file_name << "\n";
    }
    std::ofstream medit_file(output_file_name);
    CGAL::output_to_medit(medit_file, c3t3);

    std::cout << "Tetrahedral mesh has been written to " << output_file_name << "\n";
    std::cout << "\nTotal time cost: " << t_total.time() << " sec." << "\n";

    return EXIT_SUCCESS;
  }
  catch(std::exception& e)
  {
      std::cout << e.what() << "\n";
      return EXIT_FAILURE;
  }
}
