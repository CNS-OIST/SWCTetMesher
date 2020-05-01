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

#include "swc_morph.h"
#include "utility.h"

#include <iostream>
#include <math.h>

#include <boost/algorithm/string.hpp>

Segment::Segment(int start_point_idx, int end_point_type, PointR start_point, PointR end_point)
: idx (start_point_idx) 
, type(end_point_type)
, start(start_point)
, end(end_point) {
    x_dist = end.x - start.x;
    y_dist = end.y - start.y;
    z_dist = end.z - start.z;
    r_dist = end.r - start.r;
    length = sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist);
    bbox = start.getBBox() + end.getBBox();
}

void Segment::interpolate(std::vector<PointR> & interp_data, double point_distance)
{
    uint n_points = length / point_distance;

    interp_data.push_back(start);
    interp_data.push_back(end);

    if(type == 1 || n_points <= 2) {
        return;
    }
    else {
        for(auto i = 1; i < n_points; i++)
        {
            double x = start.x + static_cast<double>(i) / n_points * (x_dist);
            double y = start.y + static_cast<double>(i) / n_points * (y_dist);
            double z = start.z + static_cast<double>(i) / n_points * (z_dist);
            double r = start.r + static_cast<double>(i) / n_points * (r_dist);
            interp_data.push_back(PointR(x,y,z,r));
        }
    }
}

MorphData::MorphData(const std::string & swc_file)
{
    std::vector<std::string> data_strings;
	bool result = getFileContent(swc_file, data_strings);
	if(!result) {
        throw std::runtime_error("ERROR: Unable to read swc file\n");
    }

    for(auto& line : data_strings)
    {
        if(line.size() > 0 && line.find("#") != 0) {
			std::vector<std::string> fields;
            boost::split(fields, line, boost::is_any_of(" "));
            if (fields.size() == 7)
            {
                int idx = std::stoi(fields[0]);
                int type = std::stoi(fields[1]);
                double x = std::stod(fields[2]);
                double y = std::stod(fields[3]);
                double z = std::stod(fields[4]);
                double r = std::stod(fields[5]);
                int previous = std::stoi(fields[6]);

                r_max = std::max(r, r_max);

                PointR start;
                PointR end(x,y,z,r);
                std::shared_ptr<Segment> seg_ptr;
                if(previous == -1) {

                    seg_ptr = std::make_shared<Segment>(Segment(idx, type, start, end));
                    rootSegments.push_back(seg_ptr);
                }
                else {
                    auto previous_seg = segments.find(previous);
                    if(previous_seg != segments.end()) {
                        PointR start = previous_seg->second->getEndPoint();
                        seg_ptr = std::make_shared<Segment>(Segment(idx, type, start, end));
                    }
                    else {
                        throw std::runtime_error("ERROR: Previous segment point " + std::to_string(previous) + " not in record.\n");
                    }
                }
                segments[idx] = seg_ptr;
                bbox += seg_ptr->getBBox();
                
                if(type != 1) {
                    double length = seg_ptr->getLength();
                    if(length > 0.0) {
                        length_max = std::max(length_max, length);
                        length_min = std::min(length_min, length);
                    }
                }
            }
            else {
                throw std::runtime_error(std::string("ERROR: Unable to parser swc data line ") + line + "\n");
            }
		}
    }
}

void MorphData::interpolate(std::vector<PointR> & interp_data, double point_distance)
{
    for(auto & seg : segments) {
        seg.second->interpolate(interp_data, point_distance);
    }
}