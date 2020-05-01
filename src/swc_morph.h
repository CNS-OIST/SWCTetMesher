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

#ifndef SWCMESHER_MORPH_H
#define SWCMESHER_MORPH_H

#include <string>
#include <map>
#include <vector>
#include <limits>
#include <memory>

struct BoundingBox3 {
    double x_min {std::numeric_limits<double>::max()};
    double x_max {std::numeric_limits<double>::min()};
    double y_min {std::numeric_limits<double>::max()};
    double y_max {std::numeric_limits<double>::min()};
    double z_min {std::numeric_limits<double>::max()};
    double z_max {std::numeric_limits<double>::min()};;

    BoundingBox3& operator+=(const BoundingBox3& rhs)
    {
        x_min = std::min(x_min, rhs.x_min);
        y_min = std::min(y_min, rhs.y_min);
        z_min = std::min(z_min, rhs.z_min);
        x_max = std::max(x_max, rhs.x_max);
        y_max = std::max(y_max, rhs.y_max);
        z_max = std::max(z_max, rhs.z_max);
        return *this;
    }
    friend BoundingBox3 operator+(BoundingBox3 lhs,
                                  const BoundingBox3& rhs)
    {
        lhs += rhs;
        return lhs;
    }
};

struct PointR {
    double x,y,z,r;
    PointR() = default;
    PointR(double x_, double y_, double z_, double r_)
    : x(x_), y(y_), z(z_), r(r_)
    {}

    BoundingBox3 getBBox() {
        BoundingBox3 bbox;
        bbox.x_min = x - r;
        bbox.x_max = x + r;
        bbox.y_min = y - r;
        bbox.y_max = y + r;
        bbox.z_min = z - r;
        bbox.z_max = z + r;
        return bbox;
    }

    double distsqFromSurf(double query_x, double query_y, double query_z) {
        double x_dist = query_x - x;
        double y_dist = query_y - y;
        double z_dist = query_z - z;
        return x_dist * x_dist + y_dist * y_dist + z_dist * z_dist - r * r;
    }
};


class Segment {
public:
    Segment(int start_point_idx, int end_point_type, PointR start_point, PointR end_point);
    void interpolate(std::vector<PointR> & interp_data, double point_distance);

    void setParent(std::shared_ptr<Segment> p) { parentSeg = p; }
    void addChild(std::shared_ptr<Segment> c) { childrenSegs.push_back(c); }
    const PointR& getStartPoint() {return start;}
    const PointR& getEndPoint() {return end;}
    double getLength() {return length; }
    BoundingBox3 getBBox() {return bbox; }

private:
    int idx;
    int type;
    PointR start;
    PointR end;
    double length;
    double x_dist;
    double y_dist;
    double z_dist;
    double r_dist;
    BoundingBox3 bbox {};
    std::shared_ptr<Segment> parentSeg;
    std::vector<std::shared_ptr<Segment>> childrenSegs;
};


class MorphData {
public:
    MorphData(const std::string & swc_file);
    void interpolate(std::vector<PointR> & interp_data, double point_distance);
    BoundingBox3 getBBox() {return bbox; }
    double getMaxRadius() { return r_max; }
    double getMinSegmentLength() { return length_min; }
    double getMaxSegmentLength() { return length_max; }
    size_t getNumSegments() { return segments.size(); }
private:
    double r_max {};
    double length_min {std::numeric_limits<double>::max()};
    double length_max {};
    BoundingBox3 bbox {};
    std::vector<std::shared_ptr<Segment>> rootSegments;
    std::map<unsigned, std::shared_ptr<Segment>> segments;
};

#endif