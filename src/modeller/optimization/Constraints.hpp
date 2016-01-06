#ifndef CONSTRAINTS_HPP
#define CONSTRAINTS_HPP

#include <vector>
#include <string>

enum class geosemantic_constraints : unsigned char {
    parallel,
    orthogonal,
    collinear_axis_endpoints,
    overlapping_axis_endpoints,
    coplanar_axis_endpoints,
    coplanar_axes
};

std::string to_string(geosemantic_constraints c);
int to_int(geosemantic_constraints c);

struct geosemcon {
    geosemcon(unsigned int cp1, unsigned int cp2, const std::vector<geosemantic_constraints>& gsc);
    unsigned int component_1;
    unsigned int component_2;
    std::vector<geosemantic_constraints> constraints;
};

std::ostream& operator<<(std::ostream& out, const geosemcon& gsc);

#endif // CONSTRAINTS_HPP
