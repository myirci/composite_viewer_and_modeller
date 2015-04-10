#include "Constraints.hpp"
#include <iostream>

geosemcon::geosemcon(unsigned int cp1, unsigned int cp2, const std::vector<geosemantic_constraints>& gsc) :
    component_1(cp1), component_2(cp2), constraints(gsc) { }

std::string to_string(geosemantic_constraints c) {

    switch(c) {
    case geosemantic_constraints::parallel:
        return "parallel";
        break;
    case geosemantic_constraints::orthogonal:
        return "orthogonal";
        break;
    case geosemantic_constraints::collinear_axis_endpoints:
        return "collinear axis endpoints";
        break;
    case geosemantic_constraints::overlapping_axis_endpoints:
        return "overlapping axis endpoints";
        break;
    case geosemantic_constraints::coplanar_axis_endpoints:
        return "coplanar axis endpoints";
        break;
    case geosemantic_constraints::coplanar_axes:
        return "coplanar axes";
        break;
    default:
        std::cout << "Error: geosemantic constraint does not match!" << std::endl;
        return "";
        break;
    }
}

int to_int(geosemantic_constraints c) {
    return static_cast<int>(c);
}

std::ostream& operator<<(std::ostream& out, const geosemcon& gsc) {
    out << "<" << gsc.component_1 << ", " << gsc.component_2 << ", ";
    for(auto it = gsc.constraints.begin(); it != gsc.constraints.end(); ++it) {
        out << to_string(*it);
        if(std::next(it) != gsc.constraints.end())
            out << ", ";
    }
    out << ">" << std::endl;
}

