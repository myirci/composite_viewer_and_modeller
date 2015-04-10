#include "ModelSolver.hpp"
#include "../components/ComponentBase.hpp"
#include <iostream>
#include <algorithm>

void ModelSolver::AddComponent(ComponentBase* component) {
    m_components.push_back(component);
}

// Precondition: Geosemantic constraints must have been defined before calling this function.
void ModelSolver::Solve() {
    // step-1: construct the objective function

}

void ModelSolver::DeleteAllComponents() {

    m_components.clear();
    m_constraints.clear();
}

void ModelSolver::DeleteSelectedComponents(std::vector<int>& id_vector) {

    if(m_components.empty())
        return;

    bool found = false;
    for(int i : id_vector) {
        found = false;
        for(auto it = m_components.begin(); it != m_components.end(); ++it) {
            if((*it)->GetComponentId() == static_cast<unsigned int>(i)) {
                m_components.erase(it);
                delete_constraints_with_id(static_cast<unsigned int>(i));
                found = true;
                break;
            }
        }
        if(!found)
            std::cout << "ERROR: Component with id: " << i << " could not be found" << std::endl;
    }
}


void ModelSolver::Print() const {

    std::cout << "Num components: " << m_components.size() << std::endl;
    for(auto comp : m_components)
        comp->Print();

    std::cout << "Num constraints: " << num_constraints() << std::endl;
    for(auto c : m_constraints)
        std::cout << c << std::endl;
}

void ModelSolver::UpdateOrCreateConstraints(const unsigned int cp1, const unsigned int cp2, const std::vector<geosemantic_constraints>& gsc) {

    if(gsc.empty()) {
        for(auto it = m_constraints.begin(); it != m_constraints.end(); ++it) {
            if((it->component_1 == cp1 && it->component_2 == cp2) || (it->component_1 == cp2 && it->component_2 == cp1)) {
                m_constraints.erase(it);
                break;
            }
        }
        return;
    }
    else {
        bool updated = false;
        for(auto it = m_constraints.begin(); it != m_constraints.end(); ++it) {
            if((it->component_1 == cp1 && it->component_2 == cp2) || (it->component_1 == cp2 && it->component_2 == cp1)) {
                it->constraints.clear();
                std::copy(gsc.begin(), gsc.end(), std::back_inserter(it->constraints));
                updated = true;
                break;
            }
        }
        if(!updated)
            m_constraints.push_back(geosemcon(cp1, cp2, gsc));
    }
}

void ModelSolver::GetConstraints(const unsigned int cp1, const unsigned int cp2, std::vector<geosemantic_constraints>& gsc) const {

    for(auto c : m_constraints) {
        if((c.component_1 == cp1 && c.component_2 == cp2) || (c.component_1 == cp2 && c.component_2 == cp1)) {
            std::copy(c.constraints.begin(), c.constraints.end(), std::back_inserter(gsc));
            break;
        }
    }
}

int ModelSolver::num_constraints() const {

    int count = 0;
    for(auto c : m_constraints)
        count += c.constraints.size();
    return count;
}

void ModelSolver::delete_constraints_with_id(unsigned int comp_id) {

    if(m_constraints.empty()) { return; }
    bool found = false;
    for(auto it = m_constraints.begin(); it != m_constraints.end(); ++it) {
        if(it->component_1 == comp_id || it->component_2 == comp_id) {
            m_constraints.erase(it);
            found = true;
            break;
        }
    }
    if(found)
        delete_constraints_with_id(comp_id);
}

void ModelSolver::construct_geosemantic_constraints() {

    for(auto it1 = m_constraints.begin(); it1 != m_constraints.end(); ++it1) {
        for(auto it2 = it1->constraints.begin(); it2 != it1->constraints.end(); ++it2) {
            switch(*it2) {
            case geosemantic_constraints::parallel:
                std::cout << "Not Implemented!" << std::endl;
                break;
            case geosemantic_constraints::orthogonal:
                std::cout << "Not Implemented!" << std::endl;
                break;
            case geosemantic_constraints::collinear_axis_endpoints:
                std::cout << "Not Implemented!" << std::endl;
                break;
            case geosemantic_constraints::overlapping_axis_endpoints:
                std::cout << "Not Implemented!" << std::endl;
                break;
            case geosemantic_constraints::coplanar_axis_endpoints:
                std::cout << "Not Implemented!" << std::endl;
                break;
            case geosemantic_constraints::coplanar_axes:
                std::cout << "Not Implemented!" << std::endl;
                break;
            default:
                break;
            }
        }
    }
}

