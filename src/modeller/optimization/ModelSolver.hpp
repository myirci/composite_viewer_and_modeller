#ifndef MODEL_SOLVER_HPP
#define MODEL_SOLVER_HPP

#include "Constraints.hpp"
#include <memory>
#include <vector>
#include <list>

class ComponentBase;
class ProjectionParameters;

class ModelSolver {
public:
    ModelSolver() { }
    void Solve();
    void AddComponent(ComponentBase* component);
    void DeleteAllComponents();

    void DeleteSelectedComponents(std::vector<int>& id_vector);
    void UpdateOrCreateConstraints(const unsigned int cp1, const unsigned int cp2, const std::vector<geosemantic_constraints>& gsc);
    void GetConstraints(const unsigned int cp1, const unsigned int cp2, std::vector<geosemantic_constraints>& gsc) const;
    void Print() const;

private:
    std::list<ComponentBase*> m_components;
    std::list<geosemcon> m_constraints;
    std::shared_ptr<ProjectionParameters> m_pp;

    inline int num_constraints() const;
    void delete_constraints_with_id(unsigned int comp_id);
    void construct_geosemantic_constraints();
};

#endif // MODEL_SOLVER_HPP
