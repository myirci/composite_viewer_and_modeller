#include "ComponentRelationsDialog.hpp"
#include "../../wx/WxGuiId.hpp"
#include "../../osg/OsgWxFrame.hpp"
#include "../optimization/ModelSolver.hpp"
#include <wx/panel.h>
#include <wx/statbox.h>
#include <wx/checkbox.h>
#include <wx/sizer.h>
#include <wx/button.h>

#include <iostream>

BEGIN_EVENT_TABLE(ComponentRelationsDialog, wxDialog)
EVT_BUTTON(wxID_COMPONENT_RELATIONS_APPLY_BUTTON, ComponentRelationsDialog::OnApply)
EVT_BUTTON(wxID_COMPONENT_RELATIONS_CLOSE_BUTTON, ComponentRelationsDialog::OnClose)
END_EVENT_TABLE()

ComponentRelationsDialog::ComponentRelationsDialog(wxWindow* parent, const wxString& title) :
    wxDialog(parent, wxID_ANY, title, wxDefaultPosition, wxSize(250, 250)) {

    m_parent = dynamic_cast<OsgWxFrame*>(parent);
    if(!m_parent)
        std::cout << "ERROR: dynamic cast error" << std::endl;

    wxBoxSizer* vmainbox = new wxBoxSizer(wxVERTICAL);
    wxPanel* panel = new wxPanel(this);
    new wxStaticBox(panel, wxID_ANY, wxT("Geo-semantic constraints"), wxPoint(5, 5), wxSize(220, 180));

    std::vector<geosemantic_constraints> gsc = {
        geosemantic_constraints::parallel,
        geosemantic_constraints::orthogonal,
        geosemantic_constraints::collinear_axis_endpoints,
        geosemantic_constraints::overlapping_axis_endpoints,
        geosemantic_constraints::coplanar_axis_endpoints,
        geosemantic_constraints::coplanar_axes };

    int y = 25;
    for(auto c : gsc) {
        constraints.insert(std::pair<int, wxCheckBox*>(to_int(c), new wxCheckBox(panel, wxID_ANY, wxString(to_string(c)), wxPoint(10, y))));
        y += 25;
    }

    vmainbox->Add(panel);
    wxBoxSizer* hbox = new wxBoxSizer(wxHORIZONTAL);
    apply_button = new wxButton(this, wxID_COMPONENT_RELATIONS_APPLY_BUTTON, wxT("Apply"), wxDefaultPosition, wxSize(70, 30));
    wxButton* close_button = new wxButton(this, wxID_COMPONENT_RELATIONS_CLOSE_BUTTON, wxT("Close"), wxDefaultPosition, wxSize(70, 30));
    hbox->Add(apply_button);
    hbox->Add(close_button);
    vmainbox->Add(hbox);
    SetSizer(vmainbox);
}

void ComponentRelationsDialog::ResetConstraintsList() {

    for(auto it = constraints.begin(); it != constraints.end(); ++it)
        if(it->second->IsChecked())
            it->second->SetValue(false);
}

void ComponentRelationsDialog::UpdateConstraintsList(const std::vector<geosemantic_constraints>& gsc) {

    ResetConstraintsList();
    for(auto c : gsc) {
        auto it = constraints.find(static_cast<int>(c));
        if(it != constraints.end())
            it->second->SetValue(true);
    }
}

void ComponentRelationsDialog::OnApply(wxCommandEvent& event) {

    // get the selected components
    std::vector<unsigned int> selected_components;
    m_parent->UsrGetSelectedComponentIds(selected_components);
    if(selected_components.size() != 2) {
        std::cout << "INFO: Number of selected components must be two" << std::endl;
        return;
    }

    // get the enabled constraints
    std::vector<geosemantic_constraints> gsc;
    get_enabled_constraints(gsc);

    // update the constraints between the selected components
    ModelSolver* solver = m_parent->UsrGetModelSolver();
    solver->UpdateOrCreateConstraints(selected_components[0], selected_components[1], gsc);
}

void ComponentRelationsDialog::OnClose(wxCommandEvent& event) {

    m_parent->UsrDisplayComponentRelationsDialog();
}

void ComponentRelationsDialog::print_enabled_constraints() const {

    std::vector<geosemantic_constraints> enabled_constraints;
    get_enabled_constraints(enabled_constraints);
    for(auto c : enabled_constraints)
        std::cout << to_string(c) << std::endl;

}

void ComponentRelationsDialog::get_enabled_constraints(std::vector<geosemantic_constraints>& enabled_constraints) const {

    for(auto it = constraints.begin(); it != constraints.end(); ++it)
        if(it->second->IsChecked())
           enabled_constraints.push_back(static_cast<geosemantic_constraints>(it->first));
}
