#ifndef COMPONENT_RELATIONS_DIALOG_HPP
#define COMPONENT_RELATIONS_DIALOG_HPP

#include "../optimization/Constraints.hpp"
#include <wx/dialog.h>
#include <map>
#include <vector>

class wxCheckBox;
class OsgWxFrame;

class ComponentRelationsDialog : public wxDialog {
public:
    ComponentRelationsDialog(wxWindow* parent,
                             const wxString& title);
    void UpdateConstraintsList(
            const std::vector<geosemantic_constraints>& gsc);
    void ResetConstraintsList();
private:
    OsgWxFrame* m_parent;
    wxButton* apply_button;
    std::map<int, wxCheckBox*> constraints;
    void OnClose(wxCommandEvent& event);
    void OnApply(wxCommandEvent& event);
    void print_enabled_constraints() const;
    inline void get_enabled_constraints(
            std::vector<geosemantic_constraints>& enabled_constraints) const;
    DECLARE_EVENT_TABLE()
};

#endif // COMPONENT_RELATIONS_DIALOG_HPP
