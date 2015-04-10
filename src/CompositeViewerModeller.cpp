#include "CompositeViewerModeller.hpp"
#include "MainFrame.hpp"
#include <wx/image.h>

IMPLEMENT_APP(CompositeViewerModeller)

bool CompositeViewerModeller::OnInit() {
    wxInitAllImageHandlers(); // initialize all available image handlers
    m_main_frame = new MainFrame(NULL, wxID_ANY, wxT("Composite Viewer and Modeller"), wxPoint(50,50), wxSize(600,450));
    m_main_frame->Show();
    return true;
}
