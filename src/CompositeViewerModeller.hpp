#ifndef _COMPOSITE_VIEWER_MODELLER_HPP
#define _COMPOSITE_VIEWER_MODELLER_HPP

#include <wx/app.h>

class MainFrame;

class CompositeViewerModeller: public wxApp {

public:

    virtual bool OnInit();

private:

    MainFrame* m_main_frame;

};

DECLARE_APP(CompositeViewerModeller)

#endif
