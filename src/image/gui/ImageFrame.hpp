#ifndef _IMAGE_FRAME_HPP
#define _IMAGE_FRAME_HPP

#include <wx/frame.h>

class ImagePanel;
class MainFrame;
class wxFileDialog;
class wxToolBar;

class ImageFrame : public wxFrame {
public:
    ImageFrame(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos = wxDefaultPosition,
               const wxSize& size = wxDefaultSize, long style = wxDEFAULT_FRAME_STYLE, const wxString& name = wxFrameNameStr);
    bool UsrOpenImageFile();
    void UsrSetFrameId(int id);
private:
    ImagePanel* m_img_panel;
    MainFrame* m_parent;
    int m_id;
    wxString m_path;
    wxToolBar* m_toolbar;
    // Member functions
    inline void usrInitManubar();
    inline void usrUpdateFileTree();
    // inline void usrCreateToolBar();
    // Event Handlers:
    void OnOpen(wxCommandEvent& event);
    void OnViewOriginalSize(wxCommandEvent& event);
    void OnViewInModelWindow(wxCommandEvent& event);
    void OnViewGradientImage(wxCommandEvent& event);
    void OnRotateCW(wxCommandEvent& event);
    void OnRotateCCW(wxCommandEvent& event);
    void OnExit(wxCommandEvent& event);
    void OnClose(wxCloseEvent& event);
    void OnSaveImage(wxCommandEvent& event);
    void OnChangeOperationMode(wxCommandEvent& event);
    // void OnChangeDrawMode(wxCommandEvent& event);
    void OnChangeViewMode(wxCommandEvent& event);
    DECLARE_EVENT_TABLE()
};
#endif
