#ifndef MAIN_FRAME_HPP
#define MAIN_FRAME_HPP

#include "wx/WxUtility.hpp"
#include <wx/frame.h>
#include <wx/textctrl.h>
#include <wx/notebook.h>
#include <wx/treectrl.h>
#include <map>
#include <vector>

class ImageFrame;
class OsgWxFrame;
class wxTextCtrl;
class wxStreamToTextRedirector;
class wxNotebook;

class MainFrame : public wxFrame {

private:

    std::map<int, ImageFrame*> m_img_frames;
    std::map<int, OsgWxFrame*> m_sgraph_frames;
    std::map<char, wxTreeItemId> m_roots;

    wxNotebook* m_notebook;
    wxTextCtrl* m_textctrl;
    wxTreeCtrl* m_filetree;
    wxStreamToTextRedirector* m_tdirector;
    wxTextAttr m_default_style;
    int m_id;

public:

    static const wxString frame_text;
    MainFrame(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_FRAME_STYLE, const wxString& name = wxFrameNameStr);
    ~MainFrame();
    void UsrFrameClosedMessage(int id);
    void UsrSetLogMode(log_type ltype);
    void UsrLogErrorMessage(const std::string& str);
    void UsrAppendFileTree(char c, const wxArrayString& keyVal);
    void UsrUpdateFileTree(int m_id, const wxArrayString& keyVal);

private:

    // Member functions
    inline void UsrInitNotebook();
    inline void UsrInitLogPage();
    inline void UsrInitFilesPage();
    inline void UsrInitMenubar();
    bool UsrFindFrameInTheFileTree(const wxString& label, wxTreeItemId& itemId);
    void UsrDeleteFromFileTree(int m_id);

    // Event Handlers
    void OnOpenImage(wxCommandEvent& event);
    void OnModelFromSingleImage(wxCommandEvent& event);
    void OnOpenPointCloud(wxCommandEvent& event);
    void OnOpenModel(wxCommandEvent& event);
    void OnSaveLog(wxCommandEvent& event);
    void OnClearLog(wxCommandEvent& event);
    void OnClose(wxCloseEvent& event);
    void OnExit(wxCommandEvent& event);
    void OnNotebookPageChange(wxBookCtrlEvent& event);
    DECLARE_EVENT_TABLE()
};

#endif // MAIN_FRAME_HPP
