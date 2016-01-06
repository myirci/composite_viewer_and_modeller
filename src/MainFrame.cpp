#include "MainFrame.hpp"
#include "image/gui/ImageFrame.hpp"
#include "osg/OsgWxFrame.hpp"
#include "wx/WxGuiId.hpp"
#include <wx/statusbr.h>
#include <wx/menu.h>
#include <wx/filedlg.h>
#include <wx/datetime.h>
#include <wx/panel.h>
#include <wx/sizer.h>
#include <iostream>

BEGIN_EVENT_TABLE(MainFrame, wxFrame)
    EVT_CLOSE(MainFrame::OnClose)
    EVT_MENU(wxID_EXIT, MainFrame::OnExit)
    EVT_MENU(wxID_OPEN_IMAGE, MainFrame::OnOpenImage)
    EVT_MENU(wxID_OPEN_POINT_CLOUD, MainFrame::OnOpenPointCloud)
    EVT_MENU(wxID_OPEN_MODEL, MainFrame::OnOpenModel)
    EVT_MENU(wxID_SAVE, MainFrame::OnSaveLog)
    EVT_MENU(wxID_EDIT_CLEAR_LOG, MainFrame::OnClearLog)
    EVT_MENU(wxID_MODEL_IMAGE, MainFrame::OnModelFromSingleImage)
    EVT_NOTEBOOK_PAGE_CHANGED(wxID_ANY, MainFrame::OnNotebookPageChange)
END_EVENT_TABLE()

const wxString MainFrame::frame_text = wxT("Frame Id: ");
MainFrame::MainFrame(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style, const wxString& name) :
    wxFrame(parent, id, title, pos, size, style, name), m_id(0) {

    SetTitle(wxT("Composite Viewer and Modeller"));
    UsrInitNotebook();
    UsrInitMenubar();
    CreateStatusBar(1);
    SetStatusText(wxT("CVM"), 0);
    Centre();
}

MainFrame::~MainFrame() {

    delete m_tdirector;
}

void MainFrame::OnOpenImage(wxCommandEvent& event) {

    ImageFrame* imgFrame = new ImageFrame(this, wxID_ANY, wxT(""), wxPoint(50, 50), wxSize(600, 450));
    imgFrame->UsrSetFrameId(++m_id);
    if(imgFrame->UsrOpenImageFile()) {
        m_img_frames.insert(std::pair<int, ImageFrame*>(m_id, imgFrame));
        imgFrame->Show();
    }
    else {
        --m_id;
        imgFrame->Destroy();
    }
}

void MainFrame::OnOpenPointCloud(wxCommandEvent& event) {

    std::cout << "This menu has not been implemented yet!" << std::endl;
}

void MainFrame::OnOpenModel(wxCommandEvent& event) {

    OsgWxFrame* sgFrame = new OsgWxFrame(this, wxPoint(50, 50), wxSize(600, 450), operation_mode::displaying);
    sgFrame->UsrSetFrameId(++m_id);
    if(sgFrame->UsrOpenModelFile()) {
        m_sgraph_frames.insert(std::pair<int, OsgWxFrame*>(m_id, sgFrame));
        sgFrame->Show();
    }
    else {
        --m_id;
        sgFrame->Destroy();
    }
}

void MainFrame::OnModelFromSingleImage(wxCommandEvent& event) {

    OsgWxFrame* sgFrame = new OsgWxFrame(this, wxPoint(50, 50), wxSize(600, 450), operation_mode::modelling);
    sgFrame->UsrSetFrameId(++m_id);
    if(sgFrame->UsrOpenImageFile()) {
        m_sgraph_frames.insert(std::pair<int, OsgWxFrame*>(m_id, sgFrame));
        sgFrame->Show();
    }
    else {
        --m_id;
        sgFrame->Destroy();
    }
}

void MainFrame::OnSaveLog(wxCommandEvent& event) {

    wxFileDialog save_filedialog(this, wxT("Save log to a file"), wxT(""), wxT(""), wxT("*.txt"), wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
    if(save_filedialog.ShowModal() == wxID_CANCEL) { return; }
    m_textctrl->SaveFile(save_filedialog.GetPath());
}

void MainFrame::OnClearLog(wxCommandEvent& event) {
    m_textctrl->Clear();
}

void MainFrame::OnExit(wxCommandEvent& event) {
    Close();
}

void MainFrame::OnClose(wxCloseEvent& event) {
    if(utilityQuestionDialoag(wxT("Do you want to exit?"))) {
        Destroy();
    }
}

void MainFrame::OnNotebookPageChange(wxNotebookEvent& event) {
    std::cout << "INFO: Page Changed: " << event.GetSelection() << std::endl;
}

void MainFrame::UsrInitNotebook() {
    m_notebook = new wxNotebook(this, wxID_ANY);
    UsrInitLogPage();
    UsrInitFilesPage();
}

void MainFrame::UsrInitLogPage() {

    wxPanel* panel1 = new wxPanel(m_notebook, wxID_ANY);
    long tstyle = (wxTE_READONLY | wxTE_MULTILINE | wxTE_LEFT | wxTE_BESTWRAP | wxTE_RICH);
    m_textctrl = new wxTextCtrl(panel1, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, tstyle);
    m_tdirector = new wxStreamToTextRedirector(m_textctrl);
    m_default_style = m_textctrl->GetDefaultStyle();
    m_textctrl->SetDefaultStyle(wxTextAttr(*wxCYAN));
    wxDateTime dt = dt.Now();
    std::cout << dt.Format().char_str() << std::endl;
    m_textctrl->SetDefaultStyle(m_default_style);
    std::cout << "*****************************************************" << std::endl;
    wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);
    sizer->Add(m_textctrl, 1, wxEXPAND);
    panel1->SetSizer(sizer);
    m_notebook->AddPage(panel1, wxT("Log"));
}

void MainFrame::UsrInitFilesPage() {

    wxPanel* panel2 = new wxPanel(m_notebook, wxID_ANY);
    m_filetree = new wxTreeCtrl(panel2, wxID_ANY);
    wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);
    sizer->Add(m_filetree, 1, wxEXPAND);
    panel2->SetSizer(sizer);
    m_notebook->AddPage(panel2, wxT("Files"));
    m_roots['r'] = m_filetree->AddRoot(wxT("Files"));
    m_roots['i'] = m_filetree->AppendItem(m_roots['r'], wxT("Images"));
    m_roots['o'] = m_filetree->AppendItem(m_roots['r'], wxT("Oriented Images"));
    m_roots['m'] = m_filetree->AppendItem(m_roots['r'], wxT("Models"));
    m_roots['p'] = m_filetree->AppendItem(m_roots['r'], wxT("Point Clouds"));
}

void MainFrame::UsrInitMenubar() {

    wxMenuBar* menubar = new wxMenuBar;
    wxMenu* file = new wxMenu;
    wxMenu* open = new wxMenu;
    open->Append(wxID_OPEN_IMAGE, wxT("Image"));
    open->Append(wxID_OPEN_POINT_CLOUD, wxT("Point Cloud"));
    open->Append(wxID_OPEN_MODEL, wxT("Model"));
    file->AppendSubMenu(open, wxT("&Open"));
    file->Append(wxID_SAVE, wxT("&Save Log"));
    file->Append(wxID_EXIT, wxT("Exit"));
    menubar->Append(file, wxT("&File"));

    wxMenu* edit = new wxMenu;
    edit->Append(wxID_EDIT_CLEAR_LOG, wxT("Clear Log"));
    menubar->Append(edit, wxT("&Edit"));

    wxMenu* model = new wxMenu;
    model->Append(wxID_MODEL_IMAGE, wxT("Model from Single Image"));
    menubar->Append(model, wxT("&Model"));

    SetMenuBar(menubar);
}

void MainFrame::UsrLogErrorMessage(const std::string& str) {

    UsrSetLogMode(log_type::ERROR);
    std::cout << str << std::endl;
    UsrSetLogMode(log_type::DEFAULT);
}

void MainFrame::UsrAppendFileTree(char c, const wxArrayString& keyVal) {

    wxTreeItemId parent;
    wxTreeItemId child;
    parent = m_filetree->AppendItem(m_roots[c], keyVal.Item(0));
    for(int i = 1; i < keyVal.GetCount(); ++i) {
        child = m_filetree->AppendItem(parent, keyVal.Item(i));
        ++i;
        if(keyVal.Item(i) != wxEmptyString) {
            m_filetree->AppendItem(child, keyVal.Item(i));
        }
    }
}

void MainFrame::UsrUpdateFileTree(int m_id, const wxArrayString& keyVal) {

}

void MainFrame::UsrDeleteFromFileTree(int m_id) {

    wxString str = frame_text + wxString(std::to_string(m_id));
    wxTreeItemId itemId;
    if(UsrFindFrameInTheFileTree(str, itemId)) {
        m_filetree->Delete(itemId);
    }
    else {
        std::string estr = "\t-Frame deletion error: Frame Id-" +
                std::to_string(m_id) + " could not find!";
        UsrLogErrorMessage(estr);
    }
}

bool MainFrame::UsrFindFrameInTheFileTree(const wxString& label, wxTreeItemId& itemId) {

    wxTreeItemId child;
    wxTreeItemIdValue cookie;
    for(auto it = m_roots.begin(); it != m_roots.end(); ++it) {
        if(it->first == 'r') { continue; }
        child = m_filetree->GetFirstChild(it->second, cookie);
        if(child.IsOk()) {
            if(m_filetree->GetItemText(child) == label) {
                itemId = child;
                return true;
            }
        }
        else {
            continue;
        }
        while(true) {
            child = m_filetree->GetNextChild(it->second, cookie);
            if(child.IsOk()) {
                if(m_filetree->GetItemText(child) == label) {
                    itemId = child;
                    return true;
                }
            }
            else { break; }
        }
    }
    return false;
}

void MainFrame::UsrFrameClosedMessage(int id) {

    std::map<int, ImageFrame*>::iterator it1;
    it1 = m_img_frames.find(id);
    if(it1 != m_img_frames.end()) {
        m_img_frames.erase(it1);
    }

    std::map<int, OsgWxFrame*>::iterator it2;
    it2 = m_sgraph_frames.find(id);
    if(it2 != m_sgraph_frames.end()) {
        m_sgraph_frames.erase(it2);
    }

    UsrDeleteFromFileTree(id);
}

void MainFrame::UsrSetLogMode(log_type ltype) {

    if(ltype == log_type::DEFAULT)      m_textctrl->SetDefaultStyle(m_default_style);
    else if(ltype == log_type::ERROR)   m_textctrl->SetDefaultStyle(wxTextAttr(*wxRED));
    else if(ltype == log_type::WARNING) m_textctrl->SetDefaultStyle(wxTextAttr(*wxBLUE));
}

