#include "ImageFrame.hpp"
#include "ImagePanel.hpp"
#include "../../MainFrame.hpp"
#include "../../wx/WxUtility.hpp"
#include "../../wx/WxGuiId.hpp"
#include <wx/sizer.h>
#include <wx/menu.h>
#include <wx/filedlg.h>
#include <wx/msgdlg.h>
#include <wx/statusbr.h>
#include <wx/toolbar.h>
#include <iostream>

BEGIN_EVENT_TABLE(ImageFrame, wxFrame)
EVT_CLOSE(ImageFrame::OnClose)
EVT_MENU(wxID_EXIT, ImageFrame::OnExit)
EVT_MENU(wxID_OPEN, ImageFrame::OnOpen)
EVT_MENU(wxID_VIEW_ORIGINAL_SIZE, ImageFrame::OnViewOriginalSize)
EVT_MENU(wxID_VIEW_MODE_FIXED, ImageFrame::OnChangeViewMode)
EVT_MENU(wxID_VIEW_MODE_VARYING, ImageFrame::OnChangeViewMode)
EVT_MENU(wxID_VIEW_IN_MODEL_WINDOW, ImageFrame::OnViewInModelWindow)
EVT_MENU(wxID_EDIT_ROTATE_CW, ImageFrame::OnRotateCW)
EVT_MENU(wxID_EDIT_ROTATE_CCW, ImageFrame::OnRotateCCW)
EVT_MENU(wxID_SAVE, ImageFrame::OnSaveImage)
EVT_MENU(wxID_MODE_DEFAULT, ImageFrame::OnChangeOperationMode)
EVT_MENU(wxID_MODE_REGION_GROWING, ImageFrame::OnChangeOperationMode)
// EVT_MENU(wxID_MODE_DRAWING, ImageFrame::OnChangeOperationMode)
// EVT_TOOL(wxID_DRAWING_MODE_CUBOID, ImageFrame::OnChangeDrawMode)
// EVT_TOOL(wxID_DRAWING_MODE_CYLIDER, ImageFrame::OnChangeDrawMode)
// EVT_TOOL(wxID_DRAWING_MODE_CONE, ImageFrame::OnChangeDrawMode)
// EVT_TOOL(wxID_DRAWING_MODE_SPHERE, ImageFrame::OnChangeDrawMode)
END_EVENT_TABLE()

ImageFrame::ImageFrame(wxWindow* parent, wxWindowID id, const wxString& title,
                       const wxPoint& pos, const wxSize& size, long style, const wxString& name) :
    wxFrame(parent, id, title, pos, size, style, name), m_id(-1), m_path("") {

    m_parent = dynamic_cast<MainFrame*>(parent);
    if(!m_parent) {
        utilityShowMessageDialog(message_type::ERROR, wxT("Dynamic cast error"));
    }
    usrInitManubar();
    // usrCreateToolBar();
    CreateStatusBar(3);
    wxBoxSizer* sizer = new wxBoxSizer(wxHORIZONTAL);
    m_img_panel = new ImagePanel(this);
    sizer->Add(m_img_panel, 1, wxEXPAND);
    SetSizer(sizer);
}

// Public Member Functions:
bool ImageFrame::UsrOpenImageFile() {
    wxFileDialog open_filedialog(this, wxT("Open image"), wxT(""), wxT(""),
                                 wxT("*.jpeg;*.jpg;*.tif;*.png;"),
                                 wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if(open_filedialog.ShowModal() == wxID_CANCEL) { return false; }

    if(!m_img_panel->UsrLoadFile(open_filedialog.GetPath())) {
        std::string str("\t-File open error: ");
        str += std::string(open_filedialog.GetPath().char_str());
        m_parent->UsrLogErrorMessage(str);
        return false;
    }
    else {
        SetTitle(GetTitle() + open_filedialog.GetFilename());
        m_path = open_filedialog.GetPath();
        usrUpdateFileTree();
        std::cout << "\t-Image file: " << open_filedialog.GetPath().char_str()
                  << " is opened successfully" << std::endl;
        return true;
    }
}
void ImageFrame::UsrSetFrameId(int id) {
    m_id = id;
    wxString title("Frame Id-");
    title += wxString(std::to_string(m_id));
    title += wxT(": ");
    SetTitle(title);
}

// Private Memeber Functions:
//void ImageFrame::usrCreateToolBar() {
//    m_toolbar = CreateToolBar();
//    wxImage img(wxT("../data/icons/cube.png"), wxBITMAP_TYPE_PNG);
//    wxBitmap cube = wxBitmap(img.Rescale(64, 64, wxIMAGE_QUALITY_HIGH));
//    img.Destroy();
//    img.LoadFile(wxT("../data/icons/cylinder.png"), wxBITMAP_TYPE_PNG);
//    wxBitmap cylinder = wxBitmap(img.Rescale(64, 64, wxIMAGE_QUALITY_HIGH));
//    img.Destroy();
//    img.LoadFile(wxT("../data/icons/sphere.png"), wxBITMAP_TYPE_PNG);
//    wxBitmap sphere = wxBitmap(img.Rescale(64, 64, wxIMAGE_QUALITY_HIGH));
//    img.Destroy();
//    img.LoadFile(wxT("../data/icons/cone.png"), wxBITMAP_TYPE_PNG);
//    wxBitmap cone = wxBitmap(img.Rescale(64, 64, wxIMAGE_QUALITY_HIGH));
//    m_toolbar->AddTool(wxID_DRAWING_MODE_CUBOID, wxT("Cube"), cube,
//                       wxT("Draw Cube"), wxITEM_RADIO);
//    m_toolbar->AddTool(wxID_DRAWING_MODE_CYLIDER, wxT("Cylinder"), cylinder,
//                       wxT("Draw Cylinder"), wxITEM_RADIO);
//    m_toolbar->AddTool(wxID_DRAWING_MODE_SPHERE, wxT("Sphere"), sphere,
//                       wxT("Draw Sphere"), wxITEM_RADIO);
//    m_toolbar->AddTool(wxID_DRAWING_MODE_CONE, wxT("Cone"), cone,
//                       wxT("Draw Cone"), wxITEM_RADIO);
//    m_toolbar->Realize();
//    m_toolbar->Hide();
//}
void ImageFrame::usrUpdateFileTree() {
    wxArrayString strArr;
    wxString str = MainFrame::frame_text + wxString(std::to_string(m_id));
    strArr.Add(str);
    strArr.Add(wxT("Name"));
    strArr.Add(m_path.AfterLast('/'));
    strArr.Add(wxT("Path"));
    strArr.Add(m_path);
    m_parent->UsrAppendFileTree('i', strArr);
}
void ImageFrame::usrInitManubar() {
    wxMenuBar* menubar = new wxMenuBar;
    wxMenu* file = new wxMenu;
    file->Append(wxID_OPEN, wxT("&Open"));
    file->Append(wxID_SAVE, wxT("&Save"));
    file->Append(wxID_EXIT, wxT("Exit"));
    menubar->Append(file, wxT("&File"));

    wxMenu* view = new wxMenu;
    wxMenu* disp_mode = new wxMenu;
    wxMenuItem* fixed = new wxMenuItem(disp_mode, wxID_VIEW_MODE_FIXED, wxString(wxT("Fixed")),
                                       wxEmptyString, wxITEM_CHECK);
    disp_mode->Append(fixed);
    wxMenuItem* varying = new wxMenuItem(disp_mode, wxID_VIEW_MODE_VARYING, wxString(wxT("Varying")),
                                         wxEmptyString, wxITEM_CHECK);
    disp_mode->Append(varying);
    fixed->Check(false);
    varying->Check(true);
    varying->Enable(false);
    view->AppendSubMenu(disp_mode, wxT("View Mode"));
    view->Append(wxID_VIEW_ORIGINAL_SIZE, wxT("Original size"));
    view->AppendSeparator();

    wxMenuItem* viewInModelWindow = new wxMenuItem(
                disp_mode, wxID_VIEW_IN_MODEL_WINDOW, wxString(wxT("View in Model Window")),
                wxEmptyString, wxITEM_CHECK);
    view->Append(viewInModelWindow);
    viewInModelWindow->Check(false);
    viewInModelWindow->Enable(false);

    menubar->Append(view, wxT("&View"));

    wxMenu* edit = new wxMenu;
    edit->Append(wxID_EDIT_ROTATE_CW, wxT("Rotate 90\u00B0 CW"));
    edit->Append(wxID_EDIT_ROTATE_CCW, wxT("Rotate 90\u00B0 CCW"));
    menubar->Append(edit, wxT("&Edit"));

    wxMenu* op_mode = new wxMenu;
    op_mode->AppendRadioItem(wxID_MODE_DEFAULT, wxT("Default"));
    op_mode->AppendRadioItem(wxID_MODE_REGION_GROWING, wxT("Region Growing"));
    // op_mode->AppendRadioItem(wxID_MODE_DRAWING, wxT("Drawing"));
    menubar->Append(op_mode, wxT("Mode"));
    SetMenuBar(menubar);
}

// Event Handlers:
void ImageFrame::OnOpen(wxCommandEvent& event) {
    UsrOpenImageFile();
}
void ImageFrame::OnViewOriginalSize(wxCommandEvent& event) {
    m_img_panel->UsrViewInOriginalSize();
}
void ImageFrame::OnChangeViewMode(wxCommandEvent& event) {
    if(event.GetId() == wxID_VIEW_MODE_FIXED) {
        m_img_panel->UsrSetDisplayMode(image_display_mode::FIXED);
        GetMenuBar()->FindItem(wxID_VIEW_MODE_VARYING)->Check(false);
        GetMenuBar()->FindItem(wxID_VIEW_MODE_VARYING)->Enable(true);
        GetMenuBar()->FindItem(wxID_VIEW_MODE_FIXED)->Enable(false);
    }
    else if(event.GetId() == wxID_VIEW_MODE_VARYING) {
        m_img_panel->UsrSetDisplayMode(image_display_mode::VARYING);
        GetMenuBar()->FindItem(wxID_VIEW_MODE_FIXED)->Check(false);
        GetMenuBar()->FindItem(wxID_VIEW_MODE_FIXED)->Enable(true);
        GetMenuBar()->FindItem(wxID_VIEW_MODE_VARYING)->Enable(false);
    }
    else {
        std::cout << "\t-Id match error when setting view mode" << std::endl;
    }
}
void ImageFrame::OnChangeOperationMode(wxCommandEvent& event) {
    if(event.GetId() == wxID_MODE_DEFAULT) {
        m_img_panel->UsrSetUIOperationMode(image_operation_mode::Default);
        std::cout << "Operation mode is set to default" << std::endl;
    }
    else if(event.GetId() == wxID_MODE_REGION_GROWING) {
        m_img_panel->UsrSetUIOperationMode(image_operation_mode::Region_Growing);
        std::cout << "Operation mode is set to region growing" << std::endl;
    }
//    else if(event.GetId() == wxID_MODE_DRAWING) {
//        m_toolbar->Show();
//        m_img_panel->UsrSetUIOperationMode(image_operation_mode::modelling);
//        std::cout << "Operation mode is set to drawing" << std::endl;
//    }
    else {
        std::cout << "Error in setting operation mode" << std::endl;
    }
// if(event.GetId() != wxID_MODE_DRAWING) { m_toolbar->Hide(); }
}
//void ImageFrame::OnChangeDrawMode(wxCommandEvent& event) {
//    if(event.GetId() == wxID_DRAWING_MODE_CUBOID) {
//        m_img_panel->UsrSetDrawingMode(draw_mode::Cuboid);
//        std::cout << "Drawing mode is set to Cuboid" << std::endl;
//    }
//    else if(event.GetId() == wxID_DRAWING_MODE_CYLIDER) {
//        m_img_panel->UsrSetDrawingMode(draw_mode::Cylinder);
//        std::cout << "Drawing mode is set to Cylinder" << std::endl;
//    }
//    else if(event.GetId() == wxID_DRAWING_MODE_CONE) {
//        m_img_panel->UsrSetDrawingMode(draw_mode::Cone);
//        std::cout << "Drawing mode is set to Cone" << std::endl;
//    }
//    else if(event.GetId() == wxID_DRAWING_MODE_SPHERE) {
//        m_img_panel->UsrSetDrawingMode(draw_mode::Sphere);
//        std::cout << "Drawing mode is set to Sphere" << std::endl;
//    }
//    else {
//        std::cout << "Error in setting drawing mode" << std::endl;
//    }
//}
void ImageFrame::OnRotateCW(wxCommandEvent& event) {
    m_img_panel->UsrRotate(rotation_type::CW);
}
void ImageFrame::OnRotateCCW(wxCommandEvent& event) {
    m_img_panel->UsrRotate(rotation_type::CCW);
}
void ImageFrame::OnExit(wxCommandEvent& event) {
    Close();
}
void ImageFrame::OnClose(wxCloseEvent& event) {
    m_parent->UsrFrameClosedMessage(m_id);
    std::cout << "\t-Image File: " << GetTitle().char_str()
              << " is closed." << std::endl;
    Destroy();
}
void ImageFrame::OnSaveImage(wxCommandEvent& event) {
    wxFileDialog save_filedialog(this, wxT("Save log to a file"), wxT(""), wxT(""),
                                 wxT("*.bmp;*.jpeg;*.jpg;*.png;*.pcx;*.pnm;*.tiff;*.xpm"),
                                 wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
    if(save_filedialog.ShowModal() == wxID_CANCEL) { return; }
    if(m_img_panel->UsrSaveImage(save_filedialog.GetPath())) {
        std::cout << "\t-Image File: " << save_filedialog.GetPath().char_str()
                  << " is successfully saved." << std::endl;
    }
    else {
        std::string str("\t-Image file save error: ");
        str += save_filedialog.GetPath().char_str();
        m_parent->UsrLogErrorMessage(str);
    }
}
void ImageFrame::OnViewInModelWindow(wxCommandEvent& event) {

}
