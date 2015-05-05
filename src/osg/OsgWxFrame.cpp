#include <memory>

#include "../geometry/Circle3D.hpp"
#include "OsgWxFrame.hpp"
#include "OsgWxGLCanvas.hpp"
#include "OsgWxGraphicsWindow.hpp"
#include "OsgUtility.hpp"
#include "../MainFrame.hpp"
#include "../wx/WxUtility.hpp"
#include "../wx/WxGuiId.hpp"
#include "../modeller/ImageModeller.hpp"
#include "../modeller/CoordinateTransformations.hpp"
#include "../modeller/gui/ComponentRelationsDialog.hpp"
#include "../modeller/optimization/ModelSolver.hpp"

#include <wx/menu.h>
#include <wx/filedlg.h>

#include <osg/ValueObject>
#include <osg/Group>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osg/MatrixTransform>
#include <osgGA/TrackballManipulator>
#include <osgViewer/ViewerEventHandlers>
#include <osg/Image>

#include <OriLight/Orientation.hpp>

BEGIN_EVENT_TABLE(OsgWxFrame, wxFrame)
EVT_IDLE(OsgWxFrame::OnIdle)
EVT_CLOSE(OsgWxFrame::OnClose)
EVT_MENU(wxID_EXIT, OsgWxFrame::OnExit)
EVT_MENU(wxID_OSG_OPEN_MODEL, OsgWxFrame::OnOpenModel)
EVT_MENU(wxID_OSG_OPEN_IMAGE, OsgWxFrame::OnOpenOrientedImage)
EVT_MENU(wxID_MODES_PERSPECTIVE_PROJECTION, OsgWxFrame::OnToggleProjectionMode)
EVT_MENU(wxID_MODES_ORTHOGRAPHIC_PROJECTION, OsgWxFrame::OnToggleProjectionMode)
EVT_MENU(wxID_MODES_RENDER_MODE_POINT, OsgWxFrame::OnToggleRenderMode)
EVT_MENU(wxID_MODES_RENDER_MODE_WIREFRAME, OsgWxFrame::OnToggleRenderMode)
EVT_MENU(wxID_MODES_RENDER_MODE_FILL, OsgWxFrame::OnToggleRenderMode)
EVT_MENU(wxID_MODES_RENDER_FACE_FRONT, OsgWxFrame::OnToggleRenderFaceMode)
EVT_MENU(wxID_MODES_RENDER_FACE_BACK, OsgWxFrame::OnToggleRenderFaceMode)
EVT_MENU(wxID_MODES_RENDER_FACE_FRONT_AND_BACK, OsgWxFrame::OnToggleRenderFaceMode)
EVT_MENU(wxID_MODES_RENDER_PLANAR_SECTIONS, OsgWxFrame::OnToggleRenderType)
EVT_MENU(wxID_MODES_RENDER_PLANAR_AND_VERTICAL_SECTIONS, OsgWxFrame::OnToggleRenderType)
EVT_MENU(wxID_MODES_RENDER_TRIANGLE_STRIP, OsgWxFrame::OnToggleRenderType)
EVT_MENU(wxID_MODES_RENDER_TRIANGLE_FAN, OsgWxFrame::OnToggleRenderType)
EVT_MENU(wxID_PRINT_PROJECTION_AND_MODEL_VIEW_MATRICES, OsgWxFrame::OnPrintProjectionMatrix)
EVT_MENU(wxID_EDIT_CLEAR_VIEW, OsgWxFrame::OnDeleteModel)
EVT_MENU(wxID_MODEL_CONSTRAINTS_NO_SPINE_CONSTRAINTS, OsgWxFrame::OnToggleModellingConstraints)
EVT_MENU(wxID_MODEL_CONSTRAINTS_STRAIGHT_PLANAR_SPINE, OsgWxFrame::OnToggleModellingConstraints)
EVT_MENU(wxID_MODEL_CONSTRAINTS_PLANAR_SPINE_POINTS, OsgWxFrame::OnToggleModellingConstraints)
EVT_MENU(wxID_MODEL_CONSTRAINTS_NO_SECTION_CONSTRAINTS, OsgWxFrame::OnToggleModellingConstraints)
EVT_MENU(wxID_MODEl_CONSTRAINTS_CONSTANT_SECTIONS, OsgWxFrame::OnToggleModellingConstraints)
EVT_MENU(wxID_MODEL_CONSTRAINTS_LINEARLY_SCALED_SECTIONS, OsgWxFrame::OnToggleModellingConstraints)
EVT_MENU(wxID_MODEL_SPINE_DRAWING_MODE_CONTINUOUS, OsgWxFrame::OnToggleSpineDrawingMode)
EVT_MENU(wxID_MODEL_SPINE_DRAWING_MODE_PIECEWISE_LINEAR, OsgWxFrame::OnToggleSpineDrawingMode)
EVT_MENU(wxID_MODEL_SAVE_COMPONENT, OsgWxFrame::OnSaveLastComponent)
EVT_MENU(wxID_MODEL_SAVE_MODEL, OsgWxFrame::OnSaveModel)
EVT_MENU(wxID_MODEL_DELETE_SELECTED_COMPONENTS, OsgWxFrame::OnDeleteSelectedComponents)
EVT_MENU(wxID_MODEL_DELETE_MODEL, OsgWxFrame::OnDeleteModel)
EVT_MENU(wxID_MODEL_SOLVE, OsgWxFrame::OnSolveModel)
EVT_MENU(wxID_VIEW_DISPLAY_LOCAL_FRAMES, OsgWxFrame::OnDisplayLocalFrames)
EVT_MENU(wxID_VIEW_DISPLAY_COORDINATE_FRAME, OsgWxFrame::OnDisplayWorldCoordinateFrame)
EVT_MENU(wxID_VIEW_DISPLAY_VERTEX_NORMALS, OsgWxFrame::OnDisplayVertexNormals)
EVT_MENU(wxID_VIEW_DISPLAY_SECTION_NORMALS, OsgWxFrame::OnDisplaySectionNormals)
EVT_MENU(wxID_VIEW_DISPLAY_IMAGE, OsgWxFrame::OnToggleImageDisplay)
EVT_MENU(wxID_VIEW_DISPLAY_GRADIENT_IMAGE, OsgWxFrame::OnToggleImageDisplay)
EVT_MENU(wxID_WINDOWS_COMPONENT_RELATIONS, OsgWxFrame::OnDisplayComponentRelationsDialog)
EVT_MENU(wxID_MODES_OPERATION_MODE_DISPLAY, OsgWxFrame::OnToggleUIOperationMode)
EVT_MENU(wxID_MODES_OPERATION_MODE_MODELLING, OsgWxFrame::OnToggleUIOperationMode)
END_EVENT_TABLE()

OsgWxFrame::OsgWxFrame(wxWindow* parent, const wxPoint& pos, const wxSize& size, operation_mode md) :
    wxFrame(parent, wxID_ANY, wxEmptyString, pos, size, wxDEFAULT_FRAME_STYLE, wxEmptyString),
    m_render_mode(osg::PolygonMode::FILL), m_render_face(osg::PolygonMode::FRONT_AND_BACK),
    m_uiopmode(md), m_component_relations_win(new ComponentRelationsDialog(this, wxT("Component Relations"))),
    m_ppp(nullptr), m_bgcam(nullptr), m_bgeode(nullptr), m_model(nullptr), m_world_frame(nullptr), m_id(-1),
    m_imgdisp_mode(background_image_display_mode::image){

    m_parent = dynamic_cast<MainFrame*>(parent);
    if(!m_parent) utilityShowMessageDialog(message_type::ERROR, wxT("Dynamic cast error"));
    usrInitMenubar();
    CreateStatusBar(2);

    // Array of integers. With this parameter you can set the device context  attributes associated to this window.
    // This array is zero-terminated: it should be set up using wxGL_FLAGS constants. If a constant should be followed
    // by a value, put it in the next array position
    int* attributes = new int[7];
    attributes[0] = int(WX_GL_DOUBLEBUFFER);
    attributes[1] = WX_GL_RGBA;
    attributes[2] = WX_GL_DEPTH_SIZE;
    attributes[3] = 8;
    attributes[4] = WX_GL_STENCIL_SIZE;
    attributes[5] = 8;
    attributes[6] = 0;

    wxSize client_size = GetClientSize();
    m_canvas = new OsgWxGLCanvas(this, wxID_ANY, attributes, wxDefaultPosition, client_size, wxBORDER_SIMPLE, wxT("osgviewerWX"));
    delete[] attributes;

    m_graphics_window = new OsgWxGraphicsWindow(m_canvas);
    m_canvas->UsrSetGraphicsWindow(m_graphics_window);

    m_root = new osg::Group;
    usrSetPolygonMode(m_root.get());
    m_root->addChild(m_canvas->UsrGetSelectionBoxes());

    m_viewer = new osgViewer::Viewer;
    m_viewer->getCamera()->setGraphicsContext(m_graphics_window);
    m_viewer->getCamera()->setViewport(0, 0, client_size.GetWidth(), client_size.GetHeight());
    m_viewer->addEventHandler(new osgViewer::StatsHandler);
    // m_viewer->setThreadingModel(osgViewer::Viewer::SingleThreaded);
    m_viewer->setThreadingModel(osgViewer::Viewer::AutomaticSelection);
    m_viewer->setSceneData(m_root.get());

    if(m_uiopmode == operation_mode::displaying) {
        m_camera_manipulator = new osgGA::TrackballManipulator();
        m_viewer->setCameraManipulator(m_camera_manipulator.get());
        m_camera_manipulator->setByMatrix(osg::Matrixd::identity());
        usrEnableModellingMenus(false);
    }

    if(m_uiopmode == operation_mode::modelling)
        usrEnableModellingMenus(true);
}

// Public Member Functions:
bool OsgWxFrame::UsrOpenOrientedImageFile() {

    wxFileDialog open_filedialog(this, wxT("Open Image"), wxT(""), wxT(""), wxT("*.xml"), wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if(open_filedialog.ShowModal() == wxID_CANCEL) return false;
    if(usrLoadOrientationFile(open_filedialog.GetPath())) {
        m_path = open_filedialog.GetPath();
        usrUpdateFileTree('o');
        std::cout << "\t-Image file: " << open_filedialog.GetPath().char_str() << " is successfully opened in model window" << std::endl;
        return true;
    }
    else {
        std::stringstream ss;
        ss << "\t-File open error: " << open_filedialog.GetPath().char_str();
        UsrLogErrorMessage(ss.str());
        return false;
    }
}

bool OsgWxFrame::UsrOpenModelFile() {

    wxFileDialog open_filedialog(this, wxT("Open Model"), wxT(""), wxT(""), wxT("*.osg;*.obj;*.ply"), wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if(open_filedialog.ShowModal() == wxID_CANCEL) return false;

    if(usrLoadModelFile(open_filedialog.GetPath())) {
        SetTitle(open_filedialog.GetFilename());
        m_path = open_filedialog.GetPath();
        usrUpdateFileTree('m');
        std::cout << "\t-Model file: " << open_filedialog.GetPath().char_str() << " is opened successfully" << std::endl;
        return true;
    }
    else {
        std::stringstream ss;
        ss << "\t-File open error: " << open_filedialog.GetPath().char_str();
        UsrLogErrorMessage(ss.str());
        return false;
    }
}

bool OsgWxFrame::UsrOpenImageFile(){

    wxFileDialog open_filedialog(this, wxT("Open image"), wxT(""), wxT(""), wxT("*.jpeg;*.jpg;*.tif;*.png;"), wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if(open_filedialog.ShowModal() == wxID_CANCEL) return false;

    if(usrLoadImageFile(open_filedialog.GetPath())) {
        SetTitle(open_filedialog.GetFilename());
        m_path = open_filedialog.GetPath();
        usrUpdateFileTree('m');
        std::cout << "\t-Image File: " << open_filedialog.GetPath().char_str() << " is opened successfully" << std::endl;
        return true;
    }
    else {
        std::stringstream ss;
        ss << "\t-File open error: " << open_filedialog.GetPath().char_str();
        UsrLogErrorMessage(ss.str());
        return false;
    }
}

void OsgWxFrame::UsrSetPerspectiveProjectionMatrix(double fovy, double aspect, double near, double far) {

    m_viewer->getCamera()->setProjectionMatrixAsPerspective(fovy, aspect, near, far);
}

void OsgWxFrame::UsrSetOrthographicProjectionMatrix(double left, double right, double bottom, double top, double near, double far) {

    m_viewer->getCamera()->setProjectionMatrixAsOrtho(left, right, bottom, top, near, far);
}

void OsgWxFrame::UsrAddSelectableComponent(osg::Node* node, unsigned int component_id) {

    node->setUserValue("Selection", false);
    node->setUserValue("Selection_Box_Id", -1);
    node->setUserValue("Component_Id", static_cast<int>(component_id));
    m_model->addChild(node);
}

osg::Camera* OsgWxFrame::UsrGetMainCamera() {
    return m_viewer->getCamera();
}

osg::Camera* OsgWxFrame::UsrGetBackgroundCamera() {
    return m_bgcam.get();
}

osg::Geode* OsgWxFrame::UsrGetBackgroundNode() {
    return m_bgeode.get();
}

operation_mode OsgWxFrame::UsrGetUIOperationMode() const {
    return m_uiopmode;
}

void OsgWxFrame::UsrSetUIOperationMode(operation_mode mode) {

    m_uiopmode = mode;
    if(m_uiopmode == operation_mode::displaying) {

        // initialize camera manipulator
        if(!m_camera_manipulator.valid()) m_camera_manipulator = new osgGA::TrackballManipulator();
        m_viewer->setCameraManipulator(m_camera_manipulator.get());
        m_camera_manipulator->setByMatrix(osg::Matrixd::identity());

        // enable displaying menus, disable drawing menus
        usrEnableModellingMenus(false);
        wxMenuItem* item = GetMenuBar()->FindItem(wxID_MODES_OPERATION_MODE_DISPLAY);
        if(!item->IsChecked()) item->Check();
        std::cout << "\t-Operation mode is changed to Display" << std::endl;
    }
    else if(m_uiopmode == operation_mode::modelling) {

        // disable camera manipulator
        m_camera_manipulator->setByMatrix(osg::Matrixd::identity());
        m_viewer->updateTraversal();
        m_viewer->setCameraManipulator(nullptr);

        // enable drawing menus, disable displaying menus
        usrEnableModellingMenus(true);
        wxMenuItem* item = GetMenuBar()->FindItem(wxID_MODES_OPERATION_MODE_MODELLING);
        if(!item->IsChecked()) item->Check();
        std::cout << "\t-Operation mode is changed to Modelling" << std::endl;
    }
}

void OsgWxFrame::UsrDisplayComponentRelationsDialog() {

    m_component_relations_win->Show(!m_component_relations_win->IsShown());
    GetMenuBar()->FindItem(wxID_WINDOWS_COMPONENT_RELATIONS)->Toggle();
}

void OsgWxFrame::UsrGetSelectedComponentIds(std::vector<unsigned int>& ids) {

    bool selected = false;
    int id = 0;
    for(int i = 0; i < m_model->getNumChildren(); ++i) {
        osg::Node* child = m_model->getChild(i);
        if(child->getUserValue("Selection", selected)) {
            if(selected) {
                if(child->getUserValue("Component_Id", id))
                    ids.push_back(static_cast<unsigned int>(id));
                else
                    UsrLogErrorMessage("Component_Id could not be found!");
            }
        }
        else {
            UsrLogErrorMessage("Selection could not be found!");
        }
    }
}

ModelSolver* OsgWxFrame::UsrGetModelSolver() {
    return m_canvas->UsrGetModeller()->GetModelSolver();
}

void OsgWxFrame::UsrUpdateGeosemanticConstraints() {

    std::vector<unsigned int> selections;
    UsrGetSelectedComponentIds(selections);
    if(selections.size() == 2) {
        std::vector<geosemantic_constraints> gsc;
        UsrGetModelSolver()->GetConstraints(selections[0], selections[1], gsc);
        m_component_relations_win->UpdateConstraintsList(gsc);
    }
    else {
        m_component_relations_win->ResetConstraintsList();
    }
}

// Private Member Functions:
void OsgWxFrame::usrInitMenubar() {

    wxMenuBar* menubar = new wxMenuBar;
    wxMenu* file = new wxMenu;
    wxMenu* open = new wxMenu;
    open->Append(wxID_OSG_OPEN_MODEL, wxT("Model"));
    open->Append(wxID_OSG_OPEN_IMAGE, wxT("Oriented Image"));
    file->AppendSubMenu(open, wxT("&Open"));
    file->Append(wxID_EXIT, wxT("Exit"));
    menubar->Append(file, wxT("&File"));

    wxMenu* edit = new wxMenu;
    edit->Append(wxID_EDIT_CLEAR_VIEW, wxT("Clear Models"));
    menubar->Append(edit, wxT("&Edit"));

    wxMenu* view = new wxMenu;
    view->AppendCheckItem(wxID_VIEW_DISPLAY_COORDINATE_FRAME, wxT("Display World Coordinate Frame"));
    view->AppendCheckItem(wxID_VIEW_DISPLAY_LOCAL_FRAMES, wxT("Display Local Frames"));
    view->AppendCheckItem(wxID_VIEW_DISPLAY_SECTION_NORMALS, wxT("Display Section Normals"));
    view->AppendCheckItem(wxID_VIEW_DISPLAY_VERTEX_NORMALS, wxT("Display Vertex Normals"));
    view->AppendSeparator();
    view->AppendRadioItem(wxID_VIEW_DISPLAY_IMAGE, wxT("Image"));
    view->AppendRadioItem(wxID_VIEW_DISPLAY_GRADIENT_IMAGE, wxT("Gradient Image"));
    menubar->Append(view, wxT("&View"));

    wxMenu* modes = new wxMenu;
    wxMenu* render_mode = new wxMenu;
    render_mode->AppendRadioItem(wxID_MODES_RENDER_FACE_FRONT, wxT("Front"));
    render_mode->AppendRadioItem(wxID_MODES_RENDER_FACE_BACK, wxT("Back"));
    render_mode->AppendRadioItem(wxID_MODES_RENDER_FACE_FRONT_AND_BACK, wxT("Front and Back"));
    render_mode->AppendSeparator();
    render_mode->AppendRadioItem(wxID_MODES_RENDER_MODE_POINT, wxT("Point"));
    render_mode->AppendRadioItem(wxID_MODES_RENDER_MODE_WIREFRAME, wxT("Wireframe"));
    render_mode->AppendRadioItem(wxID_MODES_RENDER_MODE_FILL, wxT("Fill"));
    render_mode->AppendSeparator();
    render_mode->AppendRadioItem(wxID_MODES_RENDER_TRIANGLE_STRIP, wxT("Triangle Strip"));
    render_mode->AppendRadioItem(wxID_MODES_RENDER_TRIANGLE_FAN, wxT("Triangle Fan"));
    render_mode->AppendRadioItem(wxID_MODES_RENDER_PLANAR_SECTIONS, wxT("Planar Sections"));
    render_mode->AppendRadioItem(wxID_MODES_RENDER_PLANAR_AND_VERTICAL_SECTIONS, wxT("Planar and Vertical Sections"));
    modes->AppendSubMenu(render_mode, wxT("Render Mode"));

    wxMenu* op_modes = new wxMenu;
    if(m_uiopmode == operation_mode::displaying) {
        op_modes->AppendRadioItem(wxID_MODES_OPERATION_MODE_DISPLAY, wxT("Display"));
        op_modes->AppendRadioItem(wxID_MODES_OPERATION_MODE_MODELLING, wxT("Modelling"));
    }
    else {
        op_modes->AppendRadioItem(wxID_MODES_OPERATION_MODE_MODELLING, wxT("Modelling"));
        op_modes->AppendRadioItem(wxID_MODES_OPERATION_MODE_DISPLAY, wxT("Display"));
    }
    modes->AppendSubMenu(op_modes, wxT("Operation Mode"));

    wxMenu* projection_mode = new wxMenu;
    projection_mode->AppendRadioItem(wxID_MODES_PERSPECTIVE_PROJECTION, wxT("Perspective"));
    projection_mode->AppendRadioItem(wxID_MODES_ORTHOGRAPHIC_PROJECTION, wxT("Othographic"));
    modes->AppendSubMenu(projection_mode, wxT("Projection Mode"));

    menubar->Append(modes, wxT("&Modes"));

    wxMenu* model = new wxMenu;
    model->Append(wxID_MODEL_SOLVE, wxT("Solve"));
    wxMenu* model_save = new wxMenu;
    model_save->Append(wxID_MODEL_SAVE_MODEL, wxT("Save Model"));
    model_save->Append(wxID_MODEL_SAVE_COMPONENT, wxT("Save Last Component"));
    model->AppendSubMenu(model_save, wxT("Save"));

    wxMenu* model_delete = new wxMenu;
    model_delete->Append(wxID_MODEL_DELETE_SELECTED_COMPONENTS, wxT("Delete Selected Components"));
    model_delete->Append(wxID_MODEL_DELETE_MODEL, wxT("Delete Model"));
    model->AppendSubMenu(model_delete, wxT("Delete"));

    wxMenu* spncstrnts = new wxMenu;
    spncstrnts->AppendRadioItem(wxID_MODEL_CONSTRAINTS_NO_SPINE_CONSTRAINTS, wxT("None"));
    spncstrnts->AppendRadioItem(wxID_MODEL_CONSTRAINTS_PLANAR_SPINE_POINTS, wxT("Constant Depth"));
    spncstrnts->AppendRadioItem(wxID_MODEL_CONSTRAINTS_STRAIGHT_PLANAR_SPINE, wxT("Straight Planar Spine"));
    model->AppendSubMenu(spncstrnts, wxT("Spine Constraints"));

    wxMenu* sctncnstrnts = new wxMenu;
    sctncnstrnts->AppendRadioItem(wxID_MODEL_CONSTRAINTS_NO_SECTION_CONSTRAINTS, wxT("None"));
    sctncnstrnts->AppendRadioItem(wxID_MODEl_CONSTRAINTS_CONSTANT_SECTIONS, wxT("Constant Sections"));
    sctncnstrnts->AppendRadioItem(wxID_MODEL_CONSTRAINTS_LINEARLY_SCALED_SECTIONS, wxT("Linearly Scaled Sections"));
    model->AppendSubMenu(sctncnstrnts, wxT("Section Constraints"));

    wxMenu* sp_dr_mode = new wxMenu;
    sp_dr_mode->AppendRadioItem(wxID_MODEL_SPINE_DRAWING_MODE_CONTINUOUS, wxT("Continuous"));
    wxMenuItem* mi = sp_dr_mode->AppendRadioItem(wxID_MODEL_SPINE_DRAWING_MODE_PIECEWISE_LINEAR, wxT("Piecewise Linear"));
    mi->Check(true);
    model->AppendSubMenu(sp_dr_mode, wxT("Spine Drawing Mode"));
    menubar->Append(model, wxT("Model"));

    wxMenu* windows = new wxMenu;
    windows->AppendCheckItem(wxID_WINDOWS_COMPONENT_RELATIONS, wxT("Component Relations"));
    menubar->Append(windows, wxT("&Windows"));

    wxMenu* info = new wxMenu;
    info->Append(wxID_PRINT_PROJECTION_AND_MODEL_VIEW_MATRICES, wxT("Print Projection and Model-View Matrices"));
    menubar->Append(info, wxT("Info"));

    SetMenuBar(menubar);

    menubar->FindItem(wxID_MODES_RENDER_MODE_FILL)->Check(true);
    menubar->FindItem(wxID_MODES_RENDER_FACE_FRONT_AND_BACK)->Check(true);

    if(m_uiopmode == operation_mode::modelling)       usrEnableModellingMenus(true);
    else if(m_uiopmode == operation_mode::displaying) usrEnableModellingMenus(false);

}

bool OsgWxFrame::usrLoadModelFile(const wxString& fpath) {

    // load the scene.
    osg::ref_ptr<osg::Node> loadedModel = osgDB::readNodeFile(std::string(fpath.mb_str()));
    if (!loadedModel) { return false; }
    loadedModel->setUserValue("Selection", false);
    loadedModel->setUserValue("Selection_Box_Id", -1);
    m_root->addChild(loadedModel.get());
    return true;
}

bool OsgWxFrame::usrLoadOrientationFile(const wxString& fpath) {

    std::shared_ptr<Orientation> ori(new Orientation());
    if(ori->Read(fpath.ToStdString())) {
        return true;
    }
    return false;
}

bool OsgWxFrame::usrLoadImageFile(const wxString& fpath) {

    // delete the current background camera and the node
    if(m_bgcam.valid())  m_bgcam = nullptr;
    if(m_bgeode.valid()) m_bgeode = nullptr;

    // create a textured quad with the given image as texture
    wxSize img_size;
    osg::Geode* bg_image;
    if(m_imgdisp_mode == background_image_display_mode::image) {
        osg::Image* image = osgDB::readImageFile(fpath.ToStdString());
        if(!image) {
            std::cout << "Image file cannot be opened!" << std::endl;
            return false;
        }
        bg_image = create_textured_quad(image, img_size);
    }
    else if(m_imgdisp_mode == background_image_display_mode::gradient_image) {

        // if gradient image does not exist, generate the gradient image
        std::string grad_img_path = utilityInsertAfter(fpath, wxT('.'), wxT("_grad"));
        std::ifstream infile(grad_img_path);
        if(!infile.good()) {
            OtbFloatVectorImageType::Pointer img = LoadImage<OtbFloatVectorImageType>(fpath.ToStdString());
            OtbImageType::Pointer gimg = GradientMagnitudeImage(img);
            SaveImage<OtbImageType>(gimg, grad_img_path);
        }

        osg::Image* grad_image = osgDB::readImageFile(grad_img_path);
        if(!grad_image) {
            std::cout << "Image file cannot be opened!" << std::endl;
            return false;
        }
        bg_image = create_textured_quad(grad_image, img_size);
    }
    else {
        std::cerr << "Image display mode error" << std::endl;
        return false;
    }

    // create the back ground camera and and the textured quad under this camera
    m_bgcam = create_background_camera(0, img_size.x, 0, img_size.y);
    m_bgcam->addChild(bg_image);

    // crete the node for 2D user drawings on the image
    m_bgeode = new osg::Geode;
    m_bgcam->addChild(m_bgeode.get());

    // set the polygon render mode for the background camera
    usrSetCustomPolygonMode(m_bgcam.get(), osg::PolygonMode::FILL, osg::PolygonMode::FRONT_AND_BACK);

    // add background camera to the root node
    m_root->addChild(m_bgcam.get());

    // set client size to the image size
    SetClientSize(img_size);

    // perspective projection parameters of the main camera
    m_ppp = std::shared_ptr<CoordinateTransformations>(new CoordinateTransformations(45.0, img_size.x, img_size.y, 1.0, 100.0));

    // set the properties of the main camera

    m_viewer->getCamera()->setViewport(0, 0, img_size.x, img_size.y);
    m_viewer->getCamera()->setComputeNearFarMode(osg::Camera::DO_NOT_COMPUTE_NEAR_FAR);
    m_viewer->getCamera()->setViewMatrix(osg::Matrixd::identity());
    m_viewer->getCamera()->setProjectionMatrixAsPerspective(m_ppp->fovy, m_ppp->aspect, m_ppp->near, m_ppp->far);

    // initialize the modeller: this must be executed after the initialization of the m_bgeode.
    m_canvas->UsrInitializeModeller(m_ppp, fpath);

    // create the model node and add it to the root node
    m_model = new osg::Group;
    m_root->addChild(m_model.get());

    return true;
}

void OsgWxFrame::usrSetPolygonMode(osg::Node* node) {

    osg::ref_ptr<osg::PolygonMode> pm = new osg::PolygonMode;
    pm->setMode(m_render_face, m_render_mode);
    node->getOrCreateStateSet()->setAttribute(pm.get());
}

void OsgWxFrame::usrSetCustomPolygonMode(osg::Node *node, osg::PolygonMode::Mode mode, osg::PolygonMode::Face face) {

    osg::ref_ptr<osg::PolygonMode> pm = new osg::PolygonMode;
    pm->setMode(face, mode);
    node->getOrCreateStateSet()->setAttribute(pm.get());
}

void OsgWxFrame::usrUpdateFileTree(char type) {

    wxArrayString strArr;
    wxString str = MainFrame::frame_text + wxString(std::to_string(m_id));
    strArr.Add(str);
    strArr.Add(wxT("Name"));
    strArr.Add(m_path.AfterLast('/'));
    strArr.Add(wxT("Path"));
    strArr.Add(m_path);
    m_parent->UsrAppendFileTree(type, strArr);
}

void OsgWxFrame::usrEnableModellingMenus(bool flag) {

    GetMenuBar()->FindItem(wxID_MODEL_CONSTRAINTS_NO_SPINE_CONSTRAINTS)->Enable(flag);
    GetMenuBar()->FindItem(wxID_MODEL_CONSTRAINTS_PLANAR_SPINE_POINTS)->Enable(flag);
    GetMenuBar()->FindItem(wxID_MODEL_CONSTRAINTS_STRAIGHT_PLANAR_SPINE)->Enable(flag);
    GetMenuBar()->FindItem(wxID_MODEL_CONSTRAINTS_NO_SECTION_CONSTRAINTS)->Enable(flag);
    GetMenuBar()->FindItem(wxID_MODEl_CONSTRAINTS_CONSTANT_SECTIONS)->Enable(flag);
    GetMenuBar()->FindItem(wxID_MODEL_CONSTRAINTS_LINEARLY_SCALED_SECTIONS)->Enable(flag);
    GetMenuBar()->FindItem(wxID_MODEL_SPINE_DRAWING_MODE_CONTINUOUS)->Enable(flag);
    GetMenuBar()->FindItem(wxID_MODEL_SPINE_DRAWING_MODE_PIECEWISE_LINEAR)->Enable(flag);
}

// Event Handlers:
void OsgWxFrame::OnIdle(wxIdleEvent& event) {

    if(!m_viewer->isRealized()) return;
    m_viewer->frame();
    event.RequestMore();
}

void OsgWxFrame::OnClose(wxCloseEvent& event) {

    m_parent->UsrFrameClosedMessage(m_id);
    std::cout << "\t-Model File: " << GetTitle().char_str() << " is closed." << std::endl;
    Destroy();
}

void OsgWxFrame::OnExit(wxCommandEvent& event) {
    Close();
}

void OsgWxFrame::OnOpenOrientedImage(wxCommandEvent& event) {
    UsrOpenOrientedImageFile();
}

void OsgWxFrame::OnOpenModel(wxCommandEvent& event) {
    UsrOpenModelFile();
}

void OsgWxFrame::OnToggleRenderType(wxCommandEvent& event) {

    rendering_type rtype;
    switch(event.GetId()) {
    case wxID_MODES_RENDER_PLANAR_SECTIONS:
        rtype = rendering_type::planar_sections;
        std::cout << "\t-Component render type is changed to Planar Sections" << std::endl;
        break;
    case wxID_MODES_RENDER_PLANAR_AND_VERTICAL_SECTIONS:
        rtype = rendering_type::planar_and_vertical_sections;
        std::cout << "\t-Component render type is changed to Planar and Vertical Sections" << std::endl;
        break;
    case wxID_MODES_RENDER_TRIANGLE_STRIP:
        rtype = rendering_type::triangle_strip;
        std::cout << "\t-Component render type is changed to Triangle Strip" << std::endl;
        break;
    case wxID_MODES_RENDER_TRIANGLE_FAN:
        rtype = rendering_type::triangle_fan;
        std::cout << "\t-Component render type is changed to Triangle Fan" << std::endl;
        break;
    default:
        UsrLogErrorMessage("Id match error when setting component render type");
        break;
    }

    m_canvas->UsrSetRenderingType(rtype);

    for(size_t i = 0; i < m_model->getNumChildren(); ++i) {
        GeneralizedCylinder* gcyl = dynamic_cast<GeneralizedCylinder*>(m_model->getChild(i)->asGroup());
        if(gcyl) gcyl->ChangeRenderingType(rtype);
    }
    m_root->dirtyBound();
}

void OsgWxFrame::OnToggleProjectionMode(wxCommandEvent& event) {

    switch (event.GetId()) {
    case wxID_MODES_PERSPECTIVE_PROJECTION:
        m_viewer->getCamera()->setProjectionMatrixAsPerspective(m_ppp->fovy, m_ppp->aspect, m_ppp->near, m_ppp->far);
        std::cout << "\t-Projection mode is changed to perspective" << std::endl;
        break;
    case wxID_MODES_ORTHOGRAPHIC_PROJECTION:

        double right = static_cast<double>(m_ppp->width) / 2.0;
        double left = - right;
        double top = static_cast<double>(m_ppp->height) / 2.0;
        double bottom = -top;

        m_viewer->getCamera()->setProjectionMatrixAsOrtho(left, right, bottom, top, m_ppp->near, m_ppp->far);
        std::cout << "\t-Projection mode is changed to orthographic" << std::endl;
        break;
    }
}

void OsgWxFrame::OnToggleRenderMode(wxCommandEvent& event) {

    switch (event.GetId()) {
    case wxID_MODES_RENDER_MODE_POINT:
        m_render_mode = osg::PolygonMode::POINT;
        std::cout << "\t-Polygon render mode is changed to Point" << std::endl;
        break;
    case wxID_MODES_RENDER_MODE_WIREFRAME:
        m_render_mode = osg::PolygonMode::LINE;
        std::cout << "\t-Polygon render mode is changed to Wireframe" << std::endl;
        break;
    case wxID_MODES_RENDER_MODE_FILL:
        m_render_mode = osg::PolygonMode::FILL;
        std::cout << "\t-Polygon render mode is changed to Fill" << std::endl;
        break;
    default:
        UsrLogErrorMessage("Id match error when setting polygon render mode");
        break;
    }

    usrSetPolygonMode(m_viewer->getSceneData());
}

void OsgWxFrame::OnToggleRenderFaceMode(wxCommandEvent& event) {

    switch (event.GetId()) {
    case wxID_MODES_RENDER_FACE_FRONT:
        m_render_face = osg::PolygonMode::FRONT;
        std::cout << "\t-Polygon face render mode is changed to Front" << std::endl;
        break;
    case wxID_MODES_RENDER_FACE_BACK :
        m_render_face = osg::PolygonMode::BACK;
        std::cout << "\t-Polygon face render mode is changed to Back" << std::endl;
        break;
    case wxID_MODES_RENDER_FACE_FRONT_AND_BACK:
        m_render_face = osg::PolygonMode::FRONT_AND_BACK;
        std::cout << "\t-Polygon face render mode is changed to Front and Back" << std::endl;
        break;
    default:
        UsrLogErrorMessage("Id match error when setting polygon face render mode");
    }

    usrSetPolygonMode(m_viewer->getSceneData());
}

void OsgWxFrame::OnToggleUIOperationMode(wxCommandEvent& event) {

    switch (event.GetId()) {
    case wxID_MODES_OPERATION_MODE_DISPLAY:
        UsrSetUIOperationMode(operation_mode::displaying);
        break;
    case wxID_MODES_OPERATION_MODE_MODELLING:
        UsrSetUIOperationMode(operation_mode::modelling);
        break;
    default:
        UsrLogErrorMessage("Id match error when setting user interface operation mode");
        break;
    }
}

void OsgWxFrame::OnToggleModellingConstraints(wxCommandEvent& event) {

    ImageModeller* modeller = m_canvas->UsrGetModeller();
    switch (event.GetId()) {
    case wxID_MODEL_CONSTRAINTS_NO_SPINE_CONSTRAINTS:
        modeller->sp_constraints = spine_constraints::none;
        std::cout << "\t-Spine constraint is set to none" << std::endl;
        break;
    case wxID_MODEL_CONSTRAINTS_PLANAR_SPINE_POINTS:
        modeller->sp_constraints = spine_constraints::planar;
        std::cout << "\t-Spine constraint is set to planar" << std::endl;
        break;
    case wxID_MODEL_CONSTRAINTS_STRAIGHT_PLANAR_SPINE:
        modeller->sp_constraints = spine_constraints::straight_planar;
        std::cout << "\t-Spine constraint is set to straight planar" << std::endl;
        break;
    case wxID_MODEL_CONSTRAINTS_NO_SECTION_CONSTRAINTS:
        modeller->sc_constraints = section_constraints::none;
        std::cout << "\t-Section constraint is set to none" << std::endl;
        break;
    case wxID_MODEl_CONSTRAINTS_CONSTANT_SECTIONS:
        modeller->sc_constraints = section_constraints::constant;
        std::cout << "\t-Section constraint is set to constant" << std::endl;
        break;
    case wxID_MODEL_CONSTRAINTS_LINEARLY_SCALED_SECTIONS:
        modeller->sc_constraints = section_constraints::linear_scaling;
        std::cout << "\t-Section constraint is set to linear scaling" << std::endl;
        break;
    default:
        UsrLogErrorMessage("No valid constraint!");
        break;
    }
}

void OsgWxFrame::OnToggleSpineDrawingMode(wxCommandEvent& event) {

    ImageModeller* modeller = m_canvas->UsrGetModeller();
    switch (event.GetId()) {
    case wxID_MODEL_SPINE_DRAWING_MODE_CONTINUOUS :
        modeller->spd_mode = spine_drawing_mode::continuous;
        break;
    case  wxID_MODEL_SPINE_DRAWING_MODE_PIECEWISE_LINEAR:
        modeller->spd_mode = spine_drawing_mode::piecewise_linear;
        break;
    default:
        UsrLogErrorMessage("Id match error when setting spine drawing mode");
        break;
    }
}

void OsgWxFrame::OnToggleImageDisplay(wxCommandEvent& event) {

    switch (event.GetId()) {
    case wxID_VIEW_DISPLAY_IMAGE:
        if(m_imgdisp_mode == background_image_display_mode::gradient_image) {
            usrChangeBackgroundImage(background_image_display_mode::image);
            m_imgdisp_mode = background_image_display_mode::image;
        }
        break;
    case wxID_VIEW_DISPLAY_GRADIENT_IMAGE:
        if(m_imgdisp_mode == background_image_display_mode::image) {
            usrChangeBackgroundImage(background_image_display_mode::gradient_image);
            m_imgdisp_mode = background_image_display_mode::gradient_image;
        }
        break;
    }
}

void OsgWxFrame::usrChangeBackgroundImage(background_image_display_mode mode) {

    if(mode == m_imgdisp_mode) return;

    osg::Image* image;
    if(mode == background_image_display_mode::image)
        image = osgDB::readImageFile(m_path.ToStdString());
    else if(mode == background_image_display_mode::gradient_image) {
         std::string grad_img_path = utilityInsertAfter(m_path, wxT('.'), wxT("_grad"));
         std::ifstream infile(grad_img_path);
         if(!infile.good()) {
             OtbFloatVectorImageType::Pointer img = LoadImage<OtbFloatVectorImageType>(m_path.ToStdString());
             OtbImageType::Pointer gimg = GradientMagnitudeImage(img);
             SaveImage<OtbImageType>(gimg, grad_img_path);
         }
         image = osgDB::readImageFile(grad_img_path);
    }

    if(!image) {
        std::cout << "Image file cannot be opened!" << std::endl;
        return;
    }

    osg::Texture2D* texture = new osg::Texture2D();
    texture->setResizeNonPowerOfTwoHint(false);
    texture->setImage(image);
    osg::StateSet* stateset = m_bgcam->getChild(0)->asGeode()->getOrCreateStateSet();
    stateset->setTextureAttributeAndModes(0, texture, osg::StateAttribute::ON);
}

void OsgWxFrame::OnDisplayLocalFrames(wxCommandEvent& event) {
    UsrLogErrorMessage("Not implemented yet");
}

void OsgWxFrame::OnDisplayWorldCoordinateFrame(wxCommandEvent& event) {

    if(!m_world_frame.valid()) {
        m_world_frame = new osg::Switch();
        osg::Node* xyz = osgDB::readNodeFile("../data/models/coordinate_frame.osg");
        if(!xyz) {
            UsrLogErrorMessage("Coordinate_frame model is not successfully loaded!");
        }
        m_world_frame->addChild(xyz);
        m_root->addChild(m_world_frame.get());
    }
    m_world_frame->setValue(0, GetMenuBar()->FindItem(event.GetId())->IsChecked());
    // m_viewer->home();
    m_root->dirtyBound();
}

void OsgWxFrame::OnDisplayVertexNormals(wxCommandEvent& event) {

    for(size_t i = 0; i < m_model->getNumChildren(); ++i) {
        ComponentBase* component = dynamic_cast<ComponentBase*>(m_model->getChild(i)->asGroup());
        if(component) component->DisplayVertexNormals(GetMenuBar()->FindItem(event.GetId())->IsChecked());
        else UsrLogErrorMessage("Dynamic cast error");
    }
    m_root->dirtyBound();
}

void OsgWxFrame::OnDisplaySectionNormals(wxCommandEvent& event) {

    for(size_t i = 0; i < m_model->getNumChildren(); ++i) {
        GeneralizedCylinder* gcyl = dynamic_cast<GeneralizedCylinder*>(m_model->getChild(i)->asGroup());
        if(gcyl) gcyl->DisplaySectionNormals(GetMenuBar()->FindItem(event.GetId())->IsChecked());
    }
    m_root->dirtyBound();
}

void OsgWxFrame::OnDisplayComponentRelationsDialog(wxCommandEvent& event) {

    m_component_relations_win->Show(!m_component_relations_win->IsShown());
}

void OsgWxFrame::OnSaveLastComponent(wxCommandEvent& event) {

    wxFileDialog dialog(this, wxT("Save the model"), wxEmptyString, wxEmptyString, wxT("*.osg"), wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
    if(dialog.ShowModal() == wxID_CANCEL) return;
    m_canvas->UsrSaveModel(dialog.GetPath());
}

void OsgWxFrame::OnSaveModel(wxCommandEvent& event) {

    wxFileDialog dialog(this, wxT("Save the model"), wxEmptyString, wxEmptyString, wxT("*.osg"), wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
    if(dialog.ShowModal() == wxID_CANCEL) return;
    osgDB::writeNodeFile(*m_model, dialog.GetPath().ToStdString());
}

void OsgWxFrame::OnDeleteModel(wxCommandEvent& event) {

    m_model->removeChildren(0, m_model->getNumChildren());
    m_canvas->UsrGetModeller()->DeleteModel();
    osg::Switch* selection_boxes = m_canvas->UsrGetSelectionBoxes();
    selection_boxes->removeChildren(0, selection_boxes->getNumChildren());
}

void OsgWxFrame::OnDeleteSelectedComponents(wxCommandEvent& event) {

    std::vector<unsigned int> _selections;
    UsrGetSelectedComponentIds(_selections);
    if(_selections.empty()) return;

    std::vector<int> selections;
    std::for_each(_selections.begin(), _selections.end(), [&selections](unsigned int i) {
        selections.push_back(static_cast<int>(i));
    });
    int id = 0;
    for(int i = 0; i < selections.size(); ++i) {
        for(int j = 0; j < m_model->getNumChildren(); ++j) {
            osg::Node* child = m_model->getChild(j);
            if(child->getUserValue("Component_Id", id)) {
                if(id == selections[i]) {
                    m_model->removeChild(j);
                    break;
                }
            }
            else {
                UsrLogErrorMessage("Component_Id could not be found!");
            }
        }
    }

    osg::Switch* selection_boxes = m_canvas->UsrGetSelectionBoxes();
    selection_boxes->removeChildren(0, selection_boxes->getNumChildren());
    m_canvas->UsrGetModeller()->DeleteSelectedComopnents(selections);
}

void OsgWxFrame::OnSolveModel(wxCommandEvent& event) {

    m_canvas->UsrGetModeller()->GetModelSolver()->Solve();
}

void OsgWxFrame::OnPrintProjectionMatrix(wxCommandEvent& event) {

    std::cout << "*********************************" << std::endl;
    print_camera_projection_matrix(m_viewer->getCamera());
    std::cout << "---------------------------------" << std::endl;
    print_camera_frustrum(m_viewer->getCamera());
    std::cout << "---------------------------------" << std::endl;
    print_camera_modelview_matrix(m_viewer->getCamera());
    std::cout << "---------------------------------" << std::endl;
    print_camera_viewport_mapping_matrix(m_viewer->getCamera());
    std::cout << "---------------------------------" << std::endl;
    print_camera_orientation(m_viewer->getCamera());
    std::cout << "---------------------------------" << std::endl;
    print_camera_calibration_matrix(m_viewer->getCamera());
    std::cout << "---------------------------------" << std::endl;
    print_3x4_camera_projection_matrix(m_viewer->getCamera());
    std::cout << "---------------------------------" << std::endl;
    osg::Matrixd mat_vpm;
    std::cout << "constructed viewport mapping matrix:" << std::endl;
    m_ppp->construct_viewport_mapping_matrix(mat_vpm);
    print_matrix(mat_vpm, 4, 4);
    std::cout << "---------------------------------" << std::endl;
    std::cout << "perspective projection parameters:" << std::endl;
    std::cout << *m_ppp << std::endl;
    std::cout << "*********************************" << std::endl;
}

void OsgWxFrame::UsrLogErrorMessage(const std::string& str) const {
    m_parent->UsrLogErrorMessage(str);
}
