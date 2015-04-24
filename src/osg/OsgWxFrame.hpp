#ifndef _OSG_WX_FRAME_HPP
#define _OSG_WX_FRAME_HPP

#include "OsgUtility.hpp"
#include <wx/frame.h>
#include <osgViewer/Viewer>
#include <osg/PolygonMode>
#include <memory>

class MainFrame;
class OsgWxGLCanvas;
class OsgWxGraphicsWindow;
class CoordinateTransformations;
class ComponentRelationsDialog;
class ModelSolver;

enum class operation_mode : unsigned char {
    displaying,
    modelling
};

class OsgWxFrame : public wxFrame {
private:

    MainFrame* m_parent;
    OsgWxGLCanvas* m_canvas;
    OsgWxGraphicsWindow* m_graphics_window;
    wxString m_path;
    osg::ref_ptr<osgGA::CameraManipulator> m_camera_manipulator;

    osg::ref_ptr<osgViewer::Viewer> m_viewer;       // viewer
    osg::ref_ptr<osg::Group> m_root;                // root of the scene graph
    osg::ref_ptr<osg::Group> m_model;               // parent node that keeps all the model
    osg::ref_ptr<osg::Camera> m_bgcam;              // to render background image
    osg::ref_ptr<osg::Geode> m_bgeode;              // to draw on the screen

    osg::ref_ptr<osg::Switch> m_world_frame;        // display world coordinate frame

    osg::PolygonMode::Mode m_render_mode;
    osg::PolygonMode::Face m_render_face;

    int m_id;
    operation_mode m_uiopmode;
    std::shared_ptr<CoordinateTransformations> m_ppp;
    std::unique_ptr<ComponentRelationsDialog> m_component_relations_win;

public:

    OsgWxFrame(wxWindow* parent, const wxPoint& pos, const wxSize& size, operation_mode md);
    void UsrSetFrameId(int id) { m_id = id; }
    bool UsrOpenModelFile();
    bool UsrOpenOrientedImageFile();
    bool UsrOpenImageFile();
    void UsrSetProjectionMatrix(double fovy, double aspect, double near, double far);
    void UsrAddSelectableComponent(osg::Node* node, unsigned int component_id);
    osg::Camera* UsrGetMainCamera();
    osg::Camera* UsrGetBackgroundCamera();
    osg::Geode* UsrGetBackgroundNode();
    operation_mode UsrGetUIOperationMode() const;
    void UsrSetUIOperationMode(operation_mode mode);
    void UsrDisplayComponentRelationsDialog();
    void UsrGetSelectedComponentIds(std::vector<unsigned int>& ids);
    ModelSolver* UsrGetModelSolver();
    void UsrUpdateGeosemanticConstraints();
    void UsrLogErrorMessage(const std::string& str) const;
private:

    // Member functions
    bool usrLoadModelFile(const wxString& fpath);
    bool usrLoadOrientationFile(const wxString& fpath);
    bool usrLoadImageFile(const wxString& fpath);
    void usrInitMenubar();
    void usrSetPolygonMode(osg::Node* node);
    void usrSetCustomPolygonMode(osg::Node* node, osg::PolygonMode::Mode mode, osg::PolygonMode::Face face);
    void usrUpdateFileTree(char type);
    void usrEnableModellingMenus(bool flag);

    // Event handlers
    void OnIdle(wxIdleEvent& event);
    void OnClose(wxCloseEvent& event);
    void OnOpenModel(wxCommandEvent& event);
    void OnOpenOrientedImage(wxCommandEvent& event);
    void OnExit(wxCommandEvent& event);
    void OnToggleRenderMode(wxCommandEvent& event);
    void OnToggleRenderFaceMode(wxCommandEvent& event);
    void OnToggleRenderType(wxCommandEvent& event);
    void OnToggleUIOperationMode(wxCommandEvent& event);
    void OnToggleModellingConstraints(wxCommandEvent& event);
    void OnToggleSpineDrawingMode(wxCommandEvent& event);
    void OnDisplayLocalFrames(wxCommandEvent& event);
    void OnDisplayWorldCoordinateFrame(wxCommandEvent& event);
    void OnDisplayVertexNormals(wxCommandEvent& event);
    void OnDisplaySectionNormals(wxCommandEvent& event);
    void OnDisplayComponentRelationsDialog(wxCommandEvent& event);
    void OnPrintProjectionMatrix(wxCommandEvent& event);
    void OnSaveLastComponent(wxCommandEvent& event);
    void OnSaveModel(wxCommandEvent& event);
    void OnDeleteModel(wxCommandEvent& event);
    void OnDeleteSelectedComponents(wxCommandEvent& event);
    void OnSolveModel(wxCommandEvent& event);
    DECLARE_EVENT_TABLE()
};

#endif // OsgWxFrame_HPP
