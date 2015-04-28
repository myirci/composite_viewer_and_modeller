#include "OsgWxGLCanvas.hpp"
#include "OsgWxFrame.hpp"
#include "OsgSelectionHandler.hpp"
#include "../wx/WxUtility.hpp"
#include "../modeller/ImageModeller.hpp"
#include "../modeller/CoordinateTransformations.hpp"
#include "../geometry/Primitives.hpp"
#include <wx/dcclient.h>
#include <wx/image.h>

#include <osg/MatrixTransform>
#include <osg/ShapeDrawable>

BEGIN_EVENT_TABLE(OsgWxGLCanvas, wxGLCanvas)
EVT_SIZE(OsgWxGLCanvas::OnSize)
EVT_PAINT(OsgWxGLCanvas::OnPaint)
EVT_ERASE_BACKGROUND(OsgWxGLCanvas::OnEraseBackground)
EVT_CHAR(OsgWxGLCanvas::OnChar)
EVT_KEY_UP(OsgWxGLCanvas::OnKeyUp)
EVT_ENTER_WINDOW(OsgWxGLCanvas::OnMouseEnter)
EVT_LEFT_DOWN(OsgWxGLCanvas::OnMouseDown)
EVT_MIDDLE_DOWN(OsgWxGLCanvas::OnMouseDown)
EVT_RIGHT_DOWN(OsgWxGLCanvas::OnMouseDown)
EVT_LEFT_UP(OsgWxGLCanvas::OnMouseUp)
EVT_MIDDLE_UP(OsgWxGLCanvas::OnMouseUp)
EVT_RIGHT_UP(OsgWxGLCanvas::OnMouseUp)
EVT_MOTION(OsgWxGLCanvas::OnMouseMotion)
EVT_MOUSEWHEEL(OsgWxGLCanvas::OnMouseWheel)
END_EVENT_TABLE()

OsgWxGLCanvas::OsgWxGLCanvas(wxWindow *parent, wxWindowID id, int *attributes, const wxPoint& pos, const wxSize& size, long style, const wxString& name) :
    wxGLCanvas(parent, id, attributes, pos, size, style|wxFULL_REPAINT_ON_RESIZE, name),
    m_selection_handler(new OsgSelectionHandler()),
    m_modeller(nullptr) {

    // initialize the graphics context
    m_context = new wxGLContext(this);

    // set default cursor to standard
    m_oldCursor = *wxSTANDARD_CURSOR;

    m_parent = dynamic_cast<OsgWxFrame*>(parent);
    if(!m_parent) UsrLogErrorMessage("Dynamic cast_error");
}

OsgWxGLCanvas::~OsgWxGLCanvas() {

    if(m_context  != nullptr) {
        delete m_context;
        m_context = nullptr;
    }
    if(m_modeller != nullptr) {
        delete m_modeller;
        m_modeller = nullptr;
    }
}

void OsgWxGLCanvas::UsrMakeContextCurrent() {
    SetCurrent(*m_context); // Set the graphics context current
}

void OsgWxGLCanvas::UsrUseCursor(bool value) {

    if (value) SetCursor(m_oldCursor); // show the old cursor
    else {
        m_oldCursor = GetCursor(); // remember the old cursor
        SetCursor(wxStockCursor(wxCURSOR_BLANK)); // hide the cursor
    }
}

void OsgWxGLCanvas::UsrInitializeModeller(const std::shared_ptr<CoordinateTransformations>& ppp, const wxString& fpath) {

    if(m_modeller != nullptr) {
        delete m_modeller;
        m_modeller = nullptr;
    }
    m_modeller = new ImageModeller(fpath, ppp, this);
    m_modeller->Initialize2DDrawingInterface(m_parent->UsrGetBackgroundNode());
}

ImageModeller* OsgWxGLCanvas::UsrGetModeller() {
    return m_modeller;
}

osg::Switch* OsgWxGLCanvas::UsrGetSelectionBoxes() {
    return m_selection_handler->GetOrCreateSelectionBoxSwitch();
}

const osg::Camera* const OsgWxGLCanvas::UsrGetMainCamera() const {
    return m_parent->UsrGetMainCamera();
}

void OsgWxGLCanvas::UsrAddSelectableNodeToDisplay(osg::Node* node, unsigned int component_id) {
    m_parent->UsrAddSelectableComponent(node, component_id);
}

void OsgWxGLCanvas::UsrAddToBackgroundDisplay(osg::Geometry* geom) {
    m_parent->UsrGetBackgroundNode()->addDrawable(geom);
}

wxPoint OsgWxGLCanvas::usrDeviceToLogical(const wxPoint& p) const {

    wxSize size = m_parent->GetClientSize();
    return wxPoint(p.x, size.GetHeight() - p.y - 1);
}

void OsgWxGLCanvas::UsrLogicalToDevice(wxPoint& p) const {

    wxSize size = m_parent->GetClientSize();
    p.y = size.GetHeight() - p.y - 1;
}

void OsgWxGLCanvas::UsrTransformCoordinates(Point2D<int>& pt) const {

    // DeviceToLogical and LogicalToDevice tranformation is same in both ways
    wxSize size = m_parent->GetClientSize();
    pt.y = size.GetHeight() - pt.y - 1;
}

void OsgWxGLCanvas::UsrSaveModel(const wxString& path) const {
    m_modeller->SaveModel(path.ToStdString());
}

void OsgWxGLCanvas::OnPaint(wxPaintEvent& event) {
    wxPaintDC dc(this);
}

void OsgWxGLCanvas::OnSize(wxSizeEvent& event) {
    //    if (m_graphics_window.valid()) {
    //        wxSize size = m_parent->GetClientSize();
    //        m_graphics_window->getEventQueue()->windowResize(
    //                    0, 0, size.GetWidth(), size.GetHeight());
    //        m_graphics_window->resized(0, 0, size.GetWidth(), size.GetHeight());
    //    }
}

void OsgWxGLCanvas::OnEraseBackground(wxEraseEvent& event) {
    // Do nothing, to avoid flashing on MSW
}

void OsgWxGLCanvas::OnChar(wxKeyEvent &event) {

#if wxUSE_UNICODE
    int key = event.GetUnicodeKey();
#else
    int key = event.GetKeyCode();
#endif

    if(m_parent->UsrGetUIOperationMode() == operation_mode::modelling)
        if(key == WXK_ESCAPE)
            m_modeller->Reset2DDrawingInterface();

    if(key == 't' || key == 'T') {
        if(m_parent->UsrGetUIOperationMode() == operation_mode::modelling)
            m_parent->UsrSetUIOperationMode(operation_mode::displaying);
        else if(m_parent->UsrGetUIOperationMode() == operation_mode::displaying)
            m_parent->UsrSetUIOperationMode(operation_mode::modelling);
    }

    if(key == 'p' || key == 'P')
        m_modeller->DebugPrint();

    if(m_graphics_window.valid())
        m_graphics_window->getEventQueue()->keyPress(key);

    // If this key event is not processed here, we should call
    // event.Skip() to allow processing to continue.
}

void OsgWxGLCanvas::OnKeyUp(wxKeyEvent& event) {
#if wxUSE_UNICODE
    int key = event.GetUnicodeKey();
#else
    int key = event.GetKeyCode();
#endif
    if (m_graphics_window.valid())
        m_graphics_window->getEventQueue()->keyRelease(key);

    // If this key event is not processed here, we should call
    // event.Skip() to allow processing to continue.
}

void OsgWxGLCanvas::OnMouseEnter(wxMouseEvent& event) {
    // Set focus to ourselves, so keyboard events get directed to us
    SetFocus();
}

void OsgWxGLCanvas::OnMouseDown(wxMouseEvent &event) {

    operation_mode uiop_mode = m_parent->UsrGetUIOperationMode();
    if(uiop_mode == operation_mode::modelling) {

        if(m_modeller == nullptr) return;

        // clicked point
        wxPoint pt = usrDeviceToLogical(event.GetPosition());

        // Left click
        if(event.GetButton() == 1)
            m_modeller->OnLeftClick(static_cast<double>(pt.x), static_cast<double>(pt.y));

        // Right click
        if(event.GetButton() == 3)
            m_modeller->OnRightClick(static_cast<double>(pt.x), static_cast<double>(pt.y));
    }
    else {
        if(m_graphics_window.valid())
            m_graphics_window->getEventQueue()->mouseButtonPress(event.GetX(), event.GetY(), event.GetButton());
    }
}

void OsgWxGLCanvas::OnMouseUp(wxMouseEvent& event) {

    if(event.ControlDown()) {
        wxPoint pt = usrDeviceToLogical(event.GetPosition());
        m_selection_handler->HandleSelection(m_parent->UsrGetMainCamera(), pt.x, pt.y);
        m_parent->UsrUpdateGeosemanticConstraints();
    }

    if(m_graphics_window.valid()) {
        m_graphics_window->getEventQueue()->mouseButtonRelease(event.GetX(), event.GetY(), event.GetButton());
    }
}

void OsgWxGLCanvas::OnMouseMotion(wxMouseEvent& event) {

    if(m_parent->UsrGetUIOperationMode() == operation_mode::modelling) {

        if(m_modeller == nullptr) return;
        wxPoint pt = usrDeviceToLogical(event.GetPosition());
        m_modeller->OnMouseMove(static_cast<double>(pt.x), static_cast<double>(pt.y));
    }
    else if(m_parent->UsrGetUIOperationMode() == operation_mode::displaying) {
        if (m_graphics_window.valid()) {
            m_graphics_window->getEventQueue()->mouseMotion(event.GetX(), event.GetY());
        }
    }

    m_parent->SetStatusText(utilityToString(usrDeviceToLogical(event.GetPosition())), 0);
}

void OsgWxGLCanvas::OnMouseWheel(wxMouseEvent& event) {

    int delta = event.GetWheelRotation() / event.GetWheelDelta() * event.GetLinesPerAction();
    if (m_graphics_window.valid())
        m_graphics_window->getEventQueue()->mouseScroll(delta > 0 ? osgGA::GUIEventAdapter::SCROLL_UP : osgGA::GUIEventAdapter::SCROLL_DOWN);

}

void OsgWxGLCanvas::UsrLogErrorMessage(const std::string& str) const {
    m_parent->UsrLogErrorMessage(str);
}

void OsgWxGLCanvas::UsrSetRenderingType(rendering_type rtype) {
    m_modeller->SetRenderingType(rtype);
}
