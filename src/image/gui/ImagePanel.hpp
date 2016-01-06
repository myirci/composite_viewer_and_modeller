#ifndef _IMAGE_PANEL_HPP
#define _IMAGE_PANEL_HPP

#include "../../wx/WxUtility.hpp"

#include <wx/scrolwin.h>
#include <wx/image.h>
#include <wx/bitmap.h>
#include <vector>

class ImageFrame;
class DrawingLayer;

class ImagePanel : public wxScrolledWindow {
public:
    ImagePanel(ImageFrame* parent, wxString file = wxEmptyString);
    // Member functions:
    bool UsrLoadFile(const wxString& path);
    bool UsrSaveImage(const wxString& path);
    wxSize UsrGetImageSize() const;
    image_display_mode UsrGetMode() const;
    void UsrSetDisplayMode(image_display_mode mode);
    void UsrSetUIOperationMode(image_operation_mode mode);
    void UsrViewInOriginalSize();
    void UsrDisplayGraidentImage();
    void UsrRotate(rotation_type rt);
private:
    ImageFrame* m_parent;
    wxImage m_img;
    wxBitmap m_dimg;
    wxSize m_minImg;
    wxRect m_dimgRect;
    image_display_mode m_dpmode;
    image_operation_mode m_opmode;
    std::string m_path;

    // Member functions:
    void usrCalculateDisplayImageSize(double percentage);
    inline void usrUpdateDisplayImage();
    inline void usrUpdateImage();
    void usrOnOperationModeUpdated();
    // Event handlers:
    void PaintEvent(wxPaintEvent& event);
    void OnMouseMoved(wxMouseEvent& event);
    void OnLeftClick(wxMouseEvent& event);
    void OnRightClick(wxMouseEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnKeyPressed(wxKeyEvent& event);
    // some useful events
    /*
     void mouseWheelMoved(wxMouseEvent& event);
     void mouseReleased(wxMouseEvent& event);
     void mouseLeftWindow(wxMouseEvent& event);
     void keyReleased(wxKeyEvent& event);
     */
    DECLARE_EVENT_TABLE()
};

#endif
