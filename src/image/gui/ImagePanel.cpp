#include "ImagePanel.hpp"
#include "ImageFrame.hpp"
#include "../algorithms/Algorithms.hpp"
#include <wx/dcclient.h>
#include <wx/msgdlg.h>
#include <wx/string.h>
#include <iostream>
#include <cmath>

BEGIN_EVENT_TABLE(ImagePanel, wxPanel)
EVT_PAINT(ImagePanel::PaintEvent)
EVT_SIZE(ImagePanel::OnSize)
EVT_MOTION(ImagePanel::OnMouseMoved)
EVT_LEFT_DOWN(ImagePanel::OnLeftClick)
EVT_RIGHT_DOWN(ImagePanel::OnRightClick)
EVT_KEY_DOWN(ImagePanel::OnKeyPressed)

// some useful events
/*
 EVT_LEFT_UP(ImagePanel::mouseReleased)
 EVT_LEAVE_WINDOW(ImagePanel::mouseLeftWindow)
 EVT_KEY_UP(ImagePanel::keyReleased)
 EVT_MOUSEWHEEL(ImagePanel::mouseWheelMoved)
 */
END_EVENT_TABLE()

ImagePanel::ImagePanel(ImageFrame* parent, wxString file_path) : wxScrolledWindow(parent),
    m_minImg(wxSize(0,0)), m_dimgRect(), m_dpmode(image_display_mode(image_display_mode::VARYING)),
    m_opmode(image_operation_mode::Default), m_path("") {
    SetBackgroundColour(*wxWHITE);
    SetScrollRate(1, 1);
    if(file_path != wxEmptyString) {
        m_path = file_path.ToStdString();
        UsrLoadFile(file_path);
    }
    else { /* may load default image */ }
    m_parent = parent;
}

// Public Member Functions
bool ImagePanel::UsrLoadFile(const wxString& path) {

    m_img.Destroy();
    if(m_img.LoadFile(path)) {
        m_path = path.ToStdString();
        usrUpdateImage();
        return true;
    }
    return false;
}
void ImagePanel::UsrSetDisplayMode(image_display_mode mode) {

    m_dpmode = mode;
    wxSizeEvent event;
    OnSize(event);
}
void ImagePanel::UsrSetUIOperationMode(image_operation_mode mode) {

    m_opmode = mode;
}

image_display_mode ImagePanel::UsrGetMode() const {

    return m_dpmode;
}

wxSize ImagePanel::UsrGetImageSize() const {

    return wxSize(m_img.GetWidth(), m_img.GetHeight());
}

void ImagePanel::UsrDisplayGraidentImage() {

    if(!m_img.IsOk()) return;

    if(m_path.empty()) {
        std::cerr << "Image file path is empty" << std::endl;
        return;
    }

    wxImage gradImg(m_dimg.ConvertToImage());
    OtbFloatVectorImageType::Pointer img = LoadImage<OtbFloatVectorImageType>(m_path);
    GradientMagnitudeImage(img, gradImg);
    m_dimg = wxBitmap(gradImg);
    Refresh();
}

void ImagePanel::UsrViewInOriginalSize() {

    if(!m_img.IsOk()) return;
    SetVirtualSize(m_img.GetWidth(), m_img.GetHeight());
    m_dimgRect.SetPosition(wxPoint(0,0));
    m_dimgRect.SetSize(wxSize(m_img.GetWidth(), m_img.GetHeight()));
    m_parent->SetStatusText(utilityToString("disp img size: ", m_dimgRect.GetSize()), 1);
    usrUpdateDisplayImage();
    Refresh();
}

void ImagePanel::UsrRotate(rotation_type rt) {

    if(!m_img.IsOk()) return ;
    if(rt == rotation_type::CW)       m_img = m_img.Rotate90(true);
    else if(rt == rotation_type::CCW) m_img = m_img.Rotate90(false);
    usrUpdateImage();
}

bool ImagePanel::UsrSaveImage(const wxString& path) {

    wxImage img = m_dimg.ConvertToImage();
    return img.SaveFile(path);
}

// Private Member Functions
void ImagePanel::usrCalculateDisplayImageSize(double percentage) {

    int limit(0), k(0);
    wxSize screen_size = wxGetDisplaySize();
    if(screen_size.GetHeight() < screen_size.GetWidth()) {
        limit = floor(percentage * static_cast<double>(screen_size.GetHeight()));
        k = limit / m_minImg.GetHeight();
        if(k > (m_img.GetHeight() / m_minImg.GetHeight())) {
            k = m_img.GetHeight() / m_minImg.GetHeight();
            m_dimgRect.SetSize(m_minImg*k);
            return;
        }
    }
    else {
        limit = floor(percentage * static_cast<double>(screen_size.GetWidth()));
        k = limit / m_minImg.GetWidth();
        if(k > (m_img.GetWidth() / m_minImg.GetWidth())) {
            k = m_img.GetWidth() / m_minImg.GetWidth();
            m_dimgRect.SetSize(m_minImg*k);
            return;
        }
    }
    m_dimgRect.SetSize(m_img.GetSize());
}

void ImagePanel::usrUpdateDisplayImage() {

    if(!m_img.IsOk()) return;
    if((m_img.GetWidth() == m_dimgRect.GetWidth()) && (m_img.GetHeight() == m_dimgRect.GetHeight()))
        m_dimg = wxBitmap(m_img);
    else
        m_dimg = wxBitmap(m_img.Scale(m_dimgRect.GetWidth(), m_dimgRect.GetHeight(), wxIMAGE_QUALITY_HIGH));
}

void ImagePanel::usrUpdateImage() {

    m_minImg = utilitySimplify(wxSize(m_img.GetWidth(), m_img.GetHeight()));
    usrCalculateDisplayImageSize(0.6);
    m_parent->SetClientSize(m_dimgRect.GetSize());
    SetVirtualSize(m_dimgRect.GetSize());
    m_parent->SetStatusText(utilityToString("img size: ", UsrGetImageSize()), 2);
    usrUpdateDisplayImage();
    Refresh();
}

void ImagePanel::usrOnOperationModeUpdated() {

    if(m_opmode == image_operation_mode::Drawing) {

    }
}

// Event Handlers
void ImagePanel::PaintEvent(wxPaintEvent& event) {
    wxPaintDC dc(this);  // may use double-buffered dcs
    DoPrepareDC(dc);     // to prepare the device context for drawing a scrolled image
    dc.DrawBitmap(m_dimg, m_dimgRect.GetPosition());
}

void ImagePanel::OnSize(wxSizeEvent& event) {

    if(!m_img.IsOk()) return;
    wxSize diff = m_parent->GetClientSize() - m_dimgRect.GetSize();
    if(m_dpmode == image_display_mode::FIXED) {
        if(diff.GetWidth() > 0)  m_dimgRect.x = diff.GetWidth() / 2;
        else                     m_dimgRect.x = 0;
        if(diff.GetHeight() > 0) m_dimgRect.y = diff.GetHeight() / 2;
        else                     m_dimgRect.y = 0;
    }
    else if(m_dpmode == image_display_mode::VARYING) {
        int kw = diff.GetWidth() / m_minImg.GetWidth();
        int kh = diff.GetHeight() / m_minImg.GetHeight();
        int k = std::min(kw, kh);
        m_dimgRect.SetSize(m_dimgRect.GetSize() + m_minImg*k);
        SetVirtualSize(m_parent->GetClientSize());
        diff = m_parent->GetClientSize() - m_dimgRect.GetSize();
        if(diff.GetWidth() > 0)  m_dimgRect.x = diff.GetWidth() / 2;
        else                     m_dimgRect.x = 0;
        if(diff.GetHeight() > 0) m_dimgRect.y = diff.GetHeight() / 2;
        else                     m_dimgRect.y = 0;
        usrUpdateDisplayImage();
    }
    m_parent->SetStatusText(utilityToString("disp img size: ", m_dimgRect.GetSize()), 1);
    Refresh();
    event.Skip();
}

void ImagePanel::OnMouseMoved(wxMouseEvent& event) {

    if(!m_img.IsOk()) return;
    wxPoint mouse_pos = event.GetPosition();
    if(m_dimgRect.Contains(mouse_pos))
        m_parent->SetStatusText(utilityToString(mouse_pos - m_dimgRect.GetPosition() + wxPoint(GetScrollPos(wxHORIZONTAL), GetScrollPos(wxVERTICAL))), 0);
}

void ImagePanel::OnLeftClick(wxMouseEvent& event) {

    if(!m_img.IsOk()) return;
    wxPoint mouse_pos = event.GetPosition();
    if(m_opmode == image_operation_mode::RegionGrowing) {
        if(m_dimgRect.Contains(mouse_pos)) {
            wxPoint imgCoord = mouse_pos - m_dimgRect.GetPosition() + wxPoint(GetScrollPos(wxHORIZONTAL), GetScrollPos(wxVERTICAL));
            wxImage dummy(m_dimg.ConvertToImage());
            unsigned char* data = dummy.GetData();
            unsigned char* dcopy = new unsigned char[dummy.GetWidth()*dummy.GetHeight()*3];
            std::copy(data, data + dummy.GetWidth()*dummy.GetHeight()*3, dcopy);
            RegionGrowSegmentation(m_dimg.ConvertToImage(), dummy, imgCoord, 9);
            m_dimg = wxBitmap(dummy);
            Refresh();
        }
    }
}

void ImagePanel::OnRightClick(wxMouseEvent& event) {

    if(!m_img.IsOk()) return;
}

void ImagePanel::OnKeyPressed(wxKeyEvent& event) {

    if(!m_img.IsOk()) return;
}
