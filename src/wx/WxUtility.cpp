#include "WxUtility.hpp"
#include <wx/msgdlg.h>
#include <cmath>
#include <iostream>

void utilityShowMessageDialog(message_type mt, const wxString& message) {

    wxString caption;
    long style = wxOK;
    if(mt == message_type::INFO) {
        caption = wxT("Info");
    }
    else if( mt == message_type::WARNING) {
        caption = wxT("Warning");
        style |= wxICON_EXCLAMATION;
    }
    else if(mt == message_type::ERROR) {
        caption = wxT("Error");
        style |= wxICON_ERROR;
    }
    else {
        return;
    }
    wxMessageDialog* dlg = new wxMessageDialog(NULL, message, caption, style);
    dlg->ShowModal();
    dlg->Destroy();
}

bool utilityQuestionDialoag(const wxString& message) {

    return  wxMessageBox(message, wxT("Confirm"), wxICON_QUESTION | wxYES_NO) == wxYES;
}

// Euclid's greatest common divisior algorithm
int utilityGreatestCommonDivisor(int a, int b) {

    int max(a);
    int min(b);
    if(max < min) { std::swap(max,min); }
    int temp(0);
    while(min != 0) {
        temp = min;
        min = max % min;
        max = temp;
    }
    return max;
}

wxSize utilitySimplify(const wxSize& size) {

    int gcd = utilityGreatestCommonDivisor(size.GetWidth(), size.GetHeight());
    return wxSize(size.GetWidth()/gcd, size.GetHeight()/gcd);
}

void utilityPrintSize(const wxSize& size) {

    std::cout << "w: " << size.GetWidth() << " h: " << size.GetHeight() << std::endl;
}

void utilityPrintSize(const std::string& str, const wxSize& size) {

    std::cout << str << " w: " << size.GetWidth() << " h: " << size.GetHeight() << std::endl;
}

void utilityPrintPoint(const wxPoint& p) {

    std::cout << "x: " << p.x << " y: " << p.y << std::endl;
}

wxString utilityToString(const wxPoint& p) {

    std::string str = "("+ std::to_string(p.x) + ", " + std::to_string(p.y) + ")";
    return wxString(str.c_str(), wxConvUTF8);
}

wxString utilityToString(const wxSize& s) {

    std::string str = "("+ std::to_string(s.GetWidth()) + ", " + std::to_string(s.GetHeight()) + ")";
    return wxString(str.c_str(), wxConvUTF8);
}

wxString utilityToString(const std::string& str, const wxSize& s) {

    std::string st = str + "("+ std::to_string(s.GetWidth()) + ", " + std::to_string(s.GetHeight()) + ")";
    return wxString(st.c_str(), wxConvUTF8);
}

std::string utilityInsertAfter(const wxString& fpath, const wxUniChar& after, const wxString& str) {

    wxString rest;
    wxString before = fpath.BeforeLast(after, &rest);
    before = before + str + after + rest;
    return before.ToStdString();
}

void utilityPrintRect(const wxRect& rect) {

    std::cout << "position: " << "x: " << rect.x << " y: " << rect.y << std::endl;
    std::cout << "size: " << "w: " << rect.width << " h: " << rect.height << std::endl;
}

void utilityPrintDisplayMode(image_display_mode md) {

    if(md == image_display_mode::FIXED)        std::cout << "FIXED" << std::endl;
    else if(md == image_display_mode::VARYING) std::cout << "VARYING" << std::endl;

}

bool utilityIsXml(const wxString& fname) {

    wxString extension = fname.AfterLast('.');
    return (extension == wxT("xml") || extension == wxT("XML"));
}
