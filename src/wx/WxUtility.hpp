#ifndef _WX_UTILITY_HPP
#define _WX_UTILITY_HPP

#include <string>
#include <wx/unichar.h>

enum class message_type : unsigned char {
    INFO,
    WARNING,
    ERROR
};

enum class rotation_type : unsigned char {
    CW,
    CCW
};

enum class image_display_mode : unsigned char {
    FIXED,
    VARYING
};

enum class image_operation_mode : unsigned char {
    Default,
    RegionGrowing,
    Drawing
};

enum class log_type : unsigned char {
    WARNING,
    ERROR,
    DEFAULT
};

class wxString;
class wxSize;
class wxRect;
class wxPoint;

void utilityShowMessageDialog(message_type mt, const wxString& message);
bool utilityQuestionDialoag(const wxString& message);
int utilityGreatestCommonDivisor(int a, int b);
wxSize utilitySimplify(const wxSize& fraction);
void utilityPrintSize(const wxSize& size);
void utilityPrintSize(const std::string& str, const wxSize& size);
void utilityPrintRect(const wxRect& rect);
void utilityPrintPoint(const wxPoint& p);
wxString utilityToString(const wxPoint& p);
wxString utilityToString(const wxSize& s);
wxString utilityToString(const std::string& str, const wxSize& s);
std::string utilityInsertAfter(const wxString& fpath,
                               const wxUniChar& after,
                               const wxString& str);
void utilityPrintDisplayMode(image_display_mode md);
bool utilityIsXml(const wxString& fname);

#endif // WXUTILITY_HPP
