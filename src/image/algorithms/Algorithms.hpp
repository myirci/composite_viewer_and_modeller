#ifndef _ALGORITHM_HPP
#define _ALGORITHM_HPP

#include <otbImage.h>
#include <otbVectorImage.h>

class wxImage;
class wxPoint;
template <typename T> class Point2D;

// pixel type
typedef unsigned char PixelType;

// image types
typedef otb::Image<PixelType, 2>         OtbImageType;
typedef otb::VectorImage<PixelType, 2>   OtbVectorImageType;


// typedef otb::Image<float,2>         FloatImageType;
// typedef otb::VectorImage<float,2>   FloatVectorImageType;


enum class PxlValues : PixelType {
    UNKNOWN = 127,
    SELECTED = 255,
    PROPOSED = 50,
    REFUSED = 0
};

struct IndexCompare {
    bool operator()(const OtbImageType::IndexType& lhs, const OtbImageType::IndexType& rhs) const {
        if(lhs[1] != rhs[1]) { return (lhs[1] < rhs[1]); }
        return (lhs[0] < rhs[0]); }
};

OtbImageType::Pointer region_grow(const wxImage& wxImg, const wxPoint& pt, int threshold);
bool region_grow_segmentation(const wxImage& wxImg, const wxPoint& pt, int threshold);
bool region_grow_segmentation(const wxImage& wxImg, wxImage& segImg, const wxPoint& pt, int threshold);
bool ray_cast(const OtbImageType::Pointer& image, const OtbImageType::IndexType& start, const OtbImageType::IndexType& end, OtbImageType::IndexType& first_hit);
bool ray_cast(const OtbImageType::Pointer& image, const Point2D<int>& start, const Point2D<int>& end, Point2D<int>& first_hit);

#endif // ALGORITHM_HPP
