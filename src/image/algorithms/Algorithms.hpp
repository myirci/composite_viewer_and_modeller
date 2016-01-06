#ifndef _ALGORITHM_HPP
#define _ALGORITHM_HPP

#include <otbImage.h>
#include <otbVectorImage.h>
#include <otbImageFileReader.h>
#include <otbImageFileWriter.h>

class wxImage;
class wxPoint;
template <typename T> class Point2D;

// pixel type
typedef unsigned char PixelTypeUC;
typedef float         PixelTypeFL;

// image types
typedef otb::Image<PixelTypeUC, 2>         OtbImageType;
typedef otb::VectorImage<PixelTypeUC, 2>   OtbVectorImageType;
typedef otb::Image<PixelTypeFL,2>          OtbFloatImageType;
typedef otb::VectorImage<PixelTypeFL,2>    OtbFloatVectorImageType;

OtbImageType::Pointer RegionGrow(const wxImage& wxImg, const wxPoint& pt, int threshold);
bool RegionGrowSegmentation(const wxImage& wxImg, wxImage& segImg, const wxPoint& pt, int threshold);
bool BinaryImageRayCast(const OtbImageType::Pointer& image, const Point2D<int>& start, const Point2D<int>& end, Point2D<int>& first_hit);
OtbImageType::PixelType GradientImageRayCast(const OtbImageType::Pointer& image, const Point2D<int>& start, const Point2D<int>& end, Point2D<int>& hit);
void CopyToWxImageData(OtbImageType::Pointer image, unsigned char* data);
void GradientMagnitudeImage(OtbFloatVectorImageType::Pointer img, wxImage& gradImg);
OtbImageType::Pointer GradientMagnitudeImage(OtbFloatVectorImageType::Pointer image);

enum class PxlValues : PixelTypeUC {
    UNKNOWN = 127,
    SELECTED = 255,
    PROPOSED = 50,
    REFUSED = 0
};

struct IndexCompare {
    bool operator()(const OtbImageType::IndexType& lhs, const OtbImageType::IndexType& rhs) const {
        if(lhs[1] != rhs[1])
            return (lhs[1] < rhs[1]);
        return (lhs[0] < rhs[0]); }
};

template <typename ImType>
typename ImType::Pointer LoadImage(const std::string& fpath) {

    typedef otb::ImageFileReader<ImType>  ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fpath.c_str());
    reader->Update();
    return reader->GetOutput();
}

template <typename ImType>
void SaveImage(typename ImType::Pointer image, const std::string& fpath) {

    typedef otb::ImageFileWriter<ImType> ImageFileWriterType;
    typename ImageFileWriterType::Pointer writer = ImageFileWriterType::New();
    writer->SetInput(image);
    writer->SetFileName(fpath.c_str());
    writer->Update();
}


template <typename ImType, typename Functor>
void BresenhamLineAlgorithm(const typename ImType::Pointer image,
                            typename ImType::IndexType start, typename ImType::IndexType end, Functor func) {

    int delta_x = end[0] - start[0];
    int delta_y = end[1] - start[1];
    int dx_1(0), dy_1(0), dx_2(0), dy_2(0);
    if(delta_x < 0)         dx_1 = -1;
    else if (delta_x > 0)   dx_1 = 1;
    if(delta_y < 0)         dy_1 = -1;
    else if (delta_y > 0)   dy_1 = 1;
    dx_2 = dx_1;

    int longest = std::abs(delta_x);  // longest is the driving axis
    int shortest = std::abs(delta_y); // shortest is the passive axis
    if(shortest > longest) {
        std::swap(longest, shortest);
        dy_2 = dy_1;
        dx_2 = 0;
    }
    int numerator = longest >> 1;

    typename ImType::IndexType curr;
    curr[0] = start[0];
    curr[1] = start[1];

    for(int i = 0; i < longest; ++i) {
        // process the pixel with the index curr
        func(image, curr);

        // update the current segment index
        numerator += shortest;
        if(numerator >= longest) {
            numerator -= longest;
            curr[0] += dx_1;
            curr[1] += dy_1;
        } else {
            curr[0] += dx_2;
            curr[1] += dy_2;
        }
    }
}




#endif // ALGORITHM_HPP
