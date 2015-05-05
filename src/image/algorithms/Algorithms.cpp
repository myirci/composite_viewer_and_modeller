
#include <list>
#include <map>
#include <cmath>

#include "Algorithms.hpp"
#include "../../geometry/Primitives.hpp"

#include <wx/image.h>

#include <otbImageFileWriter.h>
#include <otbImportVectorImageFilter.h>
#include <otbVectorImageToImageListFilter.h>
#include <otbVectorImageToIntensityImageFilter.h>
#include <otbImageList.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkNeighborhoodIterator.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>

typedef itk::ImageRegionIterator<OtbImageType>      ItkImgIt_IteratorType;
typedef itk::ImageRegionConstIterator<OtbImageType> ItkImgIt_ConstIteratorType;
typedef itk::NeighborhoodIterator<OtbImageType>     ItkImgIt_NeighborhoodIteratorType;
typedef otb::ImageFileWriter<OtbImageType>          OtbImageWriterType;

void CopyToWxImageData(OtbImageType::Pointer image, unsigned char* data) {

    OtbImageType::IndexType index;
    OtbImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
    typedef itk::ImageRegionIterator<OtbImageType> ImageRegionIteratorType;
    ImageRegionIteratorType it(image, image->GetLargestPossibleRegion());
    for(it.GoToBegin(); !it.IsAtEnd(); ++it) {
        index = it.GetIndex();
        long pos = (index[1]* size[0] + index[0])*3;
        data[pos] = data[pos+1] = data[pos+2] = it.Get();
    }
}

bool RegionGrowSegmentation(const wxImage& wxImg, wxImage& segImg, const wxPoint& pt, int threshold) {

    // execute the region grow algorithm
    OtbImageType::Pointer img = RegionGrow(wxImg, pt, threshold);

    // convert the output to wxImg
    int numPixels = 3*wxImg.GetWidth()*wxImg.GetHeight();
    PixelTypeUC* pxdt = new PixelTypeUC[numPixels];
    CopyToWxImageData(img, pxdt);
    segImg = wxImage(wxImg.GetWidth(), wxImg.GetHeight(), pxdt, true);

    return true;
}

OtbImageType::Pointer RegionGrow(const wxImage& wxImg, const wxPoint& pt, int threshold) {

    // Step-1: Convert wxImage to OtbVectorImage
    typedef otb::ImportVectorImageFilter<OtbVectorImageType> ImporterType;
    ImporterType::Pointer importFilter = ImporterType::New();
    ImporterType::SizeType size;
    size[0] = wxImg.GetWidth();
    size[1] = wxImg.GetHeight();
    ImporterType::IndexType start;
    start.Fill(0);
    ImporterType::RegionType region;
    region.SetIndex(start);
    region.SetSize(size);
    importFilter->SetRegion(region);
    double origin[2] = {0, 0};
    importFilter->SetOrigin(origin);
    double spacing[2] = {1.0, 1.0};
    importFilter->SetSpacing(spacing);
    int numPixels(3*size[0]*size[1]);
    PixelTypeUC* dt = wxImg.GetData();
    PixelTypeUC* data = new PixelTypeUC[numPixels];
    std::copy(dt, dt + numPixels, data);
    importFilter->SetImportPointer(data, numPixels, true);
    OtbVectorImageType::Pointer vecImg = OtbVectorImageType::New();
    vecImg = dynamic_cast<OtbVectorImageType*>(importFilter->GetOutput());  

    // Step-2: Convert OtbVectorImage to OtbImage
    OtbImageType::Pointer intensity_image = OtbImageType::New();
    typedef otb::VectorImageToIntensityImageFilter<OtbVectorImageType, OtbImageType> IntensityFilter;
    IntensityFilter::Pointer filter = IntensityFilter::New();
    filter->SetInput(vecImg);
    filter->Update();
    intensity_image = filter->GetOutput();

    // Step-3: Create and initialize the labels image
    OtbImageType::Pointer labels_image = OtbImageType::New();
    labels_image->SetRegions(intensity_image->GetLargestPossibleRegion());
    labels_image->Allocate();
    ItkImgIt_IteratorType it(labels_image, labels_image->GetLargestPossibleRegion());
    for(it.GoToBegin(); !it.IsAtEnd(); ++it) {
        it.Set(static_cast<PixelTypeUC>(PxlValues::UNKNOWN));
    }

    // Step-4: Define the neighborhood iterators for the intensity and labels images
    ItkImgIt_NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);
    ItkImgIt_NeighborhoodIteratorType it1(radius, intensity_image, intensity_image->GetLargestPossibleRegion());
    ItkImgIt_NeighborhoodIteratorType it2(radius, labels_image, labels_image->GetLargestPossibleRegion());
    it1.NeedToUseBoundaryConditionOff();
    it2.NeedToUseBoundaryConditionOff();

    // Step-5: Define the containers and variables required for the analysis
    std::list<OtbImageType::IndexType> process;
    std::set<OtbImageType::IndexType, IndexCompare> processed;

    // Step-6: Push the seed pixel into the process list and m_processed set
    OtbImageType::IndexType seed;
    seed[0] = pt.x;
    seed[1] = pt.y;
    process.push_back(seed);
    processed.insert(seed);

    // Step-7: Execute the region growing process
    OtbImageType::IndexType center_index;
    it1.SetLocation(seed);
    int seed_val = static_cast<int>(it1.GetCenterPixel());
    int w = wxImg.GetWidth();
    int h = wxImg.GetHeight();
    while(!process.empty()) {
        // Update center_index
        center_index[0] = (process.front())[0];
        center_index[1] = (process.front())[1];

        // Update iterator locations by center_index
        it1.SetLocation(center_index);
        it2.SetLocation(center_index);

        // Update the value of the pixel (indicated by the center_index) with "SELECTED".
        it2.SetCenterPixel(static_cast<PixelTypeUC>(PxlValues::SELECTED));
        OtbImageType::IndexType idx = it2.GetIndex();
        if(idx[0] == 0 || idx[1] == 0 || idx[0] == w-1 || idx[1] == h-1) {
            process.pop_front();
            continue;
        }

        // Iterate over the neighbors of the center pixel and decide SELECTED or REFUSED
        for(unsigned int i = 0; i < it1.Size(); ++i) {
            if(i == 4) { continue; }
            if(processed.insert(it1.GetIndex(i)).second) {
                int pix_val = static_cast<int>(it1.GetPixel(i));
                if(std::abs( pix_val - seed_val) <= threshold) {
                    process.push_back(it1.GetIndex(i));
                }
                else {
                    it2.SetPixel(i, static_cast<PixelTypeUC>(PxlValues::REFUSED));
                }
            }
        }
        process.pop_front();
    }
    return labels_image;
}

// image is a binary image: pixel_values are either 255 or 0
// returns true if the ray hits an occupied pixel (value = 0)
// start: starting point for the ray segment
// end  : end point for the ray segment
bool BinaryImageRayCast(const OtbImageType::Pointer& image, const Point2D<int>& start, const Point2D<int>& end, Point2D<int>& first_hit) {

    int delta_x = end.x - start.x;
    int delta_y = end.y - start.y;
    int dx_1(0), dy_1(0), dx_2(0), dy_2(0);

    if(delta_x < 0)         dx_1 = -1;
    else if (delta_x > 0)   dx_1 = 1;
    if(delta_y < 0)         dy_1 = -1;
    else if (delta_y > 0)   dy_1 = 1;
    dx_2 = dx_1;

    // longest is the driving axis
    // shortest is the passive axis
    int longest = std::abs(delta_x);
    int shortest = std::abs(delta_y);

    if(shortest > longest) {
        std::swap(longest, shortest);
        dy_2 = dy_1;
        dx_2 = 0;
    }
    int numerator = longest >> 1;
    OtbImageType::IndexType pixelIndex;
    pixelIndex[0] = start.x;
    pixelIndex[1] = start.y;

    for(int i = 0; i < longest; i++) {
        if(image->GetPixel(pixelIndex) == 0) {
            first_hit.x = pixelIndex[0];
            first_hit.y = pixelIndex[1];
            return true;
        }
        numerator += shortest;
        if(numerator >= longest) {
            // increment/decrement passive and driving axis by 1
            numerator -= longest;
            pixelIndex[0] += dx_1;
            pixelIndex[1] += dy_1;
        } else {
            // increment/decrement only driving axis
            pixelIndex[0] += dx_2;
            pixelIndex[1] += dy_2;
        }

        if(pixelIndex[0] < 0 || pixelIndex[1] < 0) {
            std::cerr << "ERROR: NEGATIVE INDEX!" << std::endl;
            break;
        }
    }
    return false;
}

// input    : gradient image
// return   : maximim pixel value and corresponding pixel
// start    : starting point for the ray segment
// end      : end point for the ray segment
OtbImageType::PixelType GradientImageRayCast(const OtbImageType::Pointer& image, const Point2D<int>& start, const Point2D<int>& end, Point2D<int>& hit) {

    int delta_x = end.x - start.x;
    int delta_y = end.y - start.y;
    int dx_1(0), dy_1(0), dx_2(0), dy_2(0);

    if(delta_x < 0)         dx_1 = -1;
    else if (delta_x > 0)   dx_1 = 1;
    if(delta_y < 0)         dy_1 = -1;
    else if (delta_y > 0)   dy_1 = 1;
    dx_2 = dx_1;

    // longest is the driving axis
    // shortest is the passive axis
    int longest = std::abs(delta_x);
    int shortest = std::abs(delta_y);

    if(shortest > longest) {
        std::swap(longest, shortest);
        dy_2 = dy_1;
        dx_2 = 0;
    }

    int numerator = longest >> 1;
    OtbImageType::IndexType curr;
    curr[0] = start.x;
    curr[1] = start.y;

    OtbImageType::PixelType max_val = 0;
    OtbImageType::PixelType curr_val = 0;

    for(int i = 0; i < longest; i++) {

        curr_val = image->GetPixel(curr);
        if(curr_val > max_val) {
            curr_val = max_val;
            hit.x = curr[0];
            hit.y = curr[1];
        }

        numerator += shortest;
        if(numerator >= longest) {
            // increment/decrement passive and driving axis by 1
            numerator -= longest;
            curr[0] += dx_1;
            curr[1] += dy_1;
        } else {
            // increment/decrement only driving axis
            curr[0] += dx_2;
            curr[1] += dy_2;
        }
    }

    return max_val;
}

void GradientMagnitudeImage(OtbFloatVectorImageType::Pointer img, wxImage& gradImg) {

    OtbImageType::Pointer gImg = GradientMagnitudeImage(img);
    // convert the output to wxImg

    int numPixels = 3*gradImg.GetWidth()*gradImg.GetHeight();
    PixelTypeUC* pxdt = new PixelTypeUC[numPixels];
    CopyToWxImageData(gImg, pxdt);
    gradImg = wxImage(gradImg.GetWidth(), gradImg.GetHeight(), pxdt, true);
}

OtbImageType::Pointer GradientMagnitudeImage(OtbFloatVectorImageType::Pointer image) {

    OtbFloatImageType::Pointer float_img;
    if(image->GetNumberOfComponentsPerPixel() == 1) {
        typedef otb::ImageList<OtbFloatImageType> FloatImageListType;
        typedef otb::VectorImageToImageListFilter<OtbFloatVectorImageType, FloatImageListType> VectorImageToImageListFilterType;
        VectorImageToImageListFilterType::Pointer image_list_filter = VectorImageToImageListFilterType::New();
        image_list_filter->SetInput(image);
        image_list_filter->Update();
        float_img = image_list_filter->GetOutput()->Back();
    }
    else {
        typedef otb::VectorImageToIntensityImageFilter<OtbFloatVectorImageType, OtbFloatImageType> VectorImageToIntensityImageFilterType;
        VectorImageToIntensityImageFilterType::Pointer intensity_filter = VectorImageToIntensityImageFilterType::New();
        intensity_filter->SetInput(image);
        intensity_filter->Update();
        float_img = intensity_filter->GetOutput();
    }

    typedef itk::GradientMagnitudeImageFilter<OtbFloatImageType, OtbImageType> GradientMagnitudeImageFilterType;
    GradientMagnitudeImageFilterType::Pointer gradient_filter = GradientMagnitudeImageFilterType::New();
    gradient_filter->SetInput(float_img);
    gradient_filter->Update();

    typedef itk::RescaleIntensityImageFilter<OtbImageType, OtbImageType> RescaleIntensityImageFilterType;
    RescaleIntensityImageFilterType::Pointer rescaler = RescaleIntensityImageFilterType::New();
    rescaler->SetOutputMinimum(0);
    rescaler->SetOutputMaximum(255);
    rescaler->SetInput(gradient_filter->GetOutput());
    rescaler->Update();
    return rescaler->GetOutput();
}

