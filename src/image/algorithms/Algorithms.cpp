#include "Algorithms.hpp"
#include "../../geometry/Primitives.hpp"
#include <otbImageFileWriter.h>
#include <otbImportVectorImageFilter.h>
#include <otbVectorImageToIntensityImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkNeighborhoodIterator.h>
#include <wx/image.h>
#include <list>
#include <map>
#include <cmath>

typedef itk::ImageRegionIterator<OtbImageType>      ItkImgIt_IteratorType;
typedef itk::ImageRegionConstIterator<OtbImageType> ItkImgIt_ConstIteratorType;
typedef itk::NeighborhoodIterator<OtbImageType>     ItkImgIt_NeighborhoodIteratorType;
typedef otb::ImageFileWriter<OtbImageType>          OtbImageWriterType;

OtbImageType::Pointer region_grow(const wxImage& wxImg,
                                  const wxPoint& pt,
                                  int threshold) {

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
    PixelType* dt = wxImg.GetData();
    PixelType* data = new PixelType[numPixels];
    std::copy(dt, dt + numPixels, data);
    importFilter->SetImportPointer(data, numPixels, true);
    OtbVectorImageType::Pointer vecImg = OtbVectorImageType::New();
    vecImg = dynamic_cast<OtbVectorImageType*>(importFilter->GetOutput());

    // Step-2: Convert OtbVectorImage to OtbImage
    OtbImageType::Pointer intensity_image = OtbImageType::New();
    typedef otb::VectorImageToIntensityImageFilter<
            OtbVectorImageType, OtbImageType> IntensityFilter;
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
        it.Set(static_cast<PixelType>(PxlValues::UNKNOWN));
    }

    // Step-4: Define the neighborhood iterators for the intensity and labels images
    ItkImgIt_NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);
    ItkImgIt_NeighborhoodIteratorType it1(radius, intensity_image,
                                          intensity_image->GetLargestPossibleRegion());
    ItkImgIt_NeighborhoodIteratorType it2(radius, labels_image,
                                          labels_image->GetLargestPossibleRegion());
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
        it2.SetCenterPixel(static_cast<PixelType>(PxlValues::SELECTED));
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
                    it2.SetPixel(i, static_cast<PixelType>(PxlValues::REFUSED));
                }
            }
        }
        process.pop_front();

    }
    return labels_image;
}

bool region_grow_segmentation(const wxImage& wxImg, const wxPoint& pt, int threshold) {
    OtbImageType::Pointer img = region_grow(wxImg, pt, threshold);
    OtbImageWriterType::Pointer writer = OtbImageWriterType::New();
    writer->SetFileName("deneme.jpeg");
    writer->SetInput(img);
    writer->Update();
    return true;
}

bool region_grow_segmentation(const wxImage& wxImg,
                              wxImage& segImg,
                              const wxPoint& pt,
                              int threshold) {
    OtbImageType::Pointer img = region_grow(wxImg, pt, threshold);
    int numPixels = 3*wxImg.GetWidth()*wxImg.GetHeight();

    PixelType* pxdt = new PixelType[numPixels];
    PixelType* dit = pxdt;
    ItkImgIt_ConstIteratorType it(img, img->GetLargestPossibleRegion());
    OtbImageType::PixelType pxtype;
    for(it.GoToBegin(); !it.IsAtEnd(); ++it) {
        pxtype = it.Get();
        *dit = pxtype;
        ++dit;
        *dit = pxtype;
        ++dit;
        *dit = pxtype;
        ++dit;
    }
    segImg = wxImage(wxImg.GetWidth(), wxImg.GetHeight(), pxdt, true);
    return true;
}

// image is a binary image: pixel_values are either 255 or 0
// returns true if the ray hits an occupied pixel (value = 0)
// start: starting point for the ray segment
// end: end point for the ray segement
bool ray_cast(const OtbImageType::Pointer& image,
              const OtbImageType::IndexType& start,
              const OtbImageType::IndexType& end,
              OtbImageType::IndexType& first_hit) {
    int delta_x = end[0] - start[0];
    int delta_y = end[1] - start[1];
    int dx_1(0), dy_1(0), dx_2(0), dy_2(0);
    if(delta_x < 0) dx_1 = -1; else if (delta_x > 0) dx_1 = 1;
    if(delta_y < 0) dy_1 = -1; else if (delta_y > 0) dy_1 = 1;
    dx_2 = dx_1;

    int longest = std::abs(delta_x);
    int shortest = std::abs(delta_y);
    if(shortest > longest) {
        std::swap(longest, shortest);
        dy_2 = dy_1;
        dx_2 = 0;
    }

    int numerator = longest >> 1;
    OtbImageType::IndexType pixelIndex;
    pixelIndex[0] = start[0];
    pixelIndex[1] = start[1];

    for(int i = 0; i < longest; i++) {
        if(image->GetPixel(pixelIndex) == 0) {
            first_hit[0] = pixelIndex[0];
            first_hit[1] = pixelIndex[1];
            return true;
        }
        numerator += shortest;
        if(numerator >= longest) {
            numerator -= longest;
            pixelIndex[0] += dx_1;
            pixelIndex[1] += dy_1;
        } else {
            pixelIndex[0] += dx_2;
            pixelIndex[1] += dy_2;
        }
    }
    return false;
}

bool ray_cast(const OtbImageType::Pointer& image, const Point2D<int>& start, const Point2D<int>& end, Point2D<int>& first_hit) {

    int delta_x = end.x - start.x;
    int delta_y = end.y - start.y;
    int dx_1(0), dy_1(0), dx_2(0), dy_2(0);

    if(delta_x < 0)
        dx_1 = -1;
    else if (delta_x > 0)
        dx_1 = 1;
    if(delta_y < 0)
        dy_1 = -1;
    else if (delta_y > 0)
        dy_1 = 1;
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

