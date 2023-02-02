///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//                                              Modified:   Feng Liu
//                                              Date:       Winter 2011
//                                              Why:        Change the library file 
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
    unsigned char image_gray;
    unsigned char rgb[3];

    if (!data)
        return NULL; 

    for (int i = 0; i < width * height * 4; i+=4) {

        RGBA_To_RGB(data + i, rgb); //Removes Alpha Channel
        image_gray = (unsigned char)((float)rgb[0] * 0.299 + (float)rgb[1] * 0.587 + (float)rgb[2] * 0.114);
        data[i + RED] = image_gray;
        data[i + GREEN] = image_gray;
        data[i + BLUE] = image_gray;
    }
    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    if (!data)
        return NULL;
    unsigned char rgb[3];
    int offset;
    double alpha;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            offset = i * width * 4 + j * 4;
            RGBA_To_RGB(data + offset, rgb); //Remove alpha channel
            alpha = (double)*(data + offset + 3) / 255.0;
            *(data + offset + RED) = floor(double(rgb[RED] / 32) * 32 * alpha);
            *(data + offset + GREEN) = floor(double(rgb[GREEN] / 32) * 32 * alpha);
            *(data + offset + BLUE) = floor(double(rgb[BLUE] / 64) * 64 * alpha);
        }
    }
    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
    ClearToBlack();
    return false;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    if (!data)
        return NULL;

    To_Grayscale();

    unsigned char rgb[3];
    float threshold = 0.5;
    unsigned char white = 255;
    unsigned char black = 0;
    float normalized_pixel;

    for (int i = 0; i < width * height * 4; i += 4) {

        RGBA_To_RGB(data + i, rgb);
        normalized_pixel = rgb[0] / (float)256;

        if (normalized_pixel < threshold) {

            data[i + RED] = black;
            data[i + GREEN] = black;
            data[i + BLUE] = black;
        }
        else {

            data[i + RED] = white;
            data[i + GREEN] = white;
            data[i + BLUE] = white;
        }
    }
    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    
    if (!data)
        return NULL;
    
    To_Grayscale();

    unsigned char rgb[3];
    float threshold = 0.5;
    unsigned char white = 255;
    unsigned char black = 0;
    float normalized_pixel;
    float max = -0.2;
    float min = 0.2;

    for (int i = 0; i < width * height * 4; i += 4) {

        RGBA_To_RGB(data + i, rgb);
        float random = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
        normalized_pixel = (rgb[0] / (float)256) + random;

        if (normalized_pixel < threshold) {

            data[i + RED] = black;
            data[i + GREEN] = black;
            data[i + BLUE] = black;
        }
        else {

            data[i + RED] = white;
            data[i + GREEN] = white;
            data[i + BLUE] = white;
        }
    }
    return true;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::DitherFS_Pixel(int x, int y, double amt) {
    double v;
    if (x < 0 || x >= width)
        return false;
    if (y < 0 || y >= height)
        return false;

    v = *(data + (y * width * 4) + (x * 4));
    v += amt;
    *(data + (y * width * 4) + (x * 4) + RED) = v;
    *(data + (y * width * 4) + (x * 4) + GREEN) = v;
    *(data + (y * width * 4) + (x * 4) + BLUE) = v;
    *(data + (y * width * 4) + (x * 4) + 3) = 255;

    return true;
}

bool TargaImage::Dither_FS()
{
    if (!data)
        return NULL;

    To_Grayscale();

    int z_width = width - 1;
    int z_i = 0;
    int z_dir = 1;

    for (int j = 0; j < height; j++) {
        for (int i= z_i; (z_dir == 1) ? (i <= z_width) : (i >= z_width); i += z_dir) {
            double e;
            unsigned char bw;
            unsigned char v = *(data + (j * width * 4) + (i * 4));

            bw = (v > 128) ? 255 : 0;
            e = v - bw;

            *(data + (j * width * 4) + (i * 4) + RED) = bw;
            *(data + (j * width * 4) + (i * 4) + GREEN) = bw;
            *(data + (j * width * 4) + (i * 4) + BLUE) = bw;
            *(data + (j * width * 4) + (i * 4) + 3) = 255;

            DitherFS_Pixel(i + z_dir, j + 0, (7.0 * e) / 16.0);
            DitherFS_Pixel(i - z_dir, j + 1, (3.0 * e) / 16.0);
            DitherFS_Pixel(i + 0, j + 1, (5.0 * e) / 16.0);
            DitherFS_Pixel(i + z_dir, j + 1, (1.0 * e) / 16.0);
        }
        int k = z_width;
        z_width = z_i;
        z_i = k;
        z_dir = -z_dir;
    }
    return true;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    if (!data)
        return NULL;

    To_Grayscale();
    
    unsigned char white = 255;
    unsigned char black = 0;
    double sum=0;
    vector<unsigned char> image;

    for (int i = 0; i < width * height * 4; i += 4) {

        unsigned char rgb[3];
        RGBA_To_RGB(data + i, rgb);
        sum += rgb[0];
        image.push_back(rgb[0]);
    }

    float average_brightness = sum / ((double)width * height) / 255.0;
    float threshold_index = (1 - average_brightness) * (width * height);

    sort(image.begin(), image.end());
    float threshold_pixel = image[threshold_index];

    for (int i = 0; i < width * height * 4; i += 4) {
        unsigned char rgb[3];
        RGBA_To_RGB(data + i, rgb);
        if (rgb[0] < threshold_pixel) {

            data[i + RED] = black;
            data[i + GREEN] = black;
            data[i + BLUE] = black;
        }
        else {

            data[i + RED] = white;
            data[i + GREEN] = white;
            data[i + BLUE] = white;
        }
    }
    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    if (!data)
        return NULL;
    float cluster[4][4] = { {0.7500, 0.3750, 0.6250, 0.2500},
                            {0.0625, 1.0000, 0.8750, 0.4375},
                            {0.5000, 0.8125, 0.9375, 0.1250},
                            {0.1875, 0.5625, 0.3125, 0.6875} };
    To_Grayscale();

    unsigned char rgb[3];
    unsigned char white = 255;
    unsigned char black = 0;
    float normalized_pixel;
    int offset;

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {

            offset = ((i * height) + j) * 4;
            RGBA_To_RGB(data + offset, rgb);
            normalized_pixel = rgb[0] / (float)256;
            if (normalized_pixel < cluster[i%4][j%4]) {

                data[offset + RED] = black;
                data[offset + GREEN] = black;
                data[offset + BLUE] = black;
                data[offset + 3] = white; // Alpha channel
            }
            else {

                data[offset + RED] = white;
                data[offset + GREEN] = white;
                data[offset + BLUE] = white;
                data[offset + 3] = white; // Alpha channel
            }
        }
    }
    return true;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    ClearToBlack();
    return false;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < width * height * 4; i += 4) {
        double a = (double)data[i + 3] / 255.0;
        for (int c = 0; c < 4; c++) {
            data[i + c] += pImage->data[i + c] * (1.0 - a);
        }
    }

    return true;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < width * height * 4; i += 4) {
        double a = (double)pImage->data[i + 3] / 255.0;
        for (int c = 0; c < 4; c++) {
            data[i + c] *= a;
        }
    }
    return true;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < width * height * 4; i += 4) {
        double a = (double)pImage->data[i + 3] / 255.0;
        for (int c = 0; c < 4; c++) {
            data[i + c] *= (1.0 - a);
        }
    }

    return true;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < width * height * 4; i += 4) {
        double a_f = (double)data[i + 3] / 255.0;
        double a_g = (double)pImage->data[i + 3] / 255.0;
        for (int c = 0; c < 4; c++) {
            data[i + c] = (data[i + c] * a_g) + (pImage->data[i + c] * (1.0 - a_f));
        }
    }

    return true;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < width * height * 4; i += 4) {
        double a_f = ((double)data[i + 3] / 255.0);
        double a_g = ((double)pImage->data[i + 3] / 255.0);
        for (int c = 0; c < 4; c++) {
            data[i + c] = (data[i + c] * (1.0 - a_g)) + (pImage->data[i + c] * (1.0 - a_f));
        }
    }
    return true;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference

bool TargaImage::applyFilter(double filter[5][5])
{
    unsigned char* rgb = To_RGB();
    double sum = 0;

    // Iterate through the image pixel by pixel
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {

            // Need to perform filtering for R,G and B separately
            for (int color = 0; color < 3; color++) {

                sum = 0;
                // Traverse through the filter values
                for (int row = 0; row < 5; row++) {
                    for (int col = 0; col < 5; col++) {

                        int row_pos = i - 2 + row;
                        int col_pos = j - 2 + col;

                        if (row_pos < 0)
                            row_pos = -row_pos;
                        else if (row_pos > height)
                            row_pos = ((height * 2) - 2) - row_pos;

                        if (col_pos < 0)
                            col_pos = -col_pos;
                        else if (col_pos > width)
                            col_pos = ((width * 2) - 2) - col_pos;

                        sum += rgb[((row_pos * width + col_pos) * 3) + color] * filter[row][col];
                    }
                }
                data[((i * width + j) * 4) + color] = sum;
                if (sum > 255)
                    data[((i * width + j) * 4) + color] = 255;
                else if(sum<0)
                    data[((i * width + j) * 4) + color] = 0;
            }
        }
    }
    delete[] rgb;
    return true;
}
///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
    double filter[5][5] = { {1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0},
                       {1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0},
                       {1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0},
                       {1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0},
                       {1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0, 1.0 / 25.0} };
    applyFilter(filter);
    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    double filter[5][5] = { {1.0 / 81.0, 2.0 / 81.0, 3.0 / 81.0, 2.0 / 81.0, 1.0 / 81.0},
                          {2.0 / 81.0, 4.0 / 81.0, 6.0 / 81.0, 4.0 / 81.0, 2.0 / 81.0},
                          {3.0 / 81.0, 6.0 / 81.0, 9.0 / 81.0, 6.0 / 81.0, 3.0 / 81.0},
                          {2.0 / 81.0, 4.0 / 81.0, 6.0 / 81.0, 4.0 / 81.0, 2.0 / 81.0},
                          {1.0 / 81.0, 2.0 / 81.0, 3.0 / 81.0, 2.0 / 81.0, 1.0 / 81.0} };
    applyFilter(filter);
    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    double filter[5][5] = { {1.0 / 256.0, 4.0 / 256.0, 6.0 / 256.0, 4.0 / 256.0, 1.0 / 256.0},
                          {4.0 / 256.0, 16.0 / 256.0, 24.0 / 256.0, 16.0 / 256.0, 4.0 / 256.0},
                          {6.0 / 256.0, 24.0 / 256.0, 36.0 / 256.0, 24.0 / 256.0, 6.0 / 256.0},
                          {4.0 / 256.0, 16.0 / 256.0, 24.0 / 256.0, 16.0 / 256.0, 4.0 / 256.0},
                          {1.0 / 256.0, 4.0 / 256.0, 6.0 / 256.0, 4.0 / 256.0, 1.0 / 256.0} };
    applyFilter(filter);
    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
    ClearToBlack();
   return false;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    double filter[5][5] = { {-1.0 / 256.0, -4.0 / 256.0, -6.0 / 256.0, -4.0 / 256.0, -1.0 / 256.0},
                          {-4.0 / 256.0, -16.0 / 256.0, -24.0 / 256.0, -16.0 / 256.0, -4.0 / 256.0},
                          {-6.0 / 256.0, -24.0 / 256.0, 220.0 / 256.0, -24.0 / 256.0, -6.0 / 256.0},
                          {-4.0 / 256.0, -16.0 / 256.0, -24.0 / 256.0, -16.0 / 256.0, -4.0 / 256.0},
                          {-1.0 / 256.0, -4.0 / 256.0, -6.0 / 256.0, -4.0 / 256.0, -1.0 / 256.0} };
    applyFilter(filter);
    return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    double filter[5][5] = { {-1.0 / 256.0, -4.0 / 256.0, -6.0 / 256.0, -4.0 / 256.0, -1.0 / 256.0},
                          {-4.0 / 256.0, -16.0 / 256.0, -24.0 / 256.0, -16.0 / 256.0, -4.0 / 256.0},
                          {-6.0 / 256.0, -24.0 / 256.0, 476.0 / 256.0, -24.0 / 256.0, -6.0 / 256.0},
                          {-4.0 / 256.0, -16.0 / 256.0, -24.0 / 256.0, -16.0 / 256.0, -4.0 / 256.0},
                          {-1.0 / 256.0, -4.0 / 256.0, -6.0 / 256.0, -4.0 / 256.0, -1.0 / 256.0} };
    applyFilter(filter);
    return true;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    double filter[3][3] = {{1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0},
                           {1.0 / 8.0, 1.0 / 4.0, 1.0 / 8.0},
                           {1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0}};
    unsigned char* rgb = To_RGB();
    unsigned char* half_size = new unsigned char[width * height];
    for (int x = 0; x < height; x += 2) {
        for (int y = 0; y < width; y += 2) {
            for (int color = 0; color < 3; color++) {
                double sum = 0;
                for (int r = 0; r < 3; r++) {
                    for (int c = 0; c < 3; c++) {
                        int row_pos = x - 1 + r;
                        int col_pos = y - 1 + c;
                        if (row_pos < 0)
                            row_pos = -row_pos;
                        if (row_pos >= height)
                            row_pos = ((height * 2) - 1) - row_pos;
                        if (col_pos < 0)
                            col_pos = -col_pos;
                        if (col_pos >= width)
                            col_pos = ((width * 2) - 1) - col_pos;
                        sum += rgb[(row_pos * width + col_pos) * 3 + color] * filter[r][c];
                    }
                }
                half_size[(x / 2 * (width / 2) + y / 2) * 4 + color] = sum;
            }
            half_size[(x / 2 * (width / 2) + y / 2) * 4 + 4] = 255;
        }
    }
    delete[] data;
    data = new unsigned char[width * height];
    for (int i = 0; i < (width * height); i++) {
        data[i] = half_size[i];
    }
    height = height / 2;
    width = width / 2;
    delete[] rgb;
    return true;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    ClearToBlack();
    return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

