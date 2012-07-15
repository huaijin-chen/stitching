///////////////////////////////////////////////////////////////////////////
//
// NAME
//  BlendImages.cpp -- blend together a set of overlapping images
//
// DESCRIPTION
//  This routine takes a collection of images aligned more or less horizontally
//  and stitches together a mosaic.
//
//  The images can be blended together any way you like, but I would recommend
//  using a soft halfway blend of the kind Steve presented in the first lecture.
//
// SEE ALSO
//  BlendImages.h       longer description of parameters
//
// Copyright ?Richard Szeliski, 2001.  See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "ImageLib/ImageLib.h"
#include "BlendImages.h"
#include "WarpSpherical.h"
#include <math.h>
#include <iostream>

#ifndef M_PI 
#define M_PI    3.1415926536f
#endif // M_PI

#ifndef LoopImage
#define LoopImage(x,y,width,height)	for(int y = 0; y < height; y++)for(int x = 0; x < width; x++)
#endif

#ifndef isZero
#define isZero(x)	(abs(x) < 1e-10)
#endif

/******************* TO DO *********************
 * SetImageAlpha:
 *	INPUT:
 *		img: an image to be added to the final panorama
 *		blendRadius: radius of the blending function
 *	OUTPUT:
 *		set the alpha values of img
 *		pixels near the center of img should have higher alpha
 *		use blendRadius to determine how alpha decreases
 */
static void SetImageAlpha(CByteImage& img, float blendRadius)
{
	// *** BEGIN TODO ***
	// fill in this routine..
	CShape sh = img.Shape();
	float cx = sh.width / 2.0f, cy = sh.height / 2.0f;
	
	LoopImage(x,y,sh.width,sh.height)
	{
		float wx = 1.0f - abs(x - cx) / (sh.width / 2.0f);
		float wy = 1.0f - abs(y - cy) / (sh.height / 2.0f);
		img.Pixel(x,y,3) = (uchar)(wx * wy * 255);
	}

	// *** END TODO ***
}

/******************* TO DO *********************
 * AccumulateBlend:
 *	INPUT:
 *		img: a new image to be added to acc
 *		acc: portion of the accumulated image where img is to be added
 *		blendRadius: radius of the blending function
 *	OUTPUT:
 *		add a weighted copy of img to the subimage specified in acc
 *		the first 3 band of acc records the weighted sum of pixel colors
 *		the fourth band of acc records the sum of weight
 */
static void AccumulateBlend(CByteImage& img, CFloatImage& acc, float blendRadius)
{
	// *** BEGIN TODO ***
	// fill in this routine..
	CShape sh = img.Shape();
	if((sh.nBands != 4) || (img.alphaChannel != 3))
	{
		cout<<"AccumulateBlend: image should be converted to rgba image out of the function."<<endl;
		return;
	}

	LoopImage(x,y,sh.width,sh.height)
	{
		float alpha = img.Pixel(x,y,3) / 255.f;
		for(int c = 0; c < 3; c++)
		{
			acc.Pixel(x,y,c) = acc.Pixel(x,y,c) + img.Pixel(x,y,c) / 255.f * alpha;
		}
		acc.Pixel(x,y,c) = acc.Pixel(x,y,3) + alpha;
	}
	
	// *** END TODO ***
}

/******************* TO DO *********************
 * NormalizeBlend:
 *	INPUT:
 *		acc: input image whose alpha channel (4th channel) contains
 *		     normalizing weight values
 *		img: where output image will be stored
 *	OUTPUT:
 *		normalize r,g,b values (first 3 channels) of acc and store it into img
 */
static void NormalizeBlend(CFloatImage& acc, CByteImage& img)
{
	// *** BEGIN TODO ***
	// fill in this routine..
	CShape sh = img.Shape();
	if((sh.nBands != 4) || (img.alphaChannel != 3))
	{
		cout<<"NormalizeBlend: image should be converted to rgba image out of the function."<<endl;
		return;
	}

	LoopImage(x,y,sh.width,sh.height)
	{
		for (int c = 0; c < 3; c++)
		{
			img.Pixel(x,y,c) = (uchar)((acc.Pixel(x,y,c) / acc.Pixel(x,y,3))*255);
		}
		img.Pixel(x,y,3) = 255;
	}
	// *** END TODO ***
}

/******************* TO DO *********************
 * BlendImages:
 *	INPUT:
 *		ipv: list of input images and their global positions in the mosaic
 *		f: focal length
 *		blendRadius: radius of the blending function
 *	OUTPUT:
 *		create & return final mosaic by blending all images
 */
CByteImage BlendImages(CImagePositionV& ipv, float f, float blendRadius)
{
    // Assume all the images are of the same shape (for now)
    CByteImage& img0 = ipv[0].img;
    CShape sh        = img0.Shape();
    int width        = sh.width;
    int height       = sh.height;
    int nBands       = sh.nBands;
    int dim[2]       = {width, height};

	int nTheta = (int) (2*M_PI*f + 0.5);
	int nPhi = (int) (M_PI*f + 0.5);

    // Create a floating point accumulation image
    CShape mShape(nTheta, nPhi, nBands);
    CFloatImage accumulator(mShape);
    accumulator.ClearPixels();

    // Add in all of the images
	for (unsigned int i = 0; i < ipv.size(); i++)
    {
        // Warp the image into spherical coordinates
        CTransform3x3 M = ipv[i].position;

        CByteImage& src = ipv[i].img;
		CByteImage dst(mShape);
	    if (src.Shape().nBands != 4)
			src = ConvertToRGBA(src);
		SetImageAlpha(src, blendRadius);
		CFloatImage uv = WarpSphericalField(src.Shape(), dst.Shape(), f, M);
		WarpLocal(src, dst, uv, false, eWarpInterpLinear);

        // Perform the accumulation
        AccumulateBlend(dst, accumulator, blendRadius);
    }

    // Normalize the results
    CByteImage compImage(mShape);
    NormalizeBlend(accumulator, compImage);

    return compImage;
}
