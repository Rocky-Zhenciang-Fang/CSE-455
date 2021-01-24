#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    int size = ceil(sigma * 6);
    if (size % 2 == 0) size += 1; 
    image result = make_image(size, 1, 1); 
    for (int x = 0; x < result.w; x++) {
        float gx = x - (size - 1) / 2 ; // offset so that the mid point is 0
        float value = expf(-(gx * gx)/(2 * sigma * sigma)) / sqrtf(TWOPI * sigma);
        set_pixel(result, x, 0, 0, value);
    }
    l1_normalize(result);   // since it is a blurring filter, we need to normalize it 
    return result;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    if(1 == 0){
        image g = make_gaussian_filter(sigma);
        image s = convolve_image(im, g, 1);
        free_image(g);
        return s;
    } else {
        image hor_filter = make_1d_gaussian(sigma); 
        image ver_filter = make_image(hor_filter.h, hor_filter.w, hor_filter.c);
        for (int x = 0; x < hor_filter.w; x++) {
            set_pixel(ver_filter, 0, x, 0, get_pixel(hor_filter, x, 0, 0)); 
        }
        image result = convolve_image(im, hor_filter, 1); 
        result = convolve_image(result, ver_filter, 1); 
        return result;
    }
}

// Calculate the structure matrix of an image.
// First caculate the direvitive of the pixels, remember of burr it first
// Then store the values in the respective channels
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);
    image gx = make_gx_filter(); 
    image gy = make_gy_filter(); 
    image Ix = convolve_image(im, gx, 0);
    image Iy = convolve_image(im, gy, 0); 
    for (int x = 0; x < S.w; x++) {
        for (int y = 0; y < S.h; y++) {
            float x_val = get_pixel(Ix, x, y, 0);
            float y_val = get_pixel(Iy, x, y, 0);
            set_pixel(S, x, y, 0, x_val * x_val);
            set_pixel(S, x, y, 1, y_val * y_val);
            set_pixel(S, x, y, 2, y_val * x_val);
        }
    }
    image result = smooth_image(S, sigma); 

    free_image(S);
    free_image(Ix); 
    free_image(Iy);
    free_image(gx);
    free_image(gy);
    return result;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
    for (int x = 0; x < S.w; x++) {
        for (int y = 0; y < S.h; y++) {
            float topleft = get_pixel(S, x, y, 0);
            float buttomright = get_pixel(S, x, y, 1); 
            float buttomleft = get_pixel(S, x, y, 2);
            float topright  = buttomleft; 
            float value = topleft * buttomright - topright * buttomleft - 0.06 * powf((topleft + buttomright), 2);
            set_pixel(R, x, y, 0, value); 
        }
    }
    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    for (int x = 0; x < im.w; x++) {
        for (int y = 0; y < im.h; y++) {
            // this inner loop loop over [2w + 1] window
            float current = get_pixel(im, x, y, 0); 
            for (int nnm_x = x - w; nnm_x < x +  w + 1; nnm_x++) {
                for (int nnm_y = y - w; nnm_y < y + w + 1; nnm_y++) {
                    float compare = get_pixel(im, nnm_x, nnm_y, 0);
                    if (compare > current) { 
                        set_pixel(r, x, y, 0, -999999); 
                        break;
                    }
                }
            }
        }
    }
    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{   
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);

    int count = 0; // change this
    for (int i = 0; i < Rnms.w * Rnms.h * Rnms.c; i++) {
        if (Rnms.data[i] > thresh) {
            count += 1; 
        }
    } 
    

    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    descriptor *ptr = d; 
    for (int i = 0; i < Rnms.w * Rnms.h * Rnms.c; i++) {
        if (Rnms.data[i] > thresh) {
            d[count - 1] = describe_index(im, i); 
            count -= 1; 
        }
    } 

    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
