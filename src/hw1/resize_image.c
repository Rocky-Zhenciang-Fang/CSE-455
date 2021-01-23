#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{   
    // use build-in roundf
    return get_pixel(im, roundf(x), roundf(y), c);
}

image nn_resize(image im, int w, int h)
{  
    image result = make_image(w, h, im.c); 
    float a_w = (float)(im.w) / w; 
    float b_w = (float)(im.w - w) / (2 * w);
    float a_h = (float)(im.h) / h; 
    float b_h = (float)(im.h - h) / (2 * h);
    for (int i = 0; i < im.c; i++) { // for each channel
        for (int j = 0; j < w; j ++) { 
            for (int k = 0; k < h; k ++) {
                float old_x = a_w * j + b_w; 
                float old_y = a_h * k + b_h;
                set_pixel(result, j, k ,i, nn_interpolate(im, old_x, old_y, i)); 
            }
        }
    }
    return result;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // billnear interpolate involves four points 
    // upleft -> [floor(x), floor(y)]
    // upright -> [floor(x), ceil(y)]
    // buttomleft -> [ceil(x), floor(y)]
    // buttomright -> [ceil(x), ceil(y)]
    // Do linear interpolate twice to get the result
    float upleft = get_pixel(im, floorf(x), floorf(y), c); 
    float upright = get_pixel(im, ceilf(x), floorf(y), c); 
    float buttomleft = get_pixel(im, floorf(x), ceilf(y), c); 
    float buttomright = get_pixel(im, ceilf(x), ceilf(y), c);
    float ql = (y - floorf(y)) * buttomleft + (ceilf(y) - y) * upleft;
    float qr = (y - floorf(y)) * buttomright + (ceilf(y) - y) * upright; 
    float q = (x - floorf(x)) * qr + (ceilf(x) - x) * ql;
    return q;
}

image bilinear_resize(image im, int w, int h)
{
    // similar to nn_resize, but change nn_interpolate() to bilinear_interpolate()
    image result = make_image(w, h, im.c); 
    float a_w = (float)(im.w) / w; 
    float b_w = (float)(im.w - w) / (2 * w);
    float a_h = (float)(im.h) / h; 
    float b_h = (float)(im.h - h) / (2 * h);
    for (int i = 0; i < im.c; i++) { // for each channel
        for (int j = 0; j < w; j ++) { 
            for (int k = 0; k < h; k ++) {
                float old_x = a_w * j + b_w; 
                float old_y = a_h * k + b_h;
                set_pixel(result, j, k ,i, bilinear_interpolate(im, old_x, old_y, i)); 
            }
        }
    }
    return result;
}

