#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    /*
    normalize an image to sum to 1
    */ 
   	for (int k = 0; k < im.c; k++) {
		float channel_sum = 0;
		// first round, sum the channel
		for (int i = 0; i < im.w; i++) {
			for (int j = 0; j < im.h; j++) {
				channel_sum += get_pixel(im, i, j , k); 
			}
		}
		// second round, update the values 
		for (int i = 0; i < im.w; i++) {
			for (int j = 0; j < im.h; j++) {
				set_pixel(im, i, j, k, get_pixel(im, i, j , k) / channel_sum); 
			}
		}
	}
}

image make_box_filter(int w)
{
    /*
    Should return a square image with one channel with uniform entries that sum to 1.
    */
    image box = make_image(w, w, 1);
	for (int i = 0; i < box.w; i++) {
		for (int j = 0; j < box.h; j++) {
			for (int k = 0; k < box.c; k++) {
				set_pixel(box, i, j , k, 1); 
			}
		}
	}
    l1_normalize(box); 
    return box; 
}

image convolve_image(image im, image filter, int preserve)
{
    assert (im.c == filter.c || filter.c == 1);
    image result = make_image(im.w, im.h, im.c);

	// outter for loop updates each pixel in im
    for (int x = 0; x < im.w; x++) {
		for (int y = 0; y < im.h; y++) {
			for (int k = 0; k < im.c; k++) {
				float value = 0;	// stores the result after applying convolution on point (i, j, k)
					// inner for loop calculates the convolution
					for (int f_c = 0; f_c < filter.w; f_c+= 1) {
						for (int f_r = 0; f_r < filter.h; f_r += 1) {
                            // pitfall: the value of filter is determined by the number of its channels
                            float filter_value = (filter.c == 1) ? 
                                get_pixel(filter, f_c, f_r, 0) : get_pixel(filter, f_c, f_r, k);
							int f_x = x - (filter.w) / 2 + f_c; 
							int f_y = y - (filter.h) / 2 + f_r;
							value += filter_value * get_pixel(im, f_x, f_y, k); 
						}
					}
				set_pixel(result, x, y, k, value);
			}
		}
	}

	if (preserve != 1) {
		image sum = make_image(im.w, im.h, 1); 
		for (int i = 0; i < result.w; i++) {
			for (int j = 0; j < result.h; j++) {
				float total = 0;
				for (int k = 0; k < result.c; k++) {
					total += get_pixel(result, i, j , k); 
				}
				set_pixel(sum, i, j , 0, total);
			}
		}
        free_image(result);
		return sum; 
	} else {
		return result; 
	}
}

image make_highpass_filter()
{
    /* 
    returns a image which is [[0, -1, 0], [-1, 4, -1], [0, -1, 0]]
    */
    image filter = make_image(3, 3, 1); 
    set_pixel(filter, 0, 0, 0, 0);
    set_pixel(filter, 1, 0, 0, -1);
    set_pixel(filter, 2, 0, 0, 0);
    set_pixel(filter, 0, 1, 0, -1);
    set_pixel(filter, 1, 1, 0, 4);
    set_pixel(filter, 2, 1, 0, -1);
    set_pixel(filter, 0, 2, 0, 0);
    set_pixel(filter, 1, 2, 0, -1);
    set_pixel(filter, 2, 2, 0, 0);
    return filter;
}

image make_sharpen_filter()
{
    /* 
    returns a image which is [[0, -1, 0], [-1, 5, -1], [0, -1, 0]]
    */
    image filter = make_image(3, 3, 1); 
    set_pixel(filter, 0, 0, 0, 0);
    set_pixel(filter, 1, 0, 0, -1);
    set_pixel(filter, 2, 0, 0, 0);
    set_pixel(filter, 0, 1, 0, -1);
    set_pixel(filter, 1, 1, 0, 5);
    set_pixel(filter, 2, 1, 0, -1);
    set_pixel(filter, 0, 2, 0, 0);
    set_pixel(filter, 1, 2, 0, -1);
    set_pixel(filter, 2, 2, 0, 0);
    return filter;
}

image make_emboss_filter()
{
    /* 
    returns a image which is [[0, -1, 0], [-1, 5, -1], [0, -1, 0]]
    */
    image filter = make_image(3, 3, 1); 
    set_pixel(filter, 0, 0, 0, -2);
    set_pixel(filter, 1, 0, 0, -1);
    set_pixel(filter, 2, 0, 0, 0);
    set_pixel(filter, 0, 1, 0, -1);
    set_pixel(filter, 1, 1, 0, 1);
    set_pixel(filter, 2, 1, 0, 1);
    set_pixel(filter, 0, 2, 0, 0);
    set_pixel(filter, 1, 2, 0, 1);
    set_pixel(filter, 2, 2, 0, 2);
    return filter;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: 
//      It will not be "wrong" for using either perserve or not on all filters. They will still generate results. However, Based on the reason why the filters are used, we will have preference between using preserve or not. 
//      1. High pass filter: We should not preserve when using a high pass filter. The reason of using a high pass filter is to extract the pixels that is different from its neighbors. 
//                           We only care about the magnitude of the difference but not the soucre of difference. Thus, adding up the difference from all channel will give us a better, more contrasted visiual that is proved by 
//                           experiment.
//      2. Sharpen and emboss filter: We should preserve when using both these filters. Both of these filters is used to apply a visiual effect on the image. Thus, all channels need to be preserved otherwise the image will turn 
//                           into a black and white image.
// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: 
//      All these filters should have post-processing (such as clamp_images()) because pixel after convolution might be more then 1 or smaller than 0.

image make_gaussian_filter(float sigma)
{
    /*
    returns a image with size (ceil(sigma * 6), ceil(sigma * 6), 1), each pixel follows the formula of 2D Gaussian functions 
    */
    int size = ceil(sigma * 6);
    if (size % 2 == 0) size += 1; 
    image result = make_image(size, size, 1); 
    for (int x = 0; x < result.w; x++) {
        for (int y = 0; y < result.h; y++) {
            float gx = x - (size - 1) / 2 ; // offset so that the mid point is 0
            float gy = y - (size - 1) / 2; // offset so that the mid point is 0
            float value = expf(-(gx * gx + gy * gy)/(2 * sigma * sigma)) / (TWOPI * sigma * sigma);
            set_pixel(result, x, y, 0, value);
        }
    }
    l1_normalize(result);   // since it is a blurring filter, we need to normalize it 
    return result;
}

image add_image(image a, image b)
{
    /* 
    return a image with each pixel is the sum of the corresponding pixel from a and b
    */
    assert(a.w == b.w && a.h == b.h && a.c == b.c); 
    image result = make_image(a.w, a.h, a.c); 
    for (int x = 0; x < result.w; x++) {
        for (int y = 0; y < result.h; y++) {
            for (int c = 0; c < result.c; c++) {
                set_pixel(result, x, y, c, get_pixel(a, x, y, c) + get_pixel(b, x, y, c)); 
            }
        }
    }
    return result;
}

image sub_image(image a, image b)
{
    /* 
    return a image with each pixel is the difference of the corresponding pixel from a and b
    */
    assert(a.w == b.w && a.h == b.h && a.c == b.c); 
    image result = make_image(a.w, a.h, a.c); 
    for (int x = 0; x < result.w; x++) {
        for (int y = 0; y < result.h; y++) {
            for (int c = 0; c < result.c; c++) {
                set_pixel(result, x, y, c, get_pixel(a, x, y, c) - get_pixel(b, x, y, c)); 
            }
        }
    }
    return result;
}

image make_gx_filter()
{
    /* 
    returns a image which is [[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]]
    */
    image filter = make_image(3, 3, 1); 
    set_pixel(filter, 0, 0, 0, -1);
    set_pixel(filter, 1, 0, 0, 0);
    set_pixel(filter, 2, 0, 0, 1);
    set_pixel(filter, 0, 1, 0, -2);
    set_pixel(filter, 1, 1, 0, 0);
    set_pixel(filter, 2, 1, 0, 2);
    set_pixel(filter, 0, 2, 0, -1);
    set_pixel(filter, 1, 2, 0, 0);
    set_pixel(filter, 2, 2, 0, 1);
    return filter;
}


image make_gy_filter()
{
    /* 
    returns a image which is [[-1, -2, -1], [0, 0, 0], [1, 2, 1]]
    */
    image filter = make_image(3, 3, 1); 
    set_pixel(filter, 0, 0, 0, -1);
    set_pixel(filter, 1, 0, 0, -2);
    set_pixel(filter, 2, 0, 0, -1);
    set_pixel(filter, 0, 1, 0, 0);
    set_pixel(filter, 1, 1, 0, 0);
    set_pixel(filter, 2, 1, 0, 0);
    set_pixel(filter, 0, 2, 0, 1);
    set_pixel(filter, 1, 2, 0, 2);
    set_pixel(filter, 2, 2, 0, 1);
    return filter;
}

void feature_normalize(image im)
{
    /*
    for each channel in im, we normalize it by subtracting each pixel by the minimal value and divide by the range. If the range is 0, set all pixels fo zero
    */
    for (int channel = 0; channel < im.c; channel++) {
        float c_min = INFINITY; 
        float c_max = -INFINITY; 
        for (int x = 0; x < im.w; x++) {
            for (int y = 0; y < im.h; y++) {
                float val = get_pixel(im, x, y, channel);
                c_max = fmax(c_max, val); 
                c_min = fmin(c_min, val);
            }
        }
        for (int x = 0; x < im.w; x++) {
            for (int y = 0; y < im.h; y++) {
                if (c_min == c_max) {
                    set_pixel(im, x, y, channel, 0);
                }
                else {
                    set_pixel(im, x, y, channel, get_pixel(im, x, y, channel) / (c_max - c_min));
                }
            }
        }
    }
}

image *sobel_image(image im)
{
    /*
    Return two images, both with the size of (im.w, im.h, 1). The first one have the magnitude and the second one have the angle
    */
    image* result = calloc(2, sizeof(image));
    image x_filter = make_gx_filter(); 
    image y_filter = make_gy_filter();
    image gx = convolve_image(im, x_filter, 0); 
    image gy = convolve_image(im, y_filter, 0);
    image magnitude = make_image(im.w, im.h, 1); 
    image angle = make_image(im.w, im.h, 1); 
    for (int i = 0; i < gx.w * gx.h; i++) {
        float mag = sqrtf(powf(gx.data[i], 2) + powf(gy.data[i], 2));
        float ang = atan2(gy.data[i], gx.data[i]);
        magnitude.data[i] = mag; 
        angle.data[i] = ang; 
    }
    result[0] = magnitude; 
    result[1] = angle; 
    free_image(x_filter);
    free_image(y_filter);
    free_image(gx);
    free_image(gy);
    return result;
}

image colorize_sobel(image im)
{
    /*
    returns one image with RGB colored edge and other pixels black
    first use sobel_image() to seperate the angle and the magnitude of the edges
    then turn the magnitude into saturation and value and angle to hue
    finally transform this image to RGB and return 
    */
    image* images = sobel_image(im); 
    feature_normalize(images[0]);
    feature_normalize(images[1]);   // post-processing after applying a filter
    image result = make_image(im.w, im.h, 3);   // the result image will be stored in RGB
    for (int x = 0; x < result.w; x++) {
        for (int y = 0; y < result.h; y++) {
            set_pixel(result, x, y, 0, get_pixel(images[1], x, y, 0)); // set the hue to angle  
            set_pixel(result, x, y, 1, get_pixel(images[0], x, y, 0)); // set saturation to magnitude
            set_pixel(result, x, y, 2, get_pixel(images[0], x, y, 0)); // set value to magnitude
        }
    }
    hsv_to_rgb(result); 
    free_image(images[0]);
    free_image(images[1]);
    return result;
}
