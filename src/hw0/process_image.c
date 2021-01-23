#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    if (x < 0) {
        x = 0;
    }
    else if (x > im.w - 1) {
        x = im.w - 1;
    }
    if (y < 0) {
        y = 0;
    }
    else if (y > im.h - 1) {
        y = im.h - 1;
    }
    if (c < 0) {
        c = 0;
    }
    else if (c > 3) { 
        c = 3;
    }
    return im.data[x + (im.w * y) + (im.w * im.h * c)];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    if (0 <= x && 0 <= y && 0 <= c && x < im.w && y < im.h && c < 3) {
        im.data[x + (im.w * y) + (im.w * im.h * c)] = v;
    }
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            for (int k = 0; k < im.c; k++) {
                set_pixel(copy, i, j, k, get_pixel(im, i, j, k));
            }
        }
    }
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert (im.c == 3); 
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    for (int i = 0; i < gray.w; i++) {
        for (int j = 0; j < gray.h; j++) {
            float value = 0.299 * get_pixel(im, i, j, 0) + 0.587 * get_pixel(im, i, j, 1) + 0.114 * get_pixel(im, i, j, 2); 
            set_pixel(gray, i, j, 0, value);
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            set_pixel(im, i, j, c, get_pixel(im, i, j, c) + v);
        }
    }
}

void clamp_image(image im)
{
    // TODO Fill this in
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            for (int k = 0; k < im.c; k++) {
                if (get_pixel(im, i, j, k) < 0) {
                    set_pixel(im, i, j, k, 0); 
                }
                else if (get_pixel(im, i, j, k) > 1) {
                    set_pixel(im, i, j, k, 1);
                }
            }
        }
    }

}

// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{  
    assert (im.c == 3); 
    float r, g, b; 
    float h, s, v; 
    float h_p, c; 
    
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            r = get_pixel(im, i, j, 0); 
            g = get_pixel(im, i, j, 1);
            b = get_pixel(im, i, j, 2);
            v = three_way_max(r, g, b);
            c = v - three_way_min(r, g, b); 
            s = (v == 0) ? 0.0: c / v; 
            if (c == 0) h_p = 0;
            else if (v == r) h_p = (g - b) / c;
            else if (v == g) h_p = ((b - r) / c) + 2.0;
            else if (v == b) h_p = ((r - g) / c) + 4.0; 
            
            h = (h_p < 0) ? (h_p / 6.0) + 1.0: h_p / 6.0;  
            set_pixel(im, i, j, 0, h);
            set_pixel(im, i, j, 1, s);
            set_pixel(im, i, j, 2, v);
        }
    }
}

// Formula revised from https://www.rapidtables.com/convert/color/hsv-to-rgb.html
void hsv_to_rgb(image im)
{   
    assert (im.c == 3); 
    float h, s, v; 
    float r, g, b; 
    float c, m, x; 
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            h = get_pixel(im, i, j ,0) * 6;     // h in [0, 6)
            s = get_pixel(im, i, j, 1);     // s in [0, 1)     
            v = get_pixel(im, i, j, 2);     // s in [0, 1)    
            c = v * s; 
            x = c*(1-fabs(fmod(h, 2)-1));
            m = v - c;  
            
            if (0 <= h && h <= 1) {
                r = c; 
                g = x;
                b = 0;
            }
            else if (1 < h && h <= 2) {
                r = x; 
                g = c;
                b = 0;
            }
            else if (2 < h && h <= 3) {
                r = 0; 
                g = c;
                b = x;
            }
            else if (3 < h && h <= 4) {
                r = 0; 
                g = x;
                b = c;
            }
            else if (4 < h && h <= 5) {
                r = x; 
                g = 0;
                b = c;
            }
            else if (5 < h && h <= 6) {
                r = c; 
                g = 0;
                b = x;
            }
            set_pixel(im, i, j, 0, r + m); 
            set_pixel(im, i, j, 1, g + m); 
            set_pixel(im, i, j, 2, b + m); 
        }
    }
}

void scale_image(image im, int c, float v) {
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            set_pixel(im, i, j, c, get_pixel(im, i, j, c) * v);
        }
    }
}

void rgb_to_hcl(image im) {
    // 0. Formulas from https://www.image-engineering.de/library/technotes/958-how-to-convert-between-srgb-and-ciexyz
    // and https://observablehq.com/@mbostock/luv-and-hcl
    assert (im.c == 3); 
    image lrgb = make_image(im.w, im.h, im.c); 
    image xyz = make_image(im.w, im.h, im.c); 
    image luv = make_image(im.w, im.h, im.c); 
    
    // 1. sRGB -> Linear RGB (Gamma decomposition)
    for (int i = 0 ; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            for (int k = 0; k < im.c; k++) {
                float value = get_pixel(im, i, j, k); 
                value = (value <= 0.04045) ? value / 12.92 : powf((value + 0.055) / 1.055, 2.4); 
                set_pixel(lrgb, i, j, k, value); 
            }
        }
    }

    // 2. lRGB -> XYZ 
    float scale[3][3] = {{0.4124564, 0.3575761, 0.1804375},
                        {0.2126729, 0.7151522, 0.0721750},
                        {0.0193339, 0.1191920, 0.9503041}};
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            float lr = get_pixel(lrgb, i, j, 0); 
            float lg = get_pixel(lrgb, i, j, 1); 
            float lb = get_pixel(lrgb, i, j, 2);
            set_pixel(xyz, i, j, 0, lr * scale[0][0] + lg * scale[0][1] + lb * scale[0][2]); 
            set_pixel(xyz, i, j, 1, lr * scale[1][0] + lg * scale[1][1] + lb * scale[1][2]); 
            set_pixel(xyz, i, j, 2, lr * scale[2][0] + lg * scale[2][1] + lb * scale[2][2]); 
        }
    }

    // 3. XYZ -> CIELUV
    float un = 0.2009, vn = 0.4610, yn = 1.0; // D65 reference white, value from https://en.wikipedia.org/wiki/CIELUV
    for (int i = 0; i < xyz.w; i++) {
        for (int j = 0; j < xyz.h; j++) {
            float x = get_pixel(xyz, i, j, 0); 
            float y = get_pixel(xyz, i, j, 1); 
            float z = get_pixel(xyz, i, j, 2);
            float l_star = ((y/yn) <= powf(6.0/29, 3)) ? powf((29.0/3), 3) * y / yn : 116 * powf(y/yn, 1.0/3) - 16;
            float u_p = 4*x/(x + 15*y + 3*z); 
            float v_p = 9*y/(x + 15*y + 3*z);
            float u_star = 13 * l_star * (u_p - un); 
            float v_star = 13 * l_star * (v_p - vn); 
            set_pixel(luv, i, j, 0, l_star); 
            set_pixel(luv, i, j, 1, u_star); 
            set_pixel(luv, i, j, 2, v_star); 
        }
    }

    // 4. CIELUV -> HCL
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            float l = get_pixel(luv, i, j ,0); 
            float u = get_pixel(luv, i, j, 1); 
            float v = get_pixel(luv, i, j, 2);
            float h = atan2f(v, u);
            float c = sqrtf(u * u + v * v);
            set_pixel(im, i, j, 0, h);
            set_pixel(im, i, j, 1, c);
            set_pixel(im, i, j, 2, l);  
        }
    }
    
    // Free intermidiate images 
    free_image(lrgb);
    free_image(xyz); 
    free_image(luv); 
}

void hcl_to_rgb(image im) { 
    image luv = make_image(im.w, im.h, im.c);
    image xyz = make_image(im.w, im.h, im.c); 
    image lrgb = make_image(im.w, im.h, im.c);

    // 1. HCL -> CIEWLUV
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            float h = get_pixel(im, i, j, 0);
            float c = get_pixel(im, i, j, 1);
            float l = get_pixel(im, i, j, 2);
            float u = c * cosf(h);
            float v = c * sinf(h); 
            set_pixel(luv, i, j, 0, l); 
            set_pixel(luv, i, j, 1, u); 
            set_pixel(luv, i, j, 2, v);
        }
    }

    // 2. LUV -> xyz
    float un = 0.2009, vn = 0.4610, yn = 1.0; // D65 reference white, value from https://en.wikipedia.org/wiki/CIELUV
    for (int i = 0; i < luv.w; i++) {
        for (int j = 0; j < luv.h; j++) {
            float l = get_pixel(luv, i, j, 0);
            float u = get_pixel(luv, i, j, 1); 
            float v = get_pixel(luv, i, j, 2);
            float y = (l <= 8) ? yn * l * powf(3.0/29, 3) : yn * powf(((l + 16) / 116.0), 3); 
            float u_p = (l == 0) ? -un : u/(13*l) + un; 
            float v_p = (l == 0) ? -vn : v/(13*l) + vn;
            float x = y * (9.0 * u_p) / (4.0 * v_p); 
            float z = y * (12.0 - 3.0 * u_p - 20.0 * v_p)/(4.0 * v_p); 
            set_pixel(xyz, i, j, 0, x); 
            set_pixel(xyz, i, j, 1, y); 
            set_pixel(xyz, i, j, 2, z);
        }
    }

    // 3. xyz -> lrgb
    float transform[3][3] = {{3.2404541,    -1.5371385,    -0.4985314}, 
                             {-0.9692660,    1.8760108,     0.0415560}, 
                             {0.0556434,    -0.2040257,     1.0572252}}; 
    for (int i = 0; i < xyz.w; i++) {
        for (int j = 0; j < xyz.h; j++) {
            float x = get_pixel(xyz, i, j, 0); 
            float y = get_pixel(xyz, i, j, 1);
            float z = get_pixel(xyz, i, j, 2);
            set_pixel(lrgb, i, j, 0, x * transform[0][0] + y * transform[0][1] + z * transform[0][2]); 
            set_pixel(lrgb, i, j, 1, x * transform[1][0] + y * transform[1][1] + z * transform[1][2]); 
            set_pixel(lrgb, i, j, 2, x * transform[2][0] + y * transform[2][1] + z * transform[2][2]); 
        }
    }

    // 4. lrgb -> srgb
    for (int i = 0; i < xyz.w; i++) {
        for (int j = 0; j < xyz.h; j++) {
            for (int k = 0; k < xyz.c; k++) {
                float value = get_pixel(lrgb, i, j, k); 
                value = (value <= 0.0031308) ? value * 12.92 : 1.055 * powf(value, 1/2.4) - 0.055;
                set_pixel(im, i, j, k, value); 
            }
        }
    } 

    // Free images
    free_image(lrgb);
    free_image(xyz);
    free_image(luv);
}

