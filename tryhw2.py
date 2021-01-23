from uwimg import *
im = load_image("/Users/fangchenxiang/Documents/Rocky/Courses/CSE_455/uwimg/data/dog.jpg")
f = make_box_filter(7)
blur = convolve_image(im, f, 1)
save_image(blur, "dog-box7")

im = load_image("/Users/fangchenxiang/Documents/Rocky/Courses/CSE_455/uwimg/data/dog.jpg")
f = make_box_filter(7)
blur = convolve_image(im, f, 1)
thumb = nn_resize(blur, blur.w//7, blur.h//7)
save_image(thumb, "dogthumb")

im = load_image("/Users/fangchenxiang/Documents/Rocky/Courses/CSE_455/uwimg/data/dog.jpg")
f = make_highpass_filter()
blur = convolve_image(im, f, 0)
clamp_image(blur)
save_image(blur, "dog-high_pass_0")

im = load_image("/Users/fangchenxiang/Documents/Rocky/Courses/CSE_455/uwimg/data/dog.jpg")
f = make_sharpen_filter()
blur = convolve_image(im, f, 0)
clamp_image(blur)
save_image(blur, "dog-sharpen_0")

im = load_image("/Users/fangchenxiang/Documents/Rocky/Courses/CSE_455/uwimg/data/dog.jpg")
f = make_emboss_filter()
blur = convolve_image(im, f, 0)
clamp_image(blur)
save_image(blur, "dog-emboss_0")

im = load_image("/Users/fangchenxiang/Documents/Rocky/Courses/CSE_455/uwimg/data/dog.jpg")
f = make_gaussian_filter(2)
blur = convolve_image(im, f, 1)
save_image(blur, "dog-gauss2")

f = make_gaussian_filter(2)
lfreq = convolve_image(im, f, 1)
hfreq = im - lfreq
reconstruct = lfreq + hfreq
save_image(lfreq, "low-frequency")
save_image(hfreq, "high-frequency")
save_image(reconstruct, "reconstruct")

res = sobel_image(im)
mag = res[0]
ang = res[1]
feature_normalize(mag)
feature_normalize(ang)
save_image(mag, "magnitude")
save_image(ang, "angle")

im = load_image("/Users/fangchenxiang/Documents/Rocky/Courses/CSE_455/uwimg/data/dog.jpg")
cool = colorize_sobel(im)
feature_normalize(cool)
gussian_filter = make_gaussian_filter(1)
cool = convolve_image(cool, gussian_filter, 1)
save_image(cool, "colorize_sobel")