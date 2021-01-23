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