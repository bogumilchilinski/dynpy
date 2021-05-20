import cv2 as cv
import numpy as np 
import matplotlib.pyplot as plt 
import imutils as imutils

## ------------------------------------------------------------------Exercise 1
img_bgr = np.zeros((2, 2, 3), np.uint8)                             # Creating 2x2 matrix with 3 color channels
img_bgr[0][0] = np.array([255, 0, 0], np.uint8)                     # Creating pixel composed of red
img_bgr[0][1] = np.array([0, 255, 0], np.uint8)                     # Creating pixel composed of green
img_bgr[1][0] = np.array([255*0.7, 255*0.7, 255*0.7], np.uint8)     # Creating pixeld composed of 70% of gray
img_bgr[1][1] = np.array([0, 0, 255], np.uint8)                     # Creating pixel composed of blue

plt.figure()                                                        # Plotting
plt.title('RGB image')
plt.imshow(img_bgr, cmap="gray")                                    # Showing the result image
plt.show()


## ------------------------------------------------------------------Exercise 2
img = cv.imread('Gray21.tif')                                       # Reading the source image
plt.figure(figsize=[25,25])                                         # Creating figure with subplots
plt.subplot(1,3,1)
plt.title('Monochromic image')
imgplot = plt.imshow(img, cmap="gray") 

thereshold = np.uint(0.3*255)                                       # Threshold = 30%                                     
ret, binarized = cv.threshold(img, thereshold, 255, cv.THRESH_BINARY)       # Converting grayscale image into binary one

plt.subplot(1,3,2)
plt.title('Binary image [0.30 threshold]')
imgplot = plt.imshow(binarized, cmap="gray")

thereshold = np.uint(0.7*255)                                       # Threshold = 70% 
ret, binarized = cv.threshold(img, thereshold, 255, cv.THRESH_BINARY)       # Converting grayscale image into binary one

plt.subplot(1,3,3)
plt.title('Binary image [0.70 threshold]')
imgplot = plt.imshow(binarized, cmap="gray")
plt.show()
"""
During the converting of grayscale image into binary one level of
treshold value has to be considered. The threshold value defines at which 
level of the brightness of the gray element would be interpreted as white 
in the case of the binary image.

The all elements with higher value of brigthness of grayscale compered to 
treshold level would be stored as white color elements. 

Summing up, the binary image with higher treshold value during converting
will result with less white-color elements.
"""


## ------------------------------------------------------------------Exercise 3
img = cv.cvtColor(cv.imread('Lena.tif'), cv.COLOR_BGR2RGB)          # Reading the source image

img_cyan = img.copy()                                               # Copy an image we have already loaded
img_cyan[:, :, (2,1)] = 255                                         # Lena in cyan and white layers

plt.figure(figsize=[12,12])                                         # Plotting
plt.subplot(1,2,1)
plt.title('Original image')
imgplot = plt.imshow(img) 
plt.subplot(1,2,2)
plt.title('Original image in cyan and white shades')
plt.imshow(img_cyan) 
plt.show()


## ------------------------------------------------------------------Exercise 4
img = cv.cvtColor(cv.imread('Lena.tif'), cv.COLOR_BGR2RGB)          # Reading the source image

img_rotate = imutils.rotate(img, 30)                                # Rotating the image by 30 degrees in a clockwise direction
rotated_bound_img = imutils.rotate_bound(img, 30)                   # Adjusting the image frage using imutils library

plt.figure(figsize=[12,12])                                         # Plotting
plt.subplot(1,2,1)
plt.title('Original image')
imgplot = plt.imshow(img) 
plt.subplot(1,2,2)
plt.title('Rotated image')
imgplot = plt.imshow(rotated_bound_img) 
plt.show()


## ------------------------------------------------------------------Exercise 5
img = np.zeros((4, 4, 3), np.uint8)                                 # Creating 4x4 matrix with 3 color channels
img[1::2, 0::2] = [255, 255, 255]                                   
img[0::2, 1::2] = [255, 255, 255]                                   # Assigning the white color

img_nearest = cv.resize(img, (8,8), interpolation = cv.INTER_NEAREST)   # 8x8 pixel configuration, nearest-neighbor method
img_linear = cv.resize(img, (8, 8), interpolation=cv.INTER_LINEAR)      # 8x8 pixel configuration, bilinear method

plt.figure(figsize=[12,12])                                         # Plotting
plt.subplot(1,3,1)
plt.title('Original image')
imgplot = plt.imshow(img, cmap="gray") 

plt.subplot(1,3,2)
plt.title('Nearest interpolation')
imgplot = plt.imshow(img_nearest, cmap="gray") 

plt.subplot(1,3,3)
plt.title('Bilinear interpolation')
imgplot = plt.imshow(img_linear, cmap="gray") 
plt.show()
"""
% The 8x8 pixel configure is produced due to scaling twice the original image.

In case of Nearest-neighbor interpolation, the output pixel is assigned 
the value of the pixel that the point falls within. No other pixels are 
considered. 

In case of Bilinear interpolation, new pixel accepts the value of the
linear interpolation of pixel value from its neighborhood. The output
pixel value is a weighted average of pixels in the nearest 2 by 2
neighborhood. Due to that we can see image deviates from original in a 
sense of a gray color.
"""


## ------------------------------------------------------------------Exercise 6
img = cv.imread('Lena_mono.tif')                                    # Reading the source image
height, width = img.shape[:2]                                       # Determine height, width parameters of the image

x_shift = 100                                                       # Image is right by 100 rows
y_shift = -200                                                      # Image is moved up by 200 rows

T = np.float32([ [1,0, np.abs(x_shift)], [0,1, np.abs(y_shift)] ])  # Transformation matrix
img_translated = cv.warpAffine(img, T, 
                 ( height + np.abs(x_shift), width + np.abs(y_shift) ))       # Shifting source image regarding transformation matrix

T = np.float32([ [1,0,  x_shift*(x_shift<0) ], [0,1, y_shift*(y_shift<0)] ])  # Transformation matrix in order to adjust the image frame       
img_translated = cv.warpAffine(img_translated, T, 
                  ( height + np.abs(x_shift) , width + np.abs(y_shift) ))     # Frame scaling

plt.figure(figsize=[12,12])                                         # Plotting
plt.subplot(1,2,1)
plt.title('Original image')
imgplot = plt.imshow(img) 
plt.subplot(1,2,2)
plt.title('Translated image')
imgplot = plt.imshow(img_translated)
plt.show()


## ------------------------------------------------------------------Exercise 7
img = cv.cvtColor(cv.imread('Lena.tif'), cv.COLOR_BGR2RGB)          # Reading the source image

kernel_sharp = np.array([[0, -1, 0], [-1, 5, -1], [0, -1, 0]])      # Creating filter mask
img_sharp = cv.filter2D(img, -1, kernel_sharp)                      # Image filtering

plt.figure(figsize=[12,12])                                         # Plotting
plt.subplot(1,2,1)
plt.title('Original image')
imgplot = plt.imshow(img) 
plt.subplot(1,2,2)
plt.title('Sharpened image')
imgplot = plt.imshow(img_sharp) 
plt.show()
"""
The sharpening mask is created by combining a neutral mask and high-pass
mask. As result, the filtered image is sharp with more readable, enhanced
edges.
"""


## ------------------------------------------------------------------Exercise 8
img = cv.imread('Auto.tif')                                         # Reading the source image

thereshold = np.uint(0.1*255)                                       # Threshold = 10%                             
ret, img = cv.threshold(img, thereshold, 255, cv.THRESH_BINARY)     # Converting grayscale image into binary one

height, width = img.shape[:2]                                       # Determine height, width parameters of the image                   

SE = np.ones((3, 3), np.uint8)                                      # Creating structural element

# Outer 
grad_outer = img - cv.erode(img, SE, iterations = 1)                # Determine external half-gradient
img_outer = cv.bitwise_not(grad_outer)                              # Determine outer contour

#Inner
grad_inner = cv.dilate(img, SE, iterations = 1)                     # Determine internal half-gradient
img_inner = img + cv.bitwise_not(grad_inner)                        # Determine inner contour

# ------------------------- Plotting
plt.figure(figsize=[25,25])
plt.subplot(3,1,1)
plt.title('Original image')
plt.imshow(img) 
plt.subplot(3,1,2)
plt.title('Outer countour')
plt.imshow(img_outer)
plt.subplot(3,1,3)
plt.title('Inner countour')
plt.imshow(img_inner) 
plt.show()

plt.figure(figsize=[15,15])
plt.subplot(1,2,1)
plt.title('Marked outer contour')
plt.imshow(img) 
for x in range(len(img)):
    for y in range(len(img[x])):
        if img_outer[x,y].all() == 0:
            plt.scatter(y,x, 5, 'red')          # Plotting the red dots

plt.subplot(1,2,2)
plt.title('Marked inner contour')
plt.imshow(img) 
for x in range(len(img)):
    for y in range(len(img[x])):
        if img_inner[x,y].all() == 0:
            plt.scatter(y,x, 5, 'red')          # Plotting the red dots
plt.show()
"""
Generally, the contour of the image can be divided into inner and outer
contours, using half-gradients. Using half-gradient we are able to mark
only the chosen contour and make further operations on that area. 

Note: Because outer contour lies outside the figure and outer 
half-gradient is perfomed using erosion function, negative operation was 
performed to create coherence visualization.
"""