import cv2 as cv
import numpy as np 
import matplotlib.pyplot as plt 
import imutils as imutils



## ------------------------------------------------------------------Exercise 1
img_bgr = np.zeros((2, 2, 3), np.uint8)
img_bgr[0][0] = np.array([255, 0, 0], np.uint8)
img_bgr[0][1] = np.array([0, 255, 0], np.uint8)
img_bgr[1][0] = np.array([255*0.7, 255*0.7, 255*0.7], np.uint8)
img_bgr[1][1] = np.array([0, 0, 255], np.uint8)


plt.figure()
plt.title('RGB image')
plt.imshow(img_bgr, cmap="gray") 
plt.show()




## ------------------------------------------------------------------Exercise 2
img = cv.imread('Gray21.tif')
plt.figure()
plt.title('Monochromic image')
imgplot = plt.imshow(img, cmap="gray") 


thereshold = np.uint(0.3*255);
ret, binarized = cv.threshold(img, thereshold, 255, cv.THRESH_BINARY)
plt.figure()
plt.title('Binary image [0.30 threshold]')
imgplot = plt.imshow(binarized, cmap="gray")

thereshold = np.uint(0.7*255);
ret, binarized = cv.threshold(img, thereshold, 255, cv.THRESH_BINARY)
plt.figure()
plt.title('Binary image [0.70 threshold]')
imgplot = plt.imshow(binarized, cmap="gray")
plt.show()




## ------------------------------------------------------------------Exercise 3
img = cv.imread('Lena.tif')

img_cyan = img.copy()
img_cyan[:, :, (2,1)] = 255

plt.figure(figsize=[12,12])
plt.subplot(1,2,1)
plt.title('Original image')
imgplot = plt.imshow(img) 
plt.subplot(1,2,2)
plt.title('Original image in cyan')
plt.imshow(img_cyan) 
plt.show()
plt.savefig('Exercise3.png')
