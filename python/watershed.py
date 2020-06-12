import numpy as np
import cv2 as cv
from matplotlib import pyplot as plt
img = cv.imread('/home/jonathan/dlt2.0_cgi2016/texturesamples/PhotosculptTextures/photosculpt-pebbles-diffuse.jpg')
gray = cv.cvtColor(img,cv.COLOR_BGR2GRAY)

ret, thresh = cv.threshold(gray,0,255,cv.THRESH_BINARY_INV+cv.THRESH_OTSU)
cv.imwrite('watershed.png',thresh)
kernel = np.ones((3,3),np.uint8)
iteration = 40


for i in range(1, iteration):
    print(iteration-i)
    #cv.threshold(	src, thresh, maxval, type[, dst]	)
    erosion = cv.erode(thresh,kernel,iterations=i)
    cv.imwrite('watershed'+str(iteration-i)+'.png',erosion)

for i in range(0, iteration):
  print(i+iteration)
  #cv.threshold(	src, thresh, maxval, type[, dst]	)
  dilatation = cv.dilate(thresh,kernel,iterations=i)
  cv.imwrite('watershed'+str(i+iteration)+'.png',dilatation)
