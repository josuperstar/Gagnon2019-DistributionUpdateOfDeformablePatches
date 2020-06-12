import cv2
import numpy as np
from matplotlib import pyplot as plt

img = cv2.imread('/home/jonathan/scenes/scenes/dambreak/render/render.v001.0120.png',0)
f = np.fft.fft2(img)
fshift = np.fft.fftshift(f)
magnitude_spectrum = 20*np.log(np.abs(fshift))

imgSample = cv2.imread("/home/jonathan/dlt2.0_cgi2016/UV.png",0)
fSample = np.fft.fft2(imgSample)
fshiftSample = np.fft.fftshift(fSample)
magnitude_spectrum_sample = 20*np.log(np.abs(fshiftSample))
cv2.imwrite('spectrumsample.png',magnitude_spectrum_sample)

i=1
cv2.imwrite('spectrum'+str(i)+'.png',magnitude_spectrum)

plt.subplot(121),plt.imshow(magnitude_spectrum_sample, cmap = 'gray')
plt.title('Input Image'), plt.xticks([]), plt.yticks([])
plt.subplot(122),plt.imshow(magnitude_spectrum, cmap = 'gray')
plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])
plt.show()
