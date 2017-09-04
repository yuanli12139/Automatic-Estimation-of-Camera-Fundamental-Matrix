# Automatic Estimation of Camera Fundamental Matrix
• Implemented corner detection and used normalized correlation coefficient for feature extraction and matching with stereo images. 

• Implemented the 7-point algorithm with Sampson Error Correction in MSAC for outlier rejection. Implemented the 8-point algorithm in Direct Linear Transformation to estimate the camera’s Fundamental matrix. 

• Implemented Levenberg-Marquardt Algorithm to minimize the reprojection error with Triangulation and Optimal Correction.

[//]: # (Image References)

[fig1]: ./results/f1.png
[fig2]: ./results/f2.png
[fig3]: ./results/f3.png
[fig4]: ./results/f4.png
[fig5]: ./results/f5.png
[fig6]: ./results/f6.png

## Original images

![alt text][fig1]

## Corner detection

![alt text][fig2]

## Brute-force feature matching

![alt text][fig3]

## MSAC outlier rejection

![alt text][fig4]

## Levenberg-Marquardt algorithm optimization

![alt text][fig5]

## Final result

![alt text][fig6]
