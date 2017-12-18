Dependencies:

Eigen
Opencv3.2 or others

Input: 

* Place input images in the Input folder 
* Place the calibration matrix file in this folder
* Change the directory and image extension details in GlutSVR.cpp file

Output: Will have all the details about the vanishing points and rotation matrices
* will have both image files and txt files. 
* The text files will have the following format
vanishing points
rotation matrix (3x3)
# of lines
line probabilities, line end points 


How to run:
After compiling, just press "b". It will generate the vanishing points and rotation information for all the images. 

