//
// (c) MERL 2012-2013
//
#include "GlutSVR.h"
#include <GL/glut.h>
#include <opencv/highgui.h>
#include <opencv/cv.h>
#include <opencv/cvaux.h>
#include <Windows.h>
#include <stdio.h>
#include <stdlib.h>
#include "lsd.h"
#include "geometry.h"
#include "vp.h"
#include "rotation.h"
#include <math.h>
#include <opencv2/core/core_c.h>
#include <opencv2/core/types_c.h>
#include <opencv2/imgproc/imgproc_c.h>
#include <opencv2/imgproc.hpp>

#include <opencv2/core/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui/highgui.hpp>

#define	INPUT_DIR_NAME	"C:/Research/SceneUnderstanding/vp_detector_rotation/Input"
#define OUTPUT_DIR_NAME "C:/Research/SceneUnderstanding/vp_detector_rotation/Output"
#define CAMERA_MATRIX_FILE "C:/Research/SceneUnderstanding/vp_detector_rotation/CameraMatrixGoPro.txt"

#define	INPUT_IMAGE_EXT	".jpg"
#define MIN_LNLEN 40
#define MIN_LNLEN_VP 50
#define J_SUBSAMP 5

GlutSVR::GlutSVR()
	: IGlut(),
	  texId_(0), texWidth_(0), texHeight_(0),
	  targetImageIdx_(0)
{
	InitTexture();
	FindImageFiles(INPUT_DIR_NAME);
	LoadTargetImage();
}

GlutSVR::~GlutSVR()
{
	if (texId_)
	{
		glDeleteTextures(1, &texId_);
	}
}


void
GlutSVR::Display()
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	if (displayImage_.data)
	{
		const double winW = (double)glutGet(GLUT_WINDOW_WIDTH);
		const double winH = (double)glutGet(GLUT_WINDOW_HEIGHT);

		const double imgW = displayImage_.cols;
		const double imgH = displayImage_.rows;

		const double scaleW = winW / imgW;
		const double scaleH = winH / imgH;
		const double scale = std::min(std::min(scaleW, scaleH), 1.0);

		const double dispW = imgW * scale;
		const double dispH = imgH * scale;

		const double texW = texWidth_;
		const double texH = texHeight_;

		glViewport(0, 0, dispW, dispH);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(0.0, 1.0, 0.0, 1.0);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glColor4f(1.0, 1.0, 1.0, 1.0);

		glEnable(GL_TEXTURE_2D);
		glBegin(GL_QUADS);
		{
			glTexCoord2d(0.0, 0.0);
			glVertex2d(0.0, 1.0);

			glTexCoord2d(imgW/texW, 0.0);
			glVertex2d(1.0, 1.0);

			glTexCoord2d(imgW/texW, imgH/texH);
			glVertex2d(1.0, 0.0);

			glTexCoord2d(0.0, imgH/texH);
			glVertex2d(0.0, 0.0);
		}
		glEnd();
		glDisable(GL_TEXTURE_2D);
	}

	glutSwapBuffers();
}



void
GlutSVR::Keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'n':
		{
			targetImageIdx_ = std::min(targetImageIdx_ + 1, (int)imageFileNameVec_.size()-1);
			LoadTargetImage();
			ProcessTargetImage();
		}
		break;

	case 'p':
		{
			targetImageIdx_ = std::max(targetImageIdx_ - 1, 0);
			LoadTargetImage();
		}
		break;

	case 'a':
		{
			ProcessTargetImage();
		}
		break;

	case 'b':
		{
			targetImageIdx_ = 0;//std::min(targetImageIdx_ + 1, (int)imageFileNameVec_.size()-1);
			while(targetImageIdx_<imageFileNameVec_.size())
			{
				LoadTargetImage();
				ProcessTargetImage();
				targetImageIdx_ = targetImageIdx_ + 1;

			}
		}
		break;

	case 'g':
		{
			printf("Enter the image index:\n");
			scanf("%d",&targetImageIdx_);
			LoadTargetImage();
		}
		break;
	default:
		break;
	}

	glutPostRedisplay();
}




void
GlutSVR::InitTexture()
{
	texWidth_ = texHeight_ = 2048;

	glGenTextures(1, &texId_);

	glEnable(GL_TEXTURE_2D);
	{
		glBindTexture(GL_TEXTURE_2D, texId_);
		glTexImage2D(GL_TEXTURE_2D, 0, 4,
					 texWidth_, texHeight_,
					 0, GL_RGBA, GL_UNSIGNED_BYTE,
					 NULL);

		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE/*GL_MODULATE*/);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	}
	glDisable(GL_TEXTURE_2D);
}

void
GlutSVR::SetDisplayImage(const cv::Mat& image)
{
	image.copyTo(displayImage_);
	glEnable(GL_TEXTURE_2D);
	{
		glBindTexture(GL_TEXTURE_2D, texId_);
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0,
						displayImage_.cols, displayImage_.rows,
						GL_RGB, GL_UNSIGNED_BYTE,
						displayImage_.data);
	}
	glDisable(GL_TEXTURE_2D);
}



void
GlutSVR::FindImageFiles(const char* dirName)
{
	printf("Finding image files in %s:\n", dirName);

	WIN32_FIND_DATAA ffd;
	HANDLE hFind = INVALID_HANDLE_VALUE;
  
	std::string searchStr(dirName);
	searchStr.append("\\*");	// append "\*" to the directory name

	// find the first file in the directory
	hFind = FindFirstFileA(searchStr.c_str(), &ffd);

	if (INVALID_HANDLE_VALUE == hFind) 
	{
		printf("Error in finding image files\n");
		exit(1);
	}

	// list all the files in the directory with some info about them
	do
	{
		if (ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
		{
			printf("%s\t<DIR>\n", ffd.cFileName);
		}
		else
		{
			LARGE_INTEGER filesize;
			filesize.LowPart = ffd.nFileSizeLow;
			filesize.HighPart = ffd.nFileSizeHigh;
			printf("%s\t%ld bytes\n", ffd.cFileName, filesize.QuadPart);

			if (strstr(ffd.cFileName, INPUT_IMAGE_EXT))
			{
				imageFileNameVec_.push_back(std::string(ffd.cFileName));
			}
		}
	}
	while (FindNextFileA(hFind, &ffd) != 0);

	DWORD dwError = GetLastError();
	if (dwError != ERROR_NO_MORE_FILES) 
	{
		printf("Error in finding image files\n");
		exit(1);
	}

	FindClose(hFind);

	printf("Total of %d image files\n", imageFileNameVec_.size());
}




void
GlutSVR::LoadTargetImage()
{
	const std::string fileName = std::string(INPUT_DIR_NAME) + std::string("/") + imageFileNameVec_[targetImageIdx_].c_str();

	printf("Processing %s\n", fileName.c_str());
	targetImage_ = cv::imread(fileName.c_str(), 1);	// load as BGR	

	printf("no of channels = %d\n",targetImage_.channels());
	if (targetImage_.channels() == 1)
		cvtColor(targetImage_, targetImage_, CV_GRAY2RGB);
	else
		cvtColor(targetImage_, targetImage_, CV_BGR2RGB);

  SetDisplayImage(targetImage_);
}


Lines 
GlutSVR::extractRGBLines()
{
	const int w = targetImage_.cols;
	const int h = targetImage_.rows;

	double *image;
	Lines linesR;
	geometry geometry_;
	int noLinesR;

	image = (double *) malloc( w * h * sizeof(double) );
	
	/* Extract lines from the red component of the image */
	for(int i=0;i<(w*h);i++)
		image[i] = (double)targetImage_.data[3*(i)];
  
	/*LineSegeometry_entDetection( int * n_out,
                              double * img, int X, int Y,
                              double scale, double sigma_scale, double quant,
                              double ang_th, double log_eps, double density_th,
                              int n_bins,
                              int ** reg_img, int * reg_x, int * reg_y )*/
	linesR.lines = LineSegmentDetection(&noLinesR,image,w,h,0.8,0.6,2,22,0,0.7,1024,NULL,NULL,NULL);
	linesR.noLines=noLinesR;
	return linesR;
}



void GlutSVR::ProcessTargetImage()
{
	geometry geometry_;
	Lines linesRGB;
	vp vp_;
	cv::Mat displayImg;
	const int noCandidates=30;	
	const int w = targetImage_.cols, BD_SIZE = 5;
	const int h = targetImage_.rows;
	double lineLength;				
	int x,y,noLines;
	double *imgLines;
	Points boxep,ep;	

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* All algorithmic processing:Line Detection, vanishing point detection, line clustering, junction detection and best layout estimation.*/
	{
		linesRGB=extractRGBLines();		
		// discard short and nearby lines - 
		// while getting vp we use only long and important lines.
		// if there are two lines that are collinear we use only the longest line. 		

		noLines=linesRGB.noLines;
		imgLines=new double[7*linesRGB.noLines];

		noLines=0;
		int *validLines;
		validLines=geometry_.removeNearbyLines(linesRGB.lines,linesRGB.noLines);
		for(int i=0;i<linesRGB.noLines;i++)
		{
			lineLength=sqrt((linesRGB.lines[7*i+0]-linesRGB.lines[7*i+2])*(linesRGB.lines[7*i+0]-linesRGB.lines[7*i+2])+
						(linesRGB.lines[7*i+1]-linesRGB.lines[7*i+3])*(linesRGB.lines[7*i+1]-linesRGB.lines[7*i+3]));
			if ((lineLength>MIN_LNLEN_VP)	&&	(validLines[i]==1))
			{
				for(int j=0;j<7;j++)
					imgLines[noLines*7+j]=linesRGB.lines[7*i+j];

				noLines++;
			}
		}		
		targetImage_.copyTo(displayImg);		
		//targetImage_.copyTo(outputImg);		

		vp_=vp_.getVP(imgLines,noLines,w,h);
		
		imgLines=new double[7*linesRGB.noLines];
		
		// we discard short lines
		noLines=0;
		for(int i=0;i<linesRGB.noLines;i++)
		{
			lineLength=sqrt((linesRGB.lines[7*i+0]-linesRGB.lines[7*i+2])*(linesRGB.lines[7*i+0]-linesRGB.lines[7*i+2])+
						(linesRGB.lines[7*i+1]-linesRGB.lines[7*i+3])*(linesRGB.lines[7*i+1]-linesRGB.lines[7*i+3]));
			if ( (lineLength>MIN_LNLEN)	&&	
				 (linesRGB.lines[7*i+0] > BD_SIZE) && (linesRGB.lines[7*i+0] < w-BD_SIZE) &&
				 (linesRGB.lines[7*i+1] > BD_SIZE) && (linesRGB.lines[7*i+1] < h-BD_SIZE) &&
				 (linesRGB.lines[7*i+2] > BD_SIZE) && (linesRGB.lines[7*i+2] < w-BD_SIZE) &&
				 (linesRGB.lines[7*i+3] > BD_SIZE) && (linesRGB.lines[7*i+3] < h-BD_SIZE)  )
			{
				for(int j=0;j<7;j++)
					imgLines[noLines*7+j]=linesRGB.lines[7*i+j];

				noLines++;
			}
		}

		delete vp_.p;
		vp_.p=vp_.getVPLineProbabilities(imgLines,noLines,vp_.vps);
		vp_.lines=imgLines;//allLines.lines;
		vp_.noLines=noLines;//allLines.noLines;		
		
	}	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/* All display operations: left image- original image with lines marked with colors corresponding to the vanishing points, top scoring corners marked */
	{	
		const int winH = (int)glutGet(GLUT_WINDOW_HEIGHT); // There is an offset between OpenGL and OpenCV
		const int offset = winH - (h); // This offset is used to get the mouse position with respect to the image.
		const static int lineThickness = 2;
		const static int ptRadius = 2;
		const static cv::Scalar colorR = cv::Scalar(255, 0, 0);  
		const static cv::Scalar colorG = cv::Scalar(0, 255, 0);  
		const static cv::Scalar colorB = cv::Scalar(0, 0, 255);  
		const static cv::Scalar colorRB = cv::Scalar(255, 0, 255);
		const static cv::Scalar colorBL = cv::Scalar(0, 0, 0);
		const static cv::Scalar colorW = cv::Scalar(255, 255, 255);

		for(int i=0;i<noLines;i++)
		{
			const cv::Point s((int)(imgLines[7*i+0]+0.5),(int)(imgLines[7*i+1]+0.5));
			const cv::Point e((int)(imgLines[7*i+2]+0.5),(int)(imgLines[7*i+3]+0.5));			
			if ((vp_.p[4*i+0]>vp_.p[4*i+1])&&(vp_.p[4*i+0]>vp_.p[4*i+2])&&(vp_.p[4*i+0]>vp_.p[4*i+3]))
				cv::line(displayImg, s, e, colorR, lineThickness);
			if ((vp_.p[4*i+1]>vp_.p[4*i+0])&&(vp_.p[4*i+1]>vp_.p[4*i+2])&&(vp_.p[4*i+1]>vp_.p[4*i+3]))
				cv::line(displayImg, s, e, colorG, lineThickness);
			if ((vp_.p[4*i+2]>vp_.p[4*i+0])&&(vp_.p[4*i+2]>vp_.p[4*i+1])&&(vp_.p[4*i+2]>vp_.p[4*i+3]))
				cv::line(displayImg, s, e, colorB, lineThickness);
			if ((vp_.p[4*i+3]>vp_.p[4*i+0])&&(vp_.p[4*i+3]>vp_.p[4*i+1])&&(vp_.p[4*i+3]>vp_.p[4*i+2]))				
				cv::line(displayImg, s, e, colorRB, lineThickness);
		}
		/* Displaying the vanishing point for Z */
		const cv::Point svp(vp_.vps[4],vp_.vps[5]);
				cv::circle(displayImg, svp, 8, colorW, -1, 8);				
				cv::circle(displayImg, svp, 4, colorR, -1, 8);

		cv::Mat newDispImage(h, 2*w, CV_8UC3);
		/* On the left hand side we show the original image with lines clustered and high scoring corners */
		displayImg.copyTo(		newDispImage(cv::Range(0, h), cv::Range(0, w)) );
				
		SetDisplayImage(newDispImage);	
	}
	Eigen::Matrix3d Kmatrix;
	FILE *fp;
	fp = fopen(CAMERA_MATRIX_FILE, "r");
	for (int i = 0; i<3; i++)
		fscanf(fp, "%lf,%lf,%lf\n", &Kmatrix(i, 0), &Kmatrix(i, 1), &Kmatrix(i, 2));
	fclose(fp);

	rotation rotation_;
	rotation_.computeRotationFromVP(vp_, Kmatrix);



	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* Saving the final outputImg: image with box marked in red, txt file storing all four corners of the box */
	{
		const std::string fileName = std::string(OUTPUT_DIR_NAME) + std::string("/") + imageFileNameVec_[targetImageIdx_].c_str();		
		cvtColor(displayImg, displayImg, CV_BGR2RGB);
		cv::imwrite(fileName.c_str(),displayImg);

		/* writing the txt file storing all four corners*/
		char fileNameTXT[200];
		strcpy(fileNameTXT,OUTPUT_DIR_NAME);
		strcat(fileNameTXT,"/");
		strcat(fileNameTXT,imageFileNameVec_[targetImageIdx_].c_str());
		int filenameLength=strlen(fileNameTXT);
		printf("%s,%d\n", fileNameTXT,filenameLength);
		fileNameTXT[filenameLength-4]='\0';
		strcat(fileNameTXT,".txt");
		FILE *fp;
		fp=fopen(fileNameTXT,"w");
		fprintf(fp,"%lf,%lf,%lf,%lf,%lf,%lf\n",vp_.vps[0],vp_.vps[1],vp_.vps[2],vp_.vps[3],vp_.vps[4],vp_.vps[5]);

		for (int i = 0; i<3; i++)
			fprintf(fp, "%lf,%lf,%lf\n", rotation_.rot(i, 0), rotation_.rot(i, 1), rotation_.rot(i, 2));

		fprintf(fp,"%d\n",vp_.noLines);
		for(int i=0;i<vp_.noLines;i++)
			fprintf(fp,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",vp_.p[4*i+0],vp_.p[4*i+1],vp_.p[4*i+2],vp_.p[4*i+3],imgLines[7*i+0],imgLines[7*i+1],imgLines[7*i+2],imgLines[7*i+3]);
		fclose(fp);
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* Deleting all the pointer variables */
	{
		delete imgLines;
		delete linesRGB.lines,vp_.lines,vp_.p;
	}
}



