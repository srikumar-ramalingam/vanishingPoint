//
// (c) MERL 2012
//
#pragma once

#include "IGlut.h"
#include <opencv/cv.h>
#include <string>
#include <vector>
#include "geometry.h"

class GlutSVR : public IGlut
{
private:
	unsigned int texId_;
	int texWidth_, texHeight_;

	cv::Mat targetImage_;
	cv::Mat displayImage_;

	std::vector<std::string> imageFileNameVec_;
	int targetImageIdx_;

private:
	void InitTexture();
	void SetDisplayImage(const cv::Mat& image);
	void FindImageFiles(const char* dirName);
	void LoadTargetImage();
	void ProcessTargetImage();
	Lines extractRGBLines();

public:
	GlutSVR();

	virtual ~GlutSVR();
	virtual void Display();

	virtual void Idle() {};
	virtual void Reshape(int width, int height) {};
	virtual void Mouse(int button, int state, int x, int y){};
	virtual void Motion(int x, int y) {};
	virtual void Keyboard(unsigned char key, int x, int y);
	virtual void Special(int key, int x, int y) {};
};


