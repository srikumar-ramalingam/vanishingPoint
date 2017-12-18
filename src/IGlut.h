//
// (c) MERL 2012
//
#pragma once

// interface for glut functions
class IGlut
{
protected:
	enum MOUSE_STATE { LEFT, RIGHT, CENTER };
	MOUSE_STATE m_MouseState;
	int m_MouseX, m_MouseY;

public:
	IGlut() {};
	virtual ~IGlut() {};
	virtual void Display() = 0;
	virtual void Idle() {};
	virtual void Reshape(int width, int height) {};
	virtual void Mouse(int button, int state, int x, int y) {};
	virtual void Motion(int x, int y) {};
	virtual void Keyboard(unsigned char key, int x, int y) {};
	virtual void Special(int key, int x, int y) {};
	virtual bool needExit() { return false; }
};

