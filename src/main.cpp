//
//
#include "GlutSVR.h"
#include <GL/glut.h>


//#ifndef	_DEBUG

//#pragma comment(lib, "opencv_core320d.lib")
//#pragma comment(lib, "opencv_highgui320.lib")
//#pragma comment(lib, "opencv_imgproc320.lib")

//#else

//#pragma comment(lib, "opencv_core249d.lib")
//#pragma comment(lib, "opencv_highgui249d.lib")
//#pragma comment(lib, "opencv_imgproc249d.lib")

//#endif

static IGlut* s_IGlutClass = 0;

void display()
{
	s_IGlutClass->Display();
}

void idle()
{
	s_IGlutClass->Idle();
}

void reshape(int width, int height)
{
	s_IGlutClass->Reshape(width, height);
}

void keyboard(unsigned char key, int x, int y)
{
	// global keyboard hook
	switch (key)
	{
	case '\033':	// ESC
		{
			delete s_IGlutClass;
			exit(0);
		}
		break;

	default:
		s_IGlutClass->Keyboard(key, x, y);
		break;
	}
}

void mouse(int button, int state, int x, int y)
{
	s_IGlutClass->Mouse(button, state, x, y);
}

void motion(int x, int y)
{
	s_IGlutClass->Motion(x, y);
}

void special(int key, int x, int y)
{
	switch (key)
	{
	default:
		s_IGlutClass->Special(key, x, y);
		break;
	}
}


int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE | GLUT_ALPHA);
	glutInitWindowSize(960, 720);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("SVR");

	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutSpecialFunc(special);
	glutIdleFunc(idle);

	s_IGlutClass = new GlutSVR();

	glutMainLoop();
	return 0;
}

