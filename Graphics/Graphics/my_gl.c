
#ifdef WIN32
#include <windows.h>
#endif

#include <GL/gl.h>

// include declaration ofI_my_gl functions
#include "i_my_gl.h"
#include "inputModule.h"
// the followings are wrapper functions

// synchronize OpenGL's current matrix with i_my_gl's current matrix
void sync_matrix()
{
    GLdouble m[16];
    I_my_glGetMatrixd(m);
    glLoadMatrixd(m);
}

// switching matrix mode for both OpenGL and i_my_gl
void my_glMatrixMode(GLenum mode)
{
	if (toggle) {
		I_my_glMatrixMode(mode);
		glMatrixMode(mode);
		sync_matrix();
	}
	else {
		glMatrixMode(mode);
	}
}

// all following functions first call corresponding i_my_gl functions and then synchronize current matrix with OpenGL

void my_glLoadIdentity(void)
{
	if (toggle) {
		I_my_glLoadIdentity();
		sync_matrix();
	}
	else {
		glLoadIdentity();
	}
}

void my_glPushMatrix(void)
{
	if (toggle) {
		I_my_glPushMatrix();
		sync_matrix();
	}
	else {
		glPushMatrix();
	}
}

void my_glPopMatrix(void)
{
    if (toggle) {
		I_my_glPopMatrix();
		sync_matrix();
	}
	else {
		glPopMatrix();
	}
}

void my_glLoadMatrixf(const GLfloat *m)
{
    if (toggle) {
		I_my_glLoadMatrixf(m);
		sync_matrix();
	}
	else {
		glLoadMatrixf(m);
	}
}

void my_glLoadMatrixd(const GLdouble *m)
{
    if (toggle) {
		I_my_glLoadMatrixd(m);
		sync_matrix();    
	}
	else {
		glLoadMatrixd(m);
	}
}

void my_glTranslated(GLdouble x, GLdouble y, GLdouble z)
{
    if (toggle) {
		I_my_glTranslated(x, y, z);
		sync_matrix();
	}
	else 
		glTranslated(x, y, z);
}

void my_glTranslatef(GLfloat x, GLfloat y, GLfloat z)
{
    if (toggle) {
		I_my_glTranslatef(x, y, z);
		sync_matrix();
	}
	else 
		glTranslatef(x, y, z);
}

void my_glRotated(GLdouble angle, GLdouble x, GLdouble y, GLdouble z)
{
    if (toggle) {
		I_my_glRotated(angle, x, y, z);
		sync_matrix();
	}
	else 
		glRotated(angle, x, y, z);
}

void my_glRotatef(GLfloat angle, GLfloat x, GLfloat y, GLfloat z)
{
    if (toggle) {
		I_my_glRotatef(angle, x, y, z);
		sync_matrix();
	}
	else 
		glRotatef(angle, x, y, z);
}

void my_glScaled(GLdouble x, GLdouble y, GLdouble z)
{
    if (toggle) {
		I_my_glScaled(x, y, z);
		sync_matrix();
	}
	else 
		glScaled(x, y, z);
}

void my_glScalef(GLfloat x, GLfloat y, GLfloat z)
{
    if (toggle) {
		I_my_glScalef(x, y, z);
		sync_matrix();
	}
	else 
		glScalef(x, y, z);
}

void my_gluLookAt(GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ,
    GLdouble centerX, GLdouble centerY, GLdouble centerZ,
    GLdouble upX, GLdouble upY, GLdouble upZ)
{
    if (toggle) {
		I_my_gluLookAt(eyeX, eyeY, eyeZ, 
	        centerX, centerY, centerZ,
			upZ, upY, upZ);
		sync_matrix();
	}
	else 
		gluLookAt(eyeX, eyeY, eyeZ, 
			centerX, centerY, centerZ,
			upX, upY, upZ);
}

void my_glFrustum(GLdouble left, GLdouble right, 
    GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar)
{
    if (toggle) {
		I_my_glFrustum(left, right, bottom, top, zNear, zFar);
		sync_matrix();
	}
	else 
		glFrustum(left, right, bottom, top, zNear, zFar);
		sync_matrix();
}

void my_gluPerspective(GLdouble fovy, GLdouble aspect,
    GLdouble zNear, GLdouble zFar)
{
    if (toggle) {
		I_my_gluPerspective(fovy, aspect, zNear, zFar);
		sync_matrix();
	}
	else 
		gluPerspective(fovy, aspect, zNear, zFar);
}



