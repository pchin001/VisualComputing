
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef WIN32
#include <windows.h>
#endif

#include <GL/gl.h>

// define maximum stack capacity
#define STACK_CAP 16
// define pi for converting angles
#define PI 3.14159265359

// structure for the matrix stack, top specifies current top position on the stack, initially zero (which means one matrix in the stack)
struct matrix_stack
{
    GLdouble m[STACK_CAP][16];
    int top;
};

// define a macro for retrieving current matrix from top of current stack
#define current_matrix (current_stack->m[current_stack->top])

// identity matrix constant
const GLdouble identity[16] =
{1, 0, 0, 0,
 0, 1, 0, 0,
 0, 0, 1, 0,
 0, 0, 0, 1};

// the model view matrix stack
struct matrix_stack model_view = {{{0}}, 0};
// the projection matrix stack
struct matrix_stack projection = {{{0}}, 0};
// the current stack pointer that specifies the matrix mode
struct matrix_stack *current_stack = &model_view;

// multiply current matrix with matrix b, put the result in current matrix
// current = current * b
void matrix_multiply(const GLdouble *b)
{
    // ...
	GLdouble* m; 
	GLdouble temp[16]; 
	int cRow,cCol,i;

	m = current_matrix; 
	//do column-wise matrix multiplication
	for (cCol = 0; cCol < 4; cCol++) {
		for (cRow = 0; cRow < 4; cRow++) 
			 temp[4*cRow + cCol] = m[cCol] * b[4*cRow] //conversion from [cRow][cCol] = [4*cRow + cCol]
								+ m[4 + cCol] * b[4*cRow + 1]
								+ m[8 + cCol] * b[4*cRow + 2]
								+ m[12 + cCol] * b[4*cRow + 3];
	}
	for (i = 0; i < 16; i++) 
		m[i] = temp[i];
}

// calculating cross product of b and c, put the result in a
// a = b x c
void cross_product(GLdouble *ax, GLdouble *ay, GLdouble *az,
    GLdouble bx, GLdouble by, GLdouble bz,
    GLdouble cx, GLdouble cy, GLdouble cz)
{
    // ...
	//bxc= 	|by bz| i - |bx bz|j -	|bx by|k
	//		|cy cz|		|cx cz|		|cx cy|		
	*ax = by*cz - bz*cy;
	*ay = bz*cx - bx*cz;
	*az = bx*cy - by*cx;
}

// normalize vector (x, y, z)
void normalize(GLdouble *x, GLdouble *y, GLdouble *z)
{
    // ...
	//grab the magnitude of the vector
	//divide the magnitude from each of the coords. --> new coords = the quotient
	GLdouble magnitude; 
	magnitude = sqrt((*x) * (*x) + 
					(*y) * (*y) + 
					(*z) * (*z));
	*x = *x / magnitude; 
	*y = *y / magnitude; 
	*z = *z / magnitude;
}

// switch matrix mode by changing the current stack pointer
void I_my_glMatrixMode(GLenum mode)
{
    // ...
	switch (mode) {
		case GL_MODELVIEW: 
			current_stack = &model_view; 
			break;
		case GL_PROJECTION: 
			current_stack = &projection; 
			break;
	}
}

// overwrite current matrix with identity matrix
void I_my_glLoadIdentity(void)
{
    // ...
	GLdouble* m = current_matrix; 
	int i; 
	for (i = 0; i < 16; i++)
		m[i] = identity[i];
	glLoadMatrixd(identity);
}



// copy current matrix to m
void I_my_glGetMatrixf(GLfloat *m)
{
    // ...
	GLdouble* copy; 
	int i; 
	copy = current_matrix; 

	for (i = 0; i < 16; i++)
		m[i] = (GLfloat)copy[i]; 
}

void I_my_glGetMatrixd(GLdouble *m)
{
    // ...
	GLdouble* copy; 
	int i; 
	copy = current_matrix; 

	for (i = 0; i < 16; i++)
		m[i] = copy[i]; 
}

// overwrite currentmatrix with m
void I_my_glLoadMatrixf(const GLfloat *m)
{
    // ...
	GLdouble* temp;
	int i; 
	temp = current_matrix; 

	for (i = 0; i < 16; i++)
		temp[i] = m[i]; 
}

void I_my_glLoadMatrixd(const GLdouble *m)
{
    // ...
	GLdouble* temp;
	int i; 
	temp = current_matrix; 

	for (i = 0; i < 16; i++)
		temp[i] = m[i];
}

// push current matrix onto current stack, report error if the stack is already full
void I_my_glPushMatrix(void)
{
    // ...
	int stackTop; 
	GLfloat newMatrix[16]; 

	stackTop = current_stack -> top;

	if (stackTop == 16) 
		printf("Error: stack full");
	else {
		I_my_glGetMatrixf(newMatrix);
		(current_stack -> top)++;
		I_my_glLoadMatrixf(newMatrix);
	}
	//glPushMatrix();
}

// pop current matrix from current stack, report error if the stack has only one matrix left
void I_my_glPopMatrix(void)
{
    // ...
	int stackTop; 

	stackTop = current_stack -> top; 
	if (stackTop == 0)
		printf("Error: stack empty"); 
	else 
		(current_stack -> top)--; 
	//glPopMatrix();
}

void I_my_glTranslated(GLdouble x, GLdouble y, GLdouble z)
{
    // ...
	//multiply current_matrix with the matrix:
	// {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, x, y, z, 1} (in column-major form) 
	const GLdouble b[16]= { 1, 0, 0, 0, 
							0, 1, 0, 0, 
							0, 0, 1, 0, 
							x, y, z, 1};
	matrix_multiply(b);
}

void I_my_glTranslatef(GLfloat x, GLfloat y, GLfloat z)
{
    I_my_glTranslated((GLdouble)x, (GLdouble)y, (GLdouble)z);
}

// remember to normalize vector (x, y, z), and to convert angle from degree to radius before calling sin and cos
void I_my_glRotated(GLdouble angle, GLdouble x, GLdouble y, GLdouble z)
{
    // ...
	//rotates the current matrix of angle degrees around the vector (x, y, z) 
	GLdouble c; 
	GLdouble s; 
	GLdouble b[16];
	GLdouble radians; 

	normalize(&x, &y, &z);
	radians = (PI * angle) / 180;
	c = cos(radians);
	s = sin(radians); 

	b[0] = x*x*(1-c)+c;
	b[1] = y*x*(1-c)+z*s;
	b[2] = x*z*(1-c)-y*s;
	b[3] = 0;  
	b[4] = x*y*(1-c)-z*s;
	b[5] = y*y*(1-c)+c;
	b[6] = y*z*(1-c)+x*s;
	b[7] = 0;
	b[8] = x*z*(1-c)+y*s; 
	b[9] = y*z*(1-c)-x*s; 
	b[10] = z*z*(1-c)+c;
	b[11] = 0;
	b[12] = 0; 
	b[13] = 0; 
	b[14] = 0;
	b[15] = 1;

	matrix_multiply(b);
}

void I_my_glRotatef(GLfloat angle, GLfloat x, GLfloat y, GLfloat z)
{
    I_my_glRotated((GLdouble)angle, (GLdouble)x, (GLdouble)y, (GLdouble)z);
}

void I_my_glScaled(GLdouble x, GLdouble y, GLdouble z)
{
    // ...
	//multiply the current matrix with the scale transformation matrix. 
	const GLdouble m[16] = { x, 0, 0, 0, 
							0, y, 0, 0, 
							0, 0, z, 0, 
							0, 0, 0, 1 }; 
	matrix_multiply(m);
}

void I_my_glScalef(GLfloat x, GLfloat y, GLfloat z)
{
    I_my_glScaled((GLdouble)x, (GLdouble)y, (GLdouble)z);
}

// remember to normalize vectors
void I_my_gluLookAt(GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ, 
    GLdouble centerX, GLdouble centerY, GLdouble centerZ, 
    GLdouble upX, GLdouble upY, GLdouble upZ)
{
    // ...
	GLdouble fx, fy, fz; 
	GLdouble sx, sy, sz; 
	GLdouble ux, uy, uz; 
	GLdouble m[16]; 

	normalize (&fx, &fy, &fz);
	normalize (&upX, &upY, &upZ); 
	
	//s = f x up
	cross_product(&sx, &sy, &sz, 
				fx, fy, fz, 
				upX, upY, upZ);
	//u = s x f
	cross_product(&ux, &uy, &uz, 
				sx, sy, sz, 
				fx, fy, fz);

	m[0] = sx; 
	m[1] = ux; 
	m[2] = -fx; 
	m[3] = 0; 
	m[4] = sy; 
	m[5] = uy; 
	m[6] = -fy; 
	m[7] = 0; 
	m[8] = sz; 
	m[9] = uz; 
	m[10] = -fz; 
	m[11] = 0; 
	m[12] = 0; 
	m[13] = 0;
	m[14] = 0; 
	m[15] = 1; 

	matrix_multiply(m);
	I_my_glTranslated(-eyeX, -eyeY, -eyeZ);
}

void I_my_glFrustum(GLdouble left, GLdouble right, GLdouble bottom,
    GLdouble top, GLdouble zNear, GLdouble zFar)
{
    // ...
	GLdouble projection[16] = {1, 0, 0, 0, 
								0, 1, 0, 0, 
								0, 0, -(zFar+zNear)/(zFar-zNear), -1, 
								0, 0, -2*zFar*zNear/(zFar-zNear), 0};
	GLdouble scale[16] = {2*zNear/(right-left), 0, 0, 0, 
						0, 2*zNear/(top-bottom), 0, 0, 
						0, 0, 1, 0, 
						0, 0, 0, 1};
	GLdouble shear[16] = {1, 0, 0, 0, 
						0, 1, 0, 0, 
						(right+left)/(2*zNear), (top+bottom)/(2*zNear), 1, 0, 
						0, 0, 0, 1};
	matrix_multiply(projection);
	matrix_multiply(scale); 
	matrix_multiply(shear);
}

// remember to convert fovy from degree to radius before calling tan
void I_my_gluPerspective(GLdouble fovy, GLdouble aspect, 
    GLdouble zNear, GLdouble zFar)
{
    // ...
	GLdouble left, right, bottom, top; 
	fovy = (fovy*PI)/180;
	
	top = tan(fovy*0.5) * zNear;
	bottom = -top;
	right = aspect*top;
	left = -right;
	I_my_glFrustum(left, right, bottom, top, zNear, zFar);
}
