#pragma once

#include "GLHandles.h"

//#include <GL/glew.h>
#include <glad/glad.h>


class VertexBuffer {

public:
	VertexBuffer(GLuint index, GLint size, GLenum dataType);

	// Because we're using the VertexBufferHandle to do RAII for the buffer for us
	// and our other types are trivial or provide their own RAII
	// we don't have to provide any specialized functions here. Rule of zero
	//
	// https://en.cppreference.com/w/cpp/language/rule_of_three
	// https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#Rc-zero

	// Public interface
	void bind() const { glBindBuffer(GL_ARRAY_BUFFER, bufferID); }
	void uploadData(GLsizeiptr size, const void* data, GLenum usage);

private:
	VertexBufferHandle bufferID;
};

