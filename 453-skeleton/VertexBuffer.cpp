#include "VertexBuffer.h"

#include <utility>


VertexBuffer::VertexBuffer(GLuint index, GLint size, GLenum dataType)
	: bufferID{}
{
	bind();
	glVertexAttribPointer(index, size, dataType, GL_FALSE, 0, (void*)0);
	glEnableVertexAttribArray(index);
}


void VertexBuffer::uploadData(GLsizeiptr size, const void* data, GLenum usage) {
	bind();
	glBufferData(GL_ARRAY_BUFFER, size, data, usage);
}
