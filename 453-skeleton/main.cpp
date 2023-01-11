//#include <GL/glew.h>
//#include <GLFW/glfw3.h>
#include <glad/glad.h>

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <limits>
#include <functional>

#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Texture.h"
#include "Window.h"

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"

using namespace std;
using namespace glm;

// global variables
const vec3 red = vec3{ 1.0, 0.0, 0.0 };
const vec3 green = vec3{ 0.0, 1.0, 0.0 };
const vec3 blue = vec3{ 0.0, 0.0, 1.0 };
const vec3 black = vec3{ 0.0, 0.0, 0.0 };

const float pointDiameter = 10.0f;
const float pointHitboxMargin = 0.02f; // invisible collision hitbox around point
const int numGeneratedPointsOnCurve = 50; // no. of generated points on the curve

// camera values
const float cameraSpeed = 0.1f;
const float cameraSensitivity = 20.0f; // camera rotation
const float cameraYaw = -90.0f;
const float cameraPitch = 0.0f;

vec3* selectedPoint; // point on click
vec3 currLocation; // curr location of the mouse in 2D
vec3 last3DMousePosition; // for rotating the camera
vec3 invalidPoint = vec3{ NULL };
bool hasSelectedPoint;
bool isLeftButtonDown;

vector<vec3>::iterator selectedPointIndex;

const float PI = 3.14159265359;

// adapted from OpenGL Camera tutorial: https://github.com/JoeyDeVries/LearnOpenGL/blob/master/includes/learnopengl/camera.h
enum cameraMovement {
	FORWARD, BACKWARD, LEFT, RIGHT
};

struct Scene {
	bool hasCamera;
	bool allowToggleCurveType;
	bool allowToggleWireframe;
	bool allowSelectTensorSurfaces;
	bool showControlPoints;
	bool allowEditPoints;

	// comparison operator
	bool operator==(const Scene& other) const {
		return (hasCamera == other.hasCamera) &&
			(allowToggleCurveType == other.allowToggleCurveType) &&
			(allowToggleWireframe == other.allowToggleWireframe) &&
			(allowSelectTensorSurfaces == other.allowSelectTensorSurfaces) &&
			(showControlPoints == other.showControlPoints) &&
			(allowEditPoints == other.allowEditPoints);
	}
};

const Scene _2d = { false, true, false, false, true, true };
const Scene _3d = { true, true, false, false, true, false };
const Scene _3dRevolution = { true, false, true, false, false, false };
const Scene _3dTensor = { true, false, true, true, true, false };

static const Scene scenes[4]{ _2d, _3d, _3dRevolution, _3dTensor };
int currSceneIndex = 0;
Scene currScene = _2d;

bool isBezierCurve = true; // toggle between bezier and bspline curve
bool isWireframeView = true; // toggle between wireframe and solid view
bool isTensorSurfaceOne = true; // toggle between tensor surfaces 

// projection variables
const float fov = 45.0f; // realistic
const float zNear = 0.1f;
const float zFar = 100.0f;

const vec3 initEyePos = { 0.0f, 0.0f, 2.0f };
const vec3 initLookAtPos = { 0.0f, 0.0f, -1.0f };
const vec3 initUpVect = { 0.0f, 1.0f, 0.0f };

const vector<vector<vec3>> controlPointsSurfaceOne = {
	{{-2.0f, 0.0f, -2.0f}, {-1.0f, 0.0f, -2.0f}, {0.0f, 0.0f, -2.0f}, {1.0f, 0.0f, -2.0f}, {2.0f, 0.0f, -2.0f}},
	{{-2.0f, 0.0f, -1.0f}, {-1.0f, 1.0f, -1.0f}, {0.0f, 1.0f, -1.0f}, {1.0f, 1.0f, -1.0f}, {2.0f, 0.0f, -1.0f}},
	{{-2.0f, 0.0f, 0.0f}, {-1.0f, 1.0f, 0.0f}, {0.0f, -1.0f, 0.0f}, {1.0f, 1.0f, 0.0f}, {2.0f, 0.0f, 0.0f}},
	{{-2.0f, 0.0f, 1.0f}, {-1.0f, 1.0f, 1.0f}, {0.0f, 1.0f, 1.0f}, {1.0f, 1.0f, 1.0f}, {2.0f, 0.0f, 1.0f}},
	{{-2.0f, 0.0f, 2.0f}, {-1.0f, 0.0f, 2.0f}, {0.0f, 0.0f, 2.0f}, {1.0f, 0.0f, 2.0f}, {2.0f, 0.0f, 2.0f}},
};

const vector<vector<vec3>> controlPointsSurfaceTwo = {
	{{-1.0f, -0.4f, 0.0f}, {-0.8f, -0.4f, 0.3f}, {-0.6f, -0.4f, 0.0f}, {-0.4f, -0.4f, 0.3f}},
	{{-1.0f, -0.2f, 0.0f}, {-0.8f, -0.2f, 0.3f}, {-0.6f, -0.2f, 0.0f}, {-0.4f, -0.2f, 0.3f}},
	{{-1.0f, 0.0f, 0.0f}, {-0.8f, 0.0f, 0.3f}, {-0.6f, 0.0f, 0.0f}, {-0.4f, 0.0f, 0.3f}},
	{{-1.0f, 0.2f, 0.0f}, {-0.8f, 0.2f, 0.03}, {-0.6f, 0.2f, 0.0f}, {-0.4f, 0.2f, 0.3f}},
};

void createPointAt(CPU_Geometry& cpuGeom, vec2 location) {
	cpuGeom.verts.push_back({ location, 0.0 });
}

// check collision: AABB - AABB collision
bool isMouseOnAPoint(vec2 mouseLocation, CPU_Geometry const& pointsLocation) {
	for (const vec3& pointLocation : pointsLocation.verts) {
		bool collisionX = mouseLocation.x >= pointLocation.x - pointHitboxMargin && pointLocation.x + pointHitboxMargin >= mouseLocation.x;
		bool collisionY = mouseLocation.y >= pointLocation.y - pointHitboxMargin && pointLocation.y + pointHitboxMargin >= mouseLocation.y;
		if (collisionX && collisionY) {
			return true;
		}
	}
	return false;
}

// gets the selected point
vec3& getSelectedPoint(vec2 mouseLocation, CPU_Geometry& pointsLocation) {
	for (vec3& pointLocation : pointsLocation.verts) {
		bool collisionX = mouseLocation.x >= pointLocation.x - pointHitboxMargin && pointLocation.x + pointHitboxMargin >= mouseLocation.x;
		bool collisionY = mouseLocation.y >= pointLocation.y - pointHitboxMargin && pointLocation.y + pointHitboxMargin >= mouseLocation.y;
		if (collisionX && collisionY) {
			return pointLocation;
		}
	}
	return invalidPoint;
}

// creates the base 4 points
void initializePoints(CPU_Geometry& cpuGeom) {
	cpuGeom.verts.push_back(vec3{ -0.5, 0.5, 0 });
	cpuGeom.verts.push_back(vec3{ -0.5, -0.5, 0 });
	cpuGeom.verts.push_back(vec3{ 0.5, -0.5, 0 });
	cpuGeom.verts.push_back(vec3{ 0.5, 0.5, 0 });
}

void clearAllPointsAndCurves(CPU_Geometry& controlPointsGeom, CPU_Geometry& bezierCurveGeom, CPU_Geometry& bsplineCurveGeom) {
	controlPointsGeom.verts.clear();
	bezierCurveGeom.verts.clear();
	bsplineCurveGeom.verts.clear();

	hasSelectedPoint = false;
}

void resetPointsAndCurves(CPU_Geometry& controlPointsGeom, CPU_Geometry& bezierCurveGeom, CPU_Geometry& bsplineCurveGeom) {
	clearAllPointsAndCurves(controlPointsGeom, bezierCurveGeom, bsplineCurveGeom);

	initializePoints(controlPointsGeom);
}

class Camera {
public:
	Camera(Window& window, vec3 eyePos, vec3 lookAtPos, vec3 upVect) :
		_fov(fov), _aspect(window.getWidth() / window.getHeight()), _zNear(zNear), _zFar(zFar), 
		_eyePos(eyePos), _lookAtPos(lookAtPos), _upVect(upVect)
	{
		_worldUpVect = upVect;
		_speed = cameraSpeed;
		_sensitivity = cameraSensitivity;
		_yaw = cameraYaw;
		_pitch = cameraPitch;
		updateCameraVects();
	}

	void moveForward() {
		_eyePos += _lookAtPos * cameraSpeed;
	}

	void moveBackward() {
		_eyePos -= _lookAtPos * cameraSpeed;
	}

	void moveLeft() {
		_eyePos -= _rightVect * cameraSpeed;
		//_lookAtPos -= _rightVect * cameraSpeed;
	}

	void moveRight() {
		_eyePos += _rightVect * cameraSpeed;
		//_lookAtPos += _rightVect * cameraSpeed;
	}

	mat4 getViewMatrix() {
		if (currSceneIndex == 0) {
			return mat4(1.0f);
		} else {
			return lookAt(_eyePos, _eyePos + _lookAtPos, _upVect);
		}
	}

	mat4 getPerspectiveMatrix() {
		if (currSceneIndex == 0) {
			return mat4(1.0f);
		} else {
			return perspective(_fov, _aspect, _zNear, _zFar);
		}
	}

	void constrainPitch(float n) {
		if (_pitch > n) {
			_pitch = n;
		}
		else if (_pitch < -n) {
			_pitch = -n;
		}
	}

	// adapted from OpenGL Camera tutorial: https://github.com/JoeyDeVries/LearnOpenGL/blob/master/includes/learnopengl/camera.h
	void processMouseMovement() {
		//cout << currLocation.x << " " << currLocation.y << endl;
		float xOffset = (currLocation.x - last3DMousePosition.x) * cameraSensitivity;
		float yOffset = (last3DMousePosition.y - currLocation.y) * cameraSensitivity;
		last3DMousePosition = currLocation;

		_yaw += xOffset;
		_pitch -= yOffset;

		constrainPitch(89.0f);
		updateCameraVects();
	}

	void resetCamera() {
		_eyePos = initEyePos;
		_lookAtPos = initLookAtPos;
		_upVect = initUpVect;
		_yaw = cameraYaw;
		_pitch = cameraPitch;
	}

private:
	// to get view matrix
	float _fov, _aspect, _zNear, _zFar;

	// to get projection matrix (viewing volume)
	vec3 _eyePos, _lookAtPos, _upVect, _worldUpVect, _rightVect;

	// camera settings
	float _speed, _sensitivity;

	// euler angles
	float _yaw, _pitch;

	void updateCameraVects() {
		vec3 newLookAtPos = {
			cos(radians(_yaw)) * cos(radians(_pitch)),
			sin(radians(_pitch)),
			sin(radians(_yaw)) * cos(radians(_pitch))
		};

		_lookAtPos = normalize(newLookAtPos);
		_rightVect = normalize(cross(_lookAtPos, _worldUpVect));
		_upVect = normalize(cross(_rightVect, _lookAtPos));
	}
};


// EXAMPLE CALLBACKS
class Assignment3 : public CallbackInterface {

public:
	Assignment3(Window* window, CPU_Geometry& controlPointsGeom, CPU_Geometry& bezierCurveGeom, CPU_Geometry& bsplineCurveGeom, Camera& camera) :
		location(vec2{ 0 }),
		_controlPointsGeom(controlPointsGeom),
		_bezierCurveGeom(bezierCurveGeom),
		_bsplineCurveGeom(bsplineCurveGeom),
		_camera(camera)
	{
		_window = window;
		handCursor = glfwCreateStandardCursor(GLFW_HAND_CURSOR);
	}

	virtual void keyCallback(int key, int scancode, int action, int mods) {
		// switch scenes
		if (key == GLFW_KEY_1 && action == GLFW_PRESS) {
			currSceneIndex = 0;
		}
		else if (key == GLFW_KEY_2 && action == GLFW_PRESS) {
			currSceneIndex = 1;
		}
		else if (key == GLFW_KEY_3 && action == GLFW_PRESS) {
			currSceneIndex = 2;
		}
		else if (key == GLFW_KEY_4 && action == GLFW_PRESS) {
			currSceneIndex = 3;
		}

		// custom controls to each scene
		currScene = scenes[currSceneIndex];
		if (currScene.hasCamera) {
			if (key == GLFW_KEY_W && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
				_camera.moveForward();
			}
			else if (key == GLFW_KEY_A && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
				_camera.moveLeft();
			}
			else if (key == GLFW_KEY_S && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
				_camera.moveBackward();
			}
			else if (key == GLFW_KEY_D && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
				_camera.moveRight();
			}
			else if (key == GLFW_KEY_X && action == GLFW_PRESS) {
				_camera.resetCamera();
			}
		}
		if (currScene.allowToggleCurveType) {
			if (key == GLFW_KEY_E && action == GLFW_PRESS) {
				isBezierCurve = !isBezierCurve;
			}
		}
		if (currScene.allowToggleWireframe) {
			if (key == GLFW_KEY_Q && action == GLFW_PRESS) {
				isWireframeView = !isWireframeView;
			}
		}
		if (currScene.allowSelectTensorSurfaces) {
			if (key == GLFW_KEY_UP && action == GLFW_PRESS) {
				isTensorSurfaceOne = true;
			}
			if (key == GLFW_KEY_DOWN && action == GLFW_PRESS) {
				isTensorSurfaceOne = false;
			}
		}

		// general controls
		if (currScene.allowEditPoints) {
			if (key == GLFW_KEY_R && action == GLFW_PRESS) {
				resetPointsAndCurves(_controlPointsGeom, _bezierCurveGeom, _bsplineCurveGeom);
			}
			else if (key == GLFW_KEY_C && action == GLFW_PRESS) {
				clearAllPointsAndCurves(_controlPointsGeom, _bezierCurveGeom, _bsplineCurveGeom);
			}
		}
	
	}
	virtual void mouseButtonCallback(int button, int action, int mods) {
		currScene = scenes[currSceneIndex];

		if (button == GLFW_MOUSE_BUTTON_LEFT) {
			if (action == GLFW_PRESS) {
				if (currScene.allowEditPoints) {
					if (isMouseOnAPoint(location, _controlPointsGeom)) {
						// select existing point
						selectedPoint = &getSelectedPoint(location, _controlPointsGeom);

						hasSelectedPoint = true;
						isLeftButtonDown = true;
					}
					else {
						// create new point
						createPointAt(_controlPointsGeom, location);

						hasSelectedPoint = false;
					}
				}
				if (currScene.hasCamera) {
					// drag camera
					last3DMousePosition = { location, 0 };

					isLeftButtonDown = true;
				}
			}
			else if (action == GLFW_RELEASE) {
				isLeftButtonDown = false;
			}
		}
		if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS && hasSelectedPoint && currScene.allowEditPoints) {
			// delete selected point
			_controlPointsGeom.verts.erase(selectedPointIndex);
			hasSelectedPoint = false;
		}
	}
	virtual void cursorPosCallback(double xpos, double ypos) {
		convertLocation(xpos, ypos);

		if (isMouseOnAPoint(location, _controlPointsGeom)) {
			glfwSetCursor(_window->getWindow(), handCursor); // if cursor position is on a point, change cursor to hover
		}
		else {
			glfwSetCursor(_window->getWindow(), NULL);
		}
	}
	virtual void scrollCallback(double xoffset, double yoffset) {
	}
	virtual void windowSizeCallback(int width, int height) {
		// The CallbackInterface::windowSizeCallback will call glViewport for us
		// width/height in float. update the program for the aspect ratio. if screen is 0 height then ignore
		CallbackInterface::windowSizeCallback(width,  height);
	}
	// converts the location from a (0,0) by (1600,1600) window to the location in GL's coord system.
	void convertLocation(double x, double y) {
		vec2 s(x, y);
		s /= vec2(1600.0);
		s += vec2(-0.5);
		s = vec2(s.x, -s.y);
		s *= 2.f;
		location = s;
		currLocation = { location, 0.0 };
	}

private:
	Window* _window;
	CPU_Geometry& _controlPointsGeom, _bezierCurveGeom, _bsplineCurveGeom;
	Camera& _camera;

	vec2 location;
	GLFWcursor* handCursor;
};

// uses the de Casteljau algorithm to create a Bezier curve
void drawBezierCurve(CPU_Geometry const& controlPointsGeom, CPU_Geometry& bezierCurveGeom) {
	bezierCurveGeom.verts.clear();
	bezierCurveGeom.cols.clear();

	if (controlPointsGeom.verts.size() > 1) {
		vector<vec3> intermediatePoints = controlPointsGeom.verts; // P
		int degree = controlPointsGeom.verts.size() - 1; // d

		for (float u = 0; u < 1.0f; u += 1.0f / (float)numGeneratedPointsOnCurve) {
			for (int i = 0; i < degree; i++) {
				for (int j = 0; j < degree - i; j++) {
					intermediatePoints[j] = (1 - u) * intermediatePoints[j] + u * intermediatePoints[j + 1];
				}
			}
			bezierCurveGeom.verts.push_back(intermediatePoints[0]); // P[0]
		}
	}
}

// uses the Chaikin subdivision algorithm to draw an open quadratic B-Spline curve
void drawBSplineCurve(CPU_Geometry const& controlPointsGeom, CPU_Geometry& bsplineCurveGeom) {
	bsplineCurveGeom.verts.clear();
	bsplineCurveGeom.cols.clear();

	if (controlPointsGeom.verts.size() > 2) {
		vector<vec3> controlPoints = controlPointsGeom.verts;
		int totalIter = ::log2(numGeneratedPointsOnCurve / controlPoints.size()) * 1.5;

		for (int iter = 0; iter < totalIter; iter++) {
			bsplineCurveGeom.verts.clear();

			bsplineCurveGeom.verts.push_back(controlPoints[0]);
			bsplineCurveGeom.verts.push_back(0.5f * controlPoints[0] + 0.5f * controlPoints[1]);

			for (int i = 1; i < controlPoints.size() - 2; i++) {
				bsplineCurveGeom.verts.push_back(0.75f * controlPoints[i] + 0.25f * controlPoints[i + 1]);
				bsplineCurveGeom.verts.push_back(0.25f * controlPoints[i] + 0.75f * controlPoints[i + 1]);
			}

			bsplineCurveGeom.verts.push_back(0.5f * controlPoints[controlPoints.size() - 2] + 0.5f * controlPoints[controlPoints.size() - 1]);
			bsplineCurveGeom.verts.push_back(controlPoints[controlPoints.size() - 1]);

			controlPoints = bsplineCurveGeom.verts;
		}
	}
}

void generate3DRevolution(CPU_Geometry const& bsplineCurveGeom, CPU_Geometry& revolutionSurfaceGeom) {
	revolutionSurfaceGeom.verts.clear();
	revolutionSurfaceGeom.cols.clear();

	int numGeneratedPointsOnSurface = bsplineCurveGeom.verts.size() - 1;
	float vInc = 2.0f * PI / (float)numGeneratedPointsOnSurface;

	for (int u = 0; u < numGeneratedPointsOnSurface; u++) {
		for (float v = 0; v < 2.0f * PI; v += vInc) {
			// build quads for wireframe
			vec3 uv = {
				bsplineCurveGeom.verts[u][0] * cos(v),
				bsplineCurveGeom.verts[u][1],
				bsplineCurveGeom.verts[u][0] * sin(v)
			};

			vec3 u1v = {
				bsplineCurveGeom.verts[u + 1][0] * cos(v),
				bsplineCurveGeom.verts[u + 1][1],
				bsplineCurveGeom.verts[u + 1][0] * sin(v)
			};

			vec3 uv1 = {
				bsplineCurveGeom.verts[u][0] * cos(v + vInc),
				bsplineCurveGeom.verts[u][1],
				bsplineCurveGeom.verts[u][0] * sin(v + vInc)
			};

			vec3 u1v1 = {
				bsplineCurveGeom.verts[u + 1][0] * cos(v + vInc),
				bsplineCurveGeom.verts[u + 1][1],
				bsplineCurveGeom.verts[u + 1][0] * sin(v + vInc)
			};

			// triangle face #1 (|/)
			revolutionSurfaceGeom.verts.push_back(uv);
			revolutionSurfaceGeom.verts.push_back(uv1);
			revolutionSurfaceGeom.verts.push_back(u1v1);

			// triangle face #2 (/|)
			revolutionSurfaceGeom.verts.push_back(u1v1);
			revolutionSurfaceGeom.verts.push_back(u1v);
			revolutionSurfaceGeom.verts.push_back(uv);
		}
	}
}

void generate3DTensor(vector<vector<vec3>> controlPoints, CPU_Geometry& tensorSurfaceGeom) {
	tensorSurfaceGeom.verts.clear();
	tensorSurfaceGeom.cols.clear();

	vector<vector<vec3>> rows, cols;

	int k = controlPoints.size();

	CPU_Geometry controlPointsRow, controlPointsCol;
	CPU_Geometry bsplineCurveRow, bsplineCurveCol;

	for (int i = 0; i < k; i++) {
		controlPointsCol.verts.clear();

		for (int j = 0; j < k; j++) {
			controlPointsCol.verts.push_back(controlPoints[i][j]);
		}
		drawBSplineCurve(controlPointsCol, bsplineCurveCol);
		cols.push_back(bsplineCurveCol.verts);

		controlPointsRow.verts = controlPoints[i];
		drawBSplineCurve(controlPointsRow, bsplineCurveRow);
		rows.push_back(bsplineCurveRow.verts);
	}

	for (int i = 0; i < rows.size() - 1; i++) {
		for (int j = 0; j < rows[i].size() - 1; j++) {
			tensorSurfaceGeom.verts.push_back(rows[i][j]);
			tensorSurfaceGeom.verts.push_back(rows[i][j + 1]);
			tensorSurfaceGeom.verts.push_back(rows[i + 1][j + 1]);

			tensorSurfaceGeom.verts.push_back(rows[i + 1][j + 1]);
			tensorSurfaceGeom.verts.push_back(rows[i + 1][j]);
			tensorSurfaceGeom.verts.push_back(rows[i][j]);
		}
	}

	for (int i = 0; i < cols.size() - 1; i++) {
		for (int j = 0; j < cols[i].size() - 1; j++) {
			tensorSurfaceGeom.verts.push_back(cols[i][j]);
			tensorSurfaceGeom.verts.push_back(cols[i][j + 1]);
			tensorSurfaceGeom.verts.push_back(cols[i + 1][j + 1]);

			tensorSurfaceGeom.verts.push_back(cols[i + 1][j + 1]);
			tensorSurfaceGeom.verts.push_back(cols[i + 1][j]);
			tensorSurfaceGeom.verts.push_back(cols[i][j]);
		}
	}
}

void updateGPUGeometry(GPU_Geometry& gpuGeom, CPU_Geometry const& cpuGeom) {
	gpuGeom.bind();
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setCols(cpuGeom.cols);
}

void updateToGPU(CPU_Geometry& cpuGeom, GPU_Geometry& gpuGeom, vec3 colour) {
	// Reset the colors to draw red -> points, green -> lines
	cpuGeom.cols.clear();
	cpuGeom.cols.resize(cpuGeom.verts.size(), colour);

	updateGPUGeometry(gpuGeom, cpuGeom);
}

void updateSelectedPoint(CPU_Geometry& cpuGeom, GPU_Geometry& gpuGeom, vec3 colour) {
	if (hasSelectedPoint) {
		selectedPointIndex = find(cpuGeom.verts.begin(), cpuGeom.verts.end(), *selectedPoint);
		cpuGeom.cols.at(selectedPointIndex - cpuGeom.verts.begin()) = colour;

		updateGPUGeometry(gpuGeom, cpuGeom);
	}
}

void drawDragAndDropPoint() {
	if (isLeftButtonDown) {
		*selectedPoint = currLocation;
	}
}

void drawDragAndDropCamera(Camera& camera) {
	if (isLeftButtonDown) {
		camera.processMouseMovement();
	}
}

// update view and projection matrices
void updateShaderMatrices(ShaderProgram& shader, Camera& camera) {
	GLint uniformViewMatrix = glGetUniformLocation(shader.programId(), "viewMatrix");
	glUniformMatrix4fv(uniformViewMatrix, 1, GL_FALSE, &camera.getViewMatrix()[0][0]);

	GLint uniformProjectionMatrix = glGetUniformLocation(shader.programId(), "projectionMatrix");
	glUniformMatrix4fv(uniformProjectionMatrix, 1, GL_FALSE, &camera.getPerspectiveMatrix()[0][0]);
}

int main() {
	Log::debug("Starting main");

	// WINDOW
	glfwInit();
	Window window(1600, 1600, "CPSC 453"); // can set callbacks at construction if desired

	GLDebug::enable();

	CPU_Geometry controlPointsGeom, bezierCurveGeom, bsplineCurveGeom, revolutionSurfaceGeom, tensorSurfaceGeom, controlPointsTensorGeom;
	Camera camera(window, initEyePos, initLookAtPos, initUpVect);

	// CALLBACKS
	auto a3 = make_shared<Assignment3>(&window, controlPointsGeom, bezierCurveGeom, bsplineCurveGeom, camera);
	window.setCallbacks(a3);
	
	

	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");

	// The current CPU_Geometry and GPU_Geometry classes are defined in
	// Geometry.h/Geometry.cpp They will work for this assignment, but for some of
	// the bonuses you may have to modify them.
	initializePoints(controlPointsGeom);
	glPointSize(pointDiameter);

	GPU_Geometry pointsGPUGeom, linesGPUGeom, curveGPUGeom, surfaceGPUGeom;
	vector<vector<vec3>> controlPoints; // for tensor 

	// RENDER LOOP
	while (!window.shouldClose()) {
		glfwPollEvents();

		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_FRAMEBUFFER_SRGB);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST); // compare depth values with the z-buffer

		currScene = scenes[currSceneIndex];

		if (currScene.allowEditPoints) {
			drawDragAndDropPoint();
		}
		if (currScene.hasCamera) {
			drawDragAndDropCamera(camera);
		}

		if (currScene.showControlPoints) {
			updateToGPU(controlPointsGeom, linesGPUGeom, green);

			if (currScene == _3dTensor) {
				if (isTensorSurfaceOne) {
					controlPoints = controlPointsSurfaceOne;
				}
				else {
					controlPoints = controlPointsSurfaceTwo;
				}

				controlPointsTensorGeom.verts.clear();
				for (int i = 0; i < controlPoints.size(); i++) {
					for (int j = 0; j < controlPoints[i].size(); j++) {
						controlPointsTensorGeom.verts.push_back(controlPoints[i][j]);
					}
				}
				updateToGPU(controlPointsTensorGeom, pointsGPUGeom, red);
			}
			else {
				updateToGPU(controlPointsGeom, pointsGPUGeom, red); // colour points after lines so that other points do not get overridden by green when colouring selected point.
			}

			if (currScene.allowEditPoints) {
				updateSelectedPoint(controlPointsGeom, pointsGPUGeom, blue);
			}
		}

		if (currScene.allowToggleCurveType) {
			if (isBezierCurve) {
				drawBezierCurve(controlPointsGeom, bezierCurveGeom);
				updateToGPU(bezierCurveGeom, curveGPUGeom, black);
			}
			else {
				drawBSplineCurve(controlPointsGeom, bsplineCurveGeom);
				updateToGPU(bsplineCurveGeom, curveGPUGeom, black);
			}
		} else if (currScene == _3dRevolution) {
			drawBSplineCurve(controlPointsGeom, bsplineCurveGeom);
		}

		if (currScene.allowToggleWireframe) {
			if (isWireframeView) {
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			}
			else {
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			}
		}

		if (currScene == _3dRevolution) {
			generate3DRevolution(bsplineCurveGeom, revolutionSurfaceGeom);
			updateToGPU(revolutionSurfaceGeom, surfaceGPUGeom, blue);
		}
		else if (currScene == _3dTensor) {
			generate3DTensor(controlPoints, tensorSurfaceGeom);
			updateToGPU(tensorSurfaceGeom, surfaceGPUGeom, blue);
		}

		// apply shaders
		shader.use();
		updateShaderMatrices(shader, camera);

		// draw curve before lines so it doesn't get overlapped
		curveGPUGeom.bind();
		if (currScene.allowToggleCurveType) {
			if (isBezierCurve) {
				glDrawArrays(GL_LINE_STRIP, 0, GLsizei(bezierCurveGeom.verts.size()));
			}
			else {
				glDrawArrays(GL_LINE_STRIP, 0, GLsizei(bsplineCurveGeom.verts.size()));
			}
		}

		if (currScene.showControlPoints) {
			if (currScene == _3dTensor) {
				pointsGPUGeom.bind();
				glDrawArrays(GL_POINTS, 0, GLsizei(controlPointsTensorGeom.verts.size()));
			}
			else {
				linesGPUGeom.bind();
				glDrawArrays(GL_LINE_STRIP, 0, GLsizei(controlPointsGeom.verts.size()));

				pointsGPUGeom.bind();
				glDrawArrays(GL_POINTS, 0, GLsizei(controlPointsGeom.verts.size()));
			}
		}

		surfaceGPUGeom.bind();
		if (currScene == _3dRevolution) {
			glDrawArrays(GL_TRIANGLES, 0, GLsizei(revolutionSurfaceGeom.verts.size()));
		}
		else if (currScene == _3dTensor) {
			glDrawArrays(GL_TRIANGLES, 0, GLsizei(tensorSurfaceGeom.verts.size()));
		}

		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		window.swapBuffers();
	}

	glfwTerminate();
	return 0;
}
