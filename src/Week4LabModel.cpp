#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include "glm/glm.hpp"
#include <CanvasPoint.h>
#include <Colour.h>
#include <TextureMap.h>
#include "RedNoise.h"
#include "triangles.h"
#include <ModelTriangle.h>
#include <map>
#include <string>
#include <tuple>


#define WIDTH 320
#define HEIGHT 240

// interesting way to initialize 2d vector
std::vector<std::vector<float>> depthBuffer(HEIGHT, std::vector<float> (WIDTH, 0));

uint32_t translateColor(Colour colour) {
    uint32_t colorUint = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
    return colorUint;
}

CanvasTriangle sortCanvasPoint(CanvasPoint v0, CanvasPoint v1, CanvasPoint v2) {
    if (v1.y < v0.y) std::swap(v0, v1);
    if (v2.y < v0.y) std::swap(v0, v2);
    if (v2.y < v1.y) std::swap(v1, v2);
    CanvasTriangle triangle = CanvasTriangle(v0, v1, v2);
    return triangle;
}

CanvasPoint findPoint(CanvasPoint v0, CanvasPoint v1, CanvasPoint v2) {
    float coorX = abs(v0.x - v2.x)/abs(v0.y - v2.y) * abs(v0.y - v1.y);
    if (v0.x <= v2.x) coorX = coorX + v0.x;
    else coorX = -coorX + v0.x;
    return CanvasPoint(coorX, v1.y);
}

std::vector<ModelTriangle> readGeoFile(std::string fileName, float scale, std::map<std::string, Colour> palette) {
    std::string text;
    std::ifstream file(fileName);
    std::vector<glm::vec3> vertices;
    std::vector<ModelTriangle> connections;
    std::string name;
    while (std::getline(file, text)) {
        std::vector<std::string> values = split(text, ' ');
        if (values[0] == "usemtl") {
            name = values[1];
        } else if (values[0] == "v") {
            vertices.push_back(scale * glm::vec3 (std::stof(values[1]), std::stof(values[2]), std::stof(values[3])));
        } else if (values[0] == "f") {
//            glm::vec3 vertex1 = vertices[values[1][0] - '0' - 1]; // - '0' is how you convert char to int
            glm::vec3 vertex1 = vertices[std::stoi(split(values[1], '/')[0]) - 1]; // -1 as vertex starts from 1 in obj
            glm::vec3 vertex2 = vertices[std::stoi(split(values[2], '/')[0]) - 1];
            glm::vec3 vertex3 = vertices[std::stoi(split(values[3], '/')[0]) - 1];
            connections.push_back(ModelTriangle(vertex1, vertex2, vertex3, palette[name]));
        }
    }
    return connections;
}

std::map<std::string, Colour> readMtlFile(std::string filename) {
    std::string text;
    std::map<std::string, Colour> palette;
    std::vector<std::string> names;
    std::vector<Colour> colours;
    std::ifstream matFile(filename);
    while (std::getline(matFile, text)) {
        std::vector<std::string> values = split(text, ' ');
        if (values[0] == "newmtl") {
            names.push_back(values[1]);
        } else if (values[0] == "Kd") {
            Colour tempColour = Colour(std::stof(values[1]) * 255, std::stof(values[2]) * 255, std::stof(values[3]) * 255);
            colours.push_back(tempColour);
        }
    }
    for (unsigned int i = 0; i < names.size(); i++) {
        palette[names[i]] = colours[i];
    }
    return palette;
}

std::vector<ModelTriangle> readFiles(std::string geoFileName, std::string matFileName, float scale) {
    std::map<std::string, Colour> palette = readMtlFile(matFileName);
    std::vector<ModelTriangle> connections = readGeoFile(geoFileName, scale, palette);
    return connections;
}

// using the barycentric method to find z https://en.wikipedia.org/wiki/Barycentric_coordinate_system
float findWeights(float x, float y, CanvasTriangle connection) {
    CanvasPoint vertex1 = connection.vertices[0];
    CanvasPoint vertex2 = connection.vertices[1];
    CanvasPoint vertex3 = connection.vertices[2];

    float denominator = (vertex2.y - vertex3.y) * (vertex1.x - vertex3.x) + (vertex3.x - vertex2.x) * (vertex1.y - vertex3.y);
    float a = ((vertex2.y - vertex3.y) * (x - vertex3.x) + (vertex3.x - vertex2.x) * (y - vertex3.y)) / denominator;
    float b = ((vertex3.y - vertex1.y) * (x - vertex3.x) + (vertex1.x - vertex3.x) * (y - vertex3.y)) / denominator;
    float c = 1 - a - b;
    float z = a * vertex1.depth + b * vertex2.depth + c * vertex3.depth;
    return 1 / z;
}


void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour color, CanvasTriangle triangle) {
//    window.clearPixels();
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numOfSteps = std::max(abs(xDiff), abs(yDiff));
    float xStepSize = xDiff/numOfSteps;
    float yStepSize = yDiff/numOfSteps;
    for (float i = 0.0; i < numOfSteps; i++) {
        uint32_t colorUint = translateColor(color);
        float x = from.x + xStepSize * i;
        float y = from.y + yStepSize * i;
        float z = findWeights(x, y, triangle);
        // had x and y in reversed order, couldn't have scale factor bigger than 120
        if (depthBuffer[round(y)][round(x)] < z) {
            depthBuffer[round(y)][round(x)] = z;
            window.setPixelColour(round(x), round(y), colorUint);
        }
    }
}


void drawTriangle(DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour color) {
    CanvasTriangle triangle = CanvasTriangle(v0, v1, v2);
    drawLine(window, triangle.v0(), triangle.v1(), color, triangle);
    drawLine(window, triangle.v0(), triangle.v2(), color, triangle);
    drawLine(window, triangle.v1(), triangle.v2(), color, triangle);
}


void drawFilledTriangle(DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour color) {
    CanvasTriangle triangle = sortCanvasPoint(v0, v1, v2);
    CanvasPoint intersect = findPoint(triangle.v0(), triangle.v1(), triangle.v2());
    float numOfSteps0 = abs(triangle.v0().y - triangle.v1().y);
    float numOfSteps1 = abs((triangle.v2().y - triangle.v1().y));

    float xStepSize0 = (intersect.x - triangle.v0().x) / numOfSteps0;
    float xStepSize1 = (triangle.v2().x - intersect.x) / numOfSteps1;
    float xStepSize2 = (triangle.v1().x - triangle.v0().x) / numOfSteps0;
    float xStepSize3 = (triangle.v2().x - triangle.v1().x) / numOfSteps1;

    for (float i = 0; i < numOfSteps0; i++) {
        CanvasPoint startCanvas = CanvasPoint(round(triangle.v0().x + xStepSize0 * i), triangle.v0().y + i);
        CanvasPoint endCanvas = CanvasPoint(round(triangle.v0().x + xStepSize2 * i), triangle.v0().y + i);
        drawLine(window, startCanvas, endCanvas, color, triangle);
    }

    for (float i = 0; i < numOfSteps1; i++) {
        CanvasPoint startCanvas = CanvasPoint(round(intersect.x + xStepSize1 * i), intersect.y + i);
        CanvasPoint endCanvas = CanvasPoint(round(triangle.v1().x + xStepSize3 * i), intersect.y + i);
        drawLine(window, startCanvas, endCanvas, color, triangle);
    }

}

void drawTriangles(DrawingWindow &window, std::vector<std::tuple<Colour, CanvasTriangle>> triangles) {
    for (std::tuple<Colour, CanvasTriangle> triangleTuple : triangles) {
        Colour colour;
        CanvasTriangle triangle;
        std::tie(colour, triangle) = triangleTuple; // how you get tuple values
        drawFilledTriangle(window, triangle.v0(), triangle.v1(), triangle.v2(), colour);
    }
}

// this function is used to get the intersection point on the image plane from the actual 3d point
CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength) {
    // this is quite confusing. firstly, we are using object coordinate, this way is a bit easier
    // as towards the camera is regarded as the positive z direction
    // this is still not 100% certain in terms of x and y, we will see
    glm::vec3 cameraCoordinate = vertexPosition - cameraPosition;
    float scale_u = focalLength * ((cameraCoordinate.x) / abs(cameraCoordinate.z)) * 180; // very bad practice, but will just leave it like this, scale factor 180
    float scale_v = - 1 * focalLength * ((cameraCoordinate.y) / abs(cameraCoordinate.z)) * 180; // -1 is interesting too, pixel coordinate has reversed y, obj is constructed around center
    float image_u = scale_u + WIDTH / 2;
    float image_v = scale_v + HEIGHT / 2;
    std::cout << cameraCoordinate.z << " " << std::endl;
    CanvasPoint coordinate = CanvasPoint(image_u, image_v, abs(cameraCoordinate.z)); // storing the depth but needs to be abs
    return coordinate;
}

std::vector<std::tuple<Colour, CanvasTriangle>> convertModTriToTri(std::vector<ModelTriangle> connections) {
    std::vector<std::tuple<Colour, CanvasTriangle>> triangles;
    std::vector<CanvasPoint> points;
    std::vector<Colour> colours;
    glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 4.0);
    float focalLength = 2.0;
    for (ModelTriangle connection : connections) {
        Colour colour = connection.colour;
        for (glm::vec3 point3d : connection.vertices) {
            CanvasPoint point = getCanvasIntersectionPoint(cameraPosition, point3d, focalLength);
            points.push_back(point);
        }
        CanvasTriangle triangle = CanvasTriangle(points[0], points[1], points[2]);
        triangles.push_back({colour, triangle});
        points = std::vector<CanvasPoint>(); // this is like a neat way to erase the vector
    }
    return triangles;
}

// the colour thingy here are really bad, might need to change it later
void drawPoints(std::vector<ModelTriangle> connections, DrawingWindow &window, uint32_t colourUnit) {
    glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 4.0);
    float focalLength = 2.0;
    for (ModelTriangle connection : connections) {
        for (glm::vec3 point3d : connection.vertices) {
            CanvasPoint point = getCanvasIntersectionPoint(cameraPosition, point3d, focalLength);
            window.setPixelColour(int(round(point.x)), int(round(point.y)), colourUnit); // could write a function for scaling or include in the existing function
        }
    }
}


// deprecated lmao
// this was designed to find the interpolations of the three vertices in the object domain (3D)
void findDepth(std::vector<ModelTriangle> connections) {
    for (ModelTriangle connection : connections) {
        // this might seem a bit odd, but we want to find pairs like 0-1, 1-2, 2-0 in a single loop
        int count = 1;
        for (glm::vec3 point3d : connection.vertices) {
            float xDiff = connection.vertices[count%3].x - point3d.x;
            float yDiff = connection.vertices[count%3].y - point3d.y;
            float numOfSteps = std::max(abs(xDiff), abs(yDiff));
            float xStepSize = xDiff/numOfSteps;
            float yStepSize = yDiff/numOfSteps;
            for (float i = 0.0; i < numOfSteps; i++) {
                float x = point3d.x + xStepSize * i;
                float y = point3d.y + yStepSize * i;
                // float z = findWeights(x, y, connection); // that z is the object coordinate z
            }
            count++;
        }
    }

}


void nonOcclusionWorkflow(DrawingWindow &window) {
    Colour colour = Colour(255, 255, 255);
    uint32_t colourUnit = translateColor(colour);
    std::vector<ModelTriangle> connections = readFiles("../cornell-box.obj", "../cornell-box.mtl", 0.35);
//            drawPoints(connections, window, colourUnit);
    std::vector<std::tuple<Colour, CanvasTriangle>> triangles = convertModTriToTri(connections);
    drawTriangles(window, triangles);
    for (ModelTriangle connection : connections) {
        std::cout << connection << std::endl;
    }
}


void handleEvent(SDL_Event event, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
        else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
        else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
        else if (event.key.keysym.sym == SDLK_DOWN) {
            nonOcclusionWorkflow(window);
        }
    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
}

int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
    SDL_Event event;

    while (true) {
        if (window.pollForInputEvents(event)) handleEvent(event, window);
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }

}