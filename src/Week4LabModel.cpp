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


#define WIDTH 320
#define HEIGHT 240


uint32_t translateColor(Colour colour) {
    uint32_t colorUint = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
    return colorUint;
}

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour color) {
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
        window.setPixelColour(round(x), round(y), colorUint);
    }
}


void drawTriangle(DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour color) {
//    window.clearPixels();
    CanvasTriangle triangle = CanvasTriangle(v0, v1, v2);
    //std::cout << triangle;
    drawLine(window, triangle.v0(), triangle.v1(), color);
    drawLine(window, triangle.v0(), triangle.v2(), color);
    drawLine(window, triangle.v1(), triangle.v2(), color);
}

void drawTriangles(DrawingWindow &window, std::vector<CanvasTriangle> triangles, Colour color) {
    for (CanvasTriangle triangle : triangles) {
        drawTriangle(window, triangle.v0(), triangle.v1(), triangle.v2(), color);
    }
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
            glm::vec3 vertex1 = vertices[std::stoi(split(values[1], '/')[0]) - 1]; // - '0' is how you convert char to int
            glm::vec3 vertex2 = vertices[std::stoi(split(values[2], '/')[0]) - 1]; // - '0' is how you convert char to int
            glm::vec3 vertex3 = vertices[std::stoi(split(values[3], '/')[0]) - 1]; // - '0' is how you convert char to int
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

CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength) {
//    float offset = -1 * cameraPosition.z;
//    float z = offset + vertexPosition.z;
    glm::vec3 cameraCoordinate = vertexPosition - cameraPosition;
    float scale_u = focalLength * ((cameraCoordinate.x) / abs(cameraCoordinate.z)) * 180; // very bad practice, but will just leave it like this
    float scale_v = -1 * focalLength * ((cameraCoordinate.y) / abs(cameraCoordinate.z)) * 180;
    float image_u = scale_u + WIDTH / 2; // the width might be wrong
    float image_v = scale_v + HEIGHT / 2;
    CanvasPoint coordinate = CanvasPoint(image_u, image_v);
    return coordinate;
}

std::vector<CanvasTriangle> convertModTriToTri(std::vector<ModelTriangle> connections) {
    std::vector<CanvasTriangle> triangles;
    std::vector<CanvasPoint> points;
    glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 4.0);
    float focalLength = 2.0;
    for (ModelTriangle connection : connections) {
        for (glm::vec3 point3d : connection.vertices) {
            CanvasPoint point = getCanvasIntersectionPoint(cameraPosition, point3d, focalLength);
            points.push_back(point);
        }
        CanvasTriangle triangle = CanvasTriangle(points[0], points[1], points[2]);
        triangles.push_back(triangle);
        points = std::vector<CanvasPoint>();
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
//            std::cout << point << std::endl;
            window.setPixelColour(int(round(point.x)), int(round(point.y)), colourUnit); // could write a function for scaling or include in the existing function
        }

    }
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
        else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
        else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
        else if (event.key.keysym.sym == SDLK_DOWN) {
            Colour colour = Colour(255, 255, 255);
            uint32_t colourUnit = translateColor(colour);
            std::vector<ModelTriangle> connections = readFiles("../cornell-box.obj", "../cornell-box.mtl", 0.35);
            drawPoints(connections, window, colourUnit);
            std::vector<CanvasTriangle> triangles = convertModTriToTri(connections);
            drawTriangles(window, triangles, colour);
            for (ModelTriangle connection : connections) {
                std::cout << connection << std::endl;
            }


            std::cout << "DOWN" << std::endl;
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