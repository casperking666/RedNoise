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
#include <ModelTriangle.h>

#define WIDTH 320
#define HEIGHT 240


std::vector<ModelTriangle> readFile(std::string fileName, float scale) {
    std::string text;
    std::ifstream file(fileName);
    std::vector<glm::vec3> vertices;
    std::vector<ModelTriangle> connections;
    while (std::getline(file, text)) {
        std::vector<std::string> values = split(text, ' ');
        if (values[0] == "v") {
            vertices.push_back(scale * glm::vec3 (std::stof(values[1]), std::stof(values[2]), std::stof(values[3])));
        }
        else if (values[0] == "f") {
            glm::vec3 vertex1 = vertices[values[1][0] - '0' - 1]; // - '0' is how you convert char to int
            glm::vec3 vertex2 = vertices[values[2][0] - '0' - 1];
            glm::vec3 vertex3 = vertices[values[3][0] - '0' - 1];
            connections.push_back(ModelTriangle(vertex1, vertex2, vertex3, Colour()));
        }
    }
    return connections;
}



void handleEvent(SDL_Event event, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
        else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
        else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
        else if (event.key.keysym.sym == SDLK_DOWN) {
            std::vector<ModelTriangle> connections = readFile("../cornell-box.obj", 1);
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