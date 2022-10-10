#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include "glm/glm.hpp"
#include <CanvasPoint.h>
#include <Colour.h>

#define WIDTH 320
#define HEIGHT 240

uint32_t translateColor(Colour colour) {
    uint32_t colorUint = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
    return colorUint;
}

void drawLine(DrawingWindow &window, float fromX, float fromY, float toX, float toY) {
    window.clearPixels();
    CanvasPoint from = CanvasPoint(fromX, fromY); // useless, just for testing purposes
    CanvasPoint to = CanvasPoint(toX, toY);
    float xDiff = abs(toX - fromX);
    float yDiff = abs(toY - fromY);
    float numOfSteps = fmax(xDiff, yDiff);
    float xStepSize = xDiff/numOfSteps;
    float yStepSize = yDiff/numOfSteps;
    for (float i = 0.0; i < numOfSteps; i++) {
        Colour color = Colour(255, 255, 255);
        uint32_t colorUint = translateColor(color);
        float x = fromX + xStepSize * i;
        float y = fromY + yStepSize * i;
        window.setPixelColour(round(x), round(y), colorUint);
    }
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
        else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
        else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
        else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
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

        // testing drawLine function
        drawLine(window, 12,23,300,240);

        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }

}