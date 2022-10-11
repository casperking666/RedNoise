#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include "glm/glm.hpp"
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>

#define WIDTH 320
#define HEIGHT 240

uint32_t translateColor(Colour colour) {
    uint32_t colorUint = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
    return colorUint;
}

CanvasTriangle createTriangle(CanvasPoint v0, CanvasPoint v1, CanvasPoint v2) {
    CanvasTriangle triangle = CanvasTriangle(v0, v1, v2);
    return triangle;
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
    CanvasTriangle triangle = createTriangle(v0, v1, v2);
    //std::cout << triangle;
    drawLine(window, triangle.v0(), triangle.v1(), color);
    drawLine(window, triangle.v0(), triangle.v2(), color);
    drawLine(window, triangle.v1(), triangle.v2(), color);
}

CanvasTriangle sortCanvasPoint(CanvasPoint v0, CanvasPoint v1, CanvasPoint v2) {
    if (v1.y < v0.y) std::swap(v0, v1);
    if (v2.y < v0.y) std::swap(v0, v2);
    if (v2.y < v1.y) std::swap(v1, v2);
    CanvasTriangle triangle = createTriangle(v0, v1, v2);
    return triangle;
}

CanvasPoint findPoint(CanvasPoint v0, CanvasPoint v1, CanvasPoint v2) {
//    float diffV01 = (v1.y - v0.y);
//    float diffV02 = (v2.y - v0.y);
//    float proportion = abs(diffV02) / abs(diffV01);
//    float coorX = v0.x - (v0.x / proportion);
//    float coorX = abs(v0.x - v2.x) /
    float coorX = abs(v0.x - v2.x)/abs(v0.y - v2.y) * abs(v0.y - v1.y);
    if (v0.x <= v2.x) coorX = coorX + v0.x;
    else coorX = -coorX + v0.x;
    return CanvasPoint(coorX, v1.y);
}

void drawFilledTriangle(DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour color) {
    CanvasTriangle triangle = sortCanvasPoint(v0, v1, v2);
//    std::cout << triangle.v0() << std::endl;
//    std::cout << triangle.v1() << std::endl;
//    std::cout << triangle.v2() << std::endl;
    CanvasPoint intersect = findPoint(triangle.v0(), triangle.v1(), triangle.v2());
    //std::cout << intersect << std::endl;
    float numOfSteps0 = abs(triangle.v0().y - triangle.v1().y);
    float numOfSteps1 = abs((triangle.v2().y - triangle.v1().y));

    float xStepSize0 = (intersect.x - triangle.v0().x) / numOfSteps0;
    float xStepSize1 = (triangle.v2().x - intersect.x) / numOfSteps1;
    float xStepSize2 = (triangle.v1().x - triangle.v0().x) / numOfSteps0;
    float xStepSize3 = (triangle.v2().x - triangle.v1().x) / numOfSteps1;

    for (float i = 0; i < numOfSteps0; i++) {
        CanvasPoint startCanvas = CanvasPoint(round(triangle.v0().x + xStepSize0 * i), triangle.v0().y + i);
        CanvasPoint endCanvas = CanvasPoint(round(triangle.v0().x + xStepSize2 * i), triangle.v0().y + i);
        drawLine(window, startCanvas, endCanvas, color);
    }

    for (float i = 0; i < numOfSteps1; i++) {
//        std::cout << intersect.x << " " << v1.x << std::endl;
        CanvasPoint startCanvas = CanvasPoint(round(intersect.x + xStepSize1 * i), intersect.y + i);
        CanvasPoint endCanvas = CanvasPoint(round(triangle.v1().x + xStepSize3 * i), intersect.y + i);
        drawLine(window, startCanvas, endCanvas, color);
    }

}

void handleEvent(SDL_Event event, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
        else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
        else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
        else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
        else if (event.key.keysym.sym == SDLK_u) {

            CanvasPoint v0 = CanvasPoint(rand()%320, rand()%240); // useless, just for testing purposes
            CanvasPoint v1 = CanvasPoint(rand()%320, rand()%240);
            CanvasPoint v2 = CanvasPoint(rand()%320, rand()%240);
            Colour color = Colour(rand()%255, rand()%255, rand()%255);

            drawTriangle(window, v0, v1, v2, color);
            std::cout << "U" << std::endl;
        }

        else if (event.key.keysym.sym == SDLK_f) {
            CanvasPoint v0 = CanvasPoint(rand()%320, rand()%240); // useless, just for testing purposes
            CanvasPoint v1 = CanvasPoint(rand()%320, rand()%240);
            CanvasPoint v2 = CanvasPoint(rand()%320, rand()%240);
            Colour color = Colour(rand()%255, rand()%255, rand()%255);
            drawTriangle(window, v0, v1, v2, color);
            drawFilledTriangle(window, v0, v1, v2, color);
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

//        CanvasPoint from = CanvasPoint(240, 120); // useless, just for testing purposes
//        CanvasPoint to = CanvasPoint(20, 120);
////        testing drawLine function
//        Colour color = Colour(rand()%255, rand()%255, rand()%255);
//        drawLine(window, from, to, color);

        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }

}