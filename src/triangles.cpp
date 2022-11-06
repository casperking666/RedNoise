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

CanvasTriangle generateRandomPoints() {
    CanvasPoint v0 = CanvasPoint(rand()%320, rand()%240); // useless, just for testing purposes
    CanvasPoint v1 = CanvasPoint(rand()%320, rand()%240);
    CanvasPoint v2 = CanvasPoint(rand()%320, rand()%240);
    return CanvasTriangle(v0, v1, v2);
}

Colour generateRandomColor() {
    Colour color = Colour(rand()%255, rand()%255, rand()%255);
    return color;
}

TextureMap loadPPMImage() {
    TextureMap image = TextureMap("texture.ppm");
    return image;
}

// this is for drawing the textureMap into the screen with the same scale, basically like displaying
void drawLineTest(DrawingWindow &window, CanvasPoint from, CanvasPoint to, TextureMap image) {
//    window.clearPixels();
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numOfSteps = std::max(abs(xDiff), abs(yDiff));
    float xStepSize = xDiff/numOfSteps;
    float yStepSize = yDiff/numOfSteps;
    for (float i = 0.0; i < numOfSteps; i++) {
        float x = from.x + xStepSize * i;
        float y = from.y + yStepSize * i;
        int index = round(x + y * image.width);
        uint32_t colorUint = image.pixels[index];
//        std::cout << index << "  " << x << "  " << y << std::endl;
        window.setPixelColour(round(x), round(y), colorUint);
    }
}

// more like drawTextureLine tbh
void drawTexturePoint(DrawingWindow &window, CanvasPoint from, CanvasPoint to, TextureMap image, glm::mat3 matrix) {
//    window.clearPixels();
//    std::cout << "  " << from << "  " << to << std::endl;

    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numOfSteps = std::max(abs(xDiff), abs(yDiff));
    float xStepSize = xDiff/numOfSteps;
    float yStepSize = yDiff/numOfSteps;
    for (float i = 0.0; i < numOfSteps; i++) {
        float x = from.x + xStepSize * i;
        float y = from.y + yStepSize * i;

        glm::vec3 startVector = round(matrix * glm::vec3(x, y, 1));

        TexturePoint startTexture = TexturePoint(startVector.x, startVector.y);

        int index = round(startTexture.x + startTexture.y * image.width);
        uint32_t colorUint = image.pixels[index];
        window.setPixelColour(round(x), round(y), colorUint);
    }
}

glm::mat3 getMatrix(CanvasTriangle texture, CanvasTriangle canvas) {
    CanvasPoint c1 = canvas.v0();
    CanvasPoint c2 = canvas.v1();
    CanvasPoint c3 = canvas.v2();
    CanvasPoint t1 = texture.v0();
    CanvasPoint t2 = texture.v1();
    CanvasPoint t3 = texture.v2();

    // transform from canvas to texture; the matrix applies to the original canvas points
    glm::mat3 original(
            c1.x, c1.y, 1, // first column (not row!)
            c2.x, c2.y, 1, // second column
            c3.x, c3.y, 1);
    glm::mat3 transformed(
            t1.x, t1.y, 1,
            t2.x, t2.y, 1,
            t3.x, t3.y, 1);
    glm::mat3 inverse = glm::inverse(original);
    glm::mat3 matrix = transformed * inverse;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << matrix[i][j] << std::endl;

        }
    }
    return matrix;
}

// I just think it's silly, you gave that goddamn Colour class, now we are using uint32 again. like what?
void drawTexture(DrawingWindow &window) {
    TextureMap image = loadPPMImage();
    CanvasTriangle canvas = CanvasTriangle(CanvasPoint(160, 10), CanvasPoint(300, 230), CanvasPoint(10, 150));
    CanvasTriangle texture = CanvasTriangle(CanvasPoint(195, 5), CanvasPoint(395, 380), CanvasPoint(65, 330));
    CanvasTriangle triangle = sortCanvasPoint(canvas.v0(), canvas.v1(), canvas.v2());
    // not quite sure yet
    glm::mat3 matrix = getMatrix(texture, canvas);
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

        drawTexturePoint(window, startCanvas, endCanvas, image, matrix);
//        drawLineTest(window, startCanvas, endCanvas, image);
    }

    for (float i = 0; i < numOfSteps1; i++) {
//        std::cout << intersect.x << " " << v1.x << std::endl;
        CanvasPoint startCanvas = CanvasPoint(round(intersect.x + xStepSize1 * i), intersect.y + i);
        CanvasPoint endCanvas = CanvasPoint(round(triangle.v1().x + xStepSize3 * i), intersect.y + i);
//        drawLineTest(window, startCanvas, endCanvas, image);
        drawTexturePoint(window, startCanvas, endCanvas, image, matrix);
    }

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
//            CanvasTriangle randVertices = generateRandomPoints();
//            Colour color = generateRandomColor();
//            drawTriangle(window, randVertices.v0(), randVertices.v1(), randVertices.v2(), color);
//            std::cout << "U" << std::endl;
            drawTexture(window);
        }
        else if (event.key.keysym.sym == SDLK_f) {
            CanvasTriangle randVertices = generateRandomPoints();
            Colour color = generateRandomColor();
            drawTriangle(window, randVertices.v0(), randVertices.v1(), randVertices.v2(), color);
            drawFilledTriangle(window, randVertices.v0(), randVertices.v1(), randVertices.v2(), color);
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