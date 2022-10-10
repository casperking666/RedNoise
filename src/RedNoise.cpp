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


std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
    std::vector<float> result; // or having a namespace then no need for std
    float interval = (to - from) / (numberOfValues - 1);
    for (int i = 0; i < numberOfValues; i++) {
        result.push_back(from + interval * i);
    }
    return result;
}

std::vector<std::vector<float>> twoDimensionsInterpolation(float from, float to, int yLength, int xLength) {
    std::vector<std::vector<float>> result;
    for (int j = 0; j < yLength; j++) {
        result.push_back(interpolateSingleFloats(from, to, xLength));
    }
    return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
    std::vector<glm::vec3> results;
    std::vector<float> rVector = interpolateSingleFloats(from.r, to.r, numberOfValues);
    std::vector<float> gVector = interpolateSingleFloats(from.g, to.g, numberOfValues);
    std::vector<float> bVector = interpolateSingleFloats(from.b, to.b, numberOfValues); // ctrl g for select occurrence
    for (int i = 0; i < numberOfValues; i++) {
        results.push_back(glm::vec3(rVector[i], gVector[i], bVector[i]));
    }
    return results;
}

std::vector<std::vector<glm::vec3>> threeDimData(glm::vec3 topLeft, glm::vec3 topRight, glm::vec3 bottomLeft, glm::vec3 bottomRight, int yLength, int xLength) {
    std::vector<std::vector<glm::vec3>> results;
    glm::vec3 fraction = (bottomLeft - topLeft) / float(yLength - 1);
    for (float i = 0; i < yLength; i++) {
        results.push_back(interpolateThreeElementValues(topLeft + fraction * i, topRight + fraction * i, xLength));
    }
    return results;
}

void drawInterpolation(DrawingWindow &window, std::vector<std::vector<float>> pixels) {
    window.clearPixels();
    for (size_t y = 0; y < window.height; y++) {
        for (size_t x = 0; x < window.width; x++) {
            float greyScale = pixels[y][x];
            float red = greyScale;
            float green = greyScale;
            float blue = greyScale;
            uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
            window.setPixelColour(x, y, colour);
        }
    }
}

void drawColourfulScreen(DrawingWindow &window, std::vector<std::vector<glm::vec3>> pixels) {
    window.clearPixels();
    for (size_t y = 0; y < window.height; ++y) {
        for (int x = 0; x < window.width; ++x) {
            glm::vec3 rgb = pixels[y][x];
            float red = rgb.r;
            float green = rgb.g;
            float blue = rgb.b;
            uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
            window.setPixelColour(x, y, colour);
        }
    }
}

void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

uint32_t translateColor(Colour colour) {
    uint32_t colorUint = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
    return colorUint;
}

void drawLine(DrawingWindow &window, float fromX, float fromY, float toX, float toY) {
    window.clearPixels();
    CanvasPoint from = CanvasPoint(fromX, fromY);
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


//    test 1
//    std::vector<float> result;
//    result = interpolateSingleFloats(2.2, 8.5, 7);
//    for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
//    std::cout << std::endl;

    // test 2
    glm::vec3 from(1.0, 4.0, 9.2);
    glm::vec3 to(4.0, 1.0, 9.8);
    std::vector<glm::vec3> resultVec3 = interpolateThreeElementValues(from, to, 4);
    for(size_t i=0; i<resultVec3.size(); i++) std::cout << resultVec3[i].r << " " << resultVec3[i].g << " " << resultVec3[i].b << " ";
    std::cout << std::endl;


    while (true) {
        if (window.pollForInputEvents(event)) handleEvent(event, window);
//        std::vector<std::vector<float>> result = twoDimensionsInterpolation(255, 0, window.height, window.width);
//        drawInterpolation(window, result);
        glm::vec3 topLeft(255, 0, 0);        // red
        glm::vec3 topRight(0, 0, 255);       // blue
        glm::vec3 bottomRight(0, 255, 0);    // green
        glm::vec3 bottomLeft(255, 255, 0);   // yellow
        std::vector<std::vector<glm::vec3>> result = threeDimData(topLeft, topRight, bottomLeft, bottomRight, window.height, window.width);
//        drawColourfulScreen(window, result);

        // testing drawLine function
        drawLine(window, 12,23,300,240);

        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }

}
