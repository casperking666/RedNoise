#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include "glm/glm.hpp"
#include <glm/gtx/string_cast.hpp>
#include <CanvasPoint.h>
#include <Colour.h>
#include <TextureMap.h>
#include "RedNoise.h"
#include "triangles.h"
#include <ModelTriangle.h>
#include <map>
#include <string>
#include <tuple>
#include <cmath>
#include <RayTriangleIntersection.h>
#include <math.h>


#define WIDTH 320
#define HEIGHT 240

// interesting way to initialize 2d vector
std::vector<std::vector<float>> depthBuffer(HEIGHT, std::vector<float> (WIDTH, 0));
glm::vec3 lightSource = glm::vec3(-0.5,0.9,2.0); // not sure
//glm::vec3 lightSource = glm::vec3(0.4,0.5,3.2); // not sure
glm::mat3 cameraOrientation(
        1,0,0,
        0,1,0,
        0,0,1
);

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
            // glm::vec3 vertex1 = vertices[values[1][0] - '0' - 1]; // - '0' is how you convert char to int
            // did pretty wrong here completely ignored the '/' thingy
            glm::vec3 vertex1 = vertices[std::stoi(split(values[1], '/')[0]) - 1]; // -1 as vertex starts from 1 in obj
            glm::vec3 vertex2 = vertices[std::stoi(split(values[2], '/')[0]) - 1];
            glm::vec3 vertex3 = vertices[std::stoi(split(values[3], '/')[0]) - 1];
            ModelTriangle temp = ModelTriangle(vertex1, vertex2, vertex3, palette[name]);
            temp.normal = glm::normalize(glm::cross((vertex3 - vertex1), (vertex2 - vertex1))); // changed to - now outputs sth
            connections.push_back(temp);
//            std::cout << glm::to_string(temp.normal) << std::endl;
        }
    }
    return connections;
}

std::vector<ModelTriangle> readSphereFile(std::string fileName, float scale) {
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
            // glm::vec3 vertex1 = vertices[values[1][0] - '0' - 1]; // - '0' is how you convert char to int
            // did pretty wrong here completely ignored the '/' thingy
            glm::vec3 vertex1 = vertices[std::stoi(split(values[1], '/')[0]) - 1]; // -1 as vertex starts from 1 in obj
            glm::vec3 vertex2 = vertices[std::stoi(split(values[2], '/')[0]) - 1];
            glm::vec3 vertex3 = vertices[std::stoi(split(values[3], '/')[0]) - 1];
            ModelTriangle temp = ModelTriangle(vertex1, vertex2, vertex3, Colour(255,0,0));
            temp.normal = glm::normalize(glm::cross((vertex3 - vertex1), (vertex2 - vertex1))); // changed to - now outputs sth
            connections.push_back(temp);
//            std::cout << glm::to_string(temp.normal) << std::endl;
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
    return z;
}


void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour color, CanvasTriangle triangle) {
//    window.clearPixels();
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numOfSteps = std::max(abs(xDiff), abs(yDiff));
    float xStepSize = xDiff/numOfSteps;
    float yStepSize = yDiff/numOfSteps;
    for (float i = 0.0; i <= numOfSteps; i++) { // set to <= now the missing rasterizing points are filled
        uint32_t colorUint = translateColor(color);
        float x = from.x + xStepSize * i;
        float y = from.y + yStepSize * i;
        if (x > WIDTH - 1 or y > HEIGHT - 1 or x < 0 or y < 0) continue;
        float z = findWeights(x, y, triangle);
        // had x and y in reversed order, couldn't have scale factor bigger than 120
        if (depthBuffer[round(y)][round(x)] < z) {
            depthBuffer[round(y)][round(x)] = z;
            window.setPixelColour(round(x), round(y), colorUint);
        }
    }
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

bool checkWithinTriangle(glm::vec3 possibleSolution) {
    float t = possibleSolution[0];
    float u = possibleSolution[1];
    float v = possibleSolution[2];
    if ((u >= 0 && u <= 1) && (v >= 0 && v <= 1) && (u + v <= 1) && t > 0) {
        return true;
    } else return false;
}

RayTriangleIntersection getClosestIntersections(glm::vec3 cameraPosition, glm::vec3 rayDirection, std::vector<ModelTriangle> triangles) {
    int i = 0;
    float temp = 20; // just for testing
    RayTriangleIntersection ray;
    for (ModelTriangle triangle : triangles) {
//        std::cout << glm::to_string(triangle.normal) << std::endl;
        std::array<glm::vec3, 3> vertices = triangle.vertices;
        glm::vec3 e0 = vertices[1] - vertices[0];
        glm::vec3 e1 = vertices[2] - vertices[0];
        glm::vec3 SPVector = cameraPosition - vertices[0];
        glm::mat3 DEMatrix(-rayDirection, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector; // u, v, t

        if (checkWithinTriangle(possibleSolution)) {
            if (possibleSolution[0] < temp) {
                glm::vec3 r = vertices[0] + possibleSolution[1] * (vertices[1] - vertices[0]) + possibleSolution[2] * (vertices[2] - vertices[0]);
                ray = RayTriangleIntersection(r, possibleSolution[0], triangle, i);
                temp = possibleSolution[0];
            }
        }
        i++;
    }
//     std::cout << glm::to_string(ray.intersectedTriangle.normal) << std::endl;
    return ray;
}

glm::vec3 getReflectionVec(glm::vec3 incidenceVec, glm::vec3 normal) {
    float scale = 2.0;
    glm::vec3 reflectionVec = incidenceVec - scale * normal * (incidenceVec * normal);
    return glm::normalize(reflectionVec);
}

float specularLighting(float brightness, glm::vec3 reflectionVec, glm::vec3 rayDirection) {
    // note this should be plus, got me for a while
    float newBrightness = brightness + pow(glm::dot(reflectionVec, rayDirection), 512);
    return glm::clamp<float>(newBrightness, 0.0f, 1.0f);
}

float getBrightnessProximity(float distance) {
    float intensity = 30;
    float brightness = intensity / (4 * M_PI * distance * distance);
    return brightness;
}

float angleOfIncidenceLighting(float brightness, glm::vec3 normal, glm::vec3 vecToLight) {
    float dotProd = glm::dot(vecToLight, normal);
    float adjBrightness = brightness * dotProd;
    return adjBrightness;
}

glm::vec3 getDirection(glm::vec3 cameraPosition, CanvasPoint coordinate, float focalLength) {
    float image_u = coordinate.x;
    float image_v = coordinate.y;
    float scale_u = image_u - WIDTH / 2;
    float scale_v = image_v - HEIGHT / 2;
    float x = scale_u / (150 * focalLength); // x/z to be more specific
    float y = scale_v / (150 * focalLength);
    glm::vec3 cameraCoordinate = (glm::vec3(x, -y, -1));
//    std::cout << glm::to_string(cameraCoordinate) << std::endl;
    return cameraCoordinate;
}

glm::vec3 findNormalWeights(float x, float y, glm::vec3 vertex1, glm::vec3 vertex2, glm::vec3 vertex3) {
    float denominator = (vertex2.y - vertex3.y) * (vertex1.x - vertex3.x) + (vertex3.x - vertex2.x) * (vertex1.y - vertex3.y);
    float a = ((vertex2.y - vertex3.y) * (x - vertex3.x) + (vertex3.x - vertex2.x) * (y - vertex3.y)) / denominator;
    float b = ((vertex3.y - vertex1.y) * (x - vertex3.x) + (vertex1.x - vertex3.x) * (y - vertex3.y)) / denominator;
    float c = 1 - a - b;
//    if (a != NULL) std::cout << a << " " << b << " " << c << std::endl;
    return (glm::vec3(a, b, c));
}

glm::vec3 calculateIntersectionNormal(glm::vec3 weights, glm::vec3 vertex1Norm, glm::vec3 vertex2Norm, glm::vec3 vertex3Norm) {
    glm::vec3 normal = weights.x * vertex1Norm;
    normal += weights.y * vertex2Norm;
    normal += weights.z * vertex3Norm;
    return normal;
}

glm::vec3 findVertexNormal(glm::vec3 curVertex, std::vector<ModelTriangle> triangles) {
    float i = 0.0; // to exclude that vertex itself
    glm::vec3 vertexNormal {0,0,0};
    for (auto triangle : triangles) {
        //for (auto vertex : triangle.vertices) {
            if (triangle.vertices[0] == curVertex || triangle.vertices[1] == curVertex || triangle.vertices[2] == curVertex) {
                vertexNormal += triangle.normal;
//                std::cout << glm::to_string(vertexNormal) << std::endl;
                i += 1.0;
            }
    }
    glm::vec3 average = vertexNormal / i;
    return (average);
}

float brightnessProcessing(glm::vec3 surfaceToLightSource, glm::vec3 normal, glm::vec3 rayDirection) {
    // proximity testing
    float distance = glm::length(surfaceToLightSource);
    float brightness = getBrightnessProximity(distance);
    // spent so much time on incidence I think, cus the light source was weird, couldn't get the correct lighting, the rest was fine just without this normalize
    float angleOfIncidenceBrightness = angleOfIncidenceLighting(brightness, normal, glm::normalize(surfaceToLightSource));

    glm::vec3 reflectionVec = getReflectionVec(-glm::normalize(surfaceToLightSource), normal);
    float specularBrightness = specularLighting(angleOfIncidenceBrightness, reflectionVec, rayDirection);
    float threshold = 0.2;
    float brightnessThreshold = std::max(threshold, specularBrightness);
    return brightnessThreshold;
}

void drawGouraud(DrawingWindow &window, glm::vec3 cameraPosition, float focalLength, std::vector<ModelTriangle> triangles) {
    for (float j = 0; j < HEIGHT; j++) {
        for (float i = 0; i < WIDTH; i++) {
            glm::vec3 rayDirection = getDirection(cameraPosition, CanvasPoint(i, j), focalLength);
            RayTriangleIntersection ray = getClosestIntersections(cameraPosition, rayDirection, triangles);
            glm::vec3 surfaceToLightSource = -lightSource + ray.intersectionPoint; // think about vector calculation

            glm::vec3 vertex1Normal = findVertexNormal(ray.intersectedTriangle.vertices[0], triangles);
            glm::vec3 vertex2Normal = findVertexNormal(ray.intersectedTriangle.vertices[1], triangles);
            glm::vec3 vertex3Normal = findVertexNormal(ray.intersectedTriangle.vertices[2], triangles);

            float brightness1 = brightnessProcessing(surfaceToLightSource, vertex1Normal, rayDirection);
            float brightness2 = brightnessProcessing(surfaceToLightSource, vertex2Normal, rayDirection);
            float brightness3 = brightnessProcessing(surfaceToLightSource, vertex3Normal, rayDirection);

            glm::vec3 weights = findNormalWeights(ray.intersectionPoint.x, ray.intersectionPoint.y, ray.intersectedTriangle.vertices[0], ray.intersectedTriangle.vertices[1], ray.intersectedTriangle.vertices[2]);

            float brightnessCombine = weights.x * brightness1 + weights.y * brightness2 + weights.z * brightness3;

            RayTriangleIntersection shadowRay = getClosestIntersections(lightSource, surfaceToLightSource, triangles);
            Colour color = ray.intersectedTriangle.colour;
            Colour adjustedColor = Colour(color.red * brightnessCombine, color.green * brightnessCombine, color.blue * brightnessCombine);
            window.setPixelColour(i, j, translateColor(adjustedColor));
        }
    }
}


void drawPhong(DrawingWindow &window, glm::vec3 cameraPosition, float focalLength, std::vector<ModelTriangle> triangles) {
    for (float j = 0; j < HEIGHT; j++) {
        for (float i = 0; i < WIDTH; i++) {
            glm::vec3 rayDirection = getDirection(cameraPosition, CanvasPoint(i, j), focalLength);
            RayTriangleIntersection ray = getClosestIntersections(cameraPosition, rayDirection, triangles);
            glm::vec3 surfaceToLightSource = -lightSource + ray.intersectionPoint; // think about vector calculation

            glm::vec3 vertex1Normal = findVertexNormal(ray.intersectedTriangle.vertices[0], triangles);
            glm::vec3 vertex2Normal = findVertexNormal(ray.intersectedTriangle.vertices[1], triangles);
            glm::vec3 vertex3Normal = findVertexNormal(ray.intersectedTriangle.vertices[2], triangles);

            glm::vec3 weights = findNormalWeights(ray.intersectionPoint.x, ray.intersectionPoint.y, ray.intersectedTriangle.vertices[0], ray.intersectedTriangle.vertices[1], ray.intersectedTriangle.vertices[2]);

            glm::vec3 normal = calculateIntersectionNormal(weights, vertex1Normal, vertex2Normal, vertex3Normal);
            // had a really strange bug here where I assign normal to the ray normal and do the processing later, black dots appear when I call this new process function
            // it works fine when the function code are here, might have sth to do with reference and copy shit.
            float brightnessThreshold = brightnessProcessing(surfaceToLightSource, normal, rayDirection);

            RayTriangleIntersection shadowRay = getClosestIntersections(lightSource, surfaceToLightSource, triangles);
            Colour color = ray.intersectedTriangle.colour;
            Colour adjustedColor = Colour(color.red * brightnessThreshold, color.green * brightnessThreshold, color.blue * brightnessThreshold);
            window.setPixelColour(i, j, translateColor(adjustedColor));
        }
    }
}


void draw(DrawingWindow &window, glm::vec3 cameraPosition, float focalLength, std::vector<ModelTriangle> triangles) {
    for (float j = 0; j < HEIGHT; j++) {
        for (float i = 0; i < WIDTH; i++) {
            glm::vec3 rayDirection = getDirection(cameraPosition, CanvasPoint(i, j), focalLength);
            RayTriangleIntersection ray = getClosestIntersections(cameraPosition, rayDirection, triangles);
            glm::vec3 surfaceToLightSource = -lightSource + ray.intersectionPoint; // think about vector calculation

            float brightness = brightnessProcessing(surfaceToLightSource, ray.intersectedTriangle.normal, rayDirection);

            RayTriangleIntersection shadowRay = getClosestIntersections(lightSource, surfaceToLightSource, triangles);
            Colour color = ray.intersectedTriangle.colour;
            if (shadowRay.triangleIndex == ray.triangleIndex) {
                Colour adjustedColor = Colour(color.red * brightness, color.green * brightness, color.blue * brightness);
                // std::cout << color << std::endl;
                window.setPixelColour(i, j, translateColor(adjustedColor));
            } else window.setPixelColour(i, j, translateColor(Colour(color.red * 0.6 * brightness,color.green * 0.6 * brightness,color.blue * 0.6 * brightness)));
        }
    }
}

// this function is used to get the intersection point on the image plane from the actual 3d point
CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength) {
    // this is quite confusing. firstly, we are using object coordinate, this way is a bit easier
    // as towards the camera is regarded as the positive z direction
    // this is still not 100% certain in terms of x and y, we will see
    glm::vec3 cameraCoordinate = vertexPosition - cameraPosition;
    cameraCoordinate = cameraCoordinate * cameraOrientation; // something do with the coordinate system
    float scale_u = focalLength * ((cameraCoordinate.x) / abs(cameraCoordinate.z)) * 100; // very bad practice, but will just leave it like this, scale factor 180
    float scale_v = -1 * focalLength * ((cameraCoordinate.y) / abs(cameraCoordinate.z)) * 100; // -1 is interesting too, pixel coordinate has reversed y, obj is constructed around center
    float image_u = scale_u + WIDTH / 2;
    float image_v = scale_v + HEIGHT / 2;
//    std::cout << cameraCoordinate.z << " " << std::endl;
    CanvasPoint coordinate = CanvasPoint(image_u, image_v, abs(1 / cameraCoordinate.z)); // storing the depth but needs to be abs
    return coordinate;
}

std::vector<std::tuple<Colour, CanvasTriangle>> convertModTriToTri(std::vector<ModelTriangle> connections, glm::vec3 cameraPosition) {
    std::vector<std::tuple<Colour, CanvasTriangle>> triangles;
    std::vector<CanvasPoint> points;
    std::vector<Colour> colours;
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


glm::mat3 getRotationMatrixX(float angle) {
    glm::mat3 rotationMatrix(
            1, 0, 0,
            0, std::cos(angle), std::sin(angle),
            0, -1 * std::sin(angle), std::cos(angle));
    return rotationMatrix;
}

glm::mat3 getRotationMatrixY(float angle) {
    glm::mat3 rotationMatrix(
            std::cos(angle), 0, -1 * std::sin(angle),
            0, 1, 0,
            std::sin(angle), 0, std::cos(angle));
    return rotationMatrix;
}

float convertDegreesToRadians(float angle) {
    return angle * (M_PI / 180);
}

void changeOrientationMatrix(float angle, bool isX) {
    glm::mat3 rotationMatrix;
    if (isX) rotationMatrix = getRotationMatrixX(convertDegreesToRadians(angle));
    else rotationMatrix = getRotationMatrixY(convertDegreesToRadians(angle));
    cameraOrientation = rotationMatrix * cameraOrientation;
}

void cleanDepthBuffer() {
    for (float j = 0; j < depthBuffer.size(); j++) {
        for (float i = 0; i < depthBuffer[0].size(); i++) {
            depthBuffer[j][i] = 0.0;
        }
    }
}


void translateCamera(glm::vec3 shiftFactor, glm::vec3 *cameraPosition) {
    *cameraPosition = *cameraPosition + shiftFactor;
}

// it only works when orbitng with the y-axis
glm::mat3 getOrientationMatrixY(glm::vec3 cameraPosition) {
    glm::vec3 forward = glm::normalize(cameraPosition);
    glm::vec3 right = glm::cross(glm::vec3(0,1,0), forward);
    glm::vec3 up = glm::cross(forward, right);
    glm::mat3 orientation(
            right,up,forward
            );
    return orientation;
}

glm::mat3 getInitialOrientationMatrix() {
    glm::mat3 orientation(
            1,0,1,
            0,1,0,
            0,0,1
    );
    return orientation;
}


void nonOcclusionWorkflow(DrawingWindow &window, glm::vec3 cameraPosition) {
    std::vector<ModelTriangle> connections = readFiles("../cornell-box.obj", "../cornell-box.mtl", 0.35);
    // drawPoints(connections, window, colourUnit);
    std::vector<std::tuple<Colour, CanvasTriangle>> triangles = convertModTriToTri(connections, cameraPosition);
    drawTriangles(window, triangles);
//    for (ModelTriangle connection : connections) {
//        std::cout << connection << std::endl;
//    }
}


void handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3 *cameraPosition) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_LEFT) {
            window.clearPixels();
            cleanDepthBuffer();
            glm::vec3 shiftFactor = glm::vec3 (0.1, 0.0, 0.0);
            translateCamera(shiftFactor, cameraPosition);
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_RIGHT) {
            std::cout << "RIGHT" << std::endl;
            window.clearPixels();
            cleanDepthBuffer();
            glm::vec3 shiftFactor = glm::vec3 (-0.1, 0.0, 0.0);
            translateCamera(shiftFactor, cameraPosition);
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_UP) {
            std::cout << "UP" << std::endl;
            window.clearPixels();
            cleanDepthBuffer();
            glm::vec3 shiftFactor = glm::vec3 (0.0, -0.1, 0.0);
            translateCamera(shiftFactor, cameraPosition);
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_DOWN) {
            window.clearPixels();
            cleanDepthBuffer();
            glm::vec3 shiftFactor = glm::vec3 (0.0, 0.1, 0.0);
            translateCamera(shiftFactor, cameraPosition);
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_o) { // test for week 4 last question
            nonOcclusionWorkflow(window, *cameraPosition);
        } else if (event.key.keysym.sym == SDLK_f) {
            window.clearPixels();
            cleanDepthBuffer();
            glm::vec3 shiftFactor = glm::vec3 (0.0, 0.0, -0.1); // z is a bit different, I feel like it would change the focal length as well
            translateCamera(shiftFactor, cameraPosition);
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_b) {
            window.clearPixels();
            cleanDepthBuffer();
            glm::vec3 shiftFactor = glm::vec3 (0.0, 0.0, 0.1);
            translateCamera(shiftFactor, cameraPosition);
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_x) {
            window.clearPixels();
            cleanDepthBuffer();
            *cameraPosition = getRotationMatrixX(convertDegreesToRadians(5)) * *cameraPosition;
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_c) {
            window.clearPixels();
            cleanDepthBuffer();
            *cameraPosition = getRotationMatrixX(convertDegreesToRadians(-5)) * *cameraPosition;
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_y) {
            window.clearPixels();
            cleanDepthBuffer();
            *cameraPosition = getRotationMatrixY(convertDegreesToRadians(5)) * *cameraPosition;
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_u) {
            window.clearPixels();
            cleanDepthBuffer();
            *cameraPosition = getRotationMatrixY(convertDegreesToRadians(-5)) * *cameraPosition;
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_w) {
            window.clearPixels();
            cleanDepthBuffer();
            changeOrientationMatrix(5, true); // for those we don't change cameraPosition at all
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_s) {
            window.clearPixels();
            cleanDepthBuffer();
            changeOrientationMatrix(-5, true);
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_a) {
            window.clearPixels();
            cleanDepthBuffer();
            changeOrientationMatrix(5, false);
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_d) {
            window.clearPixels();
            cleanDepthBuffer();
            changeOrientationMatrix(-5, false);
            std::cout << glm::to_string(*cameraPosition) << std::endl;
        } else if (event.key.keysym.sym == SDLK_1) {
            window.clearPixels();
            cleanDepthBuffer();
            glm::vec3 shiftFactor = glm::vec3 (0.0, 0.1, 0.0);
            lightSource += shiftFactor;
        } else if (event.key.keysym.sym == SDLK_2) {
            window.clearPixels();
            cleanDepthBuffer();
            glm::vec3 shiftFactor = glm::vec3 (0.0, -0.1, 0.0);
            lightSource += shiftFactor;
        } else if (event.key.keysym.sym == SDLK_3) {
            window.clearPixels();
            cleanDepthBuffer();
            glm::vec3 shiftFactor = glm::vec3 (-0.1, 0.0, 0.0);
            lightSource += shiftFactor;
        } else if (event.key.keysym.sym == SDLK_4) {
            window.clearPixels();
            cleanDepthBuffer();
            glm::vec3 shiftFactor = glm::vec3 (0.1, 0.0, 0.0);
            lightSource += shiftFactor;
        } else if (event.key.keysym.sym == SDLK_5) {
            window.clearPixels();
            cleanDepthBuffer();
            glm::vec3 shiftFactor = glm::vec3 (0.0, 0.0, 0.1);
            lightSource += shiftFactor;
        } else if (event.key.keysym.sym == SDLK_6) {
            window.clearPixels();
            cleanDepthBuffer();
            glm::vec3 shiftFactor = glm::vec3 (0.0, 0.0, -0.1);
            lightSource += shiftFactor;
        }
        std::cout << glm::to_string(lightSource) << std::endl;
    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
}

void drawTrianglesLoop(DrawingWindow &window, std::vector<ModelTriangle> connections, glm::vec3 *cameraPosition) {
    window.clearPixels();
    cleanDepthBuffer();
    *cameraPosition = getRotationMatrixY(convertDegreesToRadians(-1)) * *cameraPosition;
    cameraOrientation = getOrientationMatrixY(*cameraPosition);
    std::cout << glm::to_string(*cameraPosition) << std::endl;
    std::vector<std::tuple<Colour, CanvasTriangle>> triangles = convertModTriToTri(connections, *cameraPosition);
    drawTriangles(window, triangles);
}

//void drawSphere(DrawingWindow &window)

int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
    SDL_Event event;

    // Initialization
    glm::vec3 var = glm::vec3(0.0, 0.0, 4.0);
    glm::vec3 *cameraPosition;
    cameraPosition = &var;
//    std::vector<ModelTriangle> connections = readFiles("../cornell-box.obj", "../cornell-box.mtl", 0.35);
    std::vector<ModelTriangle> connections = readSphereFile("../sphere.obj", 0.5);

    // ray tracing part, solely for testing purpose
//    glm::vec3 cameraDirection = glm::normalize(glm::vec3(0.0, 0.0, -2));
//    RayTriangleIntersection ray = getClosestIntersections(*cameraPosition, cameraDirection, connections);
//    std::cout << ray << std::endl;


    while (true) {
        if (window.pollForInputEvents(event)) handleEvent(event, window, cameraPosition);
//        std::vector<std::tuple<Colour, CanvasTriangle>> triangles = convertModTriToTri(connections, *cameraPosition); temporally no need for this line
        drawGouraud(window, *cameraPosition, 2.0, connections);

//        drawTriangles(window, triangles);
//        drawTrianglesLoop(window, connections, cameraPosition);
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }

}

// Just in case I forgot, glm matrices are fucking retarded, for the rest do m * v, for the cameraCoordinate one,
// since simon has said it that way, do v * m.

// at least I know what the problem was, fucking simon's stupid pdf slides
// the bug with cross product is true, but the main problem was it needs to be v * m instead of m * v ffs.