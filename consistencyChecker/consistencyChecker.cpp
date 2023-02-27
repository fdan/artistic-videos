// consistencyChecker
// Check consistency of forward flow via backward flow.
//
// (c) Manuel Ruder, Alexey Dosovitskiy, Thomas Brox 2016

#include <algorithm>
#include <assert.h>
#include "CTensor.h"
#include "CFilter.h"

// Which certainty value motion boundaries should get. Value between 0 (uncertain) and 255 (certain).
#define MOTION_BOUNDARIE_VALUE 0

// The amount of gaussian smoothing that sould be applied. Set 0 to disable smoothing.
#define SMOOTH_STRENGH 0.8

// readMiddlebury
bool readMiddlebury(const char* filename, CTensor<float>& flow) {
  FILE *stream = fopen(filename, "rb");
  if (stream == 0) {
    std::cout << "Could not open " << filename << std::endl;
    return false;
  }
  float help;
  int dummy;
  dummy = fread(&help,sizeof(float),1,stream);
  int aXSize,aYSize;
  dummy = fread(&aXSize,sizeof(int),1,stream);
  dummy = fread(&aYSize,sizeof(int),1,stream);
  flow.setSize(aXSize,aYSize,2);
  for (int y = 0; y < flow.ySize(); y++)
    for (int x = 0; x < flow.xSize(); x++) {
      dummy = fread(&flow(x,y,0),sizeof(float),1,stream);
      dummy = fread(&flow(x,y,1),sizeof(float),1,stream);
    }
  fclose(stream);
  return true;
}

void checkConsistency(const CTensor<float>& flow1, const CTensor<float>& flow2, CMatrix<float>& reliable,
                      int argc, char** args) {

// reliable is a white image of the same size as the motion vector

  int xSize = flow1.xSize(), ySize = flow1.ySize();
  int size = xSize * ySize;
  CTensor<float> dx(xSize,ySize,2);
  CTensor<float> dy(xSize,ySize,2);
  CDerivative<float> derivative(3);

  // perform some filtering on flow1?
  NFilter::filter(flow1,dx,derivative,1,1);
  NFilter::filter(flow1,dy,1,derivative,1);

  // create a black image of mvec dimensions
  CMatrix<float> motionEdge(xSize,ySize,0);

  // iterate through all pixels
  for (int i = 0; i < size; i++) {
    motionEdge.data()[i] += dx.data()[i]*dx.data()[i];
    motionEdge.data()[i] += dx.data()[size+i]*dx.data()[size+i];
    motionEdge.data()[i] += dy.data()[i]*dy.data()[i];
    motionEdge.data()[i] += dy.data()[size+i]*dy.data()[size+i];
  }

//loop through all pixels in the motion vector
  for (int ay = 0; ay < flow1.ySize(); ay++)
    for (int ax = 0; ax < flow1.xSize(); ax++) {
      // ax = current index in x
      // ay = current index in y
      // flow1 = the first flow to apply to the white image
      // flow2 = the second flow to apply to the white image warped by flow1

      float bx = ax+flow1(ax, ay, 0);
      // current index in x + value of first flow at current indices, first channel
      // nx in my code

      float by = ay+flow1(ax, ay, 1);
      // current index in y + value of first flow at current indices, second channel
      // ny in my code

      // bx and by are floating point numbers

      // x1, y1 are the floor values of bx and by
      int x1 = floor(bx); // nx in my code, i.e. new x to fill
      int y1 = floor(by); // ny in my code, i.e. new y to fill

      // x2, y2 are the ceiling values of bx and by
      int x2 = x1 + 1;
      int y2 = y1 + 1;

      //null conditions:
      // if new x is less than zero, or if x2 is greater than or equal to the size in x, or
      // if new y is less than zero, or if y2 is greater than or equal to the size in y
      if (x1 < 0 || x2 >= xSize || y1 < 0 || y2 >= ySize)
      { reliable(ax, ay) = 0.0f; continue; } // output is black

      // remainder from the floor function, why do we consider this the alpha?
      float alphaX = bx-x1;
      float alphaY = by-y1;

      // retrieve values from flow2
      // a over b with alpha masking?
      float a  = (1.0-alphaX) * flow2(x1, y1, 0) + alphaX * flow2(x2, y1, 0);
      float b  = (1.0-alphaX) * flow2(x1, y2, 0) + alphaX * flow2(x2, y2, 0);
      float u = (1.0-alphaY)*a+alphaY*b;
      // don't need a and b after this
      // u is flow2 1st channel value at ax+flow1(ax, ay, 0), by+flow1(ax, ay, 1)

      float a_ = (1.0-alphaX) * flow2(x1, y1, 1) + alphaX * flow2(x2, y1, 1);
      float b_ = (1.0-alphaX) * flow2(x1, y2, 1) + alphaX * flow2(x2, y2, 1);
      float v = (1.0-alphaY)*a_+alphaY*b_;
      // don't need and a and b after this
      // v is flow2 2nd channel value at ax+flow1(ax, ay, 0), by+flow1(ax, ay, 1)

      // don't need x1, x2, y1, y2 after this

//      what is the significance of cx, cy?
      float cx = bx+u; // bx is ax + flow1(ax, ay, 0)
      float cy = by+v; // by is ay + flow1(ax, ay, 1)

      // get flow1 values for current indices, first channel
      float u2 = flow1(ax, ay,0);
      // get flow1 values for current indices, last channel
      float v2 = flow1(ax, ay,1);

      // values we need for consistency check: ax, ay, cx, cy, u, u2, v, v2
      // ax and ay are current indices in the flow image
      // u2, v2 are the values for the first flow file at current indices
      // u, v, are the

//      figure 5 in the paper
      // cx-ax is bx+u-ax
      // pow(bx+u-ax, 2) + pow(by+v-ay, 2)
      // compare with u^2 + v^2 + u2^v + v2^2
      // what is 0.01 and 0.5f?: "Coefficients in inequalities (5) and (6) are taken from Sundaram et al"
      if (((cx-ax) * (cx-ax) + (cy-ay) * (cy-ay)) >= 0.01*(u2*u2 + v2*v2 + u*u + v*v) + 0.5f) {
        // Set to a negative value so that when smoothing is applied the smoothing goes "to the outside".
        // Afterwards, we clip values below 0.
        reliable(ax, ay) = -255.0f;
        continue;
      }
      if (motionEdge(ax, ay) > 0.01 * (u2*u2+v2*v2) + 0.002f) {
        reliable(ax, ay) = MOTION_BOUNDARIE_VALUE;
        continue;
      }
    }
}

int main(int argc, char** args) {
  assert(argc >= 4);

  CTensor<float> flow1,flow2;
  readMiddlebury(args[1], flow1);
  readMiddlebury(args[2], flow2);
  
  assert(flow1.xSize() == flow2.xSize());
  assert(flow1.ySize() == flow2.ySize());
  
  int xSize = flow1.xSize(), ySize = flow1.ySize();
  
  // Check consistency of forward flow via backward flow and exlucde motion boundaries
  CMatrix<float> reliable(xSize, ySize, 255.0f);
  checkConsistency(flow1, flow2, reliable, argc, args);
  
  if (SMOOTH_STRENGH > 0) {
    CSmooth<float> smooth(SMOOTH_STRENGH, 2.0f);
    NFilter::filter(reliable, smooth, smooth);
  }
  reliable.clip(0.0f, 255.0f);

  reliable.writeToPGM(args[3]);
}