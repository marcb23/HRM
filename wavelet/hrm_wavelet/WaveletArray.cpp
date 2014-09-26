#include "WaveletArray.h"

WaveletArray::WaveletArray() {
  index = 0;
  first = 0;
  length = 0;
}

void WaveletArray::appendValue(int value) {
  index %= WINDOWLENGTH;
  circArray[index] = value;
  index++;
  length++;
}

int WaveletArray::getValue(int index) {
  return circArray[(first+index)%WINDOWLENGTH];
}

void WaveletArray::getArray(unsigned int data[], int scale) {
  for(int i=0; i<WINDOWLENGTH; i++) {
    data[i] = circArray[(first+i)%WINDOWLENGTH]*scale;
  }
}

int WaveletArray::getLength() {
  return length;
}

void WaveletArray::newWindow() {
  first += WINDOWDIFF;
  first %= WINDOWLENGTH;
  length -= WINDOWDIFF;
}

void WaveletArray::reset() {
  first = 0;
  index = 0;
}
