// Implement a circular array with functionality specifically designed for the
// wavelet transform-based heart rate monitor

#include "WaveletArray.h"

WaveletArray::WaveletArray() {
  index = 0;
  first = 0;
  length = 0;
  if(PRELOAD) {
    length = WINDOWLENGTH;
  }
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

  void WaveletArray::getArray(int data[], int scale) {
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

void WaveletArray::loadArray(int preload[]) {
  for(int i=0; i<WINDOWLENGTH; i++) {
    appendValue(preload[i]);
  }
}
