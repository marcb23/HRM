#include "WaveletArray.h"
#include <MemoryFree.h>

const int RATE = 64; //Hertz
const int LEVELS = 7;
const int SCALE_FACTOR = 64; // raw input values can be scaled up by this factor for higher precision
const float LOW_CUTOFF = 0.4;
const float HIGH_CUTOFF = 8.0;

//bool sample = false;
//bool toggle1 = false;
//int pulse = 0;
int data;
//String output;
//int counter;

// dwt variables
WaveletArray rawInput;
unsigned int transformArray[WINDOWLENGTH];
unsigned int temp[WINDOWLENGTH];
float freqs[LEVELS];
float freq_widths[LEVELS];
float powers[LEVELS];


void setup() {
  //set timer1 interrupt at 1Hz
  TCCR1A = 0;// set entire TCCR1A register to 0
  TCCR1B = 0;// same for TCCR1B
  TCNT1  = 0;//initialize counter value to 0
  // set compare match register for 1hz increments
  // denominator is timer frequency
  OCR1A = 15625/RATE-1;// = (16*10^6) / (1*1024) - 1 (must be <65536)
  // turn on CTC mode
  TCCR1B |= (1 << WGM12);
  // Set CS10 and CS12 bits for 1024 prescaler
  TCCR1B |= (1 << CS12) | (1 << CS10);  
  // enable timer compare interrupt
  TIMSK1 |= (1 << OCIE1A);
  
//  data = -1;
//  pinMode(13, OUTPUT);
//  digitalWrite(13, HIGH);
//  counter = 0;
  Serial.begin(9600);
}

// data sampling handled by the interrupt. Every time 
// enough samples have accumulated for one snapshot,
// copy to the transform array and process them to find the heart rate
void loop() {
  
//    Serial.println(rawInput.getValue(rawInput.getLength() -1));
    if(rawInput.getLength() == WINDOWLENGTH) {
//      Serial.print("fm: "); Serial.println(freeMemory());
      noInterrupts();
      rawInput.getArray(transformArray);
      rawInput.newWindow();
      interrupts();
      runTransform();
    }
}

ISR(TIMER1_COMPA_vect){//timer1 interrupt 1Hz toggles pin 13 (LED)
//generates pulse wave of frequency 1Hz/2 = 0.5kHz (takes two cycles for full wave- toggle high then toggle low)
  data = analogRead(A0);
  rawInput.appendValue(data);
//  counter++;
//  if(counter == 20) {
//    if (toggle1){
//      digitalWrite(13,HIGH);
//      toggle1 = 0;
//    } else{
//      digitalWrite(13,LOW);
//      toggle1 = 1;
//    }
//    counter = 0;
//  }
}

void runTransform() {
//  Serial.print("rm: "); Serial.println(rawMean());
  wavedec();
  decodeArray();
  Serial.println(getRate());
} 

// Perform a multilevel 1D discrete wavelet transfrom based on the haar wavelet
void wavedec() {
  
  int length = WINDOWLENGTH >> 1;
//  int temp[length];
  for(int i=0; i<LEVELS; i++) {
    for(int j=0; j<length; j++) {
      int sum = transformArray[j * 2] + transformArray[j * 2 + 1];
      int diff = transformArray[j * 2] - transformArray[j * 2 + 1];
      temp[j] = sum;
      temp[j + length] = diff;
//      transformArray[j] = sum;
//      temp[j] = diff;
      
    }
    
//    for(intj=length; j < length << 1; j++) {
//      transformArray[j] = temp[j-length];
//    }
    for(int j=0; j < length << 1; j++) {
      transformArray[j] = temp[j];
    }
    length >>= 1;
  }
  
}

// Parse the wavelet transform results into arrays to be processed
void decodeArray() {
  float nyquist = RATE / 2;
  for(int i=1; i<= LEVELS; i++) {
    float lower_limit = nyquist / power(2, LEVELS-i+1);
    float upper_limit = nyquist / power(2, LEVELS-i);
    freqs[i-1] = (lower_limit + upper_limit) / 2.0;
//    Serial.print("got frequency: "); Serial.println(freqs[i-1]);
    if(i != 1) {
//      freq_widths[i-1] = (upper_limit - lower_limit) / nyquist; // remember to scale later!
      freq_widths[i-1] = freqs[i-1] - freqs[i-2];
    }
    if(freqs[i-1] <= HIGH_CUTOFF && freqs[i-1] >= LOW_CUTOFF) {
      powers[i-1] = transformMean(i);
    } else {
      powers[i-0] = 0;
    }
  }
}

// Take the centroid of the frequency vs. power plot (as if it's a solid shape)
float getRate() {
  float rate = 0;
  float numerator = 0;
  float denominator = 0;
  
  for(int i=1; i<LEVELS; i++) {
    numerator += (power(freqs[i]-freqs[i-1], 2) / 6.0) * (2.0 * powers[i] + powers[i-1]) * freq_widths[i];
    denominator += ((freqs[i]-freqs[i-1])/2.0)*(powers[i] + powers[i-1]) * freq_widths[i];
//    Serial.print("N: "); Serial.println(numerator);
//    Serial.print("D: "); Serial.println(denominator);
  }
  
  rate = numerator / denominator;
  return 60*rate;
}
  
// some mathy functions
int power(int base, int exponent) {
  int result = 1;
  for(int i=0; i<exponent; i++) {
    result *= base;
  }
  return result;
}  

float transformMean(int sectionIndex) {
  int first = WINDOWLENGTH / power(2, LEVELS-sectionIndex+1);
  int length = first;
  unsigned int sum = 0;
  for(int i=0; i<first+length; i++) {
    sum += transformArray[i];
  }
  float mean = sum / float(length);
//  Serial.print("tm: "); Serial.print(sectionIndex); Serial.print(", "); Serial.println(mean);
  return mean;
}

unsigned int rawMean() {
  long sum = 0;
  for(int i=0; i<WINDOWLENGTH; i++) {
    sum += transformArray[i];
  }
  return sum/WINDOWLENGTH;
}
