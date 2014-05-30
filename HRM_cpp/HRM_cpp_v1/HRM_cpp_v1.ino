#include "HRM.h"
HRM hrm;


long start_time;
unsigned long current_time;
unsigned int now;
byte heart_rate;

byte pulse;
boolean toggle;
boolean sample;
boolean print_bpm;

void setup() {
  //set timer1 interrupt at 1Hz
  TCCR1A = 0;// set entire TCCR1A register to 0
  TCCR1B = 0;// same for TCCR1B
  TCNT1  = 0;//initialize counter value to 0
  // set compare match register for 1hz increments
  // denominator is timer frequency
  OCR1A = 15625/50;// = (16*10^6) / (1*1024) - 1 (must be <65536)
  // turn on CTC mode
  TCCR1B |= (1 << WGM12);
  // Set CS10 and CS12 bits for 1024 prescaler
  TCCR1B |= (1 << CS12) | (1 << CS10);  
  // enable timer compare interrupt
  TIMSK1 |= (1 << OCIE1A);
  
  sample = false;
  toggle = false;
  print_bpm = false;
  pinMode(13, OUTPUT);
  
  Serial.begin(9600);
  
  start_time = millis();
}
  

void loop() {
//Serial.println("class initialized");
  if(sample) {// && counter < 200) {
    int data = map(analogRead(A0),0,1023,0,1000);
//    Serial.println(data);
    sample = false;
    current_time = millis()-start_time;
    now = current_time/10;
    
    heart_rate = hrm.updateData(now, data);
    if(print_bpm){
      Serial.print("Heart Rate: "); Serial.println(heart_rate);
      print_bpm = false;
    }
//    Serial.println(current_time);
  }
}

ISR(TIMER1_COMPA_vect){//timer1 interrupt 1Hz toggles pin 13 (LED)
  pulse++;
  if(pulse == 50) {
//    if (toggle){
//      digitalWrite(13,HIGH);
//    } else {
//      digitalWrite(13,LOW);
//    }
    pulse = 0;
//    toggle = !toggle;
    print_bpm = true;
  }
  sample = true;
}


