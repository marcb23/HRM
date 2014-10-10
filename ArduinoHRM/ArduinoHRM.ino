#include <TimerOne.h>

int inputPit = 0;
volatile int signal = 0;

void setup() {
  Serial.begin(115200);
  Timer1.initialize(1000); //read at 0.01s frequency
  Timer1.pwm(9, 512);
  Timer1.attachInterrupt(serialWriter);
}

void loop() {
  //signal = analogRead(A0);
  //Serial.println(signal);
}

void serialWriter()
{
  signal = analogRead(A0);
  Serial.println(signal);
}
