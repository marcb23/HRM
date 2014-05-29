// code structure credit to http://www.instructables.com/id/Interface-Python-and-Arduino-with-pySerial/

int count;

void setup() {
  Serial.begin(9600); // set the baud rate
//  Serial.println("Ready"); // print "Ready" once
//  findPy();
  count = 0;
}
void loop() {
  if(count == 1000)
    count = 0;
  Serial.println(count);
  count++;
  delay(1); // delay for 1/10 of a second
}

void findPy() {
  char inByte = ' ';
  while(!Serial.available()) {
    delay(1);
  }
  
  if(Serial.available()){ // only send data back if data has been sent
    char inByte = Serial.read(); // read the incoming data
    Serial.println(inByte); // send the data back in a new line so that it is not all one long line
  }
}
