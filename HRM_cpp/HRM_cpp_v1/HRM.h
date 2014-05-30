#ifndef HRM_h
#define HRM_h

#include <Arduino.h>
#include "QueueList.h"
#include "LinkedList.h"

class HRM {
  public:
    HRM();
    byte updateData(unsigned int time_in, int data_in);
  
  private:
    static const int MAX_LENGTH = 100;
    LinkedList <unsigned int> time;  // unsigned int with resolution 12.34 seconds only goes up to 10 minutes
    LinkedList <int> data;
//    QueueList <boolean> beats;

  // imma try to replace the whole beats array with this, new_beat is true if the last data point was
  // below the threshold
    boolean new_beat;
//    QueueList <int> bpm;
    LinkedList<unsigned int> time_stamps;
    
    float threshold;
//    int pos;
//    byte stamps_tail;  // hopefully don't use these cuz Queues do good things
//    byte stamps_head;
    byte buffer_pulses;
    float bpm;
    float bpm_last;
    
    void setThreshold(unsigned int current_time, int data_in);
    void findBeats(unsigned int current_time); 
    void calcBPM();
};

#endif
