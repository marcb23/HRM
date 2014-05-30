#include "HRM.h"

HRM::HRM() {
  new_beat = true;
  
  threshold = 400;
//  stamps_tail = 0;
//  stamps_head = 0;
//  pos = 0;
  buffer_pulses = 0;
  bpm = 0;
  bpm_last = 0;
}

byte HRM::updateData(unsigned int time_in, int data_in) {
  time.add(time_in);
  data.add(data_in);
  
  setThreshold(time_in, data_in);
    
  if(data_in > threshold) {
//    beats.push(true);
    if(new_beat) {
      findBeats(time_in);
    }
    new_beat = false;
  } else {
//    beats.push(false);
    new_beat = true;
  }
  calcBPM();
  return bpm;
}

///////////////// PRIVATE METHODS /////////////////////////////

void HRM::setThreshold(unsigned int current_time, int data_in) {
//  scale the threshold based on data values within the last buffer_time seconds
  int buffer_time = 2*100;
  int local_min = data_in;
  int local_max = data_in;
  
  while(time.get(time.size()-1)-time.get(0) > buffer_time) {
    time.shift();
    data.shift();
  }
  for(int i=0; i<data.size(); i++) {
    local_min = min(data.get(i), local_min);
    local_max = max(data.get(i), local_max);
  }
  
  threshold = local_max - (local_max-local_min)/2;
    
}

void HRM::findBeats(unsigned int current_time) {
//   filter the beats deque to determine which values are actually heartbeats.
//   delete false positives and add in missed beats as necessary

//   allow the heartbeat to vary by (at most) this amount each time a beat
//   is found
    float variation = 1.4;
    float diff = 0;
    byte lower_bound = 0;
    byte upper_bound = 250;
    byte missed_1 = 0;
    byte missed_2 = 0;
    float inst_bpm = 0;
    unsigned int last_stamp = time_stamps.get(time_stamps.size()-1);
    
    if(bpm_last > 0) {
      upper_bound = bpm_last*variation;
      missed_1 = bpm_last/variation;
      missed_2 = bpm_last/(2*variation);
      lower_bound = bpm_last/(3*variation);
    }
    
    if(time_stamps.size() > 0) {
      diff = (current_time - last_stamp);
      diff /= 100;
      inst_bpm = 60/diff;
    }
//    Serial.print("inst_bpm: "); Serial.println(inst_bpm);
//    Serial.print("upper_bound: "); Serial.println(upper_bound);
    
    if(inst_bpm < upper_bound) {
//      Serial.print("beat found! Threshold= "); Serial.println(threshold);
      // detect if one beat was missed
      if(bpm_last > 0 && inst_bpm < missed_1 and inst_bpm > missed_2) {
        Serial.println("add one");
        // a beat was probably missed, add it in
        time_stamps.add(diff/2 + last_stamp);
        buffer_pulses += 1;

      } else if(bpm_last > 0 && inst_bpm < missed_2 && inst_bpm > lower_bound) {
        Serial.println("add 2");
        // 2 beats were missed
        time_stamps.add(diff/3 + last_stamp);
        time_stamps.add(2*diff/3 + last_stamp);
        buffer_pulses += 2;
      }
      
      time_stamps.add(current_time);
      buffer_pulses += 1;
      //Serial.println("ignoring beat");
    } else {
      Serial.println("ignoring beat");
    }
    // else ignore the beat
}

void HRM::calcBPM() {
  unsigned int buffer_time = 10*100; // in centi-seconds <- not a real word
  byte lower_bound = 45;
  byte upper_bound = 250;
  
  // only calculate the bpm if there are at least buffer_time seconds of beats
  if(time_stamps.size() > 0 && time_stamps.get(time_stamps.size()-1)-time_stamps.get(0) > buffer_time) {
    unsigned int last_stamp = time_stamps.get(time_stamps.size()-1);
    while(time_stamps.get(time_stamps.size()-1) - time_stamps.get(0) > buffer_time) {
      buffer_pulses -= 1;
      time_stamps.shift();
    }
    float diff = (time_stamps.get(time_stamps.size()-1) - time_stamps.get(0));
    diff /= 100;
    
    float calc_bpm = 60*((buffer_pulses-1)/diff);
//    Serial.print("calc bpm: "); Serial.print(calc_bpm);
//    Serial.print(", buffer_pulses: "); Serial.print(buffer_pulses);
//    Serial.print(", diff: "); Serial.println(diff);
    
    if(calc_bpm > upper_bound or calc_bpm < lower_bound) {
      bpm = bpm_last;
    } else if(bpm_last > 0) {
      if(calc_bpm < bpm_last) {
        calc_bpm = bpm_last - 1;
      } else if(calc_bpm > bpm_last) {
        calc_bpm = bpm_last + 1;
      } 
      bpm_last = bpm;
      bpm = calc_bpm;
    } else {
      bpm_last = bpm;
      bpm = calc_bpm;
    }
//  } else {
//    Serial.print("most recent: "); 
//    Serial.print(time_stamps.get(time_stamps.size()-1));
//    Serial.print(", least recent: ");
//    Serial.print(time_stamps.get(0));
//    Serial.print(", buffer pulses: ");
//    Serial.println(buffer_pulses);
  }
}

