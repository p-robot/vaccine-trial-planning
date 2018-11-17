#pragma once
#include <queue>
#include "event.h"

class EventCompare {
public:
  bool operator()(const Event& e1, const Event& e2) {
    return e1.get_day() + e1.get_day_offset() > e2.get_day() + e2.get_day_offset();
  }
};

class EventQueue {
public:
  Event pop();
  int size() const;
  int day_of_next_event() const;
  void schedule_event(const Event &event);
private:
  std::priority_queue<Event, std::vector<Event>, EventCompare> event_queue_;
};