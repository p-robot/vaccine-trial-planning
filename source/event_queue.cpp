#include <queue>
#include "event.h"
#include "event_queue.h"

Event EventQueue::pop() {
  Event top_event = event_queue_.top();
  event_queue_.pop();
  return top_event;
}
int EventQueue::day_of_next_event() const {
  return event_queue_.top().get_day();
}
int EventQueue::size() const {
  return event_queue_.size();
}
void EventQueue::schedule_event(const Event &event) {
  event_queue_.push(event);
}

