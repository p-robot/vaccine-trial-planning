#include "calendar.h"

Calendar::Calendar() 
 : current_day_ {0}
{}

int Calendar::get_current_day() const {
  return current_day_;
}

int Calendar::increment_current_day() {
  return ++current_day_;
}
