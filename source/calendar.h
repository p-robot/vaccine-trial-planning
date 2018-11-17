#pragma once

class Calendar {
public:
  Calendar();
  int get_current_day() const;
  int increment_current_day();
private:
  int current_day_;
};