#include "date.h"
#include <iostream>

// https://howardhinnant.github.io/date/date.html#to_stream_formatting
int main()
{
  using namespace date;
  const auto curr_time = std::chrono::system_clock::now();
  //std::string time_stamp = date::format("%F %T", curr_time);
  std::string time_stamp = date::format("%Y%m%d_%H%M%S", curr_time);
  int idx = time_stamp.find(".");
  time_stamp = time_stamp.substr(0,idx);
  std::cout << time_stamp << std::endl;
}
// g++ test_date.cpp -o test_date
// ./test_date
