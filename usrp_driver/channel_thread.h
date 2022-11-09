#ifndef CHANNEL_THREAD
#define CHANNEL_THREAD

#include <mutex>
#include <condition_variable>
#include "ItemStore.hpp"

int BUFF_SIZE = 60e6;

class Channel_thread {
  ItemStore sampsIQ{BUFF_SIZE};
  bool done_producing = false;
  bool open_for_business = false;
  std::mutex m;
  std::condition_variable cv;

  friend class Producer;
  friend class Consumer;

public:
};

#endif
