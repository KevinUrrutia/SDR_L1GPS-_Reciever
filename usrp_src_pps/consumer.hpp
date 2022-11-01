#ifndef __CONSUMER__
#define __CONSUMER__

#include <string>
#include <boost/format.hpp>
#include <mutex>
#include <cstdlib>
#include <cstdio>
#include <complex>

#include "channel_thread.h"

class Consumer {
  Channel_thread* chan_;
public:
  explicit Consumer(Channel_thread* chan) : chan_(chan) {}

  void write_to_file(uhd::usrp::multi_usrp::sptr &usrp, int t_start, std::string fname);
};

void Consumer::write_to_file(uhd::usrp::multi_usrp::sptr &usrp, int t_start, std::string fname){
  std::cout << "Beginning (CONSUMER) thread" << std::endl << std::endl;
  std::cout << "File Name..." << fname << std::endl;

  size_t last_idx = 0;
  std::vector<std::complex<short>> current_items;

  // std::ofstream output_file (fname.c_str(), std::ios::out | std::ios::binary);
  FILE * pFile;
  pFile = fopen(fname.c_str(), "wb");

  while(!chan_->done_producing) {
    //Read producer comments for a more full break down of what is happening below
    //The main difference here although is that we only want the consumer to wake up when we know for a fact
    //that there is something in the shared buffer. This happens when the producer sets the open for open_for_buisness flag to true
    //Only when set to true can the consumer wake up and write to file. Then the open_for_buisness flag is immediately set to false.
    std::unique_lock<std::mutex> ul(chan_->m);
    chan_->cv.wait(ul, [&] {return (chan_->open_for_buisness) ? true : false;});
    chan_->sampsIQ.get_all(chan_->done_producing);
    chan_->open_for_buisness = false;
    ul.unlock();
    chan_->cv.notify_one();
    // output_file.write((char*)&chan_->sampsIQ.curr_items.front(), (chan_->sampsIQ.curr_last_idx)*sizeof(std::complex<short>));
    fwrite(&chan_->sampsIQ.curr_items[0], sizeof(std::complex<short>), chan_->sampsIQ.curr_last_idx, pFile);
    std::this_thread::sleep_for(std::chrono::microseconds(100));
  }
  std::cout << "(CONSUMER) has finished" << std::endl;
  // output_file.close();
  fclose(pFile);
}

#endif
