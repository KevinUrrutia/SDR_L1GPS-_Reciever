#ifndef __CONSUMER__
#define __CONSUMER__

#include <string>
#include <boost/format.hpp>
#include <mutex>
#include <cstdlib>
#include <cstdio>
#include <complex>

#include "channel_thread.h"

namespace fs = std::experimental::filesystem;

class Consumer {
  Channel_thread* chan_;
public:
  explicit Consumer(Channel_thread* chan) : chan_(chan) {}

  void write_to_file(uhd::usrp::multi_usrp::sptr &usrp, int t_start, std::string fname, std::vector<size_t> channel_nums);
};

void Consumer::write_to_file(uhd::usrp::multi_usrp::sptr &usrp, int t_start, std::string fname, std::vector<size_t> channel_nums){
  std::cout << "[CONSUMER] Starting." << std::endl;
  size_t num_chans = channel_nums.size();

  size_t last_idx = 0;
  std::vector<std::vector<std::complex<short>>> current_items;

  fs::path p(fname);
  std::string output_dir = p.parent_path();
  std::string fstem = p.stem();
  std::string ext = p.extension();

  std::vector<std::string> fnames(num_chans, "");
  FILE * pFile;
  std::vector<FILE *> pFiles(num_chans, nullptr);
  for (int i = 0; i < num_chans; i++) {
    std::stringstream fname_ss;
    fname_ss << output_dir << "/" << fstem << channel_nums[i] << ext;
    fnames[i] = fname_ss.str();
    std::cout << "Created: " << fnames[i] << std::endl;
    pFiles[i] = fopen(fnames[i].c_str(), "wb");
  }
  // std::ofstream output_file (fname.c_str(), std::ios::out | std::ios::binary);


  bool unlocked = true;
  auto sec = std::chrono::seconds(1);
  while(!chan_->done_producing) {
    //Read producer comments for a more full break down of what is happening below
    //The main difference here although is that we only want the consumer to wake up when we know for a fact
    //that there is something in the shared buffer. This happens when the producer sets the open for open_for_business flag to true
    //Only when set to true can the consumer wake up and write to file. Then the open_for_business flag is immediately set to false.
    std::unique_lock<std::mutex> ul(chan_->m);
    //chan_->cv.wait(ul, [&] {return (chan_->open_for_business) ? true : false;});
    unlocked = chan_->cv.wait_for(ul, 2*sec,[&] {return (chan_->open_for_business) ? true : false;});
    if (!unlocked) {
      std::cout << "[Consumer] Timed out!" << std::endl;
      break;
    }
    chan_->sampsIQ.get_all(chan_->done_producing);
    chan_->open_for_business = false;
    ul.unlock();
    chan_->cv.notify_one();
    // output_file.write((char*)&chan_->sampsIQ.curr_items.front(), (chan_->sampsIQ.curr_last_idx)*sizeof(std::complex<short>));
    for (int i = 0; i < num_chans; i++) {
      //std::cout << "[Consumer] Writing to file." << std::endl;
      fwrite(&chan_->sampsIQ.curr_items.at(i), sizeof(std::complex<short>), chan_->sampsIQ.curr_last_idx, pFiles[i]);
    }
    //fwrite(&chan_->sampsIQ.curr_items[0], sizeof(std::complex<short>), chan_->sampsIQ.curr_last_idx, pFile);
    std::this_thread::sleep_for(std::chrono::microseconds(100));
    // std::cout << "[CONSUMER] Number of samples saved: " << chan_->sampsIQ.curr_last_idx << std::endl;
  }
  // output_file.close();
  for (int i = 0; i < num_chans; i++) {
    fclose(pFiles[i]);
  }
  std::cout << "\n[CONSUMER] Finished." << std::endl;
}

#endif
