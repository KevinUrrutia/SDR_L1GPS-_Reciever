#ifndef __PRODUCER__
#define __PRODUCER__

//include boost libraries
#include <boost/format.hpp>

//include standard libraries
#include <string>
#include <vector>
#include <mutex>

//include custom libraries
#include "channel_thread.h"

class Producer{
  Channel_thread* chan_;

public:
  explicit Producer(Channel_thread* chan) : chan_(chan) {}

  void local_recv_num_samps(uhd::usrp::multi_usrp::sptr &usrp, int t_start, std::vector<size_t> channel_num, double num_samps, double duration);
};

void Producer::local_recv_num_samps(uhd::usrp::multi_usrp::sptr &usrp, int t_start, std::vector<size_t> channel_num, double num_samps, double duration) {
  std::cout << "[PRODUCER] Beginning producer thread" << std::endl;
  std::cout << "[PRODUCER] Num Requested Samples: " << num_samps << std::endl;

  /*
  RX a finite number of samples from the USRP
  :param num_samps: number of samps to RX
  :return: array of floating-point samples (sc16)
  */

  bool stream_now;
  size_t num_chans = channel_num.size();
  usrp->clear_command_time();


  //create a recieve streamer
  uhd::stream_args_t stream_args("sc16", "sc16");
  stream_args.channels = channel_num;
  uhd::rx_metadata_t md;
  uhd::rx_streamer::sptr rx_streamer = usrp->get_rx_stream(stream_args);


  // allocate buffers to receive with samples (one buffer per channel)
  const size_t samps_per_buff = rx_streamer->get_max_num_samps();
  std::vector<std::vector<std::complex<short>>> buffs(
      num_chans, std::vector<std::complex<short>>(samps_per_buff));

  //setup_streaming
  uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
  double t_sleep = t_start - usrp->get_time_now().get_real_secs();
  if(t_sleep > 0 || num_chans > 1) {
    stream_now = false;
  }
  else {
    stream_now = true;
  }
  stream_cmd.stream_now = stream_now;
  if(t_sleep > 0) {
    stream_cmd.time_spec = uhd::time_spec_t(static_cast<double>(t_start));
  } else if(num_chans > 1) {
    t_start = usrp->get_time_now().get_real_secs() + 0.1;
    stream_cmd.time_spec = uhd::time_spec_t(static_cast<double>(t_start));
  }

  rx_streamer->issue_stream_cmd(stream_cmd);

  // create a vector of pointers to point to each of the channel buffers
  std::vector<std::complex<short>*> buff_ptrs;
  for (size_t i = 0; i < buffs.size(); i++)
      buff_ptrs.push_back(&buffs[i].front());

  double timeout = 5.0;

  size_t num_acc_samps = 0;

  while(num_acc_samps < num_samps) {
    //recieve a single packet
    size_t num_rx_samps = rx_streamer->recv(buff_ptrs, samps_per_buff, md, timeout);
    //for subsequent packets use small timeout
    timeout = 0.5;

    //handle the error codes
    if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
      throw std::runtime_error(
        str(boost::format("[Producer]: Rx stream error: %s") % md.strerror()));
        break;
    }
    if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) {
      throw std::runtime_error(
        str(boost::format("[Producer]: Rx stream error: %s") % md.strerror()));
        break;
    }

    //The unique lock is something that gets released at the termination of the stack, therefore needs
    //To be kept alive during the entirety of the while loop.
    //Here the condition varialbe allows for the lock to bounce back and forth between the producer and comsumer
    //Here the producer has no condition to wait on, as soon as it is notified that the producer has the lock the producer can read from usb
    //The notify lets the sleeping thread know that the lock has been passed over.
    std::unique_lock<std::mutex> ul(chan_->m);
    chan_->cv.wait(ul, [] {return true;});
    chan_->sampsIQ.add(buffs);
    chan_->open_for_business = true;
    ul.unlock();
    chan_->cv.notify_one();
    num_acc_samps += num_rx_samps;
  }

  stream_cmd = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
  rx_streamer->issue_stream_cmd(stream_cmd);

  size_t samps_num = 1;
  while(samps_num > 0) {
    const size_t num_rx_samps = rx_streamer->recv(buff_ptrs, samps_per_buff, md);
    samps_num = num_rx_samps;
  }

  rx_streamer = nullptr;

  if (num_acc_samps < num_samps) {
    std::cerr << "[Producer] Receive timeout before all samples were received..." << std::endl;
  }

  std::cout << "\n[PRODUCER] Finished Streaming" << std::endl;

}



#endif
