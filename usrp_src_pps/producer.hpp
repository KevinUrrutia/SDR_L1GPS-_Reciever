#ifndef __PRODUCER__
#define __PRODUCER__

#include <string>
#include <boost/format.hpp>
#include <mutex>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <fstream>
#include <iostream>

#include "date.h"
// #include <condition_variable>

#include "channel_thread.h"

std::string now_ms() {
  // using namespace date;
  // const auto curr_time = std::chrono::system_clock::now();
  // std::string time_stamp = date::format("%Y%m%d_%H%M%S", curr_time);
  // int idx = time_stamp.find(".");
  // return time_stamp.substr(0,idx); // drop fractional seconds
  using date::operator<<;
  static char const* const fmt = "%Y-%m-%d %T";
  std::ostringstream ss;
  auto tp = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
  ss << date::format(fmt, tp);
  return ss.str();
}

class Producer {
  Channel_thread* chan_;

public:
  explicit Producer(Channel_thread* chan) : chan_(chan) {}

  void local_recv_num_samps(uhd::usrp::multi_usrp::sptr &usrp, int t_start, double num_samps, std::vector<size_t> channel_num, int duration, std::string file_name);

};

void Producer::local_recv_num_samps(uhd::usrp::multi_usrp::sptr &usrp, int t_start, double num_samps, std::vector<size_t> channel_num, int duration, std::string file_name) {
  std::cout << "Beginning (PRODUCER) thread" << std::endl << std::endl;
  std::cout << "Number of Samples..." << num_samps << std::endl;

  /*
  RX a finite number of samples for the usrp
  :param num_samps: number of samps to the RX
  :return: array of floatng-point samples (sc16)
  */
  std::string date_str  = now_ms();

  bool stream_now;
  size_t num_chans = channel_num.size();
  usrp->clear_command_time();
  usrp->set_time_unknown_pps(uhd::time_spec_t(0.0));
  std::this_thread::sleep_for(std::chrono::milliseconds(8)); //wait for external trigger pulse, 125Hz

  //open file for logging
  FILE * pFile_log;
  if (file_name != "") {
    pFile_log = fopen(file_name.c_str(), "wb");
    char log_file_header[] = "Trigger_count,System_Time,Metadata_Time,PPS_Time,Sample_Count\n";
    fwrite(log_file_header, 1, sizeof(log_file_header), pFile_log);
  }
  // std::ofstream log_file(file_name)

  //create a recieve streamer
  uhd::stream_args_t stream_args("sc16", "sc16");
  stream_args.channels = channel_num;
  uhd::rx_metadata_t md;
  uhd::rx_streamer::sptr rx_streamer = usrp->get_rx_stream(stream_args);



  uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
  stream_now = false;
  stream_cmd.stream_now = stream_now;//stream_now;

  stream_cmd.time_spec = uhd::time_spec_t(static_cast<double>(5)); //seconds in future are set to 1.5s

  rx_streamer->issue_stream_cmd(stream_cmd);

  // allocate buffers to receive with samples (one buffer per channel)
  const size_t samps_per_buff = rx_streamer->get_max_num_samps();
  std::vector<std::vector<std::complex<short>>> buffs(
      usrp->get_rx_num_channels(), std::vector<std::complex<short>>(samps_per_buff));

  // create a vector of pointers to point to each of the channel buffers
  std::vector<std::complex<short>*> buff_ptrs;
  for (size_t i = 0; i < buffs.size(); i++)
      buff_ptrs.push_back(&buffs[i].front());

  double timeout = 5;

  size_t num_acc_samps = 0; //number of accumulated samples

  auto last_pps_time = usrp->get_time_last_pps(0).get_real_secs();
  unsigned long int trig_cnt = 1;

  std::stringstream log_file;
  while (num_acc_samps < num_samps) {

    //recieve a single packet
    size_t num_rx_samps = rx_streamer->recv(buff_ptrs, samps_per_buff, md, timeout);
    auto current_pps_time = usrp->get_time_last_pps(0).get_real_secs();
    if(current_pps_time - last_pps_time > 1e-6) {
      date_str = now_ms();

      log_file << trig_cnt << ",";
      log_file << date_str << ",";
      log_file << std::setprecision(8) << md.time_spec.get_real_secs() << ",";
      log_file << std::setprecision(8) << current_pps_time << ",";
      log_file << num_acc_samps << "\n";

      if (file_name != "") {
      fwrite(log_file.str().c_str(), 1, log_file.str().size(), pFile_log);
      }
      log_file.str("");


      // std::cout << "New PPS has occurred" << std::endl;
      // std::cout << "UTC_TIME: " << date_str << std::endl;
      // std::cout << "USRP time: " << std::setprecision(15) << std::showpoint <<  md.time_spec.get_real_secs() << std::endl;
      // std::cout << "Number of Samples: " << num_acc_samps << std::endl;
      last_pps_time = current_pps_time;
      trig_cnt += 1;
    }

    // use a small timeout for subsequent packets
    timeout = 0.1;

    //handle the error codes
    if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
      break;
    }
    if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) {
      throw std::runtime_error(
        str(boost::format("Receiver error %s") % md.strerror()));
    }

    //The unique lock is something that gets released at the termination of the stack, therefore needs
    //To be kept alive during the entirety of the while loop.
    //Here the condition varialbe allows for the lock to bounce back and forth between the producer and comsumer
    //Here the producer has no condition to wait on, as soon as it is notified that the producer has the lock the producer can read from usb
    //The notify lets the sleeping thread know that the lock has been passed over.
    std::unique_lock<std::mutex> ul(chan_->m);
    chan_->cv.wait(ul, [] {return true;});
    chan_->sampsIQ.add(buffs.at(0));
    chan_->open_for_buisness = true;
    ul.unlock();
    chan_->cv.notify_one();
    num_acc_samps += num_rx_samps;
  }
  chan_->done_producing = true;
  stream_cmd = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
  rx_streamer->issue_stream_cmd(stream_cmd);

  size_t samps_num = 1;
  while(samps_num > 0) {
    const size_t num_rx_samps = rx_streamer->recv(buff_ptrs, samps_per_buff, md);
    samps_num = num_rx_samps;
  }

  rx_streamer = nullptr;

  if (num_acc_samps < num_samps) {
    std::cerr << "Recieve timeout before all samples were recieved..." << std::endl;
  }

  //Finished
  if (file_name != "") {
    fclose(pFile_log);
  }
  std::cout << std::endl << "(PRODUCER) has finished" << std::endl;

}

#endif
