//include uhd api libraries
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/thread.hpp>
#include <uhd/types/device_addr.hpp>
#include <uhd/types/metadata.hpp>
#include <uhd/types/stream_cmd.hpp>
#include <uhd/stream.hpp>

//include boost libraries
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>

//include time libraries
#include <stdio.h>
#include <stdlib.h>
#include "date.h"

//include standard libraries
#include <experimental/filesystem>
#include <string>
#include <chrono>
#include <vector>


//include custom libraries
#include "parser.cpp"
#include "channel_thread.h"
#include "producer.hpp"
#include "consumer.hpp"

//namesapces
namespace fs = std::experimental::filesystem;

struct tm * gmtime(const time_t *timer);


std::string now() {
  using namespace date;
  const auto curr_time = std::chrono::system_clock::now();
  std::string time_stamp = date::format("%Y%m%d_%H%M%S", curr_time);
  int idx = time_stamp.find(".");
  return time_stamp.substr(0,idx); // drop fractional seconds
}

//function prototypes
void setup_basic(uhd::usrp::multi_usrp::sptr &usrp, Settings conf_settings, std::vector<size_t> &channel_num);

int UHD_SAFE_MAIN(int argc, char* argv[]) {
  //parse setup parameters
  Settings conf_settings; //store configuration settings

  parse_command_line(argc, argv, &conf_settings);
  parse_config(&conf_settings);

  if(conf_settings.enable_size_map)
    std::cout << "Packet size tracking enabled - will only recieve one packet at a time!"
              << std::endl;

  //create usrp device
  std::cout << std::endl;
  std::cout << boost::format("Creating the USRP with: %s") % conf_settings.args << std::endl;
  uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(conf_settings.args);

  //setup the rest of the usrp parameters
  std::vector<size_t> channel_num;
  setup_basic(std::ref(usrp), conf_settings, std::ref(channel_num));

  //create shared class object
  std::cout << "Creating Resources." << std::endl;
  Channel_thread chan_thread;
  std::cout << "Successfully allocated resources." << std::endl;

  //setup for producer thread
  double num_samps = conf_settings.duration * conf_settings.rate;
  auto t_start = usrp->get_time_now().get_real_secs();

  //setup for consumer thread
  std::string date_str  = now();

  fs::path p(conf_settings.output_file);
  std::string output_dir = p.parent_path();
  std::string ext = p.extension();

  std::stringstream ss_file_name0;
  ss_file_name0 << output_dir << "/" << conf_settings.clock_source << "_" << conf_settings.signal_type << "_" << date_str << "_utc_" << int(conf_settings.freq) << "_" << int(conf_settings.rate) << "_" << conf_settings.channel << ext;
  std::string fname = ss_file_name0.str();
  std::cout << "IQ file name: " << fname << "\n" << std::endl;


  std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * conf_settings.setup_time)));

  std::thread t1(&Producer::local_recv_num_samps, Producer(&chan_thread), std::ref(usrp), t_start, channel_num, num_samps, conf_settings.duration);
  std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(10)));
  std::thread t2(&Consumer::write_to_file, Consumer(&chan_thread), std::ref(usrp), t_start, fname, channel_num);
  std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(10)));

  t1.join();
  t2.join();

  std::cout << "Joined threads!" << std::endl;

  return EXIT_SUCCESS;
}

void setup_basic(uhd::usrp::multi_usrp::sptr &usrp, Settings conf_settings, std::vector<size_t> &channel_num) {
  //always select the subdevice first, the channel mapping affects other settings
  if(conf_settings.subdev != "")
    usrp->set_rx_subdev_spec(conf_settings.subdev);

  std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;

  //Lock mboard clocks
  usrp->set_clock_source(conf_settings.clock_source);
  usrp->set_time_source(conf_settings.time_source);
  std::cout << "USRP clock source: " << usrp->get_clock_source(std::stoi(conf_settings.channel)) << std::endl;
  std::cout << "USRP time source: " << usrp->get_time_source(std::stoi(conf_settings.channel)) << std::endl;
  auto t0 = uhd::time_spec_t(1.0);
  usrp->set_time_next_pps(t0);
  std::this_thread::sleep_for(std::chrono::seconds(int64_t(1)));

  //Desired:
  std::cout << boost::format("Requesting RX Rate: %f Msps") % (conf_settings.rate / 1e6) << std::endl;
  std::cout << boost::format("Requesting RX Freq: %f MHz") % (conf_settings.freq / 1e6) << std::endl;
  std::cout << boost::format("Requesting RX LO Offset: %f MHz") % (conf_settings.lo_offset / 1e6) << std::endl;
  std::cout << boost::format("Requesting RX gain: %d dB") % conf_settings.gain << std::endl;
  std::cout << boost::format("Requesting RX Antenna: %s") % conf_settings.antenna << std::endl;


  //set the sampling rate
  if(conf_settings.rate <= 0.0) {
    throw std::runtime_error("Please specify a valid sampling rate!");
  }
  usrp->set_rx_rate(conf_settings.rate);

  //set the center frequency
  uhd::tune_request_t tune_request(conf_settings.freq, conf_settings.lo_offset);
  usrp->set_rx_freq(tune_request);

  //set the RF gain
  usrp->set_rx_gain(conf_settings.gain);

  //set the antenna
  usrp->set_rx_antenna(conf_settings.antenna);

  //detect which channels to use
  for(size_t ch = 0; ch < conf_settings.channel.size(); ++ch) {
    size_t chan = int(conf_settings.channel[ch]) - 48;
    if(chan >= usrp-> get_rx_num_channels()) {
      throw std::runtime_error("Invalid Channel(s) specified!");
    }
    else {
      channel_num.push_back(chan);
    }
  }

  //Acutal:
  std::cout << boost::format("Actual RX Rate: %f Msps") % (usrp->get_rx_rate() / 1e6) << std::endl;
  std::cout << boost::format("Actual RX Freq: %f MHz") % (usrp->get_rx_freq() / 1e6) << std::endl;
  std::cout << boost::format("Actual RX gain: %d dB") % usrp->get_rx_gain() << std::endl;
  std::cout << boost::format("Actual RX Antenna: %s") % usrp->get_rx_antenna() << std::endl;
}
