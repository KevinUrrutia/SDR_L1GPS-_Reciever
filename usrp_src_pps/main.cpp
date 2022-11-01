#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/thread.hpp>
#include <uhd/types/device_addr.hpp>
#include <uhd/types/metadata.hpp>
#include <uhd/types/stream_cmd.hpp>
#include <uhd/stream.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <experimental/filesystem>
#include <string>
#include <chrono>
#include <vector>

//Time libraries
//#include <time.h>
#include <stdio.h>
#include <stdlib.h>
//#include <ctime>
#include "date.h"

//custom built libraries
#include "parser.cpp"
#include "channel_thread.h"
#include "producer.hpp"
#include "consumer.hpp"
#include "find_devices.cpp"

namespace fs = std::experimental::filesystem;

struct tm * gmtime(const time_t *timer);


std::string now() {
  using namespace date;
  const auto curr_time = std::chrono::system_clock::now();
  std::string time_stamp = date::format("%Y%m%d_%H%M%S", curr_time);
  int idx = time_stamp.find(".");
  return time_stamp.substr(0,idx); // drop fractional seconds
}

// std::string now() {
//   std::time_t now = std::time(0);
//   std::tm * now_tm = std::gmtime(&now);
//   char buf[42];
//   std::strftime(buf, 42, "%Y%m%d_%H%m%OS", now_tm);
//   return buf;
// }

// std::string system_now() {
//   std::time_t t = std::time(0);
//   char buffer[9] = {0};
//   std::strftime(buffer, 9, "%H%m%OS", localtime(&t));
//   return buffer;
// }

int UHD_SAFE_MAIN(int argc, char* argv[]) {
    Settings conf_setting;

    std::string type = find_devices_types();

    parse_command_line(argc, argv, &conf_setting);
    parse_config(&conf_setting);

    if(conf_setting.enable_size_map)
      std::cout << "Packet size tracking enabled - will only recv one packet at a time!"
                << std::endl;

    //create a usrp device
    std::cout << std::endl;
    std::cout << boost::format("Creating the USRP device with: %s...") % conf_setting.args
                        << std::endl;
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(conf_setting.args);

    //Lock mboard clocks
    usrp->set_clock_source(conf_setting.clock_source);
    if(type == "B210") {
      usrp->set_time_source(conf_setting.time_source);
    }
    std::cout << "USRP clock source: " << usrp->get_clock_source(0) << std::endl;
    std::cout << "USRP time source: " << usrp->get_time_source(0) << std::endl;
    auto t0 = uhd::time_spec_t(1.0);//uhd::time_spec_t(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count() + 1 );
    usrp->set_time_next_pps(t0);
    std::this_thread::sleep_for(std::chrono::seconds(int64_t(1)));


    //always select the subdevice first, the channel mapping affects other setiings
    if(conf_setting.subdev != "") {
      usrp->set_rx_subdev_spec(conf_setting.subdev);
    }

    std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;

    //set the sampling Rate
    if(conf_setting.rate < 200.0e3) {
      std::cerr << "Please specify a valid sample rate" << std::endl;
      return ~0;
    }

    std::cout << boost::format("Setting RX Rate: %f Msps...") % (conf_setting.rate / 1e6) << std::endl;
    usrp->set_rx_rate(conf_setting.rate);
    std::cout << boost::format("Actual RX Rate: %f Msps...")
                    % (usrp->get_rx_rate() / 1e6)
              << std::endl
              << std::endl;

    //Set the Center frequency
    std::cout << boost::format("Setting RX Freq: %f MHz...") % (conf_setting.freq / 1e6)
              << std::endl;
    std::cout << boost::format("Setting RX LO offset: %f MHz...") % (conf_setting.lo_offset / 1e6)
              << std::endl << std::endl;

    uhd::tune_request_t tune_request(conf_setting.freq, conf_setting.lo_offset);
    if(conf_setting.int_n)
      tune_request.args = uhd::device_addr_t("mode_n=integer");
    usrp->set_rx_freq(tune_request);
    std::cout << boost::format("Actual RX Freq: %f MHz...")
                    % (usrp->get_rx_freq() / 1e6)
              << std::endl
              << std::endl;

    //set the RF gain
    std::cout << boost::format("Setting RX Gain: %f dB...") % conf_setting.gain << std::endl;
    usrp->set_rx_gain(conf_setting.gain);
    std::cout << boost::format("Actual RX Gain: %f dB...")
                    % usrp->get_rx_gain()
              << std::endl
              << std::endl;

    //detect which channels to use
    std::vector<std::string> channel_strings;
    std::vector<size_t> channel_num;
    boost::split(channel_strings, conf_setting.channel, boost::is_any_of("/"","));
    for(size_t ch = 0; ch < channel_strings.size(); ch++) {
      size_t chan = std::stoi(channel_strings[ch]);
      if (chan >= usrp->get_rx_num_channels()) {
        throw std::runtime_error("Invalid channel(s) specified");
      } else {
        channel_num.push_back(chan);
      }
    }

    //set the antenna
    usrp->set_rx_antenna(conf_setting.antenna);

    //Setup for producer thread
    double num_samps = conf_setting.duration * conf_setting.rate;
    auto t_start = usrp->get_time_now().get_real_secs();

    //setup for the consumer thread
    std::string date_str  = now();

    fs::path p(conf_setting.output_file);
    std::string output_dir = p.parent_path();
    std::string ext = p.extension();

    std::string system_log_name = "";
    std::stringstream ss_file_name;
    if(!conf_setting.auto_gain_tuner){
      ss_file_name << output_dir << "/" << conf_setting.clock_source << "_" << conf_setting.signal_type << "_" << date_str << "_utc_" << int(conf_setting.freq) << "_" << int(conf_setting.rate) << "_" << conf_setting.channel << ext;

      std::stringstream ss_log_file_name;
      ss_log_file_name << output_dir << "/" << conf_setting.clock_source << "_" << conf_setting.signal_type << "_" << date_str << "_utc" << "_trig_log.raw";
      system_log_name = ss_log_file_name.str();
      std::cout << system_log_name << std::endl;
    }else {
      ss_file_name << output_dir <<"/gain_test.raw";
    }
    std::string fname = ss_file_name.str();
    std::cout << fname << std::endl;



    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * conf_setting.setup_time)));

    Channel_thread chan_thread;
    std::thread t1(&Producer::local_recv_num_samps, Producer(&chan_thread), std::ref(usrp), t_start, num_samps, channel_num, conf_setting.duration, system_log_name);
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(10)));
    std::thread t2(&Consumer::write_to_file, Consumer(&chan_thread), std::ref(usrp), t_start, fname);
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(10)));

    t1.join();
    t2.join();

    std::cout << "JOINED THREADS" << std::endl;

    return EXIT_SUCCESS;
}
