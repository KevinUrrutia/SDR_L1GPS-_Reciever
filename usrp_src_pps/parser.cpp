#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <algorithm>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <experimental/filesystem>
#include <stdio.h>
#include <fcntl.h>
#include <bitset>

namespace io = boost::iostreams;
namespace po = boost::program_options;

typedef struct {
  std::string clock_source, time_source, antenna, args, config_file, output_file;
  std::string type, subdev, wirefmt, signal_type, channel;
  size_t total_num_samps, spb;
  int duration, gain;
  double freq, rate, lo_offset, bw, setup_time;
  bool bw_summary, stats, null, enable_size_map, continue_on_bad_packet, int_n, auto_gain_tuner;
} Settings;

void parse_config(Settings* conf_setting) {
  std::vector<std::string> args;

  std::cout << conf_setting->config_file << std::endl;
  int fdr = open(conf_setting->config_file.c_str(), O_RDONLY);
  if(fdr >= 0) {
    io::file_descriptor_source fdDevice(fdr, io::file_descriptor_flags::close_handle);
    io::stream <io::file_descriptor_source> in(fdDevice);
    if(fdDevice.is_open()) {
      std::string in_line;
      while(std::getline(in, in_line)){
        args.push_back(in_line);
      }
      fdDevice.close();
    }
  }

  auto it = find(args.begin(), args.end(), ""); //Delete the space between config and signal type
  args.erase(it);

  // for (unsigned i = 0; i < args.size(); ++i) {
  //   std::cout << args.at(i) << std::endl;
  // }

  for (unsigned char i = 0; i < args.size(); ++i) {
    size_t found = args.at(i).find('[');
    if(found != std::string::npos) {
      args.at(i).erase(args.at(i).begin());
      args.at(i).erase(args.at(i).end() - 1);
    }
    else {
      args.at(i) = args.at(i).substr(args.at(i).find('=') + 2, args.at(i).size()-1);
    }
  }

  // for (unsigned i = 0; i < args.size(); ++i) {
  //   std::cout << args.at(i) << std::endl;
  // }

  conf_setting->clock_source = args.at(1);
  conf_setting->time_source = args.at(2);
  conf_setting->duration = std::stod(args.at(3));
  conf_setting->signal_type = args.at(4);
  conf_setting->freq = std::stof(args.at(5));
  conf_setting->rate = std::stof(args.at(6));
  conf_setting->channel = args.at(7);
  conf_setting->antenna = args.at(8);

  // if(!conf_setting->auto_gain_tuner) {
  //   conf_setting->gain = std::stod(args.at(9));
  // }
}

void parse_command_line(int argc, char* argv[], Settings* conf_setting) {
  //setup_program options
  po::options_description desc("Allowed options");

  //clang format-off
  desc.add_options()
    ("help", "help message")
    ("args", po::value<std::string>(&conf_setting->args)->default_value(""), "multi uhd device address args")
    ("config_file,cfg", po::value<std::string>(&conf_setting->config_file)->default_value("config.ini"), "Configuration File")
    ("file", po::value<std::string>(&conf_setting->output_file)->default_value("usrp_samples.dat"), "name of the file to write binary samples to")
    ("type", po::value<std::string>(&conf_setting->type)->default_value("short"), "sample type: double, float, or short")
    ("nsamps", po::value<size_t>(&conf_setting->total_num_samps)->default_value(1000), "total number of samples to receive")
    ("gain", po::value<int>(&conf_setting->gain)->default_value(10), "Antenna Gain")
    ("spb", po::value<size_t>(&conf_setting->spb)->default_value(10000), "samples per buffer")
    ("lo-offset", po::value<double>(&conf_setting->lo_offset)->default_value(0.0),
        "Offset for frontend LO in Hz (optional)")
    ("subdev", po::value<std::string>(&conf_setting->subdev), "subdevice specification")
    ("bw", po::value<double>(&conf_setting->bw), "analog frontend filter bandwidth in Hz")
    ("wirefmt", po::value<std::string>(&conf_setting->wirefmt)->default_value("sc16"), "wire format (sc8, sc16 or s16)")
    ("setup", po::value<double>(&conf_setting->setup_time)->default_value(1.0), "seconds of setup time")
    ("progress", "periodically display short-term bandwidth")
    ("stats", "show average bandwidth on exit")
    ("sizemap", "track packet size and display breakdown on exit")
    ("null", "run without writing to file")
    ("continue", "don't abort on a bad packet")
    ("skip-lo", "skip checking LO lock status")
    ("int-n", "tune USRP with integer-N tuning")
    ("auto_gain_tuner", "Automatic Gain Tuning Set")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  //print the help message
  if(vm.count("help")) {
    std::cout << boost::format("UHD RX samples to file %s") % desc << std::endl;
    std::cout << std::endl
              << "This application streams data from a single channel of a USRP "
                 "device to a file.\n"
              << std::endl;
    abort();
  }

  conf_setting->bw_summary             = vm.count("progress") > 0;
  conf_setting->stats                  = vm.count("stats") > 0;
  conf_setting->null                   = vm.count("null") > 0;
  conf_setting->enable_size_map        = vm.count("sizemap") > 0;
  conf_setting->continue_on_bad_packet = vm.count("continue") > 0;
  conf_setting->int_n                  = vm.count("int-n") > 0;
  conf_setting->auto_gain_tuner        = vm.count("auto_gain_tuner") > 0;
}

// int main(int argc, char* argv[]) {
//   std::string conf_file = "config.ini";
//   Settings conf_setting;
//
//   parse_command_line(argc, argv, &conf_setting);
//   parse_config(&conf_setting);
//
//   std::cout << "Output File..." << conf_setting.output_file << std::endl;
//   std::cout << "Config File..." << conf_setting.config_file << std::endl;
//   std::cout << "Clock Source..." << conf_setting.clock_source << std::endl;
//   std::cout << "Time Source..." << conf_setting.time_source << std::endl;
//   std::cout << "Duration..." << conf_setting.duration << std::endl;
//   std::cout << "Gain..." << conf_setting.gain << std::endl;
//   std::cout << "Freq..." << conf_setting.freq << std::endl;
//   std::cout << "Rate..." << conf_setting.rate << std::endl;
//   std::cout << "Channel..." << conf_setting.channel << std::endl;
//   std::cout << "Antenna..." << conf_setting.antenna << std::endl;
//   std::cout << "Signal Type..."  << conf_setting.signal_type << std::endl;
//
//   return 0;
// }
