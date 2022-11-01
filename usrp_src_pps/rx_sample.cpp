#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/thread.hpp>
#include <uhd/types/device_addr.hpp>
#include <uhd/types/metadata.hpp>
#include <uhd/types/stream_cmd.hpp>
#include <uhd/stream.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <string>
#include <chrono>
#include <csignal>
#include <vector>

//Custom built libraries
#include "parser.cpp"

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

int UHD_SAFE_MAIN(int argc, char* argv[]) {
  std::string args, file, type, ant, subdev, ref, wirefmt;
  size_t channel, total_num_samps, spb;
  double rate, freq, gain, bw, total_time, setup_time, lo_offset;

  //setup the program_options
  po::options_description  desc("Allowed options");
  //clang format-off
  desc.add_options()
    ("help", "help message")
    ("args", po::value<std::string>(&args)->default_value(""), "multi uhd device address args")
    ("file", po::value<std::string>(&file)->default_value("usrp_samples.dat"), "name of the file to write binary samples to")
    ("type", po::value<std::string>(&type)->default_value("short"), "sample type: double, float, or short")
    ("nsamps", po::value<size_t>(&total_num_samps)->default_value(0), "total number of samples to receive")
    ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
    ("spb", po::value<size_t>(&spb)->default_value(10000), "samples per buffer")
    ("rate", po::value<double>(&rate)->default_value(1e6), "rate of incoming samples")
    ("freq", po::value<double>(&freq)->default_value(0.0), "RF center frequency in Hz")
    ("lo-offset", po::value<double>(&lo_offset)->default_value(0.0),
        "Offset for frontend LO in Hz (optional)")
    ("gain", po::value<double>(&gain), "gain for the RF chain")
    ("ant", po::value<std::string>(&ant), "antenna selection")
    ("subdev", po::value<std::string>(&subdev), "subdevice specification")
    ("channel", po::value<size_t>(&channel)->default_value(0), "which channel to use")
    ("bw", po::value<double>(&bw), "analog frontend filter bandwidth in Hz")
    ("ref", po::value<std::string>(&ref)->default_value("internal"), "reference source (internal, external, mimo)")
    ("wirefmt", po::value<std::string>(&wirefmt)->default_value("sc16"), "wire format (sc8, sc16 or s16)")
    ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
    ("progress", "periodically display short-term bandwidth")
    ("stats", "show average bandwidth on exit")
    ("sizemap", "track packet size and display breakdown on exit")
    ("null", "run without writing to file")
    ("continue", "don't abort on a bad packet")
    ("skip-lo", "skip checking LO lock status")
    ("int-n", "tune USRP with integer-N tuning")
  ;

  //clang format-on
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
    return ~0;
  }

  bool bw_summary             = vm.count("progress") > 0;
  bool stats                  = vm.count("stats") > 0;
  bool null                   = vm.count("null") > 0;
  bool enable_size_map        = vm.count("sizemap") > 0;
  bool continue_on_bad_packet = vm.count("continue") > 0;

  if(enable_size_map)
    std::cout << "Packet size tracking enabled - will only recv one packet at a time!"
              << std::endl;

  //create a usrp device
  std::cout << std::endl;
  std::cout << boost::format("Creating the USRP device with: %s...") % args
            << std::endl;
  uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);

  //Lock mboard clocks
  if(vm.count("ref")) {
    usrp->set_clock_source(ref);
  }

  //always select the subdevice first, the channel mapping affects other setiings
  if(vm.count("subdev")) {
    usrp->set_rx_subdev_spec(subdev);
  }

  std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;

  //set the sampling rate
  if(rate <= 0.0) {
    std::cerr << "Please specify a valid sample rate" << std::endl;
    return ~0;
  }
  std::cout << boost::format("Setting RX Rate: %f Msps...") % (rate / 1e6) << std::endl;
  usrp->set_rx_rate(rate, channel);
  std::cout << boost::format("Actual RX Rate: %f Msps...")
                  % (usrp->get_rx_rate(channel) / 1e6)
            << std::endl
            << std::endl;

  //Set the Center frequency
  if(vm.count("freq")) { //with a defualt setting of 0.0 this statement will always yeild true
    std::cout << boost::format("Setting RX Freq: %f MHz...") % (freq / 1e6)
              << std::endl;
    std::cout << boost::format("Setting RX LO offset: %f MHz...") % (lo_offset / 1e6)
              << std::endl << std::endl;

    uhd::tune_request_t tune_request(freq, lo_offset);
    if(vm.count("int-n"))
      tune_request.args = uhd::device_addr_t("mode_n=integer");
    usrp->set_rx_freq(tune_request, channel);
    std::cout << boost::format("Actual RX Freq: %f MHz...")
                    % (usrp->get_rx_freq(channel) / 1e6)
              << std::endl
              << std::endl;
  }

  //set the rf gain
  if (vm.count("gain")) {
    std::cout << boost::format("Setting RX Gain: %f dB...") % gain << std::endl;
    usrp->set_rx_gain(gain, channel);
    std::cout << boost::format("Actual RX Gain: %f dB...")
                    % usrp->get_rx_gain(channel)
              << std::endl
              << std::endl;
  }

  //set the antenna
  if(vm.count("ant"))
    usrp->set_rx_antenna(ant, channel);

  std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

  return EXIT_SUCCESS;
}
