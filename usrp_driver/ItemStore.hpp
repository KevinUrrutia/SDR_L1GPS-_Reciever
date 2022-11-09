#ifndef __ITEMSTORE__
#define __ITEMSTORE__

#include <vector>
#include <complex>

class ItemStore {
private:
  int buff_size;
  int num_chans;
  size_t last_idx = 0;
  std::vector<std::vector<std::complex<short>>> items;
  // size_t curr_last_idx = 0;
  // std::vector<std::complex<short>> curr_items;

public:
  ItemStore(int buff_size, int num_chans);
  void add(std::vector<std::vector<std::complex<short>>> & current_items);
  void get_all(bool produced = false);
  size_t get_curr_last_idx();
  std::vector<std::vector<std::complex<short>>> get_curr_items();
  std::vector<std::vector<std::complex<short>>>  curr_items;
  int curr_last_idx = 0;
};

ItemStore::ItemStore(int buff_size, int num_chans=1) {
  this->buff_size = buff_size;
  this->num_chans = num_chans;

  std::complex<short> x = (0.0, 0.0);
  std::vector<std::complex<short>> xVec(buff_size, x);
  std::vector<std::vector<std::complex<short>>> xMat(num_chans, xVec);
  items = xMat;
  curr_items = xMat;

}

void ItemStore::add(std::vector<std::vector<std::complex<short>>> & current_items) {
  int in_num_chans = current_items.size();
  int num_input_samps = current_items.at(0).size();
  if (last_idx + num_input_samps > buff_size) {
    std::cout << "I.S. Full!" << std::endl;
    return;
  }

  int new_last_idx = last_idx + num_input_samps;
  for (int i = 0; i < in_num_chans; i++) {
    std::copy(current_items.at(i).begin(), current_items.at(i).end(), items.at(i).begin() + last_idx);
  }
  last_idx = new_last_idx;
}

void ItemStore::get_all(bool produced) {
  curr_last_idx = last_idx;

  for (int i = 0; i < num_chans; i++) {
    std::copy(items.at(i).begin(), items.at(i).begin() + last_idx, curr_items.at(i).begin());
  }
  last_idx = 0;
}

size_t ItemStore::get_curr_last_idx() {
  return curr_last_idx;
}

std::vector<std::vector<std::complex<short>>> ItemStore::get_curr_items() {
  return curr_items;
}

#endif
