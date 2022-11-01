#ifndef __ITEMSTORE__
#define __ITEMSTORE__

#include <vector>
#include <complex>

class ItemStore {
private:
  int buff_size;
  size_t last_idx = 0;
  std::vector<std::complex<short>> items;
  // size_t curr_last_idx = 0;
  // std::vector<std::complex<short>> curr_items;

public:
  ItemStore(int buff_size);
  void add(std::vector<std::complex<short>>& current_items);
  void get_all(bool produced = false);
  size_t get_curr_last_idx();
  std::vector<std::complex<short>> get_curr_items();
  std::vector<std::complex<short>> curr_items;
  int curr_last_idx = 0;
};

ItemStore::ItemStore(int buff_size) {
  this->buff_size = buff_size;
  for (unsigned int i = 0; i < this->buff_size; ++i) {items.push_back(std::complex<short>(0.0, 0.0));}
  //FIX ME
  for (unsigned int i = 0; i < this->buff_size/2; ++i) {curr_items.push_back(std::complex<short>(0.0, 0.0));}
}

void ItemStore::add(std::vector<std::complex<short>>& current_items) {
  int num_input_samps = current_items.size();
  int new_last_idx = last_idx + num_input_samps;
  std::copy(current_items.begin(), current_items.end(), items.begin() + last_idx);
  last_idx = new_last_idx;
  // std::cout << "last_idx..." << last_idx << std::endl;
}

void ItemStore::get_all(bool produced) {
  curr_last_idx = last_idx;
  // std::cout << "last_idx..." << last_idx << std::endl;
  std::copy(items.begin(), items.begin() + last_idx, curr_items.begin());
  last_idx = 0;
}

size_t ItemStore::get_curr_last_idx() {
  return curr_last_idx;
}

std::vector<std::complex<short>> ItemStore::get_curr_items() {
  return curr_items;
}

#endif
