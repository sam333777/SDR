

#ifndef RDSSUPPORT_H
#define RDSSUPPORT_H

#include "dy4.h"
#include <iostream>
#include <vector>
#include <string>

// Function prototypes

// Differential decoding function
std::vector<bool> differential_decode(std::vector<real>& pre_cdr, int sps, int& symbol_state, int& bit_state, int& symb_count);

// Frame verification function
int frame_verification(const std::vector<bool>& bit_stream);

// String conversion function
std::string string_conversion(std::vector<int> &bits, unsigned short int block_type);

// Character conversion function
std::string character_conversion(const std::vector<bool>& binary);

// Frame synchronization function
void frame_sync (std::vector<bool> &bit_stream, std::vector<bool> &slot_state, bool &synced, std::string &p_service, int &decodeIdentifier_prev, int &decodeIdentifier, std::vector<std::string> &ps_name_segments, std::string &group_type, bool &ps_ready);
#endif // RDSSUPPORT_H
