#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <sstream>
#include <numeric>
#include "rdsSupport.h"

using namespace std;

// Parity matrix for RDS
vector<vector<int>> p = {
    {1,0,0,0,0,0,0,0,0,0},
    {0,1,0,0,0,0,0,0,0,0},
    {0,0,1,0,0,0,0,0,0,0},
    {0,0,0,1,0,0,0,0,0,0},
    {0,0,0,0,1,0,0,0,0,0},
    {0,0,0,0,0,1,0,0,0,0},
    {0,0,0,0,0,0,1,0,0,0},
    {0,0,0,0,0,0,0,1,0,0},
    {0,0,0,0,0,0,0,0,1,0},
    {0,0,0,0,0,0,0,0,0,1},
    {1,0,1,1,0,1,1,1,0,0},
    {0,1,0,1,1,0,1,1,1,0},
    {0,0,1,0,1,1,0,1,1,1},
    {1,0,1,0,0,0,0,1,1,1},
    {1,1,1,0,0,1,1,1,1,1},
    {1,1,0,0,0,1,0,0,1,1},
    {1,1,0,1,0,1,0,1,0,1},
    {1,1,0,1,1,1,0,1,1,0},
    {0,1,1,0,1,1,1,0,1,1},
    {1,0,0,0,0,0,0,0,0,1},
    {1,1,1,1,0,1,1,1,0,0},
    {0,1,1,1,1,0,1,1,1,0},
    {0,0,1,1,1,1,0,1,1,1},
    {1,0,1,0,1,0,0,1,1,1},
    {1,1,1,0,0,0,1,1,1,1},
    {1,1,0,0,0,1,1,0,1,1}
};

vector<int> a = {1,1,1,1,0,1,1,0,0,0};
vector<int> b = {1,1,1,1,0,1,0,1,0,0};
vector<int> c = {1,0,0,1,0,1,1,1,0,0};
vector<int> cPrime = {1,1,1,1,0,0,1,1,0,0};
vector<int> d = {1,0,0,1,0,1,1,0,0,0};

vector<bool> differential_decode(vector<real>& pre_cdr, int sps, int& symbol_state, int& bit_state, int& symb_count) {
    vector<bool> decoded_bitstream, symbols, manchester_decode;

    if (symbol_state != -1) {
        symbols.push_back(symbol_state);
    }

    while (symb_count < pre_cdr.size()) {
		if (pre_cdr[symb_count] > 0){
			symbols.push_back(true);
		} else {
			symbols.push_back(false);
		}
        symb_count += sps;
    }

    if (symbols.size() % 2 != 0) {
        symbol_state = (bool)symbols.back();
    } else {
        symbol_state = -1;
    }

    for (int i = 0; i < symbols.size() - 1; i += 2) {
		if (symbols[i] == true) {
			symbols.push_back(true);
		} else {
			symbols.push_back(false);
		}
    }

    if (!manchester_decode.empty()) {
        decoded_bitstream.push_back(bit_state ^ manchester_decode[0]);
        for (int i = 1; i < manchester_decode.size(); i++) {
            decoded_bitstream.push_back(manchester_decode[i] ^ manchester_decode[i - 1]);
        }
        bit_state = manchester_decode.back();
    }

    return decoded_bitstream;
}

int frame_verification(const vector<bool>& bit_stream) {
    vector<bool> check_word(10, false);

    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 26; j++) {
            check_word[i] = check_word[i] ^ (bit_stream[j] & p[j][i]);
        }
    }

    if (check_word == vector<bool>(a.begin(), a.end())) { // Compare to a
        return 1;
    } else if (check_word == vector<bool>(b.begin(), b.end())) { // Compare to b
        return 2;
    } else if (check_word == vector<bool>(c.begin(), c.end())) { // Compare to c
        return 3;
    } else if (check_word == vector<bool>(cPrime.begin(), cPrime.end())) { // Compare to cPrime
        return 4;
    } else if (check_word == vector<bool>(d.begin(), d.end())) { // Compare to d
        return 5;
    } else {
        return 0;
    }
    
}

string string_conversion( vector<int> &bits, unsigned short int block_type){
    string pi;
    for(int i = 0; i < 16; i++){
        pi += (bits[i] ? '1' : '0');
    }

    std::string pty = pi.substr(6, 5);

    std::map<std::string, std::string> pty_dict = {
        {"00000", "No Program"}, {"00001", "News"}, {"00010", "Information"}, {"00011", "Sports"},
        {"00100", "Talk"}, {"00101", "Rock"}, {"00110", "Classic Rock"}, {"00111", "Adult Hits"},
        {"01000", "Soft Rock"}, {"01001", "Top 40"}, {"01010", "Country"}, {"01011", "Oldies"},
        {"01100", "Soft"}, {"01101", "Nostalgia"}, {"01110", "Jazz"}, {"01111", "Classical"},
        {"10000", "Rhythm & Blues"}, {"10001", "Soft Rhythm & Blues"}, {"10010", "Foreign Language"},
        {"10011", "Religious Music"}, {"10100", "Religious Talk"}, {"10101", "Personality"},
        {"10110", "Public"}, {"10111", "College"}, {"11101", "Weather"}, {"11110", "Emergency Test"},
        {"11111", "Emergency"}
    };

    if(block_type == 1){
        // Convert binary string to integer, then to hexadecimal string
        stringstream bitStream;
        bitStream << hex << stoi(pi, nullptr, 2);
        return bitStream.str();
    } else if (block_type == 2){
        auto it = pty_dict.find(pty);
        if (it != pty_dict.end()) {
            return it->second; // Return the program type string
        } else {
            return "Unknown"; // Handle the case where the program type is not found
        }
    }
    return ""; // Add a default return statement
}

string character_conversion(const vector<bool>& binary) {
    int decimal = 0;
    for (int i = 0; i < binary.size(); ++i) {
        decimal += binary[binary.size() - 1 - i] * std::pow(2 , i);
    }
    return std::string(1, (char)(decimal));
}

void frame_sync(vector<bool> &bit_stream, vector<bool> &slot_state, bool &synced, string &p_service,
                int &decodeIdentifier_prev, int &decodeIdentifier,
                vector<string> &ps_name_segments, string &group_type, bool &ps_ready) {

    int block_pos = 0;
    vector<bool> large_bitstream(slot_state);
    large_bitstream.insert(large_bitstream.end(), bit_stream.begin(), bit_stream.end());

    ps_ready = false; // Initialize ps_ready to false at the beginning of each call

    while (block_pos + 26 < large_bitstream.size()) {
        vector<bool> slot(large_bitstream.begin() + block_pos, large_bitstream.begin() + block_pos + 26);
        int block_type = frame_verification(slot);

        if (block_type == 0) {
            block_pos += 1;
            synced = false;
        } else {
            synced = true;
            block_pos += 26;

            switch (block_type) {
                case 1:{
                    std::cerr << "a detected" << std::endl;
                    std::vector<int> slot_int(slot.begin(), slot.end());
                    std::cerr << "PI Code: " << string_conversion(slot_int, block_type) << std::endl;
                    break;
                }
                case 2: {
                    std::cerr << "b detected" << std::endl;
                    std::vector<int> slot_int(slot.begin(), slot.end());
                    std::cerr << "Program Type: " << string_conversion(slot_int, block_type) << std::endl;
                
                    decodeIdentifier = static_cast<int>(slot[14]) * 2 + static_cast<int>(slot[15]);

                    // Extract group type (5 bits)
                    stringstream gt_ss;
                    for (int i = 0; i < 5; ++i) {
                        gt_ss << (slot[i] ? '1' : '0');
                    }
                    group_type = gt_ss.str();
                    std::cerr << "Group Type: " << group_type << std::endl;

                    // Check if group_type is not in ["00000", "00001"] and any ps_name_segments are filled
                    bool any_segments_filled = false;
                    for (const auto& segment : ps_name_segments) {
                        if (!segment.empty()) {
                            any_segments_filled = true;
                            break;
                        }
                    }

                    if (group_type != "00000" && group_type != "00001" && any_segments_filled) {
                        bool all_segments_filled = true;
                        for (const auto& segment : ps_name_segments) {
                            if (segment.empty()) {
                                all_segments_filled = false;
                                break;
                            }
                        }

                        if (all_segments_filled) {
                            p_service = "";
                            for (const auto& segment : ps_name_segments) {
                                p_service += segment;
                            }
                            ps_ready = true;
                        }
                    }
                    break;
                }
                case 3:
                    std::cerr << "c detected" << std::endl;
                    break;
                case 4:
                    std::cerr << "cPrime detected" << std::endl;
                    break;
                case 5: {
                    std::cerr << "d detected" << std::endl;
                    if ((group_type == "00000" || group_type == "00001") && decodeIdentifier >= 0 && decodeIdentifier <= 3) {
                        vector<bool> c1(slot.begin(), slot.begin() + 8);
                        vector<bool> c2(slot.begin() + 8, slot.begin() + 16);
                        string char1 = character_conversion(c1);
                        string char2 = character_conversion(c2);
                        ps_name_segments[decodeIdentifier] = char1 + char2;

                        bool all_segments_filled = true;
                        for (const auto& segment : ps_name_segments) {
                            if (segment.empty()) {
                                all_segments_filled = false;
                                break;
                            }
                        }

                        if (all_segments_filled && p_service != accumulate(ps_name_segments.begin(), ps_name_segments.end(), string())) {
                            p_service = accumulate(ps_name_segments.begin(), ps_name_segments.end(), string());
                            ps_ready = true;
                        }
                    }
                    decodeIdentifier_prev = decodeIdentifier;
                    break;
                }
                default:
                    break;
            }
        }
        slot_state.assign(large_bitstream.begin() + block_pos, large_bitstream.end());
    }
}