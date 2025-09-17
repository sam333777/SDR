import numpy as np
import math

# Parity matrix for RDS
p = np.array([
    [1,0,0,0,0,0,0,0,0,0],
    [0,1,0,0,0,0,0,0,0,0],
    [0,0,1,0,0,0,0,0,0,0],
    [0,0,0,1,0,0,0,0,0,0],
    [0,0,0,0,1,0,0,0,0,0],
    [0,0,0,0,0,1,0,0,0,0],
    [0,0,0,0,0,0,1,0,0,0],
    [0,0,0,0,0,0,0,1,0,0],
    [0,0,0,0,0,0,0,0,1,0],
    [0,0,0,0,0,0,0,0,0,1],
    [1,0,1,1,0,1,1,1,0,0],
    [0,1,0,1,1,0,1,1,1,0],
    [0,0,1,0,1,1,0,1,1,1],
    [1,0,1,0,0,0,0,1,1,1],
    [1,1,1,0,0,1,1,1,1,1],
    [1,1,0,0,0,1,0,0,1,1],
    [1,1,0,1,0,1,0,1,0,1],
    [1,1,0,1,1,1,0,1,1,0],
    [0,1,1,0,1,1,1,0,1,1],
    [1,0,0,0,0,0,0,0,0,1],
    [1,1,1,1,0,1,1,1,0,0],
    [0,1,1,1,1,0,1,1,1,0],
    [0,0,1,1,1,1,0,1,1,1],
    [1,0,1,0,1,0,0,1,1,1],
    [1,1,1,0,0,0,1,1,1,1],
    [1,1,0,0,0,1,1,0,1,1]
])

a = [1,1,1,1,0,1,1,0,0,0]
b = [1,1,1,1,0,1,0,1,0,0]
c = [1,0,0,1,0,1,1,1,0,0]
cPrime = [1,1,1,1,0,0,1,1,0,0]
d = [1,0,0,1,0,1,1,0,0,0]

def differential_decode(pre_cdr, sps, symbol_state, bit_state, symb_count):
    decoded_bitstream = []
    symbols = []
    manchester_decode = []

    if symbol_state != -1:
        symbols.append(int(symbol_state))

    # Extracting symbols from post RRC filtered signal
    # symb_count is used to sample every SPS samples and maintain proper sampling between blocks
    while symb_count < len(pre_cdr):
        if pre_cdr[symb_count] > 0:
            symbols.append(True)
        else:
            symbols.append(False)
        symb_count += sps

    if len(symbols) % 2 != 0:
        symbol_state = int(symbols[-1])
    else:
        symbol_state = -1

    # Manchester decoding the symbols
    for i in range(0, len(symbols) - 1, 2):
        # Pairing based off symbol transitions
        if symbols[i]:
            manchester_decode.append(True)
        else:
            manchester_decode.append(False)

    # Differential decoding by XOR
    # Last bit of Manchester decoded bitstream is state saved for next block
    if manchester_decode:
        decoded_bitstream.append(bit_state ^ manchester_decode[0])
        for i in range(1, len(manchester_decode)):
            decoded_bitstream.append(manchester_decode[i] ^ manchester_decode[i - 1])
        bit_state = manchester_decode[-1]

    #print("Decoded Bitstream:", decoded_bitstream)

    return decoded_bitstream, symbol_state, bit_state


def frame_verification(bit_stream):
    check_word = np.zeros(10, dtype=int)

    for i in range(10):
        for j in range(26):
            check_word[i] ^= int(bit_stream[j]) * p[j][i]

    if np.array_equal(check_word, a):
        result =  1
    elif np.array_equal(check_word, b):
        result =  2
    elif np.array_equal(check_word, c):
        result =  3
    elif np.array_equal(check_word, cPrime):
        result =  4
    elif np.array_equal(check_word, d):
        result =  5
    else:
        result =  0
    return result

def string_conversion(bits, block_type):
    pi = ''.join(['1' if bit else '0' for bit in bits[:16]])
    pty = ''.join(['1' if bit else '0' for bit in bits[6:11]])

    pty_dict = {
        "00000": "No Program", "00001": "News", "00010": "Information", "00011": "Sports",
        "00100": "Talk", "00101": "Rock", "00110": "Classic Rock", "00111": "Adult Hits",
        "01000": "Soft Rock", "01001": "Top 40", "01010": "Country", "01011": "Oldies",
        "01100": "Soft", "01101": "Nostalgia", "01110": "Jazz", "01111": "Classical",
        "10000": "Rhythm & Blues", "10001": "Soft Rhythm & Blues", "10010": "Foreign Language",
        "10011": "Religious Music", "10100": "Religious Talk", "10101": "Personality",
        "10110": "Public", "10111": "College", "11101": "Weather", "11110": "Emergency Test",
        "11111": "Emergency"
    }

    if block_type == 1:
        pi_hex = hex(int(pi, 2))[2:]
        return pi_hex
    elif block_type == 2:
        return pty_dict.get(pty, "Unassigned")

def character_conversion(binary):
    decimal = sum(2**i * int(binary[7-i]) for i in range(len(binary)))
    return chr(decimal)

def frame_sync(bit_stream, slot_state, synced, p_service, di_prev, di, ps_name_segments, group_type):
    block_pos = 0
    large_bitstream = np.concatenate((slot_state, bit_stream))
    ps_ready = False
    

    while block_pos + 26 < len(large_bitstream):
        slot = large_bitstream[block_pos:block_pos + 26]
        block_type = frame_verification(slot)

        if block_type == 0:
            block_pos += 1
            synced = False
        else:
            synced = True
            block_pos += 26
            if block_type == 1:
                print("a dtetcted")
                print("PI Code:", string_conversion(slot, block_type))
            elif block_type == 2:
                # Extract Decoder Identification (DI) bits from slot (bits 14 and 15)
                # DI indicates the segment number of the Program Service (PS) name being transmitted:
                # DI = 0 -> PS characters 1 & 2
                # DI = 1 -> PS characters 3 & 4
                # DI = 2 -> PS characters 5 & 6
                # DI = 3 -> PS characters 7 & 8
                print("b dtetcted")
                print("Program Type:", string_conversion(slot, block_type))
                di = int(slot[14]) * 2 + int(slot[15])
                group_type = ''.join(['1' if bit else '0' for bit in slot[0:5]])
                print("Group Type:", group_type)
                if group_type not in ["00000", "00001"] and any(ps_name_segments):
                    
                    if all(ps_name_segments):
                        p_service = ''.join(ps_name_segments)
                        ps_ready = True


            elif block_type == 3:
                print("c dtetcted")
            elif block_type == 4:
                print("cPrime dtetcted")
            elif block_type == 5:
                print("d detected")

                if group_type in ["00000", "00001"] and 0 <= di <= 3:
                    c1 = slot[:8]
                    c2 = slot[8:16]
                    # print(f"Raw bits c1: {c1}")
                    # print(f"Raw bits c2: {c2}")
                    char1 = character_conversion(c1)
                    char2 = character_conversion(c2)
                    ps_name_segments[di] = char1 + char2
                    # print(f"Segment {di} = '{char1 + char2}'")
                    # print(f"ps_name_segments = {ps_name_segments}")

                    if all(ps_name_segments) and not p_service == ''.join(ps_name_segments):
                        p_service = ''.join(ps_name_segments)
                        ps_ready = True
                di_prev = di
    slot_state = large_bitstream[block_pos:]
    return slot_state, synced, p_service, di_prev, di, ps_name_segments, group_type, ps_ready

