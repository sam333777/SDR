# SDR

# Real-Time FM Receiver (SDR)

This project implements a real-time software-defined radio (SDR) system for receiving FM broadcasts.  
It supports mono and stereo audio decoding, as well as extracting digital data using the Radio Data System (RDS).

---
## Features

- Works with RTL-SDR dongles and Raspberry Pi. 
- FM demodulation with separate processing paths for:
  - **Mono audio** (0–15 kHz)  
  - **Stereo audio** using 19 kHz pilot tone and 38 kHz subcarrier  
  - **RDS data** (57 kHz subchannel)  
- Real-time operation with C++ backend (multi-threaded).  
- Python models provided for filter design and testing. 
