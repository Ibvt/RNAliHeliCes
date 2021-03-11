# RNAliHeliCes and RNAliHiPath
Abstract folding space analysis based on helices for aligned RNAs

Copyright (C) 2012-21 Jiabin Huang, Bjoern Voss.
Send comments/bug reports to: J.Huang <j.huang@uke.de>.
Updates: https://github.com/Ibvt/RNAliHeliCes


## Dependencies
1. (Optional) if without SUPER privilege, install libboost locally
```sh
tar xzf boost_1_58_0.tar.gz
#under python2
./bootstrap.sh --prefix=/your_path/boost_1_58_installed/
./b2 install
```

## Installation
```sh
git clone https://github.com/huang/RNAHeliCes
./configure CFLAGS="-fno-stack-protector" CPPFLAGS="-std=c++98" CXXFLAGS="-std=c++98 -fno-stack-protector"
make
sudo make install
```
Notes:
  - If there are any linking problems after installing, please check if /usr/local/lib is contained in the environmental variable LD_LIBRARY_PATH. 
  - src/libs/libRNA.a contains all Vienna package routines, this may be system dependent. If so, please reload this from Vienna RNA Package.
  - If without SUPER privilege, use 
  ./configure --with-boost-include-path=/your_path/boost_1_58_installed/include --with-boost-lib-path=/your_path/boost_1_58_installed/lib CFLAGS="-g -O2 -fno-stack-protector" CPPFLAGS="-std=c++98 -I/your_path/boost_1_58_installed/include -DNDEBUG -D_FORTIFY_SOURCE=2 -O2" CXXFLAGS="-std=c++98 -g -O2 -w -fno-stack-protector"

## Test run
```sh
RNAliHeliCes -f examples/test2.aln -k 10 -R 36.5m,41.5,27 -e -t 2 -x0
RNAliHiPath -f examples/test2.aln -k 10 -F examples/test2.ss
```

## INSTALL and COPYING
See "INSTALL"        for detailed installation instructions, and
    "COPYING"        for disclaimer and copyright.
    
## Citations
  [1] Huang J, Voß B. Simulation of Folding Kinetics for Aligned RNAs. Genes (Basel). 2021 Feb 26;12(3):347. doi: 10.3390/genes12030347. PMID: 33652983.
