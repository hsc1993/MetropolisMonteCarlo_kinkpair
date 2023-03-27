# MetropolisMonteCarlo_kinkpair
This repository uses Metropolis Monte Carlo algorithm to capture metastable state for kink-pair of BCC tungsten

## Table of Contents

- [Install](#install)
- [Usage](#usage)
- [Examples](#example)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)


## Install

The main function for simulation of this project is written in [C++](https://cplusplus.com/) and the output files are processed with [MATLAB](https://www.mathworks.com/products/matlab.html) scripts to analyze the results. 

To compile, use the following:
g++ Main.cpp simulation/Simulation.cpp dislocation/Dislocation.cpp SegmentStressCal/SegmentStressCal.cpp -o main -std=c++11

Note that the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library is used for calculation. The installation procedure of this library is available in its webpage.




## Usage

### Dislocation extraction


## Example



## Maintainers

[@AlanHe](https://github.com/hsc1993).

## Contributing


### Contributors

This project is supported by research group at UCLA, Johns Hopkins University, Hongkong City University and Pennsylvania University.


## License

[MIT](LICENSE) © Sicong He
