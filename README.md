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

"g++ Main.cpp simulation/Simulation.cpp dislocation/Dislocation.cpp SegmentStressCal/SegmentStressCal.cpp -o main -std=c++11"

Note that the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library is used for calculation. The installation procedure of this library is available in its webpage.


## Usage

Once the compilation is done, simply run the excutable with "./main"

### Metropolis Algorithm
[Simulation.cpp](simulation/Simulation.cpp) file is the where all the simulation parameters and materials constants are initialized. The code will run under a series of heights for kink-pairs, and a series of input stresses. Under each condition, a Metropolis Monte Carlo algorithm will tune the position of nodes (in 2D) on a dislocation line, trying to find the optimal shape of the line to minimize the system energy. Once the convergence criteria is met, the saddle point configuration (metastable state) will be written into 'saddle_config.txt' file.

### Data processing
[evolution.m](evolution.m) plots the convergence of each energy term of the system. And [saddleConfig.m](saddleConfig.m) plots the saddle point configuration. 

## Example
<img src="https://user-images.githubusercontent.com/56003395/228091507-fa6da8a9-1a87-474b-9746-253e77a7066f.png" width="100" height="100">

This is the saddle point configuration for zero stress applied to the system, and the dislocation line curves itself to reach lower energy state.


## Maintainers

[@AlanHe](https://github.com/hsc1993).

## Contributing


### Contributors

This project is supported by Marian group at UCLA.


## License

[MIT](LICENSE) © Sicong He
