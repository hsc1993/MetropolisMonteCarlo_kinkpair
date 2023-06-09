#ifndef __STDINC_H
#define __STDINC_H

/////////////////////////////////////////
// Standard Template Library (STL)

#include <vector>
#include <deque>
#include <list>
#include <limits>
#include <ostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <stack>
#include <complex>
#include <stdexcept>
#include <math.h>
#include <algorithm>
#include <memory>
#include <functional>
#include <cstdio>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <random>

/////////////////////////////////////////
// Standard C Library

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <locale.h>
#include <time.h>

/////////////////////////////////////////
// Boost Library

/*
#include <boost/array.hpp>
#include <boost/utility.hpp>
#include <boost/program_options.hpp>
#include <boost/pool/object_pool.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/exponential_distribution.hpp>
*/

/////////////////////////////////////////
// Our own headers

#include "Settings.h"
#include "util/Debug.h"
#include "util/linalg/Vector3.h"
#include "util/linalg/Point3.h"
#include "util/linalg/Matrix3.h"
#include "util/linalg/Tensor.h"

// Anisotropic stress calculation headers




#include <iomanip>
#include "Anisotropic_headers/ElasticL.h"
#include "Anisotropic_headers/Tensor.h"
#include "Anisotropic_headers/GreenTensor.h"
#include "Anisotropic_headers/GreenTensor0.h"
#include "Anisotropic_headers/readXY.h"


//#include "SegmentStressCal.h"


/*
#include "Anisotropic_headers/readCtensor.h"    // incorporated in Main.cpp 
#include "Anisotropic_headers/leviCivita.h"     // incorporated in SegmentStressCal.h
#include "Anisotropic_headers/cbedl.h"          // incorporated in SegmentStressCal.h
#include "Anisotropic_headers/readDtensor.h"    // no need to include 
#include "Anisotropic_headers/readBurgers.h"    
*/

/*
Giacomo's main.cpp is broken down into SegmentStressClassical.h and parts in Alan's Main.cpp 
*/

#endif	// __STDINC_H
