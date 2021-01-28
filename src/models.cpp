#include "boson/models.hpp"

#include "boson/2balls.hpp"
#include "boson/3balls.hpp"
#include "boson/4balls.hpp"


const bos::cellular_universe&
bos::select_model(unsigned nballs)
{
  static const bos::two_balls N2model {6};
  static const bos::three_balls N3model {6};
  static const bos::four_balls N4model {6};
  switch (nballs)
  {
    case 2: return N2model;
    case 3: return N3model;
    case 4: return N4model;
    default: throw std::logic_error {"undefined number of balls"};
  }
}

// namespace

