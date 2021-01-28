#ifndef BOSON_MODEL_IO_HPP
#define BOSON_MODEL_IO_HPP

#include "minuit.hpp"

#include <ostream>

namespace bos {
namespace io {
static void
write_parameters(std::ostream &out, int nparm)
{
  std::string name;
  double val, err, bnd1, bnd2;
  for (int i = 0; i < nparm; ++i)
  {
    mnpout(i+1, name, val, err, bnd1, bnd2);
    out
      << name << '\t'
      << std::hexfloat
      << val << '\t' << err << '\t' << bnd1 << '\t' << bnd2
      << std::endl;
  }
}

} // namespace bos::io
} // namespace bos

#endif
