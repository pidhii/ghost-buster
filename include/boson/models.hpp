#ifndef BOSON_MODELS_HPP
#define BOSON_MODELS_HPP

#include "boson/states.hpp"

#include <stdexcept>

namespace bos {

const cellular_universe&
select_model(unsigned nballs);

} // namespace

#endif
