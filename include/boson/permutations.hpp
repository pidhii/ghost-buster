#ifndef BOSON_PERMUTATIONS_HPP
#define BOSON_PERMUTATIONS_HPP

#include "boson/states.hpp"

#include <vector>
#include <cassert>
#include <variant>

namespace bos {

inline int
count_balls(const state &s, int start, int end)
{
  int acc = 0;
  for (int i = start; i != end; i = (i + 1) % s.n_cells())
    acc += s[i];
  return acc;
}

bool __attribute__((pure))
are_adjacent(const state &X, const state &Y) noexcept;

inline void
pick_permutations(const state_set &S, const state &X, state_set &F)
{
  for (const auto &Y : S)
  {
    if (are_adjacent(X, Y))
      F.insert(Y);
  }
}

typedef std::vector<std::vector<std::pair<int, dihedral>>> partition;
// First element in each class has identity dihedral.
void
spatial_quotient(const cellular_universe &model, partition &partition);

} // namespace bos

#endif
