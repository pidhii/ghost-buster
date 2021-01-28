#ifndef BOSON_3BALLS_HPP
#define BOSON_3BALLS_HPP

#include "boson/states.hpp"
#include "boson/permutations.hpp"

namespace bos {

void
generate_states_3(int n_cells, state_set &set)
{
  for (int i = 0; i < n_cells; ++i)
  {
    for (int j = 0; j < n_cells; ++j)
    {
      for (int k = 0; k < n_cells; ++k)
      {
        std::vector<int> mult (n_cells, 0);
        mult[i]++;
        mult[j]++;
        mult[k]++;
        set.emplace(mult);
      }
    }
  }
}

class three_balls : public cellular_universe {
  public:
  three_balls(int n_cells)
  : cellular_universe {3}
  { generate_states_3(n_cells, m_S); }

  const state_set&
  states() const override
  { return m_S; }

  protected:
  void
  list_permutations(const state &X, state_set &F) const override
  { pick_permutations(m_S, X, F); }

  private:
  state_set m_S;
};

} // namespace bos

#endif
