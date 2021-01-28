#ifndef BOSON_4BALLS_HPP
#define BOSON_4BALLS_HPP

#include "boson/states.hpp"
#include "boson/permutations.hpp"

namespace bos {

void
generate_states_4(int n_cells, state_set &set)
{
  for (int i = 0; i < n_cells; ++i)
  {
    for (int j = 0; j < n_cells; ++j)
    {
      for (int k = 0; k < n_cells; ++k)
      {
        for (int l = 0; l < n_cells; ++l)
        {
          std::vector<int> mult (n_cells, 0);
          mult[i]++;
          mult[j]++;
          mult[k]++;
          mult[l]++;
          set.emplace(mult);
        }
      }
    }
  }
}

class four_balls : public cellular_universe {
  public:
  four_balls(int n_cells)
  : cellular_universe {4}
  { generate_states_4(n_cells, m_S); }

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
