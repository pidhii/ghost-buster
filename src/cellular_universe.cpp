#include "boson/states.hpp"

const bos::state_set&
bos::cellular_universe::permutations(const state &x) const noexcept
{
  auto it = m_permmap.find(x);
  if (it == m_permmap.end())
  {
    it = m_permmap.emplace(x, state_set{}).first;
    list_permutations(x, it->second);
  }
  return it->second;
}

const std::vector<bos::state>&
bos::cellular_universe::ordered_states() const noexcept
{
  if (m_vec.empty())
  {
    m_vec.resize(states().size());
    std::copy(states().begin(), states().end(), m_vec.begin());
    std::sort(m_vec.begin(), m_vec.end());
  }
  return m_vec;
}

const std::unordered_map<bos::state, int>&
bos::cellular_universe::index_map() const noexcept
{
  if (m_idxmap.empty())
  {
    for (int xidx = 0; xidx < states().size(); ++xidx)
      m_idxmap[ordered_states()[xidx]] = xidx;
  }
  return m_idxmap;
}

const std::vector<std::pair<int, int>>&
bos::cellular_universe::multiplicity_pairs() const noexcept
{
  if (m_multpairs.empty())
  {
    m_multpairs.reserve(n_balls()*n_balls());
    for (int nu = 0; nu < n_balls(); ++nu)
    {
      for (int nv = 0; nv < n_balls(); ++nv)
        m_multpairs.emplace_back(nu, nv);
    }
  }
  return m_multpairs;
}

const std::vector<int>&
bos::cellular_universe::final_states_with_multiplicities(int xidx,
    cell_pair cells, mult_pair mults) const noexcept
{
  mult_table &multmap = m_multmap[cells];
  if (multmap.empty())
  { // for each X
    for (int xidx = 0; xidx < states().size(); ++xidx)
    { // partition non-vanishing transitions according to final multiplicities
      auto &xmultmap = multmap[xidx];
      for (const state &y : permutations(xidx))
      {
        const int nu = y[cells.first];
        const int nv = y[cells.second];
        const int yidx = state_index(y);
        xmultmap[{nu, nv}].push_back(yidx);
      }
    }
  }
  return multmap[xidx][mults];
}

