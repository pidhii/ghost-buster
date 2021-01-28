#ifndef BOSON_STATES_HPP
#define BOSON_STATES_HPP

#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <ostream>
#include <algorithm>
#include <iostream>

#include "dihedral.hpp"

namespace bos {

class state {
  public:
  state() = default;

  state(const std::vector<int>& mult)
  {
    m_mult.resize(mult.size());
    std::copy(mult.begin(), mult.end(), m_mult.begin());
  }

  template <size_t N, typename T>
  state(const T mult[N])
  {
    m_mult.resize(mult.size());
    std::copy(mult, mult + N, m_mult.begin());
  }

  bool
  operator == (const state &other) const noexcept
  { return m_mult == other.m_mult; }

  bool
  operator < (const state &other) const noexcept;

  unsigned __attribute__((pure))
  operator [] (size_t i) const noexcept
  { return m_mult[i]; }

  state
  operator () (const dihedral &g) const noexcept
  {
    std::vector<int> mult (n_cells());
    for (int i = 0; i < n_cells(); ++i)
      mult[g(i)] = m_mult[i];
    return {mult};
  }

  size_t __attribute__((pure))
  n_cells() const noexcept
  { return m_mult.size(); }

  private:
  std::vector<unsigned> m_mult;
};
}

namespace std {
template <>
struct hash<bos::state> {
  size_t
  operator () (const bos::state &s) const noexcept
  {
    size_t acc = 0;
    for (int i = 0; i < s.n_cells(); ++i)
      acc = acc*10 + s[i];
    return acc;
  }
};
} // namespace std

namespace bos {

using state_set = std::unordered_set<state>;

class cellular_universe {
  public:
  cellular_universe(int n_balls) : m_nballs {n_balls} { }

  virtual const state_set&
  states() const = 0;

  int
  n_balls() const noexcept
  { return m_nballs; }

  int
  n_cells() const noexcept
  { return ordered_states()[0].n_cells(); }

  const std::vector<std::pair<int, int>>&
  multiplicity_pairs() const noexcept;

  const state_set&
  permutations(const state &x) const noexcept;

  const state_set&
  permutations(int xidx) const noexcept
  { return permutations(ordered_states()[xidx]); }

  const std::vector<state>&
  ordered_states() const noexcept;

  const std::unordered_map<state, int>&
  index_map() const noexcept;

  struct mult_pair {
    int first, second;
    mult_pair(int _1, int _2) : first {_1}, second {_2} { }
    mult_pair(const std::pair<int, int> &pair)
    : first {pair.first}, second {pair.second} { }
    bool operator < (mult_pair other) const noexcept;
    bool operator == (mult_pair other) const noexcept
    { return first == other.first and second == other.second; }
  };

  struct cell_pair {
    int first, second;
    cell_pair(int _1, int _2) : first {_1}, second {_2} { }
    cell_pair(const std::pair<int, int> &pair)
    : first {pair.first}, second {pair.second} { }
    bool operator < (cell_pair other) const noexcept;
    bool operator == (cell_pair other) const noexcept
    { return first == other.first and second == other.second; }
  };

  const std::vector<int>&
  final_states_with_multiplicities(int xidx, cell_pair cells, mult_pair mults)
    const noexcept;

  int
  state_index(const state &x) const noexcept
  { return index_map().at(x); }

  protected:
  virtual void
  list_permutations(const state &X, state_set &F) const = 0;

  typedef std::vector<std::unordered_map<mult_pair, std::vector<int>>>
    mult_table;

  private:
  int m_nballs;
  mutable std::unordered_map<state, state_set> m_permmap;
  mutable std::vector<state> m_vec;
  mutable std::unordered_map<state, int> m_idxmap;
  mutable std::vector<std::pair<int, int>> m_multpairs;
  mutable std::map<cell_pair, mult_table> m_multmap;
};

} // namespace bos

namespace std {

  template <>
  struct hash<bos::cellular_universe::mult_pair> {
    hash() = default;
    size_t
    operator () (const bos::cellular_universe::mult_pair &pair) const noexcept
    { return pair.first * 10 + pair.second; }
  }; // struct std::hash<bos::cellular_universe::mult_pair>

  template <>
  struct hash<bos::cellular_universe::cell_pair> {
    hash() = default;
    size_t
    operator () (const bos::cellular_universe::cell_pair &pair) const noexcept
    { return pair.first * 10 + pair.second; }
  }; // struct std::hash<bos::cellular_universe::cell_pair>

} // namespace std

inline bool
bos::cellular_universe::mult_pair::operator < (mult_pair other) const noexcept
{
  std::hash<mult_pair> hash;
  return hash(*this) < hash(other);
}

inline bool
bos::cellular_universe::cell_pair::operator < (cell_pair other) const noexcept
{
  std::hash<cell_pair> hash;
  return hash(*this) < hash(other);
}

inline std::ostream&
operator << (std::ostream &os, const bos::state &x)
{
  for (int i = 0; i < x.n_cells(); ++i)
    os << x[i];
  return os;
}

inline bool
bos::state::operator < (const bos::state &other) const noexcept
{
  std::hash<bos::state> hash;
  return hash(*this) < hash(other);
}



#endif
