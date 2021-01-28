#include "boson/permutations.hpp"

bool
bos::are_adjacent(const state &X, const state &Y) noexcept
{
  assert(X.n_cells() == Y.n_cells());
  const int V = X.n_cells();

  bool XY = true, YX = true;
  for (int dc = 0; dc < V - 2; ++dc)
  {
    for (int c = 0; c < V; ++c)
    {
      {
        const int nX = count_balls(X, c, (c + dc) % V);
        const int nY = count_balls(Y, (c + V-1) % V, (c + dc + 1) % V);
        XY = (XY and nY >= nX);
      }
      {
        const int nY = count_balls(Y, c, (c + dc) % V);
        const int nX = count_balls(X, (c + V-1) % V, (c + dc + 1) % V);
        YX = (YX and nX >= nY);
      }
    }
  }
  return XY and YX;
}

void
bos::spatial_quotient(const cellular_universe &model, partition &partition)
{
  const std::vector<state> &states = model.ordered_states();
  const int ncells = states[0].n_cells();
  const int nstates = states.size();
  const dihedral e {ncells};
  const dihedral s = e * dihedral::s;

  std::vector<bool> memory (nstates, false);
  std::unordered_map<int, dihedral> buf;
  for (int xidx = 0; xidx < nstates; ++xidx)
  {
    if (memory[xidx])
      continue;

    buf.clear();
    // r, r^2, r^3, ..., r^n, s, sr, sr^2, sr^3, ..., sr^n
    for (int r = 0; r < ncells; ++r)
    {
      const state &x = states[xidx];
      const dihedral gy = e*r;
      const dihedral gz = s*r;
      const int yidx = model.state_index(x(gy));
      const int zidx = model.state_index(x(gz));
      buf.emplace(yidx, gy);
      buf.emplace(zidx, gz);
      memory[yidx] = true;
      memory[zidx] = true;
    }
    partition.emplace_back(buf.begin(), buf.end());

    // sort so that first element had the identity dihedral
    std::sort(partition.back().begin(), partition.back().end(),
      [] (const auto &lhs, const auto &rhs) -> bool {
        return lhs.second.hash() < rhs.second.hash();
      });
  }
}

