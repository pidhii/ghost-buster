#ifndef PERM_MATRIX_HPP
#define PERM_MATRIX_HPP

#include <vector>

class perm_matrix {
  public:
  perm_matrix(int dim, const std::vector<std::pair<int, int>> &perms)
  : m_rows (dim)
  {
    // make identity
    for (int i = 0; i < dim; ++i)
      m_rows[i] = i;
    // implement permutations (`perms` are 2-cycle-decomposition of permutation)
    for (const auto [i, f] : perms)
    {
      m_rows[f] = i;
      m_rows[i] = f;
    }
  }

  class row {
    double
    operator [] (int j) const noexcept
    { return !!(m_nonzeroidx == j); }

    int
    nonzero_index() const noexcept
    { return m_nonzeroidx; }

    private:
    friend class perm_matrix;
    row(int j) : m_nonzeroidx {j} { }

    int m_nonzeroidx;
  };

  row
  operator [] (int i) const noexcept
  { return row {i}; }

  perm_matrix
  transpose() const
  {
    std::vector<int> transp (m_rows.size());
    for (int i = 0; i < m_rows.size(); ++i)
      transp[m_rows[i]] = i;
    return {std::move(transp)};
  }

  private:
  perm_matrix(const std::vector<int> &rows) : m_rows {rows} { }

  std::vector<int> m_rows;
};


#endif
