#ifndef BOSON_MATRIX_HPP
#define BOSON_MATRIX_HPP

#include "boson/states.hpp"
#include "boson/measure.hpp"
#include "boson/permutations.hpp"
#include "perm_matrix.hpp"

#include <vector>
#include <cmath>
#include <cfloat>

namespace bos {

enum class sym_def {
  sum_of_squares,
  sqrt_sum_of_squares,
  sum_of_modules,
  gauss_of_sum_of_squares,
  cup_of_sum_of_squares,
};

struct row_entry {
  row_entry() = default;
  row_entry(size_t _Yidx, size_t _Xcnt = SIZE_MAX, double _val = DBL_MAX)
  : state_idx {_Yidx}, sym_ent_cnt {_Xcnt}, val {_val}
  { }

  size_t state_idx, sym_ent_cnt;
  double val;
};

typedef std::vector<row_entry> matrix_row;

struct mult_prob_marker {
  int nfin;
  bos::state substate;

  mult_prob_marker(int _mult, const bos::state &_substate)
  : nfin {_mult}, substate {_substate}
  { }

  bool
  operator == (const mult_prob_marker &other) const noexcept
  { return nfin == other.nfin and substate == other.substate; }
};

typedef std::vector<std::pair<int, double>> mult_probs;
typedef std::unordered_map<mult_prob_marker, mult_probs> mult_prob_table;

struct matrix_base {
  matrix_base(const cellular_universe &model)
  : m_model {model},
    m_states {model.ordered_states()}
  { }

  virtual int
  n_rows() const = 0;

  virtual int
  n_ext_parm() const = 0;

  virtual void
  update(const double *extparm) = 0;

  virtual void
  calc_rxs(const bos::measure &fu, const bos::measure &fv,
      std::vector<std::pair<double, double>> &out) const = 0;

  struct symmetry {
    sym_def def;
    double sigma;

    static symmetry
    sum_of_squares() noexcept
    { return symmetry { .def = sym_def::sum_of_squares, .sigma = NAN }; }

    static symmetry
    sqrt_sum_of_squares() noexcept
    { return symmetry { .def = sym_def::sqrt_sum_of_squares, .sigma = NAN }; }

    static symmetry
    sum_of_modules() noexcept
    { return symmetry { .def = sym_def::sum_of_modules, .sigma = NAN }; }

    static symmetry
    gauss_of_sum_of_squares(double sigma) noexcept
    { return symmetry{.def = sym_def::gauss_of_sum_of_squares, .sigma = sigma};}

    static symmetry
    cup_of_sum_of_squares(double sigma) noexcept
    { return symmetry{.def = sym_def::cup_of_sum_of_squares, .sigma = sigma};}
  };

  virtual void
  calc_derivatives_rx(const measure &fu, const measure &fv, const int *iswelldef,
      symmetry sym, const double *extparm, double *rxgrad, double *symgrad)
    const = 0;

  virtual void
  calc_corrfunc(const std::pair<int, int> &cells,
      std::vector<std::pair<double, double>> &out) const = 0;

  virtual void
  calc_mult_probs(const int cell, mult_prob_table &out) const = 0;

  virtual double
  symmetry_term() const = 0;

  virtual double
  symmetry_term_abs() const = 0;

  virtual void
  write(FILE *out) const = 0;

  virtual void
  read(std::istream &in) = 0;

  virtual void
  to_ext_parms(std::vector<double> &extparms) const = 0;

  const cellular_universe&
  model() const noexcept
  { return m_model; }

  const std::vector<state>&
  states() const noexcept
  { return m_states; }

  protected:
  const cellular_universe &m_model;
  const std::vector<state> &m_states;
};

class matrix : public matrix_base {
  public:
  matrix(const cellular_universe &model);

  int
  n_rows() const noexcept override
  { return m_rows.size(); }

  int
  row_length(int Xidx) const noexcept
  { return m_rows[Xidx].size(); }

  int
  n_ext_parm(int Xidx) const noexcept
  { return m_rows[Xidx].size() - 1; }

  int
  n_ext_parm() const noexcept override;

  int
  row_entry_to_state_index(int rowidx, int entcnt) const noexcept
  { return m_rows[rowidx][entcnt].state_idx; }

  int
  state_index_to_row_entry(int xidx, int yidx) const noexcept
  {
    const matrix_row &xrow = m_rows[xidx];
    for (int iy = 0; iy < xrow.size(); ++iy)
    {
      if (xrow[iy].state_idx == yidx)
        return iy;
    }
    abort();
  }

  double
  at(int xidx, int yidx) const noexcept
  { return m_rows[xidx][state_index_to_row_entry(xidx, yidx)].val; }

  double
  operator () (int xidx, int yidx) const noexcept
  { return this->at(xidx, yidx); }

  double
  symmetric_entry(int xidx, int ycnt) const noexcept
  {
    const int yidx = m_rows[xidx][ycnt].state_idx;
    const int xcnt = m_rows[xidx][ycnt].sym_ent_cnt;
    return m_rows[yidx][xcnt].val;
  }

  void
  update(const double *extparm) override;

  void
  calc_rxs(const measure &fu, const measure &fv,
      std::vector<std::pair<double, double>> &out) const override;

  void
  calc_corrfunc(const std::pair<int, int> &cells,
      std::vector<std::pair<double, double>> &out) const;

  void
  calc_mult_probs(const int cell, mult_prob_table &out) const;

  void
  calc_derivatives_rx(const measure &fu, const measure &fv, const int *iswelldef,
      symmetry sim, const double *extparm, double *rxgrad, double *symgrad)
    const override;

  double
  symmetry_term() const override;

  double
  symmetry_term_abs() const override;

  void
  write(FILE *out) const override;

  void
  read(std::istream &in) override;

  void
  to_ext_parms(std::vector<double> &extparms) const override { abort(); }

  private:
  std::vector<matrix_row> m_rows;

  friend class spatsym_matrix;
  friend class sym_matrix;
}; // class bos::matrix

class spatsym_matrix : public matrix_base {
  public:
  spatsym_matrix(const cellular_universe &model);

  int
  n_rows() const noexcept override
  { return m_base.n_rows(); }

  int
  n_ext_parm() const noexcept override
  {
    int n = 0;
    for (const auto &c : m_part)
      n += m_base.m_rows[c[0].first].size();
    return n;
  }

  int
  n_ext_parm(int iclass) const noexcept
  { return m_base.m_rows[m_part[iclass][0].first].size() - 1; }

  void
  update(const double *extparm) override;

  void
  calc_rxs(const bos::measure &fu, const bos::measure &fv,
      std::vector<std::pair<double, double>> &out) const override
  { return m_base.calc_rxs(fu, fv, out); }

  void
  calc_derivatives_rx(const measure &fu, const measure &fv, const int *iswelldef,
      symmetry sym, const double *extparm, double *rxgrad, double *symgrad)
    const override;

  void
  calc_corrfunc(const std::pair<int, int> &cells,
      std::vector<std::pair<double, double>> &out) const
  { return m_base.calc_corrfunc(cells, out); }

  void
  calc_mult_probs(const int cell, mult_prob_table &out) const
  { return m_base.calc_mult_probs(cell, out); }

  double
  symmetry_term() const override
  { return m_base.symmetry_term(); }

  double
  symmetry_term_abs() const override
  { return m_base.symmetry_term_abs(); }

  void
  write(FILE *out) const override
  { m_base.write(out); }

  void
  read(std::istream &in) override
  { m_base.read(in); }

  void
  to_ext_parms(std::vector<double> &extparms) const override;

  private:
  matrix m_base;
  partition m_part;
  std::vector<std::vector<int>> m_gymap;
  const std::unordered_map<state, int> &m_idxmap;
};

} // namespace bos

namespace std {
template <>
struct hash<bos::mult_prob_marker> {
  size_t
  operator () (const bos::mult_prob_marker &m) const noexcept
  { return m_statehash(m.substate) + m.nfin*pow(10, m.substate.n_cells()+1); }
  private:
  std::hash<bos::state> m_statehash;
};
}

#endif
