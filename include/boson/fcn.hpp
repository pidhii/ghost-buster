#ifndef BOSON_FCN_HPP
#define BOSON_FCN_HPP

#include "boson/measure.hpp"
#include "boson/matrix.hpp"
#include "gnuplot.hpp"

#include <variant>

#define FCN_EPSILON 1E-10

namespace bos {

struct loss_def {
  enum class def { rx_sum, corr_func, gauss_corr_func, cup_corr_func };

  static loss_def
  rx_sum()
  { return {.def = def::rx_sum, .sigma = 0 }; }

  static loss_def
  corr_func()
  { return {.def = def::corr_func, .sigma = 0 }; }

  static loss_def
  gauss_corr_func(double sigma)
  { return {.def = def::gauss_corr_func, .sigma = sigma }; }

  static loss_def
  cup_corr_func(double sigma)
  { return {.def = def::cup_corr_func, .sigma = sigma }; }

  loss_def::def def;
  double sigma;
};

class fcn {
  public:
  fcn(bos::matrix_base &matrix, loss_def lossdef = loss_def::rx_sum(),
      matrix_base::symmetry sym = matrix_base::symmetry::sum_of_squares());

  double
  min_loss() const noexcept
  { return m_minloss; }

  void
  attach(gnuplot &gp) noexcept
  { m_gp = &gp; }

  void
  add_measures(const bos::measure &fu, const bos::measure &fv)
  { m_measures.push_back({&fu, &fv}); }

  void
  add_corr_cells(int first, int second)
  { m_cellpairs.emplace_back(first, second); }

  const bos::matrix_base&
  matrix() const noexcept
  { return m_matrix; }

  void
  print() const;

  void
  plot() const;

  void
  set_gscale(double scale)
  { m_gscale = scale; }

  void
  set_sscale(double scale)
  { m_sscale = scale; }

  void
  set_yscale(double scale)
  { m_yscale = scale; }

  void
  update_matrix(const double *extparm)
  { m_matrix.update(extparm); }

  double
  operator () (const double *extparm = nullptr) const
  {
    const double gloss = ghosts_loss();
    const double yloss = yhwh_loss();
    const double sloss = calc_symmetry_term();
    const double loss = gloss*m_gscale + yloss*m_yscale + sloss*m_sscale;

    if (m_minloss - loss > FCN_EPSILON)
    {
      m_minloss = loss;
      print();
      if (m_drawcnt++ % 10 == 0)
        plot();
    }

    return loss;
  }

  void
  grad(const double *extparm, double *grad)
  {
    //switch (m_lossdef.def)
    //{
      //case loss_def::def::rx_sum:
        //return grad_rx_sum(extparm, grad);
      //case loss_def::def::corr_func:
      //case loss_def::def::gauss_corr_func:
        //break;
      //default:
        //abort();
    //}
  }

  std::pair<int, double>
  find_rx_max(int imeas) const;

  int
  n_measures() const noexcept
  { return m_measures.size(); }

  std::pair<const measure&, const measure&>
  fu_fv(int imeas) const noexcept
  { return {*m_measures[imeas].first, *m_measures[imeas].second}; }

  const std::vector<int>&
  is_rx_well_defined(int imeas) const;

  bool
  is_rx_well_defined(int imeas, int xidx) const
  { return is_rx_well_defined(imeas)[xidx]; }

  struct rx_stats { double rx_sum, rx_max, rx_mean, rx_stddev; };
  struct rx_stats_tot { double rx_max, rx_mean, rx_stddev; };
  void
  calc_rx_stats(rx_stats_tot &totstats, std::vector<rx_stats> &stats) const;

  private:
  double
  calc_symmetry_term() const
  {
    switch (m_sym.def)
    {
      case sym_def::sum_of_squares:
        return m_matrix.symmetry_term();
      case sym_def::sqrt_sum_of_squares:
        return sqrt(m_matrix.symmetry_term());
      case sym_def::sum_of_modules:
        return m_matrix.symmetry_term_abs();
      case sym_def::gauss_of_sum_of_squares:
      {
        const double st = m_matrix.symmetry_term();
        return 1-exp(-0.5*(st*st)/(m_sym.sigma*m_sym.sigma));
      }
      case sym_def::cup_of_sum_of_squares:
      {
        const double st = m_matrix.symmetry_term();
        return exp(st/m_sym.sigma) - 1;
      }
      default:
        abort();
    }
  }

  double
  rx_sum_loss() const;

  double
  corrfunc_loss() const;

  double
  gauss_corrfunc_loss() const;

  double
  cup_corrfunc_loss() const;

  double
  ghosts_loss() const
  {
    switch (m_lossdef.def)
    {
      case loss_def::def::rx_sum: return rx_sum_loss();
      case loss_def::def::corr_func: return corrfunc_loss();
      case loss_def::def::gauss_corr_func: return gauss_corrfunc_loss();
      case loss_def::def::cup_corr_func: return cup_corrfunc_loss();
      default: abort();
    }
  }

  double
  yhwh_loss() const;

  void
  grad_rx_sum(const double *extparm, double *grad);

  void
  plot_matrix() const;

  void
  plot_rx() const;

  private:
  bos::matrix_base &m_matrix;
  std::vector<std::pair<const bos::measure*, const bos::measure*>> m_measures;
  std::vector<std::pair<int, int>> m_cellpairs;

  gnuplot *m_gp;
  mutable size_t m_drawcnt;

  matrix_base::symmetry m_sym;
  loss_def m_lossdef;

  mutable double m_minloss;
  mutable std::vector<std::vector<int>> m_isrxwelldef;

  double m_sscale;
  double m_gscale;
  double m_yscale;
};

} // namespace bos

#endif
