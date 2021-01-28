#ifndef BOSON_FCN_HPP

#include "boson/fcn.hpp"

#include <iostream>
#include <cassert>

#define PUPV_CUTOFF 1e-30

bos::fcn::fcn(bos::matrix_base &matrix, loss_def lossdef,
    matrix_base::symmetry sym)
: m_matrix {matrix},
  m_gp {nullptr},
  m_drawcnt {0},
  m_sym {sym},
  m_lossdef {lossdef},
  m_minloss {1e7},
  m_sscale {1},
  m_gscale {1},
  m_yscale {1}
{ }

void
bos::fcn::print() const
{
  const double sym_sqares = m_matrix.symmetry_term();
  const double sym_abs = m_matrix.symmetry_term_abs();
  const double ghst = ghosts_loss();
  const double yhwh = yhwh_loss();
  const double sym = calc_symmetry_term();
  const double loss = this->operator()();
  printf("loss=%1.2E, sym(B^2)=%1.2E, sym(|B|)=%1.2E\n"
         "GM=(%1.2E x %1.2E), YHWH=(%1.2E x %1.2E), SM=(%1.2E x %1.2E)\n",
      loss, sym_sqares, sym_abs, ghst, m_gscale, yhwh, m_yscale, sym, m_sscale);
  for (int imeas = 0; imeas < m_measures.size(); ++imeas)
  {
    const measure &fu = *m_measures[imeas].first;
    const measure &fv = *m_measures[imeas].second;
    std::vector<std::pair<double, double>> rxs;
    m_matrix.calc_rxs(fu, fv, rxs);
    double rx_sum = 0;
    for (int xidx = 0; xidx < m_matrix.n_rows(); ++xidx)
    {
      const auto [top, bot] = rxs[xidx];
      if (not is_rx_well_defined(imeas, xidx) or std::isnan(bot) or bot == 0)
        continue;

      const double rx = top/bot;
      rx_sum += std::fabs(rx);
    }
    const double maxrx = find_rx_max(imeas).second;
    printf("- measure #%d: sum(Rx)=%1.2E, max(Rx)=%1.2E\n",
        imeas, rx_sum, maxrx);
  }
  for (int icells = 0; icells < m_cellpairs.size(); ++icells)
  {
    std::vector<std::pair<double, double>> probs;
    m_matrix.calc_corrfunc(m_cellpairs[icells], probs);
    double meanabscorr = 0;
    double maxabscorr = DBL_MIN;
    double loss = 0;
    for (const auto [puv, pupv] : probs)
    {
      if (std::isnan(pupv) or pupv < PUPV_CUTOFF)
        continue;
      const double abscorr = std::fabs(puv/pupv - 1);
      meanabscorr += abscorr;
      maxabscorr = std::max(maxabscorr, abscorr);
      loss += abscorr*abscorr;
    }
    meanabscorr /= probs.size();
    printf("- cells (%d,%d): <|Puv/PuPv-1|>=%1.2E, max{|Puv/PuPv-1|}=%1.2E, sum((Puv/PuPv-1)^2)=%1.2E\n",
        m_cellpairs[icells].first, m_cellpairs[icells].second, meanabscorr,
        maxabscorr, loss);
  }
}

void
bos::fcn::plot() const
{
  if (m_gp)
  {
    for (int i = 0; i < 2; ++i) // XXX
    {
      m_gp->command("set multiplot layout 1,2");
      plot_matrix();
      plot_rx();
      fflush(*m_gp);
    }
  }
}

void
bos::fcn::plot_matrix() const
{
  if (m_gp)
  {
    const int nstates = m_matrix.n_rows();
    char buf[256];

    sprintf(buf, "set xrange [-1:%d]", nstates);
    m_gp->command(buf);
    sprintf(buf, "set yrange [-1:%d]", nstates);
    m_gp->command(buf);
    m_gp->command("plot '-' matrix with image");
    m_matrix.write(*m_gp);
    m_gp->command("e");
    m_gp->command("unset xrange");
    m_gp->command("unset yrange");
  }
}

double
bos::fcn::rx_sum_loss() const
{
  const int n_states = m_matrix.n_rows();

  double loss = 0;
  for (int imeas = 0; imeas < m_measures.size(); ++imeas)
  {
    const measure &fu = *m_measures[imeas].first;
    const measure &fv = *m_measures[imeas].second;

    std::vector<std::pair<double, double>> rxs;
    m_matrix.calc_rxs(fu, fv, rxs);
    for (int xidx = 0; xidx < m_matrix.n_rows(); ++xidx)
    {
      const auto [top, bot] = rxs[xidx];
      if (not is_rx_well_defined(imeas, xidx) or std::isnan(bot) or bot == 0)
        continue;

      const double rx = top/bot;
      loss += (top*top)/(bot*bot);
    }
  }
  return loss;
}

double
bos::fcn::corrfunc_loss() const
{
  const int n_states = m_matrix.n_rows();

  double loss = 0;
  for (int icells = 0; icells < m_cellpairs.size(); ++icells)
  {
    std::vector<std::pair<double, double>> probs;
    m_matrix.calc_corrfunc(m_cellpairs[icells], probs);
    for (const auto [puv, pupv] : probs)
    {
      if (std::isnan(pupv) or pupv < PUPV_CUTOFF)
        continue;
      const double corr = puv/pupv - 1;
      //std::cout << "corr = " << corr << std::endl;
      loss += corr*corr;
    }
  }
  return loss;
}

double
bos::fcn::gauss_corrfunc_loss() const
{
  const int n_states = m_matrix.n_rows();

  double loss = 0;
  for (int icells = 0; icells < m_cellpairs.size(); ++icells)
  {
    std::vector<std::pair<double, double>> probs;
    m_matrix.calc_corrfunc(m_cellpairs[icells], probs);
    for (const auto [puv, pupv] : probs)
    {
      if (std::isnan(pupv) or pupv < PUPV_CUTOFF)
        continue;
      const double corr = puv/pupv - 1;
      loss += 1-exp(-0.5*(corr*corr)/(m_lossdef.sigma*m_lossdef.sigma));
    }
  }
  return loss;
}


double
bos::fcn::cup_corrfunc_loss() const
{
  const int n_states = m_matrix.n_rows();

  double loss = 0;
  for (int icells = 0; icells < m_cellpairs.size(); ++icells)
  {
    std::vector<std::pair<double, double>> probs;
    m_matrix.calc_corrfunc(m_cellpairs[icells], probs);
    for (const auto [puv, pupv] : probs)
    {
      if (std::isnan(pupv) or pupv < PUPV_CUTOFF)
        continue;
      const double corr = puv/pupv - 1;
      loss += exp(corr*corr/m_lossdef.sigma*m_lossdef.sigma);
    }
  }
  return loss;
}


static double
fold_mult_probs(const bos::mult_prob_table &mptab)
{
  double ret = 0;
  for (const auto &[mark, probs] : mptab) // for each substate and n^fin
  {
    // all corresponding probabilites must equal

    for (int i = 0, j = 1; j < probs.size(); ++i, ++j)
    {
      const double pxi = probs[i].second;
      const double pxj = probs[j].second;
      ret += (pxi/pxj-1)*(pxi/pxj-1);
    }

    //for (int i = 0; i < probs.size(); ++i)
    //{
      //for (int j = i+1; j < probs.size(); ++j)
      //{
        //const double pxi = probs[i].second;
        //const double pxj = probs[j].second;
        //ret += (pxi/pxj-1)*(pxi/pxj-1);
      //}
    //}
  }
  return ret;
}

double
bos::fcn::yhwh_loss() const
{
  double loss = 0;
  mult_prob_table mptab;
  for (int icells = 0; icells < m_cellpairs.size(); ++icells)
  {
    mptab.clear();
    m_matrix.calc_mult_probs(m_cellpairs[icells].first, mptab);
    loss += fold_mult_probs(mptab);

    mptab.clear();
    m_matrix.calc_mult_probs(m_cellpairs[icells].second, mptab);
    loss += fold_mult_probs(mptab);
  }
  return loss;
}

std::pair<int, double>
bos::fcn::find_rx_max(int imeas) const
{
  const int n_states = m_matrix.n_rows();
  const auto &[fu, fv] = fu_fv(imeas);

  std::vector<std::pair<double, double>> rxs;
  m_matrix.calc_rxs(fu, fv, rxs);

  int xidxmax = -1;
  double rxmax = -1;
  for (int xidx = 0; xidx < n_states; ++xidx)
  {
    const auto [top, bot] = rxs[xidx];
    if (not is_rx_well_defined(imeas, xidx) or std::isnan(bot) or bot == 0)
      continue;

    double rx = std::fabs(top/bot);
    if (rx > rxmax)
    {
      rxmax = rx;
      xidxmax = xidx;
    }
  }
  return {xidxmax, rxmax};
}

void
bos::fcn::plot_rx() const
{
  std::vector<std::vector<double>> rxs (m_measures.size());
  for (int imeas = 0; imeas < m_measures.size(); ++imeas)
  {
    const auto &[fu, fv] = fu_fv(imeas);

    std::vector<std::pair<double, double>> topbot;
    m_matrix.calc_rxs(fu, fv, topbot);

    for (int xidx = 0; xidx < m_matrix.n_rows(); ++xidx)
    {
      const auto fufv = m_measures[imeas];
      const auto [top, bot] = topbot[xidx];

      double rx;
      if (not is_rx_well_defined(imeas, xidx) or std::isnan(bot) or bot == 0)
        rx = 1e-20;
      else
        rx = std::fabs(top/bot);

      rxs[imeas].push_back(rx);
    }
  }

  m_gp->command("set logscale y");
  m_gp->command("set yrange [1e-15:]");
  //m_gp->command("set yrange [1e-20:1]");
  m_gp->command("set format y \"%E\"");
  std::string cmd = "plot";
  for (int imeas = 0; imeas < m_measures.size(); ++imeas)
  {
    if (imeas > 0)
      cmd += ",";
    std::string measname = "f" + std::to_string(imeas+1);
    cmd += " '-' title '"+measname+"'";
  }
  m_gp->command(cmd);
  for (int imeas = 0; imeas < m_measures.size(); ++imeas)
  {
    for (int xidx = 0; xidx < m_matrix.n_rows(); ++xidx)
      *m_gp << xidx << '\t' << rxs[imeas][xidx] << std::endl;
    *m_gp << 'e' << std::endl;
  }
  m_gp->command("unset logscale y");
  m_gp->command("unset yrange");
  m_gp->command("unset format y");
}

void
bos::fcn::grad_rx_sum(const double *extparm, double *grad)
{
  const int nstates = m_matrix.n_rows();
  const int nextparm = m_matrix.n_ext_parm();

  double rxgrad[nextparm];
  double symgrad[nextparm];
  std::fill(rxgrad, rxgrad + nextparm, 0);
  std::fill(symgrad, symgrad + nextparm, 0);

  for (int imeas = 0; imeas < m_measures.size(); ++imeas)
  {
    const auto &[fu, fv] = fu_fv(imeas);
    m_matrix.calc_derivatives_rx(fu, fv, is_rx_well_defined(imeas).data(),
        m_sym, extparm, rxgrad, symgrad);
  }

  for (int iext = 0; iext < nextparm; ++iext)
    grad[iext] = rxgrad[iext]*m_gscale + symgrad[iext]*m_sscale;
}

const std::vector<int>&
bos::fcn::is_rx_well_defined(int imeas) const
{
  if (m_isrxwelldef.empty())
  {
    m_isrxwelldef.resize(m_measures.size());
    for (int imeas = 0; imeas < m_measures.size(); ++imeas)
    {
      std::vector<int> &isdef = m_isrxwelldef[imeas];
      isdef.resize(m_matrix.n_rows());
      const auto &[fu, fv] = fu_fv(imeas);
      for (int xidx = 0; xidx < m_matrix.n_rows(); ++xidx)
      {
        // check if this Rx is undefined
        const bos::state &x = m_matrix.states()[xidx];
        const bos::state_set &xperms = m_matrix.model().permutations(x);
        const double refval_u = fu(xidx, xidx);
        const double refval_v = fv(xidx, xidx);
        bool same_u = true;
        bool same_v = true;
        for (const bos::state &y : xperms)
        {
          const int yidx = m_matrix.model().state_index(y);
          same_u = same_u and fu(xidx, yidx) == refval_u;
          same_v = same_v and fv(xidx, yidx) == refval_v;
        }
        isdef[xidx] = (not same_u) and (not same_v);
      }
    }
  }
  return m_isrxwelldef[imeas];
}

void
bos::fcn::calc_rx_stats(rx_stats_tot &totstats, std::vector<rx_stats> &stats)
  const
{
  stats.resize(m_measures.size());
  int ndef_tot = 0;
  double rx_max_tot = DBL_MIN;
  double rx_sum_tot = 0;
  std::vector<std::vector<std::pair<double, double>>> rxss;
  rxss.resize(m_measures.size());
  for (int imeas = 0; imeas < m_measures.size(); ++imeas)
  {
    const auto [fu, fv] = fu_fv(imeas);
    std::vector<std::pair<double, double>> &rxs = rxss[imeas];
    m_matrix.calc_rxs(fu, fv, rxs);

    double rx_sum = 0;
    double rx_max = -1;
    int ndef = 0;
    for (int xidx = 0; xidx < m_matrix.n_rows(); ++xidx)
    {
      const auto [top, bot] = rxs[xidx];
      if (not is_rx_well_defined(imeas, xidx) or std::isnan(bot) or bot == 0)
        continue;

      const double rx = std::abs(top/bot);

      ndef ++;
      rx_max = std::max(rx_max, rx);
      rx_sum += rx;


      ndef_tot ++;
      rx_max_tot = std::max(rx_max_tot, rx);
      rx_sum_tot += rx;
    }

    const double rx_mean = rx_sum / ndef;
    double rx_stddev = 0;
    for (int xidx = 0; xidx < m_matrix.n_rows(); ++xidx)
    {
      const auto [top, bot] = rxs[xidx];
      if (not is_rx_well_defined(imeas, xidx) or std::isnan(bot) or bot == 0)
        continue;

      const double rx = std::abs(top/bot);
      rx_stddev += (rx - rx_mean)*(rx - rx_mean);
    }
    rx_stddev = std::sqrt(rx_stddev / ndef);

    stats[imeas].rx_sum = rx_sum;
    stats[imeas].rx_max = rx_max;
    stats[imeas].rx_mean = rx_mean;
    stats[imeas].rx_stddev = rx_stddev;
  }

  const double rx_mean_tot = rx_sum_tot / ndef_tot;
  double rx_stddev_tot = 0;
  for (int imeas = 0; imeas < m_measures.size(); ++imeas)
  {
    for (int xidx = 0; xidx < m_matrix.n_rows(); ++xidx)
    {
      const auto [top, bot] = rxss[imeas][xidx];
      if (not is_rx_well_defined(imeas, xidx) or std::isnan(bot) or bot == 0)
        continue;

      const double rx = std::abs(top/bot);
      rx_stddev_tot += (rx - rx_mean_tot)*(rx - rx_mean_tot);
    }
  }
  rx_stddev_tot = std::sqrt(rx_stddev_tot / ndef_tot);

  totstats.rx_max = rx_max_tot;
  totstats.rx_mean = rx_mean_tot;
  totstats.rx_stddev = rx_stddev_tot;
}

#endif
