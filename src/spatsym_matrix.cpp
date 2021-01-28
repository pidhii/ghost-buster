#include "boson/matrix.hpp"


bos::spatsym_matrix::spatsym_matrix(const cellular_universe &model)
: matrix_base {model},
  m_base {model},
  m_idxmap {model.index_map()}
{
  spatial_quotient(model, m_part);

  m_gymap.resize(n_rows());
  for (int iclass = 0; iclass < m_part.size(); ++iclass)
  {
    const int repridx = m_part[iclass][0].first;
    matrix_row &reprrow = m_base.m_rows[repridx];

    for (const auto [xidx, g] : m_part[iclass])
    {
      matrix_row &xrow = m_base.m_rows[xidx];
      m_gymap[xidx].resize(xrow.size());
      for (int irepry = 0; irepry < xrow.size(); ++irepry)
      {
        const state &repry = m_states[reprrow[irepry].state_idx];
        const state &y = repry(g);
        const int yidx = m_idxmap.at(y);
        // now we need to find offset of this `y` in xidx-row
        const auto it = std::find_if(xrow.begin(), xrow.end(),
          [=] (const auto &rowent) -> bool {
            return rowent.state_idx == yidx;
          });
        const int iy = it - xrow.begin();
        m_gymap[xidx][iy] = irepry;
      }
    }
  }
}

void
bos::spatsym_matrix::update(const double *extparm)
{
  for (int iclass = 0; iclass < m_part.size(); ++iclass)
  {
    // representative of the class
    const int repridx = m_part[iclass][0].first;
    matrix_row &reprrow = m_base.m_rows[repridx];

    double s = 1;
    for (int iy = 0; iy < reprrow.size() - 1; s += extparm[iy++]);
    for (int iy = 0, iext = 0; iy < reprrow.size(); ++iy)
    {
      if (reprrow[iy].state_idx == repridx)
        reprrow[iy].val = 1/s;
      else
        reprrow[iy].val = extparm[iext++]/s;
    }

    for (const auto &[xidx, g] : m_part[iclass])
    {
      matrix_row &row = m_base.m_rows[xidx];
      assert(row.size() == reprrow.size());
      for (int iy = 0; iy < row.size(); ++iy)
        row[iy].val = reprrow[m_gymap[xidx][iy]].val;
    }

    // shift parameters
    extparm += reprrow.size() - 1;
  }
}

void
bos::spatsym_matrix::calc_derivatives_rx(const measure &fu, const measure &fv,
    const int *iswelldef, symmetry sym, const double *extparm,
    double *rxgrad, double *symgrad) const
{
  const int nstates = n_rows();

  double sym_squares = NAN;
  double sym_gauss_prefix = NAN;
  double sym_modules = NAN;
  switch (sym.def)
  {
    case sym_def::sum_of_squares:
    case sym_def::sqrt_sum_of_squares:
      sym_squares = symmetry_term();
      break;
    case sym_def::sum_of_modules:
      sym_modules = symmetry_term_abs();
      break;
    default:
      abort();
  }

  std::vector<std::pair<double, double>> rxs;
  calc_rxs(fu, fv, rxs);

  for (int iclass = 0;
      iclass < m_part.size();
      extparm += n_ext_parm(iclass),
      rxgrad += n_ext_parm(iclass),
      symgrad += n_ext_parm(iclass),
      ++iclass)
  {
    // representative of the class
    const int repridx = m_part[iclass][0].first;
    const matrix_row &reprrow = m_base.m_rows[repridx];

    const int nextparm = reprrow.size() - 1;
    const int ny = nextparm + 1;

    // convert external parameters to probabilites
    double s = 0;
    double reprb[ny];
    for (int iext = 0; iext < nextparm; s += extparm[iext++]);
    for (int iy = 0, iext = 0; iy < ny; ++iy)
    {
      const int yidx = reprrow[iy].state_idx;
      if (yidx == repridx)
        reprb[iy] = 1/(1+s);
      else
        reprb[iy] = extparm[iext++]/(1+s);
    }

    // calculate derivatives of B wrt external parameters
    double dreprbdext[ny][nextparm];
    for (int iext = 0; iext < nextparm; ++iext)
    {
      // diagonal B is 1/(1+s)
      int iy = 0;
      for (; iy < ny; ++iy)
      {
        const int yidx = reprrow[iy].state_idx;
        if (yidx == repridx)
          break;
        dreprbdext[iy][iext] = -extparm[iext]/((1+s)*(1+s));
        if (iy == iext)
          dreprbdext[iy][iext] += 1/(1+s);
      }
      dreprbdext[iy++][iext] = -1/((1+s)*(1+s));
      for (; iy < ny; ++iy)
      {
        dreprbdext[iy][iext] = -extparm[iext]/((1+s)*(1+s));
        if (iy-1 == iext)
          dreprbdext[iy][iext] += 1/(1+s);
      }
    }

    for (const auto [xidx, g] : m_part[iclass])
    {
      const auto [top, bot] = rxs[xidx];
      if (not iswelldef[xidx] or std::isnan(bot) or bot == 0)
        continue;

      const matrix_row &xrow = m_base.m_rows[xidx];

      // map B-entries of representative to the entries of the current row
      int gymap[ny];
      std::copy(m_gymap[xidx].begin(), m_gymap[xidx].end(), gymap);

      // calculate derivatives of "top" and "bot" wrt B-entries
      double dtopdb[ny], dbotdb[ny];
      double sumbfu = 0, sumbfv = 0, sumbfu2 = 0, sumbfv2 = 0;
      for (int iy = 0; iy < ny; ++iy)
      {
        const int yidx = xrow[iy].state_idx;
        const double fuxy = fu(xidx,yidx), fvxy = fv(xidx,yidx);
        sumbfu += reprb[gymap[iy]]*fuxy;
        sumbfv += reprb[gymap[iy]]*fvxy;
        sumbfu2 += reprb[gymap[iy]]*fuxy*fuxy;
        sumbfv2 += reprb[gymap[iy]]*fvxy*fvxy;
      }
      for (int iy = 0; iy < ny; ++iy)
      {
        const int yidx = xrow[iy].state_idx;
        const double fuxy = fu(xidx,yidx), fvxy = fv(xidx,yidx);
        // top
        dtopdb[iy] = fuxy*fvxy - fuxy*sumbfv - fvxy*sumbfu;
        // bot
        const double du2 = sumbfu2 - sumbfu*sumbfu;
        const double dv2 = sumbfv2 - sumbfv*sumbfv;
        const double uterm = (fuxy*fuxy - 2*fuxy*sumbfu)*dv2;
        const double vterm = (fvxy*fvxy - 2*fvxy*sumbfv)*du2;
        dbotdb[iy] = 0.5*(uterm + vterm)/bot;
      }

      // calculate derivatives of "top" and "bot" wrt external parameters
      double dtopdext[nextparm], dbotdext[nextparm];
      for (int iext = 0; iext < nextparm; ++iext)
      {
        double sum_top = 0, sum_bot = 0;
        for (int iy = 0; iy < ny; ++iy)
        {
          sum_top += dtopdb[iy]*dreprbdext[gymap[iy]][iext];
          sum_bot += dbotdb[iy]*dreprbdext[gymap[iy]][iext];
        }
        dtopdext[iext] = sum_top;
        dbotdext[iext] = sum_bot;
      }

      // final calculations
      for (int iext = 0; iext < nextparm; ++iext)
      {
        // Rx term
        const double topterm = dtopdext[iext]/bot;
        const double botterm = (top*dbotdext[iext])/(bot*bot);
        const double rxterm = 2.*(top/bot)*(topterm - botterm);

        // symmetric term
        double symterm = 0;
        for (int iy = 0; iy < ny; ++iy)
        {
          const int yidx = xrow[iy].state_idx;
          const double sign = xidx > yidx ? -1 : 1;
          const double bxy = reprb[gymap[iy]];
          const double byx = m_base.symmetric_entry(xidx, iy);
          const double dbdext = dreprbdext[gymap[iy]][iext];
          switch (sym.def)
          {
            case sym_def::sum_of_squares:
              symterm += sign*2*(byx-bxy)*dbdext;
              break;
            case sym_def::sqrt_sum_of_squares:
              symterm += 1/sqrt(sym_squares)*sign*2*(byx-bxy)*dbdext;
              break;
            case sym_def::sum_of_modules:
              symterm += ((byx-bxy) > 0 ? 1 : -1)*sign*dbdext;
              break;
            case sym_def::gauss_of_sum_of_squares:
              symterm += sym_gauss_prefix*sign*2*(byx - bxy)*dbdext;
              break;
            default:
              abort();
          }
        }

        rxgrad[iext] += rxterm;
        symgrad[iext] += symterm;
      }
    }
  }
}

void
bos::spatsym_matrix::to_ext_parms(std::vector<double> &extparms) const
{
  extparms.clear();
  extparms.reserve(n_ext_parm());
  for (int iclass = 0; iclass < m_part.size(); ++iclass)
  {
    // representative of the class
    const int repridx = m_part[iclass][0].first;
    const matrix_row &reprrow = m_base.m_rows[repridx];

    // find diagonal entry
    int iy_diag = 0;
    for (int iy = 0; iy < reprrow.size(); ++iy)
    {
      if (reprrow[iy].state_idx == repridx)
      {
        iy_diag = iy;
        break;
      }
    }

    // find `s`
    const double s = 1/reprrow[iy_diag].val - 1;

    // calc ext parameters
    for (int iext = 0, iy = 0; iext < reprrow.size() - 1; ++iext, ++iy)
    {
      if (iy == iy_diag)
        iy++;
      extparms.push_back(reprrow[iy].val*(1 + s));
    }
  }
}

