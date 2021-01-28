#include <iostream>

#include "boson/matrix.hpp"
#include "boson/permutations.hpp"

bos::matrix::matrix(const cellular_universe &model)
: matrix_base(model)
{
  const int n = m_states.size();
  m_rows.reserve(n);
  for (int i = 0; i < n; ++i)
  {
    m_rows.emplace_back();
    matrix_row &row = m_rows.back();

    for (int j = 0; j < n; ++j)
    {
      if (are_adjacent(m_states[i], m_states[j]))
        row.emplace_back(j);
    }
  }

  for (int Xidx = 0; Xidx < n; ++Xidx)
  {
    matrix_row &Xrow = m_rows[Xidx];
    for (int iy = 0; iy < Xrow.size(); ++iy)
    {
      const matrix_row &Yrow = m_rows[Xrow[iy].state_idx];
      for (int Zcnt = 0; Zcnt < Yrow.size(); ++Zcnt)
      {
        if (Yrow[Zcnt].state_idx == Xidx)
        {
          Xrow[iy].sym_ent_cnt = Zcnt;
          break;
        }
      }
    }
  }
}

int
bos::matrix::n_ext_parm() const noexcept
{
  int n = 0;
  for (int i = 0; i < m_rows.size(); ++i)
    n += m_rows[i].size() - 1;
  return n;
}

double
bos::matrix::symmetry_term() const
{
  double sum = 0;
  for (int Xidx = 0; Xidx < m_rows.size(); ++Xidx)
  {
    for (int iy = 0; iy < m_rows[Xidx].size(); ++iy)
    {
      const int i = Xidx;
      const int j = m_rows[Xidx][iy].state_idx;
      if (i > j) continue;

      const double Bxy = m_rows[i][iy].val;
      const double Byx = m_rows[j][m_rows[Xidx][iy].sym_ent_cnt].val;
      sum += (Bxy - Byx)*(Bxy - Byx);
    }
  }
  return sum;
}

double
bos::matrix::symmetry_term_abs() const
{
  double sum = 0;
  for (int Xidx = 0; Xidx < m_rows.size(); ++Xidx)
  {
    for (int iy = 0; iy < m_rows[Xidx].size(); ++iy)
    {
      const int i = Xidx;
      const int j = m_rows[Xidx][iy].state_idx;
      if (i > j) continue;

      const double Bxy = m_rows[i][iy].val;
      const double Byx = m_rows[j][m_rows[Xidx][iy].sym_ent_cnt].val;
      sum += std::fabs(Bxy - Byx);
    }
  }
  return sum;
}

void
bos::matrix::write(FILE *out) const
{
//#define SCALE(x) log10(x)
//#define SCALE(x) ((x) == 0 ? 0 : log10(x))
#define SCALE(x) (x)

  for (int Xidx = 0; Xidx < m_rows.size(); ++Xidx)
  {
    int yidx = 0;
    for (int iy = 0; iy < m_rows[Xidx].size(); ++iy)
    {
      for (; yidx < m_rows[Xidx][iy].state_idx; ++yidx)
        fprintf(out, " %a", SCALE(0.));
      fprintf(out, " %a", SCALE(m_rows[Xidx][iy].val));
      yidx++;
    }
    for (; yidx < m_rows.size(); ++yidx)
      fprintf(out, " %a", SCALE(0.));
    fputc('\n', out);
  }
}

void
bos::matrix::update(const double *extparm)
{
  for (int xidx = 0; xidx < n_rows(); extparm += n_ext_parm(xidx), ++xidx)
  {
    matrix_row &row = m_rows[xidx];

    double s = 1;
    for (int iy = 0; iy < row.size() - 1; s += extparm[iy++]);
    for (int iy = 0, iext = 0; iy < row.size(); ++iy)
    {
      if (row[iy].state_idx == xidx)
        row[iy].val = 1/s;
      else
        row[iy].val = extparm[iext++]/s;
    }
  }
}

void
bos::matrix::calc_rxs(const bos::measure &fu, const bos::measure &fv,
    std::vector<std::pair<double, double>> &out) const
{
  out.reserve(n_rows());
  for (int xidx = 0; xidx < n_rows(); ++xidx)
  {
    const matrix_row &row = m_rows[xidx];

    // calculate mean values
    double meanU = 0;
    double meanU2 = 0;
    double meanV = 0;
    double meanV2 = 0;
    double meanUV = 0;
    for (int iy = 0; iy < row.size(); ++iy)
    {
      const int yidx = row[iy].state_idx;
      const double By = row[iy].val;
      const double U = fu(xidx, yidx);
      const double V = fv(xidx, yidx);
      meanU += By*U;
      meanU2 += By*U*U;
      meanV += By*V;
      meanV2 += By*V*V;
      meanUV += By*U*V;
    }

    // calculate numerator & denomenator of Rx
    const double top = meanUV - meanU*meanV;
    const double bot = sqrt((meanU2 - meanU*meanU) * (meanV2 - meanV*meanV));
    out.emplace_back(top, bot);
  }
}

void
bos::matrix::calc_corrfunc(const std::pair<int, int> &cells,
    std::vector<std::pair<double, double>> &out) const
{
  const int nstates = model().states().size();
  const int nballs = model().n_balls();

  // for each X
  for (int xidx = 0; xidx < nstates; ++xidx)
  { // for each pair of multiplicities, n_u, n_v,  calculate Px(n_u, n_v) and
    // Px(n_u) Px(n_v)
    const matrix_row &xrow = m_rows[xidx];

    double px_u[nballs+1];
    double px_v[nballs+1];
    double px_uv[nballs+1][nballs+1];
    std::fill(px_u, px_u+nballs+1, 0);
    std::fill(px_v, px_v+nballs+1, 0);
    for (int i = 0; i <= nballs; ++i)
      std::fill(px_uv[i], px_uv[i]+nballs+1, 0);

    for (int iy = 0; iy < xrow.size(); ++iy)
    {
      // this probability will participate in:
      // - Px(n_u^Y, n_v^Y | X);
      // - Px(n_u^Y | X);
      // - Px(n_v^Y | X).
      const int yidx = xrow[iy].state_idx;
      const state &y = model().ordered_states()[yidx];
      const int nu = y[cells.first];
      const int nv = y[cells.second];
      const double bxy = xrow[iy].val;
      px_u[nu] += bxy;
      px_v[nv] += bxy;
      px_uv[nu][nv] += bxy;
    }

    for (int nu = 0; nu <= nballs; ++nu)
    {
      for (int nv = 0; nv <= nballs; ++nv)
        out.emplace_back(px_uv[nu][nv], px_u[nu]*px_v[nv]);
    }
  }
}

void
bos::matrix::calc_mult_probs(const int cell, mult_prob_table &out) const
{
  const int ncells = model().n_cells();
  const std::vector<int> adjcells = {
    (cell + ncells-1) % ncells, cell, (cell + 1) % ncells
  };
  std::vector<int> buf (3);

  for (int xidx = 0; xidx < n_rows(); ++xidx) {
    const state &x = states()[xidx];
    const matrix_row &xrow = m_rows[xidx];

    // extract substate arround `cell`
    for (int i = 0; i < buf.size(); ++i)
      buf[i] = x[adjcells[i]];
    state substate {buf};

    // calculate probabilites to get n^fin balls in final state given X
    double pnfin[model().n_balls()+1];
    std::fill(pnfin, pnfin+model().n_balls()+1, 0);
    for (int iy = 0; iy < xrow.size(); ++iy)
    {
      const state &y = states()[xrow[iy].state_idx];
      pnfin[y[cell]] += xrow[iy].val;
    }

    // record this probabilites
    mult_prob_marker mark {0, substate};
    for (int nfin = 0; nfin <= model().n_balls(); ++nfin)
    {
      if (pnfin[nfin] != 0)
      {
        mark.nfin = nfin;
        out[mark].emplace_back(xidx, pnfin[nfin]);
      }
    }
  }
}

void
bos::matrix::calc_derivatives_rx(const measure &fu, const measure &fv,
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
    case sym_def::gauss_of_sum_of_squares:
    case sym_def::cup_of_sum_of_squares:
      abort();
  }

  std::vector<std::pair<double, double>> rxs;
  calc_rxs(fu, fv, rxs);

  for (int xidx = 0; xidx < nstates;
      extparm += n_ext_parm(xidx),
      rxgrad += n_ext_parm(xidx),
      symgrad += n_ext_parm(xidx),
      ++xidx)
  {
    const auto [top, bot] = rxs[xidx];
    if (not iswelldef[xidx] or std::isnan(bot) or bot == 0)
      continue;

    const int nextparm_x = n_ext_parm(xidx);
    const int ny = row_length(xidx);

    // convert external parameters to probabilites
    double s = 0;
    double b[ny];
    for (int iext = 0; iext < nextparm_x; s += extparm[iext++]);
    for (int iy = 0, iext = 0; iy < ny; ++iy)
    {
      const int yidx = row_entry_to_state_index(xidx, iy);
      if (yidx == xidx)
        b[iy] = 1/(1+s);
      else
        b[iy] = extparm[iext++]/(1+s);
    }

    // calculate derivatives of "top" and "bot" wrt B-entries
    double dtopdb[ny], dbotdb[ny];
    double sumbfu = 0, sumbfv = 0, sumbfu2 = 0, sumbfv2 = 0;
    for (int iy = 0; iy < ny; ++iy)
    {
      const int yidx = row_entry_to_state_index(xidx, iy);
      const double fuxy = fu(xidx,yidx), fvxy = fv(xidx,yidx);
      sumbfu += b[iy]*fuxy;
      sumbfv += b[iy]*fvxy;
      sumbfu2 += b[iy]*fuxy*fuxy;
      sumbfv2 += b[iy]*fvxy*fvxy;
    }
    for (int iy = 0; iy < ny; ++iy)
    {
      const int yidx = row_entry_to_state_index(xidx, iy);
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

    // calculate derivatives of B wrt external parameters
    double dbdext[ny][nextparm_x];
    for (int iext = 0; iext < nextparm_x; ++iext)
    {
      // diagonal B is 1/(1+s)
      int iy = 0;
      for (; iy < ny; ++iy)
      {
        const int yidx = row_entry_to_state_index(xidx, iy);
        if (yidx == xidx)
          break;
        dbdext[iy][iext] = -extparm[iext]/((1+s)*(1+s));
        if (iy == iext)
          dbdext[iy][iext] += 1/(1+s);
      }
      dbdext[iy++][iext] = -1/((1+s)*(1+s));
      for (; iy < ny; ++iy)
      {
        dbdext[iy][iext] = -extparm[iext]/((1+s)*(1+s));
        if (iy-1 == iext)
          dbdext[iy][iext] += 1/(1+s);
      }
    }

    // calculate derivatives of "top" and "bot" wrt external parameters
    double dtopdext[nextparm_x], dbotdext[nextparm_x];
    for (int iext = 0; iext < nextparm_x; ++iext)
    {
      double sum_top = 0, sum_bot = 0;
      for (int iy = 0; iy < ny; ++iy)
      {
        sum_top += dtopdb[iy]*dbdext[iy][iext];
        sum_bot += dbotdb[iy]*dbdext[iy][iext];
      }
      dtopdext[iext] = sum_top;
      dbotdext[iext] = sum_bot;
    }

    // final calculations
    for (int iext = 0; iext < nextparm_x; ++iext)
    {
      // Rx term
      const double topterm = dtopdext[iext]/bot;
      const double botterm = (top*dbotdext[iext])/(bot*bot);
      const double rxterm = 2.*(top/bot)*(topterm - botterm);

      // symmetric term
      double symterm = 0;
      for (int iy = 0; iy < ny; ++iy)
      {
        const int yidx = row_entry_to_state_index(xidx, iy);
        const double biy_t = symmetric_entry(xidx, iy);
        const double sign = xidx > yidx ? -1 : 1;
        switch (sym.def)
        {
          case sym_def::sum_of_squares:
            symterm += sign*2*(biy_t - b[iy])*dbdext[iy][iext];
            break;
          case sym_def::sqrt_sum_of_squares:
            symterm +=
              1/sqrt(sym_squares)*sign*2*(biy_t - b[iy])*dbdext[iy][iext];
            break;
          case sym_def::sum_of_modules:
            symterm += ((biy_t - b[iy]) > 0 ? 1 : -1)*sign*dbdext[iy][iext];
            break;
          case sym_def::gauss_of_sum_of_squares:
            symterm += sym_gauss_prefix*sign*2*(biy_t - b[iy])*dbdext[iy][iext];
            break;
          case sym_def::cup_of_sum_of_squares:
            abort();
        }
      }

      rxgrad[iext] += rxterm;
      symgrad[iext] += symterm;
    }
  }
}

void
bos::matrix::read(std::istream &in)
{
  std::string buf;

  for (int xidx = 0; xidx < m_rows.size(); ++xidx)
  {
    int yidx = 0;
    for (int iy = 0; iy < m_rows[xidx].size(); ++iy)
    {
      for (; yidx < m_rows[xidx][iy].state_idx; ++yidx)
        in >> buf; // zero
      in >> buf;
      m_rows[xidx][iy].val = std::strtod(buf.c_str(), nullptr);
      yidx++;
    }
    for (; yidx < m_rows.size(); ++yidx)
      in >> buf; // zero

    double sum = 0;
    for (int iy = 0; iy < m_rows[xidx].size(); ++iy)
      sum += m_rows[xidx][iy].val;
    printf("sum(row[%d])=%E\n", xidx, sum);
    //std::cout << "sum(row[" << xidx << "]) = " << sum << std::endl;
  }
}

