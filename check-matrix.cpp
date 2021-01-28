#include <iostream>
#include <optional>
#include <fstream>
#include <cstring>
#include <cassert>
#include <vector>
#include <iomanip>

#include <getopt.h>

#include "boson/models.hpp"
#include "boson/matrix.hpp"
#include "gnuplot.hpp"

typedef long double REAL;

template <typename T>
inline void
copy_matrix(T *dst, int m, int n, const T *src)
{ std::memcpy(dst, src, sizeof(T)*m*n); }

template <typename T>
void
matmult(T *C, int m, int n, const T *A, int o, int p, const T *B)
{
  assert(n == o);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < p; ++j) {
      T sum = 0;
      for (int k = 0; k < n; ++k) {
        const T aik = A[n*i+k];
        const T bkj = B[p*k+j];
        sum += aik*bkj;
      }
      C[i*p+j] = sum;
    }
  }
}

int
main(int argc, char **argv)
{
  option longopts[] = {
    {0,0,0,0}
  };

  int nevolve = 100;
  std::optional<int> nballs;

  int opt;
  while ((opt = getopt_long(argc, argv, "e:N:", longopts, NULL)) > 0) {
    switch (opt) {
      case 'e':
        nevolve = atoi(optarg);
        break;

      case 'N':
        nballs = atoi(optarg);
        break;

      default:
        std::cout << "> error: undefined option, " << argv[optind-1] << std::endl;
        return EXIT_FAILURE;
    }
  }

  if (not nballs.has_value()) {
    std::cerr << "> error: number of balls not specified" << std::endl;
    return EXIT_FAILURE;
  }
  if (optind == argc) {
    std::cerr << "> error: no matrix provided" << std::endl;
    return EXIT_FAILURE;
  }
  const std::string mpath = argv[optind];

  const bos::cellular_universe &model = bos::select_model(nballs.value());

  //////////////////////////////////////////////////////////////////////////////
  //                             READ MATRIX
  //
  const int n = model.states().size();
  REAL B[n][n];
  std::cout << "> reading matrix" << std::endl;
  std::ifstream matfile {mpath};
  std::string buf;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      matfile >> buf;
      B[i][j] = std::strtod(buf.c_str(), nullptr);
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  //                                 GOD
  //
  std::vector<std::vector<std::tuple<int, int, REAL, REAL>>> precepts;
  for (int xidx = 0; xidx < n; ++xidx) {
    const bos::state &x = model.ordered_states()[xidx];
    precepts.emplace_back();
    auto &xprec = precepts.back();

    for (int yidx = 0; yidx < n; ++yidx) {
      const bos::state &y = model.ordered_states()[yidx];

      if (x[0]!=y[0] or x[1]!=y[1] or x[2]!=y[2])
        continue;

      for (int nfin = 0; nfin <= model.n_balls(); ++nfin) {
        REAL px = 0, py = 0;
        for (int zidx = 0; zidx < n; ++zidx) {
          const bos::state &z = model.ordered_states()[zidx];
          if (z[1] != nfin)
            continue;
          px += B[xidx][zidx];
          py += B[yidx][zidx];
        }
        xprec.emplace_back(yidx, nfin, px, py);
      }
    }
  }

  std::cout << "> searching for a God..." << std::endl;
  for (int xidx = 0; xidx < n; ++xidx) {
    REAL maxdiff = DBL_MIN;
    int maxstate = -1;
    int maxnfin = -1;
    for (const auto [yidx, nfin, px, py] : precepts[xidx]) {
      //const REAL cmp = std::fabs(px/py-1);
      const REAL cmp = std::fabs(px-py);
      if (cmp > maxdiff) {
        maxdiff = cmp;
        maxstate = yidx;
        maxnfin = nfin;
      }
    }
    if (maxstate >= 0)
      std::cout << "  Xi=" << model.ordered_states()[xidx] << " -> " << maxdiff
        << "(Xj=" << model.ordered_states()[maxstate] << ", n^{fin}=" << maxnfin << ")"
        << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////////
  //                              EVOLUTION
  //
  REAL buf1[n][n], buf2[n][n];
  copy_matrix(&buf1[0][0], n, n, &B[0][0]);
  REAL *oldacc = &buf1[0][0];
  REAL *newacc = &buf2[0][0];
  std::cout << "> running evolution...";
  std::cout.flush();
  for (int i = 1; i <= nevolve; ++i) {

    //matmult(newacc, n, n, oldacc, n, n, &B[0][0]);
    matmult(newacc, n, n, oldacc, n, n, oldacc);
    std::swap(newacc, oldacc);

    std::cout << "\e[2K\r";
    std::streamsize ss = std::cout.precision();
    std::cout << "> running evolution... "
      << std::fixed << std::setprecision(1)
      << (REAL(i*100)/nevolve) << "%";
    std::cout.unsetf(std::ios::fixed);
    std::cout.unsetf(std::ios_base::floatfield);
    std::cout << std::setprecision(ss);
    std::defaultfloat(std::cout);
    //std::cout.unsetf(std::fixed);
    std::cout.flush();
  }
  std::cout << "\e[2K\r";
  std::cout << "> running evolution... done" << std::endl;

  //////////////////////////////////////////////////////////////////////////////
  //                               STD. DEV.
  //
  {
    const REAL mean_theo = 1./n;
    REAL mean_real = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        const REAL bij = newacc[i*n+j];
        mean_real += bij;
      }
    }
    mean_real = mean_real/(n*n);
    REAL stddev_theo = 0, stddev_real = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        const REAL bij = newacc[i*n+j];
        stddev_theo += (mean_theo - bij)*(mean_theo - bij);
        stddev_real += (mean_real - bij)*(mean_real - bij);
      }
    }
    stddev_theo = sqrt(stddev_theo/n);
    stddev_real = sqrt(stddev_real/n);
    std::cout << "> theo: mean=" << mean_theo << ", std.dev.=" << stddev_theo << std::endl;
    std::cout << "> real: mean=" << mean_real << ", std.dev.=" << stddev_real << std::endl;
    std::cout << "> theo - real: " << std::fabs(mean_theo - mean_real) << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////////
  //                                 PLOT
  //
  gnuplot gp;
  gp.init();

  char cmdbuf[256];
  sprintf(cmdbuf, "set xrange [-1:%d]", n);
  gp.command(cmdbuf);
  sprintf(cmdbuf, "set yrange [-1:%d]", n);
  gp.command(cmdbuf);
  gp.command("plot '-' matrix with image");
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
      gp << ' ' << newacc[i*n+j];
    gp << '\n';
  }
  gp.command("e");
  gp.command("unset xrange");
  gp.command("unset yrange");
}
