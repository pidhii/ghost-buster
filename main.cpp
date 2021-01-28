#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <type_traits>
#include <array>
#include <cfloat>
#include <random>
#include <fstream>
#include <getopt.h>
#include <cstring>
#include <optional>

#include "boson/states.hpp"
#include "boson/permutations.hpp"
#include "boson/matrix.hpp"
#include "boson/measure.hpp"
#include "boson/models.hpp"
#include "boson/fcn.hpp"
#include "gnuplot.hpp"
#include "minuit.hpp"


static void
write_parameters(std::ostream &out, int nparm);

static void
read_parameters(std::istream &in, int nparm, double errscale = 1,
    std::optional<double> newbnd1 = std::nullopt,
    std::optional<double> newbnd2 = std::nullopt);

static void
fcn_proxy(int*, REAL *grad, REAL *fval, REAL *xval, int *iflag, void *futil);

static bos::matrix_base&
select_matrix(const std::string &matrixdef, const bos::cellular_universe &model)
{
  static bos::matrix defaultmatrix {model};
  static bos::spatsym_matrix spatsymmatrix {model};
  if (matrixdef == "default")
    return defaultmatrix;
  else if (matrixdef == "spatsym")
    return spatsymmatrix;
  else
  {
    std::cerr << "> error: undefined value of matrix definition (" << matrixdef
      << ")" << std::endl;
    exit(EXIT_FAILURE);
  }
}

int
main(int argc, char **argv)
{
  option longopts[] = {
    {"help", false, nullptr, 'h'},
    {"input", true, nullptr, 'i'},
    {"output", true, nullptr, 'o'},
    {"update", true, nullptr, 'u'},
    {"cycles", true, nullptr, 'c'},
    {"interactive", false, nullptr, 'I'},
    {"loss-scale", true, nullptr, 'R'},
    {"sym-scale", true, nullptr, 'S'},
    {"yhwh-scale", true, nullptr, 'Y'},
    {"sym-def", true, nullptr, 's'},
    {"error-scale", true, nullptr, 'E'},
    {"measure", true, nullptr, 'm'},
    {"cells-pair", true, nullptr, 'p'},
    {"all-cells", false, nullptr, 0xC1},
    {"cmd", true, nullptr, 0xC2},
    {"bnd1", true, nullptr, 0xC3},
    {"bnd2", true, nullptr, 0xC4},
    {"matrix-def", true, nullptr, 0xC5},
    {"confirm-save", false, nullptr, 0xC6},
    {"write-matrix", true, nullptr, 0xC7},
    {"load-matrix", true, nullptr, 0xC8},
    {"sval", true, nullptr, 0xC9},
    {"loss-def", true, nullptr, 0xCA},
    {0,0,0,0}
  };

  std::string ipath, opath;
  int cycles = 2;
  bool interactive = false;
  std::optional<double> gscale, sscale, yscale;
  double errscale = 1;
  int nballs = 3;
  std::vector<std::string> measlist;
  std::vector<std::string> extracommands;
  std::optional<double> bnd1, bnd2;
  bool confirm_save = false;
  std::string mpath;
  std::vector<std::pair<int, int>> meascells = {{2, 5}};
  //std::vector<std::pair<int, int>> meascells = {{0,3}, {1,4}, {2,5}};
  bos::matrix_base::symmetry sym = bos::matrix_base::symmetry::sum_of_squares();
  std::string matrixdef = "spatsym";
  bos::loss_def lossdef = bos::loss_def::corr_func();
  std::string lpath;
  std::optional<std::pair<double, double>> sval;
  std::vector<std::pair<int, int>> cellpairs;

  int opt;
  while ((opt = getopt_long(argc, argv, "hN:i:o:u:c:IR:S:Y:m:p:E:s:", longopts, nullptr)) > 0)
  {
    switch (opt)
    {
      case 'h':
        printf("usage: %s -N <balls> [OPTIONS]\n", argv[0]);
        printf("\n"
               "\e[1mRead/write matrix parameterization:\e[0m\n"
               "  --input  -i <path>   Read parameters from specified file.\n"
               "  --output -o <path>   Write parameters to specified file after the fit.\n"
               "  --update -u <path>   Combination of the two options above.\n"
               "\n"
               "\e[1mLoss function definition:\e[0m\n"
               "  --loss-def         <type>    Specify a strategy to minimize correlations. Values must be one of:\n"
               "                               - rx-sum -- Minimize Rx for supplied definitions of measures.\n"
               "                               - corr-func -- Minimize factorization of probabilities directly (default).\n"
               "                               - gauss(<sigma>)\n"
               "                               - cup(<sigma>)\n"
               "  --sym-def       -s <type>    Specify implementation of symmetry constrain. Value must be one of:\n"
               "                               - sum-of-squares\n"
               "                               - sqrt-sum-of-squares\n"
               "                               - sum-of-modules\n"
               "                               - gauss(<sigma>)\n"
               "                               - cup(<sigma>)\n"
               "  --matrix-def=spatsym         Specify matrix type: 'default' or 'spatsym'.\n"
               "  --measure       -m <type>[x<n>]  Add a measure specified by \e[4mtype\e[0m in the loss function.\n"
               "                               Available values are:\n"
               "                               - f1 -- From paper: 1/(1 + n_u^I + n_u^F)\n"
               "                               - f2 -- From paper: 1/(1 + n_u^I) + 1/(1 + n_u^F)\n"
               "                               - full -- Generate random graph of function f(X, n_u).\n"
               "                               - mult -- Generate random graph of function f(n_u, n_u).\n"
               "                               In the last two cases separate function is generated for each cell to be measured.\n"
               "  --cells-pair    -p <first>,<second>|'all'\n"
               "  --all-cells                  Measure correlations for all separated pairs of cells.\n"
               "  --rxterm-scale  -R <scale>   Scale Rx-term of the loss function.\n"
               "  --symterm-scale -S <scale>   Scale symmetry-term of the loss function.\n"
               "\n"
               "\e[1mFit configuration:\e[0m\n"
               "  --cycles      -c <n=2>\n"
               "  --interactive -I\n"
               "  --cmd <command>\n"
               "\n"
               "\e[1mMiscelenious:\e[0m\n"
               "  --error-scale -E <x>\n"
               "  --bnd1 <x>\n"
               "  --bnd2 <x>\n"
               "  --sval <val>,<err>\n"
               "  --write-matrix   <path>\n"
               "  --load-matrix    <path>\n"
               "\n"
            );
        return EXIT_SUCCESS;

      case 'N':
        nballs = atoi(optarg);
        break;

      case 'i':
        ipath = optarg;
        break;

      case 'o':
        opath = optarg;
        break;

      case 'u':
      {
        // use as input path if file exists
        char buf[0x100];
        std::sprintf(buf, "test -e %s", optarg);
        if (std::system(buf) == 0)
          ipath = optarg;
        // use as output path
        opath = optarg;
        break;
      }

      case 'c':
        cycles = atoi(optarg);
        break;

      case 'I':
        interactive = true;
        break;

      case 'R':
        gscale = strtod(optarg, nullptr);
        break;

      case 'S':
        sscale = strtod(optarg, nullptr);
        break;

      case 'Y':
        yscale = strtod(optarg, nullptr);
        break;

      case 's':
        if (std::strcmp(optarg, "sum-of-squares") == 0)
          sym = bos::matrix_base::symmetry::sum_of_squares();
        else if (std::strcmp(optarg, "sqrt-sum-of-squares") == 0)
          sym = bos::matrix_base::symmetry::sqrt_sum_of_squares();
        else if (std::strcmp(optarg, "sum-of-modules") == 0)
          sym = bos::matrix_base::symmetry::sum_of_modules();
        else if (std::strncmp(optarg, "gauss", 5) == 0)
        {
          double sigma;
          sscanf(optarg, "gauss(%lf)", &sigma);
          sym = bos::matrix_base::symmetry::gauss_of_sum_of_squares(sigma);
        }
        else if (std::strncmp(optarg, "cup", 3) == 0)
        {
          double sigma;
          sscanf(optarg, "cup(%lf)", &sigma);
          sym = bos::matrix_base::symmetry::cup_of_sum_of_squares(sigma);
        }
        else
        {
          std::cerr << "> error: undefined value for symmetric term definition"
            << std::endl;
          return EXIT_FAILURE;
        }
        break;

      case 'E':
        errscale = std::strtod(optarg, nullptr);
        break;

      case 'm':
        if (char *p = std::strchr(optarg, 'x'))
        {
          int n;
          std::sscanf(p, "x%d", &n);
          const std::string meas (optarg, p);
          while (n--)
            measlist.emplace_back(meas);
        }
        else
          measlist.emplace_back(optarg);
        break;

      case 'p':
      {
        if (std::strcmp(optarg, "all") == 0)
        {
          for (int i = 0; i < 3; ++i)
            cellpairs.emplace_back(i, i+3);
        }
        else
        {
          int first, second;
          sscanf(optarg, "%d,%d", &first, &second);
          cellpairs.emplace_back(first, second);
        }
        break;
      }

      case 0xC1:
        meascells = {{0,3}, {1,4}, {2,5}};
        break;

      case 0xC2:
        extracommands.emplace_back(optarg);
        break;

      case 0xC3:
        bnd1 = std::strtod(optarg, nullptr);
        break;

      case 0xC4:
        bnd2 = std::strtod(optarg, nullptr);
        break;

      case 0xC5:
        matrixdef = optarg;
        break;

      case 0xC6:
        confirm_save = true;
        break;

      case 0xC7:
        mpath = optarg;
        break;

      case 0xC8:
        lpath = optarg;
        break;

      case 0xC9:
      {
        double val, err;
        sscanf(optarg, "%lf,%lf", &val, &err);
        sval = {val, err};
        break;
      }

      case 0xCA:
        if (std::strcmp(optarg, "rx-sum") == 0)
          lossdef = bos::loss_def::rx_sum();
        else if (std::strcmp(optarg, "corr-func") == 0)
          lossdef = bos::loss_def::corr_func();
        else if (std::strncmp(optarg, "gauss", 5) == 0)
        {
          double sigma;
          sscanf(optarg, "gauss(%lf)", &sigma);
          lossdef = bos::loss_def::gauss_corr_func(sigma);
        }
        else if (std::strncmp(optarg, "cup", 3) == 0)
        {
          double sigma;
          sscanf(optarg, "cup(%lf)", &sigma);
          lossdef = bos::loss_def::cup_corr_func(sigma);
        }
        else
        {
          std::cerr << "> error: undefined value for loss definition"
            << std::endl;
          return EXIT_FAILURE;
        }
        break;

      default:
        std::cout << "> error: undefined option, " << argv[optind-1] << std::endl;
        return EXIT_FAILURE;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  //                            CELL MODEL
  const bos::cellular_universe &model = bos::select_model(nballs);
  bos::matrix_base &matrix = select_matrix(matrixdef, model);

  //////////////////////////////////////////////////////////////////////////////
  //                           RANDOM GENERATOR
  std::random_device rnddev;
  std::mt19937 rndgen {rnddev()};


  //////////////////////////////////////////////////////////////////////////////
  //                           INITIALIZE MINUIT
  mninit();
  mncomd(fcn_proxy, "SET PRINT 0", nullptr);
  mncomd(fcn_proxy, "SET STRATEGY 2", nullptr);
  if (not lpath.empty())
  {
    std::cout << "> deduce parameters from matrix " << lpath << std::endl;
    std::ifstream matfile {lpath};
    matrix.read(matfile);
    std::vector<double> extparm;
    matrix.to_ext_parms(extparm);
    char parnam[42];
    for (int i = 0; i < matrix.n_ext_parm(); ++i)
    {
      sprintf(parnam, "%d", i);
      const double val = extparm[i];
      const double err = val*1e-2;
      mnparm(i+1, parnam, val, err, bnd1.value_or(0), bnd2.value_or(0));
    }
  }
  else if (not ipath.empty())
  {
    std::cout << "> read parameters from " << ipath << std::endl;
    std::ifstream parmfile {ipath};
    read_parameters(parmfile, matrix.n_ext_parm(), errscale, bnd1, bnd2);
  }
  else if (sval.has_value())
  {
    char parnam[42];
    int parmcnt = 1;
    const double val = sval->first;
    const double err = sval->second;
    std::cout << "> set starting values for all parameters to " << val << "±"
      << err << std::endl;
    for (int i = 0; i < matrix.n_ext_parm(); ++i)
    {
      sprintf(parnam, "%d", i);
      mnparm(parmcnt++, parnam, val, err, bnd1.value_or(0), bnd2.value_or(0));
    }
  }
  else
  { // Generate random initial parameters.
    char parnam[42];
    int parmcnt = 1;
    std::uniform_real_distribution<double> rnddis {1e-2, 1};
    std::cout << "> generate random parameters" << std::endl;
    for (int i = 0; i < matrix.n_ext_parm(); ++i)
    {
      sprintf(parnam, "%d", i);
      const double val = rnddis(rndgen);
      const double err = val*1e-2;
      mnparm(parmcnt++, parnam, val, err, bnd1.value_or(0), bnd2.value_or(0));
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  //                           SET UP FCN
  bos::fcn fcn {matrix, lossdef, sym};
  for (const auto [first, second] : cellpairs)
    fcn.add_corr_cells(first, second);
  for (const auto [u, v] : meascells)
  {
    for (const std::string &meastype : measlist)
    {
      if (meastype == "full")
      {
        auto *fu_full = new bos::measure_graph_full {model, u};
        auto *fv_full = new bos::measure_graph_full {model, v};
        fu_full->fill_random(rndgen, -10, 10);
        fv_full->fill_random(rndgen, -10, 10);
        fcn.add_measures(*fu_full, *fv_full);
      }
      else if (meastype == "mult")
      {
        auto *fu_mult = new bos::measure_graph_mult {model, u};
        auto *fv_mult = new bos::measure_graph_mult {model, v};
        fu_mult->fill_random(rndgen, -10, 10);
        fv_mult->fill_random(rndgen, -10, 10);
        fcn.add_measures(*fu_mult, *fv_mult);
      }
      else if (meastype == "f1")
      {
        auto *fu_1 = new bos::measure_paper_1 {model, u};
        auto *fv_1 = new bos::measure_paper_1 {model, v};
        fcn.add_measures(*fu_1, *fv_1);
      }
      else if (meastype == "f2")
      {
        auto *fu_2 = new bos::measure_paper_2 {model, u};
        auto *fv_2 = new bos::measure_paper_2 {model, v};
        fcn.add_measures(*fu_2, *fv_2);
      }
      else
      {
        std::cerr << "> error: undefined mesure type (" << meastype << ")"
          << std::endl;
        return EXIT_FAILURE;
      }
    }
  }

  // set loss scales
  if (gscale.has_value())
    fcn.set_gscale(gscale.value());
  if (sscale.has_value())
    fcn.set_sscale(sscale.value());
  if (yscale.has_value())
    fcn.set_yscale(yscale.value());


  //////////////////////////////////////////////////////////////////////////////
  //                              PLOTTING
  gnuplot gp;
  if (not getenv("DISABLE_GNUPLOT"))
  {
    gp.init();
    //gp.command("set terminal qt size 1600,800");
    gp.command("set terminal qt");
    fcn.attach(gp);
  }

  mncomd(fcn_proxy, "CALL FCN", &fcn);
  fcn.print();
  fcn.plot();


  //////////////////////////////////////////////////////////////////////////////
  //                                FIT
  if (interactive)
  {
    std::cout << "> starting interactive mode" << std::endl;
    mnintr(fcn_proxy, &fcn);
    std::cout << "> left interactive mode" << std::endl;
    mncomd(fcn_proxy, "CALL FCN", &fcn);
    fcn.print();
    fcn.plot();
  }
  else
  {
    for (const std::string &cmd : extracommands)
      mncomd(fcn_proxy, cmd.c_str(), &fcn);

    std::cout << "> begin fit-loop" << std::endl;
    for (int cyc = 0; cyc < cycles; ++cyc)
    {
      mncomd(fcn_proxy, "MIGRAD 100000", &fcn);
      mncomd(fcn_proxy, "CALL FCN", &fcn);
      fcn.print();
      fcn.plot();
    }
    std::cout << "> done" << std::endl;
  }

  if (not mpath.empty())
  {
    std::cout << "> write matrix to " << mpath << std::endl;
    FILE *file = fopen(mpath.c_str(), "w");
    fcn.matrix().write(file);
    fclose(file);
  }

  if (not opath.empty())
  {
    if (confirm_save)
    {
      std::cout << "> confirm saving parameters to " << opath
        << "?\n\e[1m>\e[0m ";
      std::cout.flush();
      std::string ans;
      std::cin >> ans;
      if (std::toupper(ans.at(0)) != 'Y')
      {
        std::cout << "> wont save current parameters" << std::endl;
        goto skip_save;
      }
    }
    else
      std::cout << "> save parameters to " << opath << std::endl;

    std::ofstream parmfile {opath};
    write_parameters(parmfile, matrix.n_ext_parm());
  }
skip_save:

#ifdef DUMP_STATS
  bos::fcn::rx_stats_tot totstats;
  std::vector<bos::fcn::rx_stats> stats;
  fcn.calc_rx_stats(totstats, stats);
  std::cout << std::hexfloat;
  std::cout << "@STATS BEGIN@" << std::endl;
  std::cout << '{' << std::endl;
  std::cout << "\tsym_squares=" << fcn.matrix().symmetry_term() << ',' << std::endl;
  std::cout << "\tsym_abs=" << fcn.matrix().symmetry_term_abs() << ',' << std::endl;
  std::cout
    << "\ttot_stats={\n"
    << "\t\trx_max=" << totstats.rx_max << ",\n"
    << "\t\trx_mean=" << totstats.rx_mean << ",\n"
    << "\t\trx_stddev=" << totstats.rx_stddev << ",\n"
    << "\t}," << std::endl;
  std::cout << "\tstat_list=[" << std::endl;
  for (int i = 0; i < stats.size(); ++i)
  {
    std::cout << "\t\t{"
      << "id=" << i << ','
      << "rx_sum=" << stats[i].rx_sum << ','
      << "rx_max=" << stats[i].rx_max << ','
      << "rx_mean=" << stats[i].rx_mean << ','
      << "rx_stddev=" << stats[i].rx_stddev
      << "}," << std::endl;
  }
  std::cout << "\t]" << std::endl;
  std::cout << "}" << std::endl;
  std::cout << "@STATS END@" << std::endl;
#endif

  gp.command("pause mouse close");
  return EXIT_SUCCESS;
}



static void
write_parameters(std::ostream &out, int nparm)
{
  std::string name;
  double val, err, bnd1, bnd2;
  for (int i = 0; i < nparm; ++i)
  {
    mnpout(i+1, name, val, err, bnd1, bnd2);
    out
      << name << '\t'
      << std::hexfloat
      << val << '\t' << err << '\t' << bnd1 << '\t' << bnd2
      << std::endl;
  }
}

static void
read_parameters(std::istream &in, int nparm, double errscale,
    std::optional<double> newbnd1, std::optional<double> newbnd2)
{
  std::string name, valbuf, errbuf, bnd1buf, bnd2buf;
  for (int i = 0; i < nparm; ++i)
  {
    in >> name >> valbuf >> errbuf >> bnd1buf >> bnd2buf;
    const double val = std::strtod(valbuf.c_str(), nullptr);
    double err = std::strtod(errbuf.c_str(), nullptr);
    const double bnd1 = newbnd1.value_or(std::strtod(bnd1buf.c_str(), nullptr));
    const double bnd2 = newbnd2.value_or(std::strtod(bnd2buf.c_str(), nullptr));
    // "Fix" parameter error
    //if (val - err < bnd1 or err < 1e-7)
    if (val - err < bnd1)
    {
      std::cout << "> parameter exceeds a limit: " << val << " - " << err
        << " < " << bnd1 << std::endl;
      err = (val - bnd1)*1e-1;
      std::cout << "> error is changed to " << err << std::endl;
    }
    else
      err *= errscale;
    // Initialize MINUIT parameter
    mnparm(i+1, name, val, err, bnd1, bnd2);
  }
}

static void
fcn_proxy(int*, REAL *grad, REAL *fval, REAL *xval, int *iflag,
    void *futil)
{
  bos::fcn &fcn = *static_cast<bos::fcn*>(futil);

  fcn.update_matrix(xval);

  if (*iflag == 2)
    fcn.grad(xval, grad);

  *fval = fcn(xval);
}

