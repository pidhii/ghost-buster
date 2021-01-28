#ifndef BOSON_MEASURE_HPP
#define BOSON_MEASURE_HPP

#include "boson/states.hpp"

#include <random>
#include <array>
#include <iostream>

namespace bos {

struct measure {
  virtual double
  operator () (int Xidx, int Yidx) const = 0;
};

class measure_graph_full : public measure {
  public:
  measure_graph_full(const cellular_universe &model, int cellno)
  : m_model {model},
    m_cellno {cellno},
    m_graph (model.states().size() * (model.n_balls()+1), 0)
  { }

  template <typename RndGen>
  void
  fill_random(RndGen &rndgen, double min, double max)
  {
    std::uniform_real_distribution<double> rnddis {min, max};
    for (double &fval : m_graph)
      fval = rnddis(rndgen);
  }

  double
  operator () (int Xidx, int Yidx) const noexcept override
  {
    const int row_start = Xidx*(m_model.n_balls()+1);
    const int col = m_model.ordered_states()[Yidx][m_cellno];
    return m_graph[row_start + col];
  }

  private:
  const cellular_universe &m_model;
  int m_cellno;
  std::vector<double> m_graph;
};

class measure_graph_mult : public measure {
  public:
  measure_graph_mult(const cellular_universe &model, int cellno)
  : m_model {model},
    m_cellno {cellno},
    m_graph ((model.n_balls()+1) * (model.n_balls()+1), 0)
  { }

  template <typename RndGen>
  void
  fill_random(RndGen &rndgen, double min, double max)
  {
    std::uniform_real_distribution<double> rnddis {min, max};
    for (double &fval : m_graph)
      fval = rnddis(rndgen);
  }

  double
  operator () (int Xidx, int Yidx) const noexcept override
  {
    const int row_start =
      m_model.ordered_states()[Xidx][m_cellno]*(m_model.n_balls()+1);
    const int col = m_model.ordered_states()[Yidx][m_cellno];
    return m_graph[row_start + col];
  }

  private:
  const cellular_universe &m_model;
  int m_cellno;
  std::vector<double> m_graph;
};

class measure_paper_1 : public measure {
  public:
  measure_paper_1(const cellular_universe &model, int cellno)
  : m_model {model},
    m_cellno {cellno}
  { }

  double
  operator () (int Xidx, int Yidx) const noexcept override
  {
    const double ni = m_model.ordered_states()[Xidx][m_cellno];
    const double nf = m_model.ordered_states()[Yidx][m_cellno];
    return 1. / (1 + ni + nf);
  }

  private:
  const cellular_universe &m_model;
  int m_cellno;
};

class measure_paper_2 : public measure {
  public:
  measure_paper_2(const cellular_universe &model, int cellno)
  : m_model {model},
    m_cellno {cellno}
  { }

  double
  operator () (int Xidx, int Yidx) const noexcept override
  {
    const double ni = m_model.ordered_states()[Xidx][m_cellno];
    const double nf = m_model.ordered_states()[Yidx][m_cellno];
    return 1./(1. + ni) + 1./(1. + nf);
  }

  private:
  const cellular_universe &m_model;
  int m_cellno;
};

} // namespace bos

#endif
