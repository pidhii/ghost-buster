#ifndef GNUPLOT_HPP
#define GNUPLOT_HPP

#include <iostream>
#include <cstdio>

class gnuplot : private std::streambuf, public std::ostream {
  public:
  gnuplot() : std::ostream {this} { }

  ~gnuplot()
  {
    if (m_gppipe)
      pclose(m_gppipe);
  }

  void
  init()
  { m_gppipe = popen("gnuplot -p --slow", "w"); }

  int
  overflow(int c) noexcept override
  {
    if (m_gppipe)
      fputc(c, m_gppipe);
    return 0;
  }

  operator FILE* () noexcept
  { return m_gppipe; }

  void
  command(const std::string &cmd)
  {
    if (m_gppipe)
    {
      fputs(cmd.c_str(), m_gppipe);
      fputc('\n', m_gppipe);
      fflush(m_gppipe);
    }
  }

  private:
  FILE *m_gppipe;
};

#endif
