#ifndef MINUIT_HPP
#define MINUIT_HPP

#include <string>
#include <stdexcept>

typedef double REAL;

typedef void (*fcn_t)(int *npar, REAL *grad, REAL *fval, REAL *xval, int *iflag,
    void *futil);

extern "C" void
mninit_(int *ird, int *irw, int *isav);

extern "C" void
mnparm_(int *num, const char *chnam, REAL *stval, REAL *step, REAL *bnd1,
    REAL *bnd2, int *ierflg, int chnam_len);

extern "C" void
mncomd_(fcn_t fcn, const char *chstr, int *icondn, void *futil, int chstr_len);

extern "C" void
mnintr_(fcn_t fcn, void *futil);

extern "C" void
mnpout_(int *num, char *chnam, REAL *val, REAL *error, REAL *bnd1, REAL *bnd2,
    int *ivarbl, int chnam_len);

struct minuit_error : std::runtime_error {
  using std::runtime_error::runtime_error;
};

void
mninit(int ird = 5, int iwr = 6, int isav = 7)
{ mninit_(&ird, &iwr, &isav); }

int
mnparm(int num, const std::string &chnam, REAL stval, REAL step,
    REAL bnd1 = 0, REAL bnd2 = 0)
{
  int ierflag = 0;
  mnparm_(&num, chnam.c_str(), &stval, &step,&bnd1,&bnd2,&ierflag,chnam.size());
  return ierflag;
}

int
mncomd(fcn_t fcn, const std::string &chstr, void *futil = nullptr)
{
  int icondn = 0;
  mncomd_(fcn, chstr.c_str(), &icondn, futil, chstr.size());
  return icondn;
}

void
mnintr(fcn_t fcn, void *futil = nullptr)
{ mnintr_(fcn, futil); }

void
mnpout(int num, std::string &name, REAL &val, REAL &err, REAL &bnd1, REAL &bnd2)
{
  int ivarbl;
  name.resize(42);
  mnpout_(&num, (char*)name.c_str(), &val, &err, &bnd1, &bnd2, &ivarbl,
      name.size());
}

void
mnpout(int num, std::string &name, REAL &val, REAL &err)
{
  REAL bnd1, bnd2;
  name.resize(42);
  mnpout(num, name, val, err, bnd1, bnd2);
}

void
mnpout(int num, std::string &name, REAL &val)
{
  REAL err;
  name.resize(42);
  mnpout(num, name, val, err);
}

#endif
