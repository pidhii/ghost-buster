#ifndef DIHEDRAL_HPP
#define DIHEDRAL_HPP

#include <functional>

// Dihedral group.
// Elements are accumulated in form:
//   g = s r^n
struct dihedral {
  struct s_type { };
  static constexpr s_type s {};

  dihedral(int _n) : n {_n}, refl {false}, rot {0} { }

  dihedral
  operator * (int r) const noexcept
  {
    if (r < 0)
      r = n - (-r % n);
    return {n, refl, (rot + r) % n};
  }

  dihedral
  operator * (s_type _) const noexcept
  { return dihedral {n, !refl, n - rot}; }

  int
  operator () (int k) const noexcept
  {
    if (k < 0)
      k = n - (-k % n);

    k = (k + rot) % n; // rotation
    if (refl) k = (n - k) % n; // reflection
    return k;
  }

  bool
  operator == (const dihedral &other) const noexcept
  { return refl == other.refl and rot == other.rot; }

  size_t
  hash() const noexcept
  { return (!!refl)*n + rot; }

  private:
  dihedral(int _n, bool _s, int _r) : n {_n}, refl {_s}, rot {_r} { }

  int n;
  bool refl;
  int rot;
};

namespace std {
template <>
struct hash<dihedral> {
  size_t
  operator () (const dihedral &g) const noexcept
  { return g.hash(); }
};
} // namespace std

#endif
