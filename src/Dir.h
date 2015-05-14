#pragma once

/// A really simple class to make it easy to get the next direction
/// modulo the dimension.

#include "SAMRAI/tbox/Dimension.h"

namespace Gamra
{
class Dir
{
public:
  SAMRAI::tbox::Dimension::dir_t d;

  Dir(const SAMRAI::tbox::Dimension::dir_t &D): d(D) {}
  
  operator SAMRAI::tbox::Dimension::dir_t() const
  {
    return d;
  }

  Dir next(const SAMRAI::tbox::Dimension::dir_t &dim) const
  {
    Dir dd(*this);
    ++dd;
    dd.d%=dim;
    return dd;
  }

  Dir & operator++()
  {
    ++d;
    return *this;
  }

  static Dir from_int(const int &D)
  {
    return Dir(static_cast<SAMRAI::tbox::Dimension::dir_t>(D));
  }
};
}

namespace std {
  template<> class numeric_limits<Gamra::Dir> {
  public:
    static Gamra::Dir max()
    {
      return std::numeric_limits<SAMRAI::tbox::Dimension::dir_t>::max();
    }
  };
} 
