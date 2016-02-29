#pragma once

/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

/// Interpolate up to quadratic order from three coarse points (C) to
/// two fine points (f) with the setup

///    C--|--|--f--C--f--|--|--C


inline void quad_offset_interpolate(const double &plus, const double &center,
                                    const double &minus,
                                    double &fine_plus, double &fine_minus)
{
  const double d_plus=plus-center;
  const double d_minus=minus-center;

  fine_plus=center + (5*d_plus - 3*d_minus)/32;
  fine_minus=center + (5*d_minus - 3*d_plus)/32;
}

/// Interpolate to only the fine point closest to the "plus" coarse point.

inline double quad_offset_interpolate(const double &plus, const double &center,
                                      const double &minus)
{
  const double d_plus=plus-center;
  const double d_minus=minus-center;

  return center + (5*d_plus - 3*d_minus)/32;
}


