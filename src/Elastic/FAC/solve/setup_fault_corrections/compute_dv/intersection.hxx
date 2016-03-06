#pragma once

/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <FTensor.hpp>

int intersection(const FTensor::Tensor1<double,3> &ntt,
                 const FTensor::Tensor2<double,3,3> &rot,
                 const FTensor::Tensor1<double,3> &dx,
                 const double fault[],
                 const int &dim);

