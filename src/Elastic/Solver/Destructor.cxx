/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Solver.hxx"

Elastic::Solver::~Solver()
{
  deallocateSolverState();
}
