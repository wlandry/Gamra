/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACSolver.hxx"

/// Do not allow FAC preconditioner and Elastic FAC operators to be
/// set from database, as that may cause them to be inconsistent with
/// this object if user does not coordinate the inputs correctly.

/// SAMRAI::tbox::Database::isBool() and friends are not const, so we
/// can not pass a const database :(
void Elastic::FACSolver::getFromInput (SAMRAI::tbox::Database &database)
{
  if (database.isBool("enable_logging"))
    { operators->logging=database.getBool("enable_logging"); }
  if (database.isDouble("coarse_solver_tolerance"))
    { setCoarsestLevelSolverTolerance(database.getDouble
                                      ("coarse_solver_tolerance")); }
  if (database.isInteger("coarse_solver_max_iterations"))
    { setCoarsestLevelSolverMaxIterations(database.getInteger
                                          ("coarse_solver_max_iterations")); }
}
