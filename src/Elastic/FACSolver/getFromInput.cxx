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
    enableLogging(database.getBool("enable_logging"));
  if (database.isString("coarse_fine_discretization"))
    setCoarseFineDiscretization(database.getString
                                ("coarse_fine_discretization"));
  if (database.isString("v_prolongation_method"))
    set_V_ProlongationMethod(database.getString("v_prolongation_method"));
  if (database.isDouble("coarse_solver_tolerance"))
    setCoarsestLevelSolverTolerance(database.getDouble
                                    ("coarse_solver_tolerance"));
  if (database.isInteger("coarse_solver_max_iterations"))
    setCoarsestLevelSolverMaxIterations(database.getInteger
                                        ("coarse_solver_max_iterations"));
}
