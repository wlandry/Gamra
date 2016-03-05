#pragma once

/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>
#include <SAMRAI/appu/VisDerivedDataStrategy.h>
#include <SAMRAI/appu/VisItDataWriter.h>

#include "Elastic/FACSolver.hxx"
#include <FTensor.hpp>

void compute_intersections_2D(const FTensor::Tensor1<double,3> &ntt,
                              const FTensor::Tensor1<double,3> &xyz,
                              const FTensor::Tensor2<double,3,3> &rot,
                              const FTensor::Tensor1<double,3> dx[],
                              const double fault[],
                              const int &dim,
                              const int &ix,
                              int &intersect_diagonal,
                              int intersect_mixed[2]);

void compute_intersections_3D(const FTensor::Tensor1<double,3> &ntt,
                              const FTensor::Tensor1<double,3> &xyz,
                              const FTensor::Tensor2<double,3,3> &rot,
                              const FTensor::Tensor1<double,3> dx[],
                              const double fault[],
                              const int &dim,
                              const int &ix,
                              int &intersect_diagonal,
                              int intersect_mixed[4],
                              int intersect_corner[4]);

namespace Elastic
{
  class FAC: public SAMRAI::mesh::StandardTagAndInitStrategy,
             public SAMRAI::appu::VisDerivedDataStrategy
  {
  public:
    FAC(const SAMRAI::tbox::Dimension& dim, SAMRAI::tbox::Database &database);
    virtual ~FAC() {}

    virtual void
    initializeLevelData
    (const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
     const int level_number,
     const double init_data_time,
     const bool can_be_refined,
     const bool initial_time,
     const boost::shared_ptr<SAMRAI::hier::PatchLevel>& old_level,
     const bool allocate_data);

    virtual void
    resetHierarchyConfiguration
    (const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& new_hierarchy,
     int coarsest_level,
     int finest_level);

    virtual void
    applyGradientDetector
    (const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
     const int level_number,
     const double error_data_time,
     const int tag_index,
     const bool initial_time,
     const bool uses_richardson_extrapolation);

    virtual bool
    packDerivedDataIntoDoubleBuffer(double* buffer,
                                    const SAMRAI::hier::Patch& patch,
                                    const SAMRAI::hier::Box& region,
                                    const std::string& variable_name,
                                    int depth, double) const;
    bool solve();

    void
    setupPlotter(SAMRAI::appu::VisItDataWriter& plotter) const;

  private:
    void fix_moduli();
    const SAMRAI::tbox::Dimension d_dim;

    boost::shared_ptr<SAMRAI::hier::PatchHierarchy> d_hierarchy;

  public:
    Boundary_Conditions d_boundary_conditions;
  private:
    Elastic::FACSolver d_elastic_fac_solver;
    boost::shared_ptr<SAMRAI::hier::VariableContext> d_context;
    double d_adaption_threshold;
    int min_full_refinement_level;
  public:
    int cell_moduli_id, edge_moduli_id, v_id, v_rhs_id, dv_diagonal_id,
      dv_mixed_id, level_set_id;

    // FIXME: This should be a std::array or std::vector, but we have
    // to wait for C++11.
    Input_Expression lambda, mu, v_rhs[3], v_initial[3], level_set;

    /// Whether to offset the vector when outputing.  This is useful
    /// for convergence tests, because it removes interpolation
    /// errors.  This is especially important near faults where there
    /// is a jump in values
    bool offset_vector_on_output;

    std::vector<double> faults;
    std::vector<double> refinement_points;

    bool have_embedded_boundary() const
    {
      return level_set.is_valid;
    }
    template<class T> void setup_fault_corrections();
  };
}
