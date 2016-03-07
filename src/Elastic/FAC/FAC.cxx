/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Constants.hxx"
#include "Elastic/FAC.hxx"

#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/pdat/SideVariable.h>
#include <SAMRAI/pdat/NodeVariable.h>
#include <SAMRAI/pdat/EdgeVariable.h>

Elastic::FAC::FAC(const SAMRAI::tbox::Dimension& dimension,
                  SAMRAI::tbox::Database &Database):
  database(Database),
  dimension(dimension),
  cell_moduli_id(invalid_id),
  edge_moduli_id(invalid_id),
  v_id(invalid_id),
  v_rhs_id(invalid_id),
  dv_diagonal_id(invalid_id),
  dv_mixed_id(invalid_id),
  level_set_id(invalid_id),
  lambda("lambda",database,dimension),
  mu("mu",database,dimension),
  level_set("level_set",database,dimension),
  offset_vector_on_output(database.getBoolWithDefault
                          ("offset_vector_on_output",false))
{
  const int dim(dimension.getValue());

  std::string xyz("xyz");
  for(int d=0; d<dim; ++d)
    {
      v_rhs[d]=Input_Expression(std::string("v_rhs_")+xyz[d],database,
                                dimension,true);
      v_initial[d]=Input_Expression(std::string("v_initial_")+xyz[d],database,
                                    dimension);
    }

  SAMRAI::hier::VariableDatabase* vdb =
    SAMRAI::hier::VariableDatabase::getDatabase();

  /// Get a unique context for variables owned by this object.
  context = vdb->getContext("Elastic::FAC:Context");

  /// Register variables with SAMRAI::hier::VariableDatabase and get
  /// the descriptor indices for those variables.  Ghost cells width
  /// are 1 just in case it is needed.

  int depth=2;
  boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >
    cell_moduli_ptr(new SAMRAI::pdat::CellVariable<double>
                    (dimension,"Elastic::FAC:cell_moduli",depth));
  cell_moduli_id = vdb->registerVariableAndContext
    (cell_moduli_ptr,context,SAMRAI::hier::IntVector::getOne(dimension));

  if(database.keyExists("faults"))
    { faults=database.getDoubleVector("faults"); }
  if(faults.size()%9!=0)
    TBOX_ERROR("The number of points in faults must be "
               "divisible by 9.  Read "
               << faults.size() << " points");
  if(!faults.empty())
    {
      boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >
        dv_diagonal_ptr(new SAMRAI::pdat::CellVariable<double>
                        (dimension,"Elastic::FAC:dv_diagonal",dim));
      dv_diagonal_id = vdb->registerVariableAndContext
        (dv_diagonal_ptr,context,SAMRAI::hier::IntVector::getOne(dimension));

      boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
        dv_mixed_ptr(new SAMRAI::pdat::SideVariable<double>
                     (dimension,"Elastic::FAC:dv_mixed",
                      SAMRAI::hier::IntVector::getOne(dimension),
                      dim==2 ? 2 : 8));
      dv_mixed_id = vdb->registerVariableAndContext
        (dv_mixed_ptr,context,SAMRAI::hier::IntVector::getOne(dimension));
    }
  if(database.keyExists("refinement_points"))
    { refinement_points=database.getDoubleVector("refinement_points"); }
  if(refinement_points.size()%3!=0)
    { TBOX_ERROR("The number of points in refinement_points must be "
                 "divisible by 3.  Read "
                 << refinement_points.size() << " points"); }

  if(have_embedded_boundary())
    {
      boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
        level_set_ptr(new SAMRAI::pdat::SideVariable<double>
                      (dimension,"Elastic::FAC:level_set",
                       SAMRAI::hier::IntVector::getOne(dimension),depth));
      level_set_id = vdb->registerVariableAndContext
        (level_set_ptr,context,SAMRAI::hier::IntVector::getOne(dimension));
    }

  if(dim==2)
    {
      boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> >
        edge_moduli_ptr(new SAMRAI::pdat::NodeVariable<double>
                        (dimension,"Elastic::FAC:edge_moduli",depth));
      edge_moduli_id = vdb->registerVariableAndContext
        (edge_moduli_ptr,context,SAMRAI::hier::IntVector::getOne(dimension));
    }
  else if(dim==3)
    {
      boost::shared_ptr<SAMRAI::pdat::EdgeVariable<double> >
        edge_moduli_ptr(new SAMRAI::pdat::EdgeVariable<double>
                        (dimension,"Elastic::FAC:edge_moduli",depth));
      edge_moduli_id = vdb->registerVariableAndContext
        (edge_moduli_ptr,context,SAMRAI::hier::IntVector::getOne(dimension));
    }

  boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
    v_ptr(new SAMRAI::pdat::SideVariable<double>
          (dimension, "Elastic::FAC:v",
           SAMRAI::hier::IntVector::getOne(dimension), 1));
  v_id = vdb->registerVariableAndContext
    (v_ptr,context,SAMRAI::hier::IntVector::getOne(dimension));

  boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
    v_rhs_ptr(new SAMRAI::pdat::SideVariable<double>
              (dimension,"Elastic::FAC:v right hand side",
               SAMRAI::hier::IntVector::getOne(dimension)));
  v_rhs_id =
    vdb->registerVariableAndContext(v_rhs_ptr,context,
                                    SAMRAI::hier::IntVector::getOne(dimension));

  d_adaption_threshold=database.getDoubleWithDefault("adaption_threshold",
                                                      1.0e-15);
  min_full_refinement_level
    =database.getIntegerWithDefault("min_full_refinement_level",0);
}
