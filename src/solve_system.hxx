#pragma once

/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/tbox/Database.h>
#include <SAMRAI/tbox/InputDatabase.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#include <SAMRAI/mesh/BergerRigoutsos.h>
#include <SAMRAI/mesh/TreeLoadBalancer.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>

template<class T>
bool solve_system
(T &fac,
 boost::shared_ptr<SAMRAI::tbox::Database> &main_db,
 boost::shared_ptr<SAMRAI::tbox::InputDatabase> &input_db,
 boost::shared_ptr<SAMRAI::hier::PatchHierarchy> &patch_hierarchy,
 const std::string &base_name,
 const SAMRAI::tbox::Dimension &dim)
{
  boost::shared_ptr<SAMRAI::mesh::StandardTagAndInitialize>
    tag_and_initializer(new SAMRAI::mesh::StandardTagAndInitialize
                        ("CellTaggingMethod",&fac,
                         input_db->getDatabase("StandardTagAndInitialize")));

  boost::shared_ptr<SAMRAI::mesh::BergerRigoutsos>
    box_generator(new SAMRAI::mesh::BergerRigoutsos(dim));
  boost::shared_ptr<SAMRAI::mesh::TreeLoadBalancer>
    load_balancer(new SAMRAI::mesh::TreeLoadBalancer
                  (dim,
                   "load balancer",
                   boost::shared_ptr<SAMRAI::tbox::Database>()));
  load_balancer->setSAMRAI_MPI(SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());

  SAMRAI::mesh::GriddingAlgorithm
    gridding_algorithm(patch_hierarchy,
                       "Gridding Algorithm",
                       input_db->getDatabase("GriddingAlgorithm"),
                       tag_and_initializer,
                       box_generator,
                       load_balancer);
  gridding_algorithm.makeCoarsestLevel(0.0);

  bool use_visit=true;
  if (main_db->keyExists("vis_writer"))
    {
      std::vector<std::string> vis_writer;
      vis_writer = main_db->getStringVector("vis_writer");
      use_visit=(std::find(vis_writer.begin(),vis_writer.end(),"VisIt")
                 !=vis_writer.end());
    }

  bool intermediate_output(main_db->getBoolWithDefault("intermediate_output",
                                                       false));

  boost::shared_ptr<SAMRAI::appu::VisItDataWriter> visit_writer;
  std::string vis_filename(main_db->getStringWithDefault("vis_filename",
                                                         base_name));
  if (use_visit)
    {
      visit_writer = boost::make_shared<SAMRAI::appu::VisItDataWriter>
        (dim,"Visit Writer",vis_filename + ".visit");
      fac.setupPlotter(*visit_writer);
    }

  if(main_db->getBoolWithDefault("print_input_file",true))
    input_db->printClassData(SAMRAI::tbox::plog);

  int lnum = 0;
  bool converged(fac.solve());
  for (;patch_hierarchy->levelCanBeRefined(lnum) && converged; ++lnum)
    {
      if (use_visit && intermediate_output)
        visit_writer->writePlotData(patch_hierarchy, lnum);
      std::vector<int> tag_buffer(patch_hierarchy->getMaxNumberOfLevels(),1);
      gridding_algorithm.regridAllFinerLevels(0,tag_buffer,0,0.0);
      SAMRAI::tbox::plog << "Newly adapted hierarchy\n";

      converged=fac.solve();
      SAMRAI::tbox::TimerManager::getManager()->print(SAMRAI::tbox::plog);
    }
  if (use_visit)
    visit_writer->writePlotData(patch_hierarchy, lnum);

  return converged;
}


