#ifndef GAMRA_SOLVE_SYSTEM_H
#define GAMRA_SOLVE_SYSTEM_H

template<class T>
void solve_system(T &fac,
                  boost::shared_ptr<SAMRAI::tbox::Database> &main_db,
                  boost::shared_ptr<SAMRAI::tbox::InputDatabase> &input_db,
                  boost::shared_ptr<SAMRAI::hier::PatchHierarchy> &patch_hierarchy,
                  const string &base_name,
                  const SAMRAI::tbox::Dimension &dim)
{

  /*
   * Create the tag-and-initializer, box-generator and load-balancer
   * object references required by the gridding_algorithm object.
   */
  boost::shared_ptr<SAMRAI::mesh::StandardTagAndInitialize>
    tag_and_initializer(new SAMRAI::mesh::StandardTagAndInitialize
                        (dim,"CellTaggingMethod",&fac,
                         input_db->getDatabase("StandardTagAndInitialize")));

  boost::shared_ptr<SAMRAI::mesh::BergerRigoutsos>
    box_generator(new SAMRAI::mesh::BergerRigoutsos(dim));
  boost::shared_ptr<SAMRAI::mesh::TreeLoadBalancer>
    load_balancer(new SAMRAI::mesh::TreeLoadBalancer
                  (dim,
                   "load balancer",
                   boost::shared_ptr<SAMRAI::tbox::Database>()));
  load_balancer->setSAMRAI_MPI(SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());

  /*
   * Create the gridding algorithm used to generate the SAMR grid
   * and create the grid.
   */
  boost::shared_ptr<SAMRAI::mesh::GriddingAlgorithm>
    gridding_algorithm(new SAMRAI::mesh::GriddingAlgorithm(patch_hierarchy,
                                                           "Gridding Algorithm",
                                                           input_db->getDatabase("GriddingAlgorithm"),
                                                           tag_and_initializer,
                                                           box_generator,
                                                           load_balancer));
  // SAMRAI::tbox::plog << "Gridding algorithm:" << endl;
  // gridding_algorithm->printClassData(SAMRAI::tbox::plog);

  /*
   * Make the coarsest patch level where we will be solving.
   */
  gridding_algorithm->makeCoarsestLevel(0.0);
  // bool done = false;
  // for (int lnum = 0;
  //      patch_hierarchy->levelCanBeRefined(lnum) && !done; lnum++) {
  //   SAMRAI::tbox::plog << "Adding finner levels with lnum = " << lnum << endl;
  //   gridding_algorithm->makeFinerLevel(0.0,true,0);
  //   SAMRAI::tbox::plog << "Just added finer levels with lnum = " << lnum << endl;
  //   done = !(patch_hierarchy->finerLevelExists(lnum));
  // }

  /*
   * Set up the plotter for the hierarchy just created.  The FAC
   * object handles the data and has the function
   * setupExternalPlotter to register its data with the plotter.
   */
  SAMRAI::tbox::Array<string> vis_writer(1);
  vis_writer[0] = "Visit";
  if (main_db->keyExists("vis_writer")) {
    vis_writer = main_db->getStringArray("vis_writer");
  }
  bool use_visit = false;
  for (int i = 0; i < vis_writer.getSize(); i++) {
    if (vis_writer[i] == "VisIt") use_visit = true;
  }
#ifdef HAVE_HDF5
  boost::shared_ptr<SAMRAI::appu::VisItDataWriter> visit_writer;
  string vis_filename =
    main_db->getStringWithDefault("vis_filename", base_name);
  if (use_visit) {
    visit_writer = boost::make_shared<SAMRAI::appu::VisItDataWriter>
      (dim,"Visit Writer",vis_filename + ".visit");
    fac.setupPlotter(*visit_writer);
  }
#endif

  /*
   * After creating all objects and initializing their state,
   * we print the input database and variable database contents
   * to the log file.
   */
  SAMRAI::tbox::plog << "\nCheck input data and variables before simulation:"
                     << endl;
  SAMRAI::tbox::plog << "Input database..." << endl;
  input_db->printClassData(SAMRAI::tbox::plog);

  /*
   * Solve.
   */
  fac.solve();


  bool done(false);
  for (int lnum = 0;
       patch_hierarchy->levelCanBeRefined(lnum) && !done; lnum++)
    {
      SAMRAI::tbox::Array<int> tag_buffer(patch_hierarchy->getMaxNumberOfLevels());
      for (int ln = 0; ln < tag_buffer.getSize(); ++ln) {
        tag_buffer[ln] = 1;
      }
      gridding_algorithm->regridAllFinerLevels(
                                               0,
                                               0.0,
                                               tag_buffer);
      SAMRAI::tbox::plog << "Newly adapted hierarchy\n";


      // SAMRAI::tbox::plog << "Adding finner levels with lnum = " << lnum << endl;
      // gridding_algorithm->makeFinerLevel(0.0,true,0);
      // SAMRAI::tbox::plog << "Just added finer levels with lnum = " << lnum << endl;
      done = !(patch_hierarchy->finerLevelExists(lnum));
      fac.solve();
    }

  // {
  //     SAMRAI::tbox::Array<int> tag_buffer(patch_hierarchy->getMaxNumberOfLevels());
  //     for (int ln = 0; ln < tag_buffer.getSize(); ++ln) {
  //       tag_buffer[ln] = 1;
  //     }
  //     gridding_algorithm->regridAllFinerLevels(
  //                                              0,
  //                                              0.0,
  //                                              tag_buffer);
  // }
#ifdef HAVE_HDF5
  /*
   * Plot.
   */
  if (use_visit) {
    visit_writer->writePlotData(patch_hierarchy, 0);
  }
#endif

  SAMRAI::tbox::TimerManager::getManager()->print(SAMRAI::tbox::plog);
}

#endif
