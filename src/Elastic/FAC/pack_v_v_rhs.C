/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Elastic solver 
 *
 ************************************************************************/
#include "Elastic/FAC.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/solv/SimpleCellRobinBcCoefs.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

/*
*************************************************************************
* Write derived data to the given stream.                               *
*************************************************************************
*/
bool
Elastic::FAC::pack_v_v_rhs(double* buffer,
                           const SAMRAI::hier::Patch& patch,
                           const SAMRAI::hier::Box& region,
                           const std::string& variable_name,
                           const int &depth) const
{
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr;
  if (variable_name == "Displacement") {
    v_ptr = boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
      (patch.getPatchData(v_id));
  }
  else if ("Fault Correction" == variable_name)
    {
      v_ptr = boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (patch.getPatchData(v_rhs_id));
    }
  else
    {
      // Did not register this name.
      TBOX_ERROR("Unregistered variable name '" << variable_name << "' in\n"
                 << "Elastic::FAC::packDerivedDataIntoDoubleBuffer");
    }

  SAMRAI::pdat::SideData<double>& v = *v_ptr;
  if(d_dim.getValue()==2)
    {
      const SAMRAI::hier::Index ip(1,0), jp(0,1);

      SAMRAI::pdat::CellData<double>::iterator iend(region,false);
      for (SAMRAI::pdat::CellData<double>::iterator icell(region,true);
           icell!=iend; ++icell) {

        SAMRAI::pdat::CellIndex center(*icell);
        const SAMRAI::pdat::SideIndex
          x(center,0,SAMRAI::pdat::SideIndex::Lower),
          y(center,1,SAMRAI::pdat::SideIndex::Lower);
	
        double vx=(v(x+ip) + v(x))/2.;
        double vy=(v(y+jp) + v(y))/2.;

        // if(variable_name=="Displacement")
        //   SAMRAI::tbox::plog << variable_name << " "
        //                      << depth << " "
        //                      << center << " "
        //                      << v(x) << " "
        //                      << v(x+ip) << " "
        //                      << v(y) << " "
        //                      << v(y+jp) << " "
        //                      << "\n";

        if (0==depth)
          {
            *buffer = vx;
          }
        else
          {
            *buffer = vy;
          }
        buffer = buffer + 1;
      }
    }
  else
    {
      const SAMRAI::hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
      SAMRAI::pdat::CellData<double>::iterator iend(region,false);
      for (SAMRAI::pdat::CellData<double>::iterator icell(region,true);
           icell!=iend; ++icell) {

        SAMRAI::pdat::CellIndex center(*icell);
        const SAMRAI::pdat::SideIndex
          x(center,0,SAMRAI::pdat::SideIndex::Lower),
          y(center,1,SAMRAI::pdat::SideIndex::Lower),
          z(center,2,SAMRAI::pdat::SideIndex::Lower);
	
        double vx=(v(x+ip) + v(x))/2.;
        double vy=(v(y+jp) + v(y))/2.;
        double vz=(v(z+kp) + v(z))/2.;

        // if(variable_name=="Displacement")
        //   SAMRAI::tbox::plog << variable_name
        //                      << center << " "
        //                      << v(x) << " "
        //                      << v(x+ip) << " "
        //                      << v(y) << " "
        //                      << v(y+jp) << " "
        //                      << v(z) << " "
        //                      << v(z+kp) << " "
        //                      << "\n";

        if (0==depth)
          {
            *buffer = vx;
          }
        else if (1==depth)
          {
            *buffer = vy;
          }
        else
          {
            *buffer = vz;
          }
        buffer = buffer + 1;
      }
    }
  // Return true if this patch has derived data on it.
  // False otherwise.
  return true;
}
