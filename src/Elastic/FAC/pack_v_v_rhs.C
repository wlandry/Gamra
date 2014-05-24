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

/// Write derived data to the given stream.

void
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
  else if ("Fault Correction + RHS" == variable_name)
    {
      v_ptr = boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (patch.getPatchData(v_rhs_id));
    }
  else
    {
      TBOX_ERROR("Unregistered variable name '" << variable_name << "' in\n"
                 << "Elastic::FAC::packDerivedDataIntoDoubleBuffer");
    }

  SAMRAI::pdat::SideData<double>& v = *v_ptr;
  if(d_dim.getValue()==2)
    {
      const SAMRAI::hier::Index ip(1,0), jp(0,1);

      SAMRAI::pdat::CellIterator iend(SAMRAI::pdat::CellGeometry::end(region));
      for (SAMRAI::pdat::CellIterator
             icell(SAMRAI::pdat::CellGeometry::begin(region));
           icell!=iend; ++icell)
        {

          const SAMRAI::pdat::CellIndex &center(*icell);
          const SAMRAI::pdat::SideIndex
            x(center,0,SAMRAI::pdat::SideIndex::Lower),
            y(center,1,SAMRAI::pdat::SideIndex::Lower);
	
          double vx=v(x);
          double vy=v(y);

          if(!offset_vector_on_output)
            {
              vx=(v(x+ip) + vx)/2;
              vy=(v(y+jp) + vy)/2;
            }

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
      SAMRAI::pdat::CellIterator iend(SAMRAI::pdat::CellGeometry::end(region));
      for (SAMRAI::pdat::CellIterator
             icell(SAMRAI::pdat::CellGeometry::begin(region));
           icell!=iend; ++icell)
        {

          const SAMRAI::pdat::CellIndex &center(*icell);
          const SAMRAI::pdat::SideIndex
            x(center,0,SAMRAI::pdat::SideIndex::Lower),
            y(center,1,SAMRAI::pdat::SideIndex::Lower),
            z(center,2,SAMRAI::pdat::SideIndex::Lower);
	
          double vx=v(x);
          double vy=v(y);
          double vz=v(z);

          if(!offset_vector_on_output)
            {
              vx=(v(x+ip) + vx)/2;
              vy=(v(y+jp) + vy)/2;
              vz=(v(z+kp) + vz)/2;
            }

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
          ++buffer;
        }
    }
}
