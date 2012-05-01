/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Stokes solver 
 *
 ************************************************************************/
#include "FACStokes.h"

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

namespace SAMRAI {

  /*
*************************************************************************
* Write derived data to the given stream.                               *
*************************************************************************
*/
  bool FACStokes::packDerivedDataIntoDoubleBuffer(
                                                  double* buffer,
                                                  const hier::Patch& patch,
                                                  const hier::Box& region,
                                                  const std::string&
                                                  variable_name,
                                                  int depth_id) const
  {
	pdat::CellData<double>::Iterator icell(region);

	tbox::Pointer<pdat::SideData<double> > v_ptr;
	if (variable_name == "Displacement") {
		v_ptr = patch.getPatchData(v_id);
	}
	else if ("Equivalent body force" == variable_name)
	{
		v_ptr = patch.getPatchData(v_rhs_id);
	}
	else
	{
		// Did not register this name.
		TBOX_ERROR(
		"Unregistered variable name '" << variable_name << "' in\n"
		<<
		"FACStokesX::packDerivedDataIntoDoubleBuffer");
	}

	pdat::SideData<double>& v = *v_ptr;
        if(d_dim.getValue()==2)
          {
            const hier::Index ip(1,0), jp(0,1);
            for ( ; icell; icell++) {

              pdat::CellIndex center(*icell);
              const pdat::SideIndex
		x(center,0,pdat::SideIndex::Lower),
		y(center,1,pdat::SideIndex::Lower);
	
              double vx=(v(x+ip) + v(x))/2.;
              double vy=(v(y+jp) + v(y))/2.;

              if (0==depth_id)
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
            const hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
            for ( ; icell; icell++) {

              pdat::CellIndex center(*icell);
              const pdat::SideIndex
		x(center,0,pdat::SideIndex::Lower),
		y(center,1,pdat::SideIndex::Lower),
		z(center,2,pdat::SideIndex::Lower);
	
              double vx=(v(x+ip) + v(x))/2.;
              double vy=(v(y+jp) + v(y))/2.;
              double vz=(v(z+kp) + v(z))/2.;

              if (0==depth_id)
		{
                  *buffer = vx;
		}
              else if (1==depth_id)
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

}
