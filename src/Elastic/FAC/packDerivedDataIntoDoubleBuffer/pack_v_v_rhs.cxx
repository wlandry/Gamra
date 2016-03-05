/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FAC.hxx"

void pack_v_v_rhs(double* buffer,
                  const SAMRAI::hier::Patch& patch,
                  const SAMRAI::hier::Box& region,
                  const std::string& variable_name,
                  const int &depth,
                  const SAMRAI::tbox::Dimension &d_dim,
                  const bool &offset_vector_on_output,
                  const int &v_id,
                  const int &v_rhs_id)
{
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr;
  if (variable_name == "Displacement")
    {
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
