/*
 *  helix_rebin_ssrb.c
 *
 *  Sinlge slice rebinnnig method for the reconstruction of cone-beam tomography data collected 
 *  with a helical source trajectory.
 *
 * in
 *	cg			ct_geom()
 *	ig			image_geom()
 *	proj	[ns nt na]	cone-beam projection views (line integrals)
 *  short    1 for short-scan fanbeam; 0 for 360 fanbeam 
 *
 *
 * out
 *	img	[nx ny nz]	reconstructed image
 *
 *
 *  Created by Gregory Handy on 7/15/10.
 *  Copyright 2010. All rights reserved.
 *
 */

#include <math.h>
#include "helix_rebin_ssrb.h"

sof helix_rebin_ssrb(
		     float *sino,
		     float *z_orbit_start,
		     cint na_rebin,
		     cint   nz,
		     cfloat dz,
		     cfloat offset_z,
		     cfloat dso,
		     cfloat dsd,
		     cfloat dfs,
		     cint ns,
		     cint nt,
		     cint na,
		     cfloat ds, // horizontal ray spacing
		     cfloat dt, // vertical ray spacing)
		     cfloat offset_s, // channel offset [pixels]
		     cfloat offset_t, // vertical offset on detector [pixels]
		     cfloat pitch,
		     cfloat source_z0,    // first z position of source
		     cfloat orbit,          // total rotation angle [degrees]
		     cfloat orbit_start,    // starting angle [degrees]
		     cfloat  *proj)          // [ns nt na]   I: projection
{
  double num_turns = orbit/360; 
  int num_betas = floor(na/num_turns);
  double myPitch = pitch * nt * (dso / dsd) * dt;
  double currentZ;
  double upper_limit, lower_limit, deltaZ, spoints[ns], tpoints; 
  double scale, alpha; 
  double t_index, y_0, y_1;
  int ia, iz, is, zindex, x_0, x_1;
  double betas[2];
  double z_slices[nz]; 
  double *fan_beam_proj;
  double new_orbit;
  double change_betas;
  int orbit_found;
  double zloc;
  int angle, nextAngle;
  double smax=0, rmax;
  double image_z0;

  Mem0(fan_beam_proj, ns*num_betas*nz);

  for (is = 0; is < ns; ++is) 
    {
      spoints[is]= ds *( is - ((ns-1)/2+offset_s));
      
      if (smax < abs(spoints[is])) 
	{
	  smax = abs(spoints[is]);
	}
    }
  
  if (dfs == 0) 
    {
      rmax = dso * sin(smax/dsd);
    }
  else 
    {
      rmax = dso * sin(atan(smax / dsd));
    }
  
  double delta = asin(rmax/dso);
  double dist = .5 * myPitch * (M_PI + 2 * delta)/(2*M_PI);
  int newNA = na_rebin;
 
  // used to calculate change_betas farther down
  betas[0] = orbit_start * (M_PI/180);
  betas[1] = (orbit_start + orbit * 1/na)*(M_PI/180);
 
  // z locations of the z-slices
  image_z0 = (0-((nz-1)/2+offset_z))*dz;
  for (iz=0; iz < nz; ++iz) 
    {
      z_slices[iz] = image_z0 + iz*dz;
    }
  
  // loop over the different view angles
  for (ia = 0; ia < na; ++ia) 
    {
      currentZ = source_z0 + ia * (myPitch/na) * num_turns;
      upper_limit = currentZ + dist;
      lower_limit = currentZ - dist;
      
      zindex = 0;
      
      // enter into the acceptable range for the z-slices
      while ((zindex < nz) && (z_slices[zindex] < lower_limit))
	{
	  zindex = zindex + 1;
	}
      
      // loop over the acceptable z-slices for the current view angle
      while ((zindex < nz) && (z_slices[zindex] < upper_limit))
	{
	  deltaZ = z_slices[zindex]-currentZ;
	  
	  // Calculate the values of t to be used for each s
	  
	  for (is = 0; is<ns; ++is) 
	    {
	      tpoints = ((spoints[is]*spoints[is]+dsd*dsd)/ (dso*dsd)) * deltaZ;  
	      t_index = (tpoints / dt) + ((nt-1)/ 2+ offset_t);
	      x_0 = floor(t_index);
	      x_1 = 1 + x_0;
  
	      // Scaling factor due to different ray lenghts
	      scale = sqrt((spoints[is]*spoints[is]+dsd*dsd))/sqrt(spoints[is]*spoints[is]+tpoints*tpoints+dsd*dsd);
	      
	      alpha = t_index - x_0;
	      y_1 = proj[is+ns*x_1+ns*nt*ia];
	      y_0 = proj[is+ns*x_0+ns*nt*ia];
	      fan_beam_proj[is+ns*(ia%num_betas)+ns*num_betas*zindex] = scale*(alpha*y_1+(1-alpha)*y_0);   
	    }

	  zindex = zindex + 1;
	}
    }
  
  // Next portion of the code prepares the sinogram for FB2
  
  // decide on beginning orbit for each z-slice
  new_orbit = orbit_start*(M_PI/180);
  
  zloc = source_z0;
  
  change_betas = betas[1]-betas[0];
  for (iz = 0; iz < nz; ++iz)
    {
      orbit_found = 0;
      while (orbit_found == 0)
	{
	  if (z_slices[iz] < zloc + dist)
	    {
	      z_orbit_start[iz] = new_orbit * (180/M_PI);
	      orbit_found = 1;
	    }
	  else
	    {
	      new_orbit = new_orbit + change_betas;
	      zloc = zloc + (myPitch/na) * num_turns;
	    }
	}
    }
  
  // deletes the empty rows of the sinograms and reorganizes the rows 
  for(iz = 0; iz<nz; ++iz)
    {
      angle = round(((z_orbit_start[iz]-orbit_start)*(M_PI/180)) / change_betas);
      for (ia = 0; ia<newNA; ++ia)
	{
	  nextAngle = angle + ia;
	  for (is=0; is<ns; ++is) 
	    {
	      sino[is+ns*ia+ns*newNA*iz] = fan_beam_proj[is+ns*(nextAngle%num_betas)+ns*num_betas*iz];
	    }
	}
    }
  
  // should apply Parker Weight to each sinogram
 
  Free0(fan_beam_proj);
	Ok;
}
