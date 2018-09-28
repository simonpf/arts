/* Copyright (C) 2002-2012
   Sreerekha T.R. <rekha@uni-bremen.de>
   Claudia Emde <claudia.emde@dlr.de>
   Cory Davies <cory@met.ed.ac.uk>
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

/*!
  \file   m_optproperties.cc
  \author Sreerekha T.R. <rekha@uni-bremen.de>, 
          Claudia Emde <claudia.emde@dlr.de>
          Cory Davies <cory@met.ed.ac.uk>
  \date   Mon Jun 10 11:19:11 2002 
  \brief  This filecontains workspace methods for calculating the optical 
  properties for the radiative transfer. 

  Optical properties are the extinction matrix, absorption vector and
  scattering vector.  The optical properties for the gases can be
  calculated with or without Zeeman effect.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cfloat>
#include <cmath>
#include "arts.h"
#include "exceptions.h"
#include "array.h"
#include "matpackIII.h"
#include "matpackVII.h"
#include "logic.h"
#include "interpolation.h"
#include "messages.h"
#include "xml_io.h"
#include "montecarlo.h"
#include "optproperties.h"
#include "math_funcs.h"
#include "sorting.h"
#include "check_input.h"
#include "auto_md.h" 

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;

#define PART_TYPE scat_data[i_ss][i_se].ptype
#define F_DATAGRID scat_data[i_ss][i_se].f_grid
#define T_DATAGRID scat_data[i_ss][i_se].T_grid
#define ZA_DATAGRID scat_data[i_ss][i_se].za_grid
#define AA_DATAGRID scat_data[i_ss][i_se].aa_grid
#define PHA_MAT_DATA scat_data[i_ss][i_se].pha_mat_data
#define EXT_MAT_DATA scat_data[i_ss][i_se].ext_mat_data
#define ABS_VEC_DATA scat_data[i_ss][i_se].abs_vec_data

#define PND_LIMIT 1e-12 // If particle number density is below this value, 
                        // no transformations will be performed


/* Workspace method: Doxygen documentation will be auto-generated */
void pha_mat_sptFromData( // Output:
                         Tensor5& pha_mat_spt,
                         // Input:
                         const ArrayOfArrayOfSingleScatteringData& scat_data,
                         const Vector& scat_za_grid,
                         const Vector& scat_aa_grid,
                         const Index& scat_za_index, // propagation directions
                         const Index& scat_aa_index,
                         const Index& f_index,
                         const Vector& f_grid,
                         const Numeric& rtp_temperature,
                         const Tensor4& pnd_field, 
                         const Index& scat_p_index,
                         const Index& scat_lat_index,
                         const Index& scat_lon_index,
                         const Verbosity& verbosity
                         )
{
  const Index stokes_dim = pha_mat_spt.ncols();
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                        "must be 1,2,3 or 4");
  }

  // Determine total number of scattering elements
  const Index N_se_total = TotalNumberOfElements(scat_data);
  if( N_se_total != pnd_field.nbooks() )
    {
      ostringstream os;
      os << "Total number of scattering elements in scat_data "
         << "inconsistent with size of pnd_field.";
      throw runtime_error(os.str());
    }
  // as pha_mat_spt is typically initiallized from pnd_field, this theoretically
  // checks the same as the runtime_error above. Still, we keep it to be on the
  // save side.
  assert( pha_mat_spt.nshelves() == N_se_total );
  
  // Check that we don't have scat_data_mono here. Only checking the first
  // scat element, assuming the other elements have been processed in the same
  // manner. That's save against having mono data, but not against having
  // individual elements produced with only a single frequency. This, however,
  // has been checked by scat_data_raw reading routines (ScatSpecies/Element*Add/Read).
  // Unsafe, however, remain when ReadXML is used directly or if scat_data(_raw) is
  // (partly) produced from scat_data_singleTmatrix.
  if( scat_data[0][0].f_grid.nelem() < 2 )
  {
      ostringstream os;
      os << "Scattering data seems to be *scat_data_mono* (1 freq point only),\n"
         << "but frequency interpolable data (*scat_data* with >=2 freq points) "
         << "is expected here.";
      throw runtime_error( os.str() );
  }
          
  const Index N_ss = scat_data.nelem();

  // Phase matrix in laboratory coordinate system. Dimensions:
  // [frequency, za_inc, aa_inc, stokes_dim, stokes_dim]
  Tensor5 pha_mat_data_int;


  Index i_se_flat = 0;
  // Loop over scattering species
  for (Index i_ss = 0; i_ss < N_ss; i_ss++)
  {
      const Index N_se = scat_data[i_ss].nelem();

      // Loop over the included scattering elements
      for (Index i_se = 0; i_se < N_se; i_se++)
      {
          // If the particle number density at a specific point in the
          // atmosphere for the i_se scattering element is zero, we don't need
          // to do the transfromation!
          if (pnd_field(i_se_flat, scat_p_index, scat_lat_index, scat_lon_index)
              > PND_LIMIT)
          {

              // First we have to transform the data from the coordinate system
              // used in the database (depending on the kind of ptype) to the
              // laboratory coordinate system.

              // Frequency and temperature interpolation:

              // Container for data at one frequency and one temperature.
              pha_mat_data_int.resize(PHA_MAT_DATA.nshelves(),
                                      PHA_MAT_DATA.nbooks(),
                                      PHA_MAT_DATA.npages(),
                                      PHA_MAT_DATA.nrows(),
                                      PHA_MAT_DATA.ncols());


              // Gridpositions:
              GridPos freq_gp;
              gridpos(freq_gp, F_DATAGRID, f_grid[f_index]);
              GridPos t_gp;
              Vector itw;

              Index ti=-1;

              if( PHA_MAT_DATA.nvitrines() == 1 ) // just 1 T_grid element
              {
                  ti=0;
              }
              else if( rtp_temperature < 0. ) // coding for 'not interpolate, but
                                              // pick one temperature'
              {
                if( rtp_temperature > -10. )      // lowest T-point
                {
                  ti = 0;
                }
                else if( rtp_temperature > -20. ) // highest T-point
                {
                  ti = T_DATAGRID.nelem()-1;
                }
                else                              // median T-point
                {
                  ti = T_DATAGRID.nelem()/2;
                }
              }

              if( ti<0 ) // temperature interpolation
              {
                  ostringstream os;
                  os << "In pha_mat_sptFromData.\n"
                     << "The temperature grid of the scattering data does not\n"
                     << "cover the atmospheric temperature at cloud location.\n"
                     << "The data should include the value T = "
                     << rtp_temperature << " K.";
                  chk_interpolation_grids( os.str(), T_DATAGRID, 
                                           rtp_temperature );

                  gridpos(t_gp, T_DATAGRID, rtp_temperature);

                  // Interpolation weights:
                  itw.resize(4);
                  interpweights(itw, freq_gp, t_gp);

                  for (Index i_za_sca = 0;
                       i_za_sca < PHA_MAT_DATA.nshelves(); i_za_sca++)
                    for (Index i_aa_sca = 0;
                         i_aa_sca < PHA_MAT_DATA.nbooks(); i_aa_sca++)
                      for (Index i_za_inc = 0;
                           i_za_inc < PHA_MAT_DATA.npages(); i_za_inc++)
                        for (Index i_aa_inc = 0;
                             i_aa_inc < PHA_MAT_DATA.nrows(); i_aa_inc++)
                          for (Index i = 0; i < PHA_MAT_DATA.ncols(); i++)
                            // Interpolation of phase matrix:
                            pha_mat_data_int(i_za_sca, i_aa_sca,
                                             i_za_inc, i_aa_inc, i) =
                              interp(itw,
                                     PHA_MAT_DATA(joker, joker,
                                     i_za_sca, i_aa_sca, i_za_inc, i_aa_inc, i),
                                     freq_gp, t_gp);
              }
              else
              {
                  // Interpolation weights:
                  itw.resize(2);
                  interpweights(itw, freq_gp);
                  for (Index i_za_sca = 0;
                       i_za_sca < PHA_MAT_DATA.nshelves(); i_za_sca++)
                    for (Index i_aa_sca = 0;
                         i_aa_sca < PHA_MAT_DATA.nbooks(); i_aa_sca++)
                      for (Index i_za_inc = 0;
                           i_za_inc < PHA_MAT_DATA.npages(); i_za_inc++)
                        for (Index i_aa_inc = 0;
                             i_aa_inc < PHA_MAT_DATA.nrows(); i_aa_inc++)
                          for (Index i = 0; i < PHA_MAT_DATA.ncols(); i++)
                            // Interpolation of phase matrix:
                            pha_mat_data_int(i_za_sca, i_aa_sca,
                                             i_za_inc, i_aa_inc, i) =
                              interp(itw,
                                     PHA_MAT_DATA(joker, ti,
                                     i_za_sca, i_aa_sca, i_za_inc, i_aa_inc, i),
                                     freq_gp);
              }

              // Do the transformation into the laboratory coordinate system.
              for (Index za_inc_idx = 0; za_inc_idx < scat_za_grid.nelem();
                   za_inc_idx ++)
              {
                  for (Index aa_inc_idx = 0; aa_inc_idx < scat_aa_grid.nelem();
                       aa_inc_idx ++)
                  {
                      pha_matTransform(pha_mat_spt(i_se_flat,
                                                   za_inc_idx, aa_inc_idx,
                                                   joker, joker),
                                       pha_mat_data_int,
                                       ZA_DATAGRID, AA_DATAGRID,
                                       PART_TYPE, scat_za_index, scat_aa_index,
                                       za_inc_idx, 
                                       aa_inc_idx, scat_za_grid, scat_aa_grid,
                                       verbosity);
                  }
              }
          }
          i_se_flat++;
      }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pha_mat_sptFromDataDOITOpt(// Output:
                                Tensor5& pha_mat_spt,
                                // Input:
                                const ArrayOfTensor7& pha_mat_sptDOITOpt,
                                const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                                const Index& doit_za_grid_size,
                                const Vector& scat_aa_grid,
                                const Index& scat_za_index, // propagation directions
                                const Index& scat_aa_index,
                                const Numeric& rtp_temperature,
                                const Tensor4&  pnd_field, 
                                const Index& scat_p_index,
                                const Index&  scat_lat_index,
                                const Index& scat_lon_index,
                                const Verbosity&)
{
  const Index N_se_total = TotalNumberOfElements(scat_data_mono);

  if( N_se_total != pnd_field.nbooks() )
    {
      ostringstream os;
      os << "Total number of scattering elements in scat_data_mono "
         << "inconsistent with size of pnd_field.";
      throw runtime_error(os.str());
    }
  // as pha_mat_spt is typically initiallized from pnd_field, this theoretically
  // checks the same as the runtime_error above. Still, we keep it to be on the
  // save side.
  assert( pha_mat_spt.nshelves() == N_se_total );

  // atmosphere_dim = 3
  if (pnd_field.ncols() > 1)
    {
      assert(pha_mat_sptDOITOpt.nelem() == N_se_total);
      // Assuming that if the size is o.k. for one scattering element, it will
      // also be o.k. for the other scattering elements. 
      assert(pha_mat_sptDOITOpt[0].nlibraries() ==
             scat_data_mono[0][0].T_grid.nelem());
      assert(pha_mat_sptDOITOpt[0].nvitrines() == doit_za_grid_size);
      assert(pha_mat_sptDOITOpt[0].nshelves() == scat_aa_grid.nelem() );
      assert(pha_mat_sptDOITOpt[0].nbooks() == doit_za_grid_size);
      assert(pha_mat_sptDOITOpt[0].npages() == scat_aa_grid.nelem()); 
    }
  
  // atmosphere_dim = 1, only zenith angle required for scattered directions. 
  else if ( pnd_field.ncols() == 1 )
    {
      //assert(is_size(scat_theta, doit_za_grid_size, 1,
      //                doit_za_grid_size, scat_aa_grid.nelem()));
      
      assert(pha_mat_sptDOITOpt.nelem() == TotalNumberOfElements(scat_data_mono));
      // Assuming that if the size is o.k. for one scattering element, it will
      // also be o.k. for the other scattering elements. 
      assert(pha_mat_sptDOITOpt[0].nlibraries() == scat_data_mono[0][0].T_grid.nelem());
      assert(pha_mat_sptDOITOpt[0].nvitrines() == doit_za_grid_size);
      assert(pha_mat_sptDOITOpt[0].nshelves() == 1);
      assert(pha_mat_sptDOITOpt[0].nbooks() == doit_za_grid_size);
      assert(pha_mat_sptDOITOpt[0].npages() == scat_aa_grid.nelem()); 
    }
  
  assert(doit_za_grid_size > 0);
  
  // Check that we do indeed have scat_data_mono here. Only checking the first
  // scat element, assuming the other elements have been processed in the same
  // manner. That's save against having scat_data here if that originated from
  // scat_data_raw reading routines (ScatSpecies/Element*Add/Read), it's not safe
  // against data read by ReadXML directly or if scat_data(_raw) has been (partly)
  // produced from scat_data_singleTmatrix. That would be too costly here,
  // though.
  // Also, we can't check here whether data is at the correct frequency since we
  // don't know f_grid and f_index here (we could pass it in, though).
  if( scat_data_mono[0][0].f_grid.nelem() > 1 )
  {
      ostringstream os;
      os << "Scattering data seems to be *scat_data* (several freq points),\n"
         << "but *scat_data_mono* (1 freq point only) is expected here.";
      throw runtime_error( os.str() );
  }
  
  // Create equidistant zenith angle grid
  Vector za_grid;
  nlinspace(za_grid, 0, 180, doit_za_grid_size);  
  
  const Index N_ss = scat_data_mono.nelem();
  const Index stokes_dim = pha_mat_spt.ncols();
  
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                        "must be 1,2,3 or 4");
  }
  
  GridPos T_gp;
  Vector itw(2);
  
  // Initialisation
  pha_mat_spt = 0.;

  Index i_se_flat = 0;

  // Do the transformation into the laboratory coordinate system.
  for (Index i_ss = 0; i_ss < N_ss; i_ss++)
  {
      const Index N_se = scat_data_mono[i_ss].nelem();

      for (Index i_se = 0; i_se < N_se; i_se++)
      {
          // If the particle number density at a specific point in the
          // atmosphere for the i_se scattering element is zero, we don't need
          // to do the transformation!
          if (pnd_field(i_se_flat, scat_p_index, scat_lat_index, scat_lon_index)
              > PND_LIMIT) //TRS
          {

              Index nT = scat_data_mono[i_ss][i_se].pha_mat_data.nvitrines();
              Index ti=-1;

              if( nT == 1 ) // just 1 T_grid element
              {
                  ti=0;
              }
              else if( rtp_temperature < 0. ) // coding for 'not interpolate, but
                                              // pick one temperature'
              {
                if( rtp_temperature > -10. )      // lowest T-point
                {
                  ti = 0;
                }
                else if( rtp_temperature > -20. ) // highest T-point
                {
                  ti = nT-1;
                }
                else                              // median T-point
                {
                  ti = nT/2;
                }
              }
              else
              {
                  ostringstream os;
                  os << "In pha_mat_sptFromDataDOITOpt.\n"
                     << "The temperature grid of the scattering data does not\n"
                     << "cover the atmospheric temperature at cloud location.\n"
                     << "The data should include the value T = "
                     << rtp_temperature << " K.";
                  chk_interpolation_grids( os.str(), 
                                           scat_data_mono[i_ss][i_se].T_grid, 
                                           rtp_temperature);

                  // Gridpositions:
                  gridpos( T_gp, scat_data_mono[i_ss][i_se].T_grid, 
                           rtp_temperature);
                  // Interpolation weights:
                  interpweights(itw, T_gp);
              }



              for (Index za_inc_idx = 0; za_inc_idx < doit_za_grid_size;
                   za_inc_idx ++)
              {
                  for (Index aa_inc_idx = 0; aa_inc_idx < scat_aa_grid.nelem();
                       aa_inc_idx ++)
                  {
                      if( ti<0 ) // Temperature interpolation
                      {
                          for (Index i = 0; i< stokes_dim; i++)
                          {
                              for (Index j = 0; j< stokes_dim; j++)
                              {
                                  pha_mat_spt(i_se_flat, za_inc_idx, aa_inc_idx,
                                              i, j)=
                                  interp(itw,pha_mat_sptDOITOpt[i_se_flat]
                                         (joker, scat_za_index,
                                          scat_aa_index, za_inc_idx,
                                          aa_inc_idx, i, j) , T_gp);
                              }
                          }
                      }
                      else
                      {
                          pha_mat_spt(i_se_flat, za_inc_idx, aa_inc_idx,
                                      joker, joker) =
                          pha_mat_sptDOITOpt[i_se_flat](ti, scat_za_index,
                                                        scat_aa_index, za_inc_idx,
                                                        aa_inc_idx, joker, joker);
                      }
                  }
              }
          }// TRS

          i_se_flat++;
      }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void opt_prop_sptFromData(// Output and Input:
                          ArrayOfPropagationMatrix& ext_mat_spt,
                          ArrayOfStokesVector& abs_vec_spt,
                          // Input:
                          const ArrayOfArrayOfSingleScatteringData& scat_data,
                          const Vector& scat_za_grid,
                          const Vector& scat_aa_grid,
                          const Index& scat_za_index, // propagation directions
                          const Index& scat_aa_index,
                          const Index& f_index,
                          const Vector& f_grid,
                          const Numeric& rtp_temperature, 
                          const Tensor4& pnd_field, 
                          const Index& scat_p_index,
                          const Index& scat_lat_index,
                          const Index& scat_lon_index,
                          const Verbosity& verbosity)
{
  
  const Index N_ss = scat_data.nelem();
  const Numeric za_sca = scat_za_grid[scat_za_index];
  const Numeric aa_sca = scat_aa_grid[scat_aa_index];

  DEBUG_ONLY(const Index N_se_total = TotalNumberOfElements(scat_data);
  if(N_ss)
  {
    assert( ext_mat_spt[0].NumberOfFrequencies() == N_se_total );
    assert( abs_vec_spt[0].NumberOfFrequencies() == N_se_total );
  }
  );
  
  // Check that we don't have scat_data_mono here. Only checking the first
  // scat element, assuming the other elements have been processed in the same
  // manner. That's save against having mono data, but not against having
  // individual elements produced with only a single frequency. This, however,
  // has been checked by scat_data_raw reading routines (ScatSpecies/Element*Add/Read).
  // Unsafe, however, remain when ReadXML is used directly or if scat_data(_raw) is
  // (partly) produced from scat_data_singleTmatrix.
  if( scat_data[0][0].f_grid.nelem() < 2 )
  {
      ostringstream os;
      os << "Scattering data seems to be *scat_data_mono* (1 freq point only),\n"
         << "but frequency interpolable data (*scat_data* with >=2 freq points) "
         << "is expected here.";
      throw runtime_error( os.str() );
  }
          
  // Phase matrix in laboratory coordinate system. Dimensions:
  // [frequency, za_inc, aa_inc, stokes_dim, stokes_dim]
  Tensor3 ext_mat_data_int;
  Tensor3 abs_vec_data_int;
  
  // Initialisation
  ext_mat_spt = 0.;
  abs_vec_spt = 0.;


  Index i_se_flat = 0;
  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss < N_ss; i_ss++)
  {
      const Index N_se = scat_data[i_ss].nelem();

      // Loop over the included scattering elements
      for (Index i_se = 0; i_se < N_se; i_se++)
      {
          // If the particle number density at a specific point in the
          // atmosphere for the i_se scattering element is zero, we don't need
          // to do the transformation

          if (pnd_field(i_se_flat, scat_p_index, scat_lat_index, scat_lon_index)
              > PND_LIMIT)
          {
              // First we have to transform the data from the coordinate system
              // used in the database (depending on the kind of ptype) to the
              // laboratory coordinate system.

              // Frequency interpolation:

              // The data is interpolated on one frequency.
              //
              // Resize the variables for the interpolated data:
              //
              ext_mat_data_int.resize(EXT_MAT_DATA.npages(),
                                      EXT_MAT_DATA.nrows(),
                                      EXT_MAT_DATA.ncols());
              //
              abs_vec_data_int.resize(ABS_VEC_DATA.npages(),
                                      ABS_VEC_DATA.nrows(),
                                      ABS_VEC_DATA.ncols());


              // Gridpositions:
              GridPos freq_gp;
              gridpos(freq_gp, F_DATAGRID, f_grid[f_index]);
              GridPos t_gp;
              Vector itw;

              if ( T_DATAGRID.nelem() > 1)
              {
                  ostringstream os;
                  os << "In opt_prop_sptFromData.\n"
                     << "The temperature grid of the scattering data does not\n"
                     << "cover the atmospheric temperature at cloud location.\n"
                     << "The data should include the value T = "
                     << rtp_temperature << " K.";
                  chk_interpolation_grids( os.str(), T_DATAGRID, 
                                           rtp_temperature );

                  gridpos(t_gp, T_DATAGRID, rtp_temperature);

                  // Interpolation weights:
                  itw.resize(4);
                  interpweights(itw, freq_gp, t_gp);

                  for (Index i_za_sca = 0; i_za_sca < EXT_MAT_DATA.npages();
                       i_za_sca++)
                  {
                      for(Index i_aa_sca = 0; i_aa_sca < EXT_MAT_DATA.nrows();
                          i_aa_sca++)
                      {
                          //
                          // Interpolation of extinction matrix:
                          //
                          for (Index i = 0; i < EXT_MAT_DATA.ncols(); i++)
                          {
                              ext_mat_data_int(i_za_sca, i_aa_sca, i) =
                              interp(itw, EXT_MAT_DATA(joker, joker,
                                                           i_za_sca, i_aa_sca, i),
                                     freq_gp, t_gp);
                          }
                      }
                  }

                  for (Index i_za_sca = 0; i_za_sca < ABS_VEC_DATA.npages();
                       i_za_sca++)
                  {
                      for(Index i_aa_sca = 0; i_aa_sca < ABS_VEC_DATA.nrows();
                          i_aa_sca++)
                      {
                          //
                          // Interpolation of absorption vector:
                          //
                          for (Index i = 0; i < ABS_VEC_DATA.ncols(); i++)
                          {
                              abs_vec_data_int(i_za_sca, i_aa_sca, i) =
                              interp(itw, ABS_VEC_DATA(joker, joker, i_za_sca,
                                                           i_aa_sca, i),
                                     freq_gp, t_gp);
                          }
                      }
                  }
              }
              else
              {
                  // Interpolation weights:
                  itw.resize(2);
                  interpweights(itw, freq_gp);

                  for (Index i_za_sca = 0; i_za_sca < EXT_MAT_DATA.npages();
                       i_za_sca++)
                  {
                      for(Index i_aa_sca = 0; i_aa_sca < EXT_MAT_DATA.nrows();
                          i_aa_sca++)
                      {
                          //
                          // Interpolation of extinction matrix:
                          //
                          for (Index i = 0; i < EXT_MAT_DATA.ncols(); i++)
                          {
                              ext_mat_data_int(i_za_sca, i_aa_sca, i) =
                              interp(itw, EXT_MAT_DATA(joker, 0,
                                                           i_za_sca, i_aa_sca, i),
                                     freq_gp);
                          }
                      }
                  }

                  for (Index i_za_sca = 0; i_za_sca < ABS_VEC_DATA.npages();
                       i_za_sca++)
                  {
                      for(Index i_aa_sca = 0; i_aa_sca < ABS_VEC_DATA.nrows();
                          i_aa_sca++)
                      {
                          //
                          // Interpolation of absorption vector:
                          //
                          for (Index i = 0; i < ABS_VEC_DATA.ncols(); i++)
                          {
                              abs_vec_data_int(i_za_sca, i_aa_sca, i) =
                              interp(itw, ABS_VEC_DATA(joker, 0, i_za_sca,
                                                           i_aa_sca, i),
                                     freq_gp);
                          }
                      }
                  }
              }


              //
              // Do the transformation into the laboratory coordinate system.
              //
              // Extinction matrix:
              //
              ext_matTransform(ext_mat_spt[i_se_flat],
                              ext_mat_data_int,
                              ZA_DATAGRID, AA_DATAGRID, PART_TYPE,
                              za_sca, aa_sca,
                              verbosity);
              //
              // Absorption vector:
              //
              abs_vecTransform(abs_vec_spt[i_se_flat],
                              abs_vec_data_int,
                              ZA_DATAGRID, AA_DATAGRID, PART_TYPE,
                              za_sca, aa_sca, verbosity);
          }

          i_se_flat++;
      }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void opt_prop_sptFromScat_data(// Output and Input:
                          ArrayOfPropagationMatrix& ext_mat_spt,
                          ArrayOfStokesVector& abs_vec_spt,
                          // Input:
                          const ArrayOfArrayOfSingleScatteringData& scat_data,
                          const Index& scat_data_checked,
                          const Vector& scat_za_grid,
                          const Vector& scat_aa_grid,
                          const Index& scat_za_index, // propagation directions
                          const Index& scat_aa_index,
                          const Index& f_index,
                          const Numeric& rtp_temperature, 
                          const Tensor4& pnd_field, 
                          const Index& scat_p_index,
                          const Index& scat_lat_index,
                          const Index& scat_lon_index,
                          const Verbosity& verbosity)
{
  
  if( scat_data_checked != 1 )
    throw runtime_error( "The scattering data must be flagged to have "
                         "passed a consistency check (scat_data_checked=1)." );

  const Index N_ss = scat_data.nelem();
  const Index stokes_dim = ext_mat_spt[0].StokesDimensions();
  const Numeric za_sca = scat_za_grid[scat_za_index];
  const Numeric aa_sca = scat_aa_grid[scat_aa_index];
  
  if (stokes_dim > 4 || stokes_dim < 1)
  {
    throw runtime_error("The dimension of the stokes vector \n"
                        "must be 1,2,3 or 4");
  }

  DEBUG_ONLY(const Index N_se_total = TotalNumberOfElements(scat_data);)
  assert( ext_mat_spt.nelem() == N_se_total );
  assert( abs_vec_spt.nelem() == N_se_total );

  // Phase matrix in laboratory coordinate system. Dimensions:
  // [frequency, za_inc, aa_inc, stokes_dim, stokes_dim]
  Tensor3 ext_mat_data_int;
  Tensor3 abs_vec_data_int;
  
  // Initialisation
  for(auto& pm : ext_mat_spt)
    pm.SetZero();
  for(auto& sv : abs_vec_spt)
    sv.SetZero();

  Index this_f_index;

  Index i_se_flat = 0;
  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss < N_ss; i_ss++)
  {
      const Index N_se = scat_data[i_ss].nelem();

      // Loop over the included scattering elements
      for (Index i_se = 0; i_se < N_se; i_se++)
      {
          // If the particle number density at a specific point in the
          // atmosphere for the i_se scattering element is zero, we don't need
          // to do the transformation

          if (pnd_field(i_se_flat, scat_p_index, scat_lat_index, scat_lon_index)
              > PND_LIMIT)
          {
              // First we have to transform the data from the coordinate system
              // used in the database (depending on the kind of ptype) to the
              // laboratory coordinate system.

              // Resize the variables for the interpolated data (1freq, 1T):
              ext_mat_data_int.resize(EXT_MAT_DATA.npages(),
                                      EXT_MAT_DATA.nrows(),
                                      EXT_MAT_DATA.ncols());
              abs_vec_data_int.resize(ABS_VEC_DATA.npages(),
                                      ABS_VEC_DATA.nrows(),
                                      ABS_VEC_DATA.ncols());

              // Gridpositions and interpolation weights;
              GridPos t_gp;
              Vector itw;
              if ( EXT_MAT_DATA.nbooks()>1 || ABS_VEC_DATA.nbooks()>1 )
                {
                  ostringstream os;
                  os << "In opt_prop_sptFromScat_data.\n"
                     << "The temperature grid of the scattering data does not\n"
                     << "cover the atmospheric temperature at cloud location.\n"
                     << "The data should include the value T = "
                     << rtp_temperature << " K.";
                  chk_interpolation_grids( os.str(), T_DATAGRID, 
                                           rtp_temperature );

                  gridpos(t_gp, T_DATAGRID, rtp_temperature);

                  // Interpolation weights:
                  itw.resize(2);
                  interpweights(itw, t_gp);
                }

              // Frequency extraction and temperature interpolation

              if( EXT_MAT_DATA.nshelves()==1 )
                this_f_index = 0;
              else
                this_f_index = f_index;

              if ( EXT_MAT_DATA.nbooks() > 1)
                {
                  // Interpolation of extinction matrix:
                  for (Index i_za_sca = 0; i_za_sca < EXT_MAT_DATA.npages();
                       i_za_sca++)
                  {
                      for(Index i_aa_sca = 0; i_aa_sca < EXT_MAT_DATA.nrows();
                          i_aa_sca++)
                      {
                          for (Index i = 0; i < EXT_MAT_DATA.ncols(); i++)
                          {
                              ext_mat_data_int(i_za_sca, i_aa_sca, i) =
                              interp(itw, EXT_MAT_DATA(this_f_index, joker,
                                                           i_za_sca, i_aa_sca, i),
                                     t_gp);
                          }
                      }
                  }
                }
              else
                {
                  ext_mat_data_int = EXT_MAT_DATA(this_f_index, 0,
                                                      joker, joker, joker);
                  /*
                  for (Index i_za_sca = 0; i_za_sca < EXT_MAT_DATA.npages();
                       i_za_sca++)
                  {
                      for(Index i_aa_sca = 0; i_aa_sca < EXT_MAT_DATA.nrows();
                          i_aa_sca++)
                      {
                          for (Index i = 0; i < EXT_MAT_DATA.ncols(); i++)
                          {
                              ext_mat_data_int(i_za_sca, i_aa_sca, i) =
                                EXT_MAT_DATA(this_f_index, 0,
                                                 i_za_sca, i_aa_sca, i);
                          }
                      }
                  } */
                }

              if( ABS_VEC_DATA.nshelves()==1 )
                this_f_index = 0;
              else
                this_f_index = f_index;

              if ( ABS_VEC_DATA.nbooks() > 1)
                {
                  // Interpolation of absorption vector:
                  for (Index i_za_sca = 0; i_za_sca < ABS_VEC_DATA.npages();
                       i_za_sca++)
                  {
                      for(Index i_aa_sca = 0; i_aa_sca < ABS_VEC_DATA.nrows();
                          i_aa_sca++)
                      {
                          for (Index i = 0; i < ABS_VEC_DATA.ncols(); i++)
                          {
                              abs_vec_data_int(i_za_sca, i_aa_sca, i) =
                              interp(itw, ABS_VEC_DATA(this_f_index, joker,
                                                           i_za_sca, i_aa_sca, i),
                                     t_gp);
                          }
                      }
                  }
                }
              else
                {
                  abs_vec_data_int = ABS_VEC_DATA(this_f_index, 0,
                                                      joker, joker, joker);
                  /*
                  for (Index i_za_sca = 0; i_za_sca < ABS_VEC_DATA.npages();
                       i_za_sca++)
                  {
                      for(Index i_aa_sca = 0; i_aa_sca < ABS_VEC_DATA.nrows();
                          i_aa_sca++)
                      {
                          for (Index i = 0; i < ABS_VEC_DATA.ncols(); i++)
                          {
                              abs_vec_data_int(i_za_sca, i_aa_sca, i) =
                                ABS_VEC_DATA(this_f_index, 0,
                                                 i_za_sca, i_aa_sca, i);
                          }
                      }
                  } */
                }

              //
              // Do the transformation into the laboratory coordinate system.
              //
              // Extinction matrix:
              ext_matTransform(ext_mat_spt[i_se_flat],
                               ext_mat_data_int,
                               ZA_DATAGRID, AA_DATAGRID, PART_TYPE,
                               za_sca, aa_sca,
                               verbosity);
              // Absorption vector:
              abs_vecTransform(abs_vec_spt[i_se_flat],
                               abs_vec_data_int,
                               ZA_DATAGRID, AA_DATAGRID, PART_TYPE,
                               za_sca, aa_sca, verbosity);
          }
          i_se_flat++;
      }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void opt_prop_bulkCalc(// Output and Input:
                       PropagationMatrix& ext_mat,
                       StokesVector&  abs_vec,
                       // Input:
                 const ArrayOfPropagationMatrix& ext_mat_spt,
                 const ArrayOfStokesVector& abs_vec_spt,
                 const Tensor4& pnd_field,
                 const Index&   scat_p_index,
                 const Index&   scat_lat_index,
                 const Index&   scat_lon_index,
                 const Verbosity&)
{
  Index N_se = abs_vec_spt.nelem();
  //assert( ext_mat_spt.npages()==N_se )
  if( ext_mat_spt.nelem() not_eq N_se )
    {
      ostringstream os;
      os << "Number of scattering elements in *abs_vec_spt* and *ext_mat_spt*\n"
         << "does not agree.";
      throw runtime_error( os.str() );
    }
  
  Index stokes_dim = abs_vec_spt[0].StokesDimensions();
  //assert( ext_mat_spt.ncols()==stokes_dim && ext_mat_spt.nrows()==stokes_dim )
  if( ext_mat_spt[0].StokesDimensions() not_eq stokes_dim )
  {
    ostringstream os;
    os << "*stokes_dim* of *abs_vec_spt* and *ext_mat_spt* does not agree.";
    throw runtime_error( os.str() );
  }
  if (stokes_dim > 4 || stokes_dim < 1){
      ostringstream os;
      os << "The dimension of stokes vector can only be 1, 2, 3, or 4.";
      throw runtime_error( os.str() );
  }

  ext_mat = PropagationMatrix( 1, stokes_dim );
  ext_mat.SetZero();                  // Initialize to zero!
  abs_vec = StokesVector( 1, stokes_dim );
  abs_vec.SetZero();                  // Initialize to zero!

  PropagationMatrix ext_mat_part(1, stokes_dim);
  ext_mat_part.SetZero();
  StokesVector abs_vec_part(1, stokes_dim);
  abs_vec_part.SetZero();

  // this is the loop over the different scattering elements
  for (Index l = 0; l < N_se; l++)
  { 
    abs_vec_part.MultiplyAndAdd(pnd_field(l, scat_p_index, scat_lat_index, scat_lon_index), abs_vec_spt[l]);
    ext_mat_part.MultiplyAndAdd(pnd_field(l, scat_p_index, scat_lat_index, scat_lon_index), ext_mat_spt[l]);
  }

  //Add absorption due single scattering element.
  abs_vec += abs_vec_part;
  //Add extinction matrix due single scattering element to *ext_mat*.
  ext_mat += ext_mat_part;
} 


/* Workspace method: Doxygen documentation will be auto-generated */
void ext_matAddGas(PropagationMatrix& ext_mat,
                   const ArrayOfPropagationMatrix& propmat_clearsky,
                   const Verbosity&)
{
  // Number of Stokes parameters:
  const Index stokes_dim = ext_mat.StokesDimensions();

  // The second dimension of ext_mat must also match the number of
  // Stokes parameters:
  if ( stokes_dim != propmat_clearsky[0].StokesDimensions() )
    throw runtime_error("Col dimension of propmat_clearsky "
                        "inconsistent with col dimension in ext_mat.");

  // Number of frequencies:
  const Index f_dim = ext_mat.NumberOfFrequencies();

  // This must be consistent with the second dimension of
  // propmat_clearsky. Check this:
  if ( f_dim != propmat_clearsky[0].NumberOfFrequencies() )
    throw runtime_error("Frequency dimension of ext_mat and propmat_clearsky\n"
                        "are inconsistent in ext_matAddGas.");

  for(auto& pm : propmat_clearsky)
    ext_mat += pm;
      
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_vecAddGas(StokesVector& abs_vec,
                   const ArrayOfPropagationMatrix& propmat_clearsky,
                   const Verbosity&)
{
  // Number of frequencies:
  const Index f_dim = abs_vec.NumberOfFrequencies();
  const Index stokes_dim = abs_vec.StokesDimensions();
  
  // This must be consistent with the second dimension of
  // propmat_clearsky. Check this:
  if ( f_dim != propmat_clearsky[0].NumberOfFrequencies() )
    throw runtime_error("Frequency dimension of abs_vec and propmat_clearsky\n"
                        "are inconsistent in abs_vecAddGas.");
  if ( stokes_dim != propmat_clearsky[0].StokesDimensions() )
    throw runtime_error("Stokes dimension of abs_vec and propmat_clearsky\n"
                        "are inconsistent in abs_vecAddGas.");
    
  // Loop all frequencies. Of course this includes the special case
  // that there is only one frequency.
  for(auto& pm : propmat_clearsky)
    abs_vec += pm; // Defined to only add to the 

  // We don't have to do anything about higher elements of abs_vec,
  // since scalar gas absorption only influences the first element.
}


/* Workspace method: Doxygen documentation will be auto-generated */
/*
void ext_matAddGasZeeman( Tensor3&      ext_mat,
                          const Tensor3&  ext_mat_zee,
                          const Verbosity&)
{
  // Number of Stokes parameters:
  const Index stokes_dim = ext_mat.ncols();

  // The second dimension of ext_mat must also match the number of
  // Stokes parameters:
  if ( stokes_dim != ext_mat.nrows() )
    throw runtime_error("Row dimension of ext_mat inconsistent with "
                        "column dimension."); 

  for ( Index i=0; i<stokes_dim; ++i )
    {
      for ( Index j=0; j<stokes_dim; ++j )
        {
          // Add the zeeman extinction to extinction matrix.
          ext_mat(joker,i,j) += ext_mat_zee(joker, i, j);
        }
      
    }
}
*/


/* Workspace method: Doxygen documentation will be auto-generated */
/*
void abs_vecAddGasZeeman( Matrix&      abs_vec,
                          const Matrix& abs_vec_zee,
                          const Verbosity&)
{
  // Number of Stokes parameters:
  const Index stokes_dim = abs_vec_zee.ncols();
  // that there is only one frequency.
  for ( Index j=0; j<stokes_dim; ++j )
    {
      abs_vec(joker,j) += abs_vec_zee(joker,j);
    }
}
*/


/* Workspace method: Doxygen documentation will be auto-generated */
void pha_matCalc(Tensor4& pha_mat,
                 const Tensor5& pha_mat_spt,
                 const Tensor4& pnd_field,
                 const Index& atmosphere_dim,
                 const Index& scat_p_index,
                 const Index& scat_lat_index,
                 const Index& scat_lon_index,
                 const Verbosity&)
{

  Index N_se = pha_mat_spt.nshelves();
  Index Nza = pha_mat_spt.nbooks();
  Index Naa = pha_mat_spt.npages();
  Index stokes_dim = pha_mat_spt.nrows();
 
  pha_mat.resize(Nza, Naa, stokes_dim, stokes_dim);

  // Initialisation
  pha_mat = 0.0;

  Index ilat=0;
  Index ilon=0;
  if (atmosphere_dim > 1)
    ilat = scat_lat_index;
  if (atmosphere_dim > 2)
    ilon = scat_lon_index;

  // this is a loop over the different scattering elements
  for (Index pt_index = 0; pt_index < N_se; ++ pt_index)
    // these are loops over zenith angle and azimuth angle
    for (Index za_index = 0; za_index < Nza; ++ za_index)
      for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
        // now the last two loops over the stokes dimension.
        for (Index stokes_index_1 = 0; stokes_index_1 < stokes_dim; 
             ++  stokes_index_1)
          for (Index stokes_index_2 = 0; stokes_index_2 < stokes_dim;
               ++ stokes_index_2)
            //summation of the product of pnd_field and 
            //pha_mat_spt.
            pha_mat(za_index, aa_index, stokes_index_1, stokes_index_2) += 
              ( pha_mat_spt(pt_index, za_index, aa_index,  
                            stokes_index_1, stokes_index_2) * 
                pnd_field(pt_index,scat_p_index, ilat, ilon) );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void scat_dataCheck( //Input:
                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                     const String& check_type,
                     const Numeric& threshold,
                     const Verbosity& verbosity )
{
    CREATE_OUT0;
    CREATE_OUT1;
    CREATE_OUT2;
    //CREATE_OUT3;

    // FIXME:
    // so far, this works for both scat_data and scat_data_raw. Needs to be
    // adapted, though, once we have WSM that can create Z/K/a with different
    // f/T dimensions than scat_data_single.f/T_grid.

    const Index N_ss = scat_data.nelem();

    // 1) any negative values in Z11, K11, or a1? K11>=a1?
    // 2) scat_data containing any NaN?
    // 3) sca_mat norm sufficiently good (int(Z11)~=K11-a1?)

    // Loop over the included scattering species
    out2 << " checking for negative values in Z11, K11, and a1, and for K11<a1\n";
    for (Index i_ss = 0; i_ss < N_ss; i_ss++)
    {
      const Index N_se = scat_data[i_ss].nelem();

      // Loop over the included scattering elements
      for (Index i_se = 0; i_se < N_se; i_se++)
      {
        for (Index f = 0; f < F_DATAGRID.nelem(); f++)
        {
          for (Index zai=0; zai<ABS_VEC_DATA.npages(); zai++)
            for (Index aai=0; aai<ABS_VEC_DATA.nrows(); aai++)
            {
              for (Index t = 0; t < T_DATAGRID.nelem(); t++)
              {
                if( EXT_MAT_DATA(f,t,zai,aai,0)<0 ||
                    ABS_VEC_DATA(f,t,zai,aai,0)<0 )
                  {
                    ostringstream os;
                    os << "Scatt. species #" << i_ss << " element #" << i_se
                       << " contains negative K11 or a1 at f#"
                       << f << ", T#" << t << ", za#" << zai << ", aa#" << aai
                       << "\n";
                    throw runtime_error( os.str() );
                  }
                if( EXT_MAT_DATA(f,t,zai,aai,0) <
                    ABS_VEC_DATA(f,t,zai,aai,0) )
                  {
                    ostringstream os;
                    os << "Scatt. species #" << i_ss << " element #" << i_se
                       << " has K11<a1 at f#"
                       << f << ", T#" << t << ", za#" << zai << ", aa#" << aai
                       << "\n";
                    throw runtime_error( os.str() );
                  }
              }
              // Since allowing pha_mat to have a single T entry only (while
              // T_grid, ext_mat, abs_vec have more), we need a separate T loop
              // for pha_mat
              Index nTpha = PHA_MAT_DATA.nvitrines();
              for (Index t = 0; t < nTpha; t++)
              {
                for (Index zas=0; zas<PHA_MAT_DATA.nshelves(); zas++)
                  for (Index aas=0; aas<PHA_MAT_DATA.nbooks(); aas++)
                    if( PHA_MAT_DATA(f,t,zas,aas,zai,aai,0)<0 )
                    {
                      ostringstream os;
                      os << "Scatt. species #" << i_ss << " element #" << i_se
                         << " contains negative Z11 at f#" << f
                         << ", T#" << t << " (of " << nTpha << "), za_sca#"
                         << zas << ", aa_sca#" << aas << ", za_inc#" << zai
                         << ", aa_inc#" << aai << "\n";
                      throw runtime_error( os.str() );
                    }
              }
            }
        }
      }
    }

    // Loop over the included scattering species
    out2 << " checking for NaN anywhere in Z, K, and a\n";
    for (Index i_ss = 0; i_ss < N_ss; i_ss++)
    {
      const Index N_se = scat_data[i_ss].nelem();

      // Loop over the included scattering elements
      for (Index i_se = 0; i_se < N_se; i_se++)
      {
        for (Index f = 0; f < F_DATAGRID.nelem(); f++)
        {
          for (Index zai=0; zai<ABS_VEC_DATA.npages(); zai++)
            for (Index aai=0; aai<ABS_VEC_DATA.nrows(); aai++)
            {
              for (Index t = 0; t < T_DATAGRID.nelem(); t++)
              {
                for (Index st=0; st<ABS_VEC_DATA.ncols(); st++)
                  if( std::isnan(ABS_VEC_DATA(f,t,zai,aai,st)) )
                  {
                    ostringstream os;
                    os << "Scatt. species #" << i_ss << " element #" << i_se
                       << " contains NaN in abs_vec at f#" << f << ", T#"
                       << t << ", za#" << zai << ", aa#" << aai << ", stokes #"
                       << st << "\n";
                    throw runtime_error( os.str() );
                  }
                for (Index st=0; st<EXT_MAT_DATA.ncols(); st++)
                  if( std::isnan(EXT_MAT_DATA(f,t,zai,aai,st)) )
                  {
                    ostringstream os;
                    os << "Scatt. species #" << i_ss << " element #" << i_se
                       << " contains NaN in ext_mat at f#" << f << ", T#"
                       << t << ", za#" << zai << ", aa#" << aai << ", stokes #"
                       << st << "\n";
                    throw runtime_error( os.str() );
                  }
              }
              Index nTpha = PHA_MAT_DATA.nvitrines();
              for (Index t = 0; t < nTpha; t++)
              {
                for (Index zas=0; zas<PHA_MAT_DATA.nshelves(); zas++)
                  for (Index aas=0; aas<PHA_MAT_DATA.nbooks(); aas++)
                    for (Index st=0; st<PHA_MAT_DATA.ncols(); st++)
                      if( std::isnan(PHA_MAT_DATA(f,t,zas,aas,zai,aai,st)) )
                      {
                        ostringstream os;
                        os << "Scatt. species #" << i_ss << " element #" << i_se
                           << " contains NaN in pha_mat at f#" << f << ", T#"
                           << t << " (of " << nTpha << "), za_sca#" << zas
                           << ", aa_sca#" << aas << ", za_inc#" << zai
                           << ", aa_inc#" << aai << ", stokes #" << "\n";
                        throw runtime_error( os.str() );
                      }
              }
            }
        }
      }
    }

    if( check_type.toupper() == "ALL" )
    {
      // Loop over the included scattering species
      out2 << " checking normalization of scattering matrix\n";
      for (Index i_ss = 0; i_ss < N_ss; i_ss++)
      {
        const Index N_se = scat_data[i_ss].nelem();

        // Loop over the included scattering elements
        for (Index i_se = 0; i_se < N_se; i_se++)
          if( T_DATAGRID.nelem() == PHA_MAT_DATA.nvitrines() )
            switch (PART_TYPE)
            {
              case PTYPE_TOTAL_RND:
              {
                for (Index f = 0; f < F_DATAGRID.nelem(); f++)
                {
                  for (Index t = 0; t < T_DATAGRID.nelem(); t++)
                  {
                    Numeric Csca = AngIntegrate_trapezoid(
                                     PHA_MAT_DATA(f, t, joker, 0, 0, 0, 0),
                                     ZA_DATAGRID);
                    Numeric Cext_data = EXT_MAT_DATA(f,t,0,0,0);
                    //Numeric Cabs = Cext_data - Csca;
                    Numeric Cabs_data = ABS_VEC_DATA(f,t,0,0,0);
                    Numeric Csca_data = Cext_data - Cabs_data;

                    /*
                    out3 << "  Coefficients in database: "
                         << "Cext: " << Cext_data << " Cabs: " << Cabs_data
                         << " Csca: " << Csca_data << "\n"
                         << "  Calculated coefficients: "
                         << "Cabs calc: " << Cabs
                         << " Csca calc: " << Csca << "\n"
                         << "  Deviations "
                         << "Cabs: " << 1e2*Cabs/Cabs_data-1e2
                         << "% Csca: " << 1e2*Csca/Csca_data-1e2
                         << "% Alb: " << (Csca-Csca_data)/Cext_data << "\n";
                    */
                  
                    //if (abs(Csca/Csca_data-1.)*Csca_data/Cext_data > threshold)
                    // below equivalent to the above
                    // (it's actually the (absolute) albedo deviation!)
                    if (abs(Csca-Csca_data)/Cext_data > threshold)
                    {
                      ostringstream os;
                      os << "  Deviations in scat_data too large:\n"
                         << "  scat dev [%] " << 1e2*Csca/Csca_data-1e2
                         << " at nominal (actual) albedo of "
                         << Csca_data/Cext_data << " ("
                         << Csca/Cext_data << ").\n"
                         << "  Check entry for scattering element " << i_se
                         << " of scattering species " << i_ss << " at "
                         << f << ".frequency and " << t << ".temperature!\n";
                      throw runtime_error( os.str() );
                    }
                  }
                }
                break;
              }
                    
              case PTYPE_AZIMUTH_RND:
              {
                for (Index f = 0; f < F_DATAGRID.nelem(); f++)
                {
                  for (Index t = 0; t < T_DATAGRID.nelem(); t++)
                  {
                    for (Index iza = 0; iza < ABS_VEC_DATA.npages(); iza++)
                    {
                      Numeric Csca = 2 * AngIntegrate_trapezoid(
                                       PHA_MAT_DATA(f, t, joker, joker, iza, 0, 0),
                                       ZA_DATAGRID, AA_DATAGRID );
                      Numeric Cext_data = EXT_MAT_DATA(f,t,iza,0,0);
                      //Numeric Cabs = Cext_data - Csca;
                      Numeric Cabs_data = ABS_VEC_DATA(f,t,iza,0,0);
                      Numeric Csca_data = Cext_data - Cabs_data;

                      /*
                      out3 << "  Coefficients in database: "
                           << "Cext: " << Cext_data << " Cabs: " << Cabs_data
                           << " Csca: " << Csca_data << "\n"
                           << "  Calculated coefficients: "
                           << "Cabs calc: " << Cabs
                           << " Csca calc: " << Csca << "\n"
                           << "  Deviations "
                           << "Cabs: " << 1e2*Cabs/Cabs_data-1e2
                           << "% Csca: " << 1e2*Csca/Csca_data-1e2
                           << "% Alb: " << (Csca-Csca_data)/Cext_data << "\n";
                      */

                      //if (abs(Csca/Csca_data-1.)*Csca_data/Cext_data > threshold)
                      // below equivalent to the above
                      // (it's actually the (absolute) albedo deviation!)
                      if (abs(Csca-Csca_data)/Cext_data > threshold)
                      {
                        ostringstream os;
                        os << "  Deviations in scat_data too large:\n"
                           << "  scat dev [%] " << 1e2*Csca/Csca_data-1e2
                           << " at nominal (actual) albedo of "
                           << Csca_data/Cext_data << " ("
                           << Csca/Cext_data << ").\n"
                           << "  Check entry for scattering element " << i_se
                           << " of scattering species " << i_ss << " at "
                           << f << ". frequency, " << t << ". temperature, and "
                           << iza << ". incident polar angle!\n";
                        throw runtime_error( os.str() );
                      }
                    }
                  }
                }
                break;
              }

              default:
              {
                out0 << "  WARNING:\n"
                     << "  scat_data consistency check not implemented (yet?!) for\n"
                     << "  ptype " << PART_TYPE << "!\n";
              }
            }
          else
            out2 << "  WARNING:\n"
                 << "  scat_data norm check can not be performed for pha_mat-only"
                 << " T-reduced scattering elements\n"
                 << "  as found in scatt element #" << i_se
                 << " of scatt species #" << i_ss << "!\n";
      }
    }
    else if (check_type.toupper() == "SANE")
    {
      out1 << "  WARNING:\n"
           << "  Normalization check on pha_mat switched off.\n"
           << "  Scattering solution might be wrong.\n";
    }
    else
    {
        ostringstream os;
        os << "Invalid value for argument *check_type*: '" << check_type << "'.\n";
        os << "Valid values are 'all' or 'none'.";
        throw runtime_error( os.str() );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void DoitScatteringDataPrepare(Workspace& ws,//Output:
                               ArrayOfTensor7& pha_mat_sptDOITOpt,
                               ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                               Tensor7& pha_mat_doit,
                               //Input:
                               const Index& doit_za_grid_size,
                               const Vector& scat_aa_grid,
                               const ArrayOfArrayOfSingleScatteringData& scat_data,
                               const Index& scat_data_checked,
                               const Index& f_index,
                               const Index& atmosphere_dim,
                               const Index& stokes_dim,
                               const Tensor3& t_field,
                               const ArrayOfIndex& cloudbox_limits,
                               const Tensor4& pnd_field,
                               const Agenda& pha_mat_spt_agenda,
                               const Verbosity& verbosity)
{
  if( scat_data_checked != 1 )
    throw runtime_error( "The scattering data must be flagged to have "
                         "passed a consistency check (scat_data_checked=1)." );

  // Number of azimuth angles.
  const Index Naa = scat_aa_grid.nelem();
  Vector grid_stepsize(2);
  grid_stepsize[0] = 180./(Numeric)(doit_za_grid_size - 1);
  grid_stepsize[1] = 360./(Numeric)(Naa - 1);

  // Initialize variable *pha_mat_spt*
  Tensor5 pha_mat_spt_local(pnd_field.nbooks(), doit_za_grid_size,
                            scat_aa_grid.nelem(), stokes_dim, stokes_dim, 0.);
  Tensor4 pha_mat_local(doit_za_grid_size, Naa,
                        stokes_dim, stokes_dim, 0.);
  Tensor6 pha_mat_local_out(cloudbox_limits[1]-cloudbox_limits[0]+1,doit_za_grid_size,
                            doit_za_grid_size, Naa, stokes_dim, stokes_dim, 0.);

  // Interpolate all the data in frequency
  scat_data_monoExtract(scat_data_mono, scat_data, f_index, verbosity);
  
  // For 1D calculation the scat_aa dimension is not required:
  Index N_aa_sca;
  if  ( atmosphere_dim == 1 )
    N_aa_sca = 1;
  else
    N_aa_sca = scat_aa_grid.nelem();
  
  Vector za_grid;
  nlinspace(za_grid, 0, 180, doit_za_grid_size);

  assert( scat_data.nelem() == scat_data_mono.nelem() );
  
  const Index N_ss = scat_data.nelem();
  // FIXME: We need this still as a workspace variable because pha_mat_spt_agenda
  // contains a WS method that requires it as input
  pha_mat_sptDOITOpt.resize(TotalNumberOfElements(scat_data));

  Index i_se_flat = 0;
  for (Index i_ss = 0; i_ss < N_ss; i_ss++)
  {
      const Index N_se = scat_data[i_ss].nelem();

      for (Index i_se = 0; i_se < N_se; i_se++)
      {
          Index N_T = scat_data_mono[i_ss][i_se].T_grid.nelem();
          pha_mat_sptDOITOpt[i_se_flat].resize(N_T, doit_za_grid_size, N_aa_sca,
                                          doit_za_grid_size, scat_aa_grid.nelem(),
                                          stokes_dim, stokes_dim);

          //    Initialize:
          pha_mat_sptDOITOpt[i_se_flat]= 0.;

          // Calculate all scattering angles for all combinations of incoming
          // and scattered directions and interpolation.
          for (Index t_idx = 0; t_idx < N_T; t_idx ++)
          {
              // These are the scattered directions as called in *scat_field_calc*
              for (Index za_sca_idx = 0; za_sca_idx < doit_za_grid_size; za_sca_idx ++)
              {
                  for (Index aa_sca_idx = 0; aa_sca_idx < N_aa_sca; aa_sca_idx ++)
                  {
                      // Integration is performed over all incoming directions
                      for (Index za_inc_idx = 0; za_inc_idx < doit_za_grid_size;
                           za_inc_idx ++)
                      {
                          for (Index aa_inc_idx = 0; aa_inc_idx <
                               scat_aa_grid.nelem();
                               aa_inc_idx ++)
                          {
                              pha_matTransform(pha_mat_sptDOITOpt[i_se_flat]
                                               (t_idx, za_sca_idx,
                                                aa_sca_idx, za_inc_idx, aa_inc_idx,
                                                joker, joker),
                                               scat_data_mono[i_ss][i_se].
                                               pha_mat_data
                                               (0,t_idx,joker,joker,joker,
                                                joker,joker),
                                               scat_data_mono[i_ss][i_se].za_grid,
                                               scat_data_mono[i_ss][i_se].aa_grid,
                                               scat_data_mono[i_ss][i_se].ptype,
                                               za_sca_idx,
                                               aa_sca_idx,
                                               za_inc_idx,
                                               aa_inc_idx,
                                               za_grid,
                                               scat_aa_grid,
                                               verbosity);
                          }
                      }
                  }
              }
          }
          
          i_se_flat++;
      }
  }
    // Interpolate phase matrix to current grid
    pha_mat_doit.resize(cloudbox_limits[1] - cloudbox_limits[0] + 1,
                        doit_za_grid_size, N_aa_sca, doit_za_grid_size,
                        Naa, stokes_dim, stokes_dim);
    pha_mat_doit = 0;

    if ( atmosphere_dim == 1)
    {
        Index scat_aa_index_local = 0;
        
        // Get pha_mat at the grid positions
        // Since atmosphere_dim = 1, there is no loop over lat and lon grids
        for (Index p_index = 0; p_index<=cloudbox_limits[1]-cloudbox_limits[0] ;
             p_index++)
        {
            Numeric rtp_temperature_local =
            t_field(p_index + cloudbox_limits[0], 0, 0);
            //There is only loop over zenith angle grid ; no azimuth angle grid.
            for (Index scat_za_index_local = 0;
                 scat_za_index_local < doit_za_grid_size; scat_za_index_local ++)
            {
                // Dummy index
                Index index_zero = 0;
                
                // Calculate the phase matrix of individual scattering elements
                pha_mat_spt_agendaExecute(ws, pha_mat_spt_local,
                                          scat_za_index_local,
                                          index_zero,
                                          index_zero,
                                          p_index,
                                          scat_aa_index_local,
                                          rtp_temperature_local,
                                          pha_mat_spt_agenda);
                
                // Sum over all scattering elements
                pha_matCalc(pha_mat_local, pha_mat_spt_local, pnd_field,
                            atmosphere_dim, p_index, 0,
                            0, verbosity);
                pha_mat_doit(p_index ,scat_za_index_local, 0,
                                  joker ,joker ,joker ,joker)
                = pha_mat_local;
            }
        }
    }   
}


/* Workspace method: Doxygen documentation will be auto-generated */
void scat_dataCalc(ArrayOfArrayOfSingleScatteringData& scat_data,
             const ArrayOfArrayOfSingleScatteringData& scat_data_raw,
             const Vector& f_grid,
             const Index& interp_order,
             const Verbosity&)
// FIXME: when we allow K, a, Z to be on different f and T grids, their use in
// the scatt solvers needs to be reviewed again and adaptedto this!
{
  //Extrapolation factor:
  //const Numeric extpolfac = 0.5;

  Index nf=f_grid.nelem();
  
  // Check, whether single scattering data contains the right frequencies:
  // The check was changed to allow extrapolation at the boundaries of the 
  // frequency grid.
  const String which_interpolation="scat_data_raw.f_grid to f_grid";
  for (Index i_ss = 0; i_ss<scat_data_raw.nelem(); i_ss++)
  {
    for (Index i_se = 0; i_se<scat_data_raw[i_ss].nelem(); i_se++)
    {
      // Check for the special case that ssd.f_grid f_grid have only one
      // element. If identical, that's  fine. If not, throw error.
      if (scat_data_raw[i_ss][i_se].f_grid.nelem()==1 && nf==1)
        if ( !is_same_within_epsilon(scat_data_raw[i_ss][i_se].f_grid[0],
                                     f_grid[0],2*DBL_EPSILON) )
        {
          ostringstream os;
          os << "There is a problem with the grids for the following "
             << "interpolation:\n" << which_interpolation << "\n"
             << "If original grid has only 1 element, the new grid must also have\n"
             << "only a single element and hold the same value as the original grid.";
          throw runtime_error( os.str() );
        }

      // check with extrapolation
      chk_interpolation_grids(which_interpolation,
                              scat_data_raw[i_ss][i_se].f_grid, f_grid,
                              interp_order);
    }
  }

  //Initialise scat_data
  scat_data.resize(scat_data_raw.nelem());

  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss<scat_data_raw.nelem(); i_ss++)
  {
    const Index N_se = scat_data_raw[i_ss].nelem();

    //Initialise scat_data
    scat_data[i_ss].resize(N_se);

    // Loop over the included scattering elements
    for (Index i_se = 0; i_se < N_se; i_se++)
    {
      //Stuff that doesn't need interpolating
      PART_TYPE   = scat_data_raw[i_ss][i_se].ptype;
      F_DATAGRID  = f_grid;
      T_DATAGRID  = scat_data_raw[i_ss][i_se].T_grid;
      ZA_DATAGRID = scat_data_raw[i_ss][i_se].za_grid;
      AA_DATAGRID = scat_data_raw[i_ss][i_se].aa_grid;

      //Sizing SSD data containers
      PHA_MAT_DATA.resize(nf,
                              scat_data_raw[i_ss][i_se].pha_mat_data.nvitrines(),
                              scat_data_raw[i_ss][i_se].pha_mat_data.nshelves(),
                              scat_data_raw[i_ss][i_se].pha_mat_data.nbooks(),
                              scat_data_raw[i_ss][i_se].pha_mat_data.npages(),
                              scat_data_raw[i_ss][i_se].pha_mat_data.nrows(),
                              scat_data_raw[i_ss][i_se].pha_mat_data.ncols());
      EXT_MAT_DATA.resize(nf,
                              scat_data_raw[i_ss][i_se].ext_mat_data.nbooks(),
                              scat_data_raw[i_ss][i_se].ext_mat_data.npages(),
                              scat_data_raw[i_ss][i_se].ext_mat_data.nrows(),
                              scat_data_raw[i_ss][i_se].ext_mat_data.ncols());
      ABS_VEC_DATA.resize(nf,
                              scat_data_raw[i_ss][i_se].abs_vec_data.nbooks(),
                              scat_data_raw[i_ss][i_se].abs_vec_data.npages(),
                              scat_data_raw[i_ss][i_se].abs_vec_data.nrows(),
                              scat_data_raw[i_ss][i_se].abs_vec_data.ncols());

      const bool single_se_fgrid=(scat_data_raw[i_ss][i_se].f_grid.nelem()==1);
      if( !single_se_fgrid )
      {
        // Gridpositions:
        ArrayOfGridPosPoly freq_gp( nf );;
        gridpos_poly(freq_gp, scat_data_raw[i_ss][i_se].f_grid, f_grid, interp_order);

        // Interpolation weights:
        Matrix itw(nf, interp_order+1);
        interpweights(itw, freq_gp);

        //Phase matrix data
        for (Index t_index = 0; t_index < scat_data_raw[i_ss][i_se].pha_mat_data.nvitrines(); t_index ++)
        {
          for (Index i_za_sca = 0; i_za_sca < scat_data_raw[i_ss][i_se].pha_mat_data.nshelves(); i_za_sca++)
          {
            for (Index i_aa_sca = 0; i_aa_sca < scat_data_raw[i_ss][i_se].pha_mat_data.nbooks(); i_aa_sca++)
            {
              for (Index i_za_inc = 0; i_za_inc < scat_data_raw[i_ss][i_se].pha_mat_data.npages(); i_za_inc++)
              {
                for (Index i_aa_inc = 0; i_aa_inc < scat_data_raw[i_ss][i_se].pha_mat_data.nrows(); i_aa_inc++)
                {
                  for (Index i = 0; i < scat_data_raw[i_ss][i_se].pha_mat_data.ncols(); i++)
                  {
                    interp(scat_data[i_ss][i_se].pha_mat_data(joker, t_index,
                                                              i_za_sca, i_aa_sca,
                                                              i_za_inc, i_aa_inc, i),
                           itw,
                           scat_data_raw[i_ss][i_se].pha_mat_data(joker, t_index,
                                                                  i_za_sca, i_aa_sca,
                                                                  i_za_inc, i_aa_inc, i),
                           freq_gp);
                  }
                }
              }
            }
          }
        }

        //Extinction matrix data
        for (Index t_index = 0; t_index < scat_data_raw[i_ss][i_se].ext_mat_data.nbooks(); t_index ++)
        {
          for (Index i_za_sca = 0; i_za_sca < scat_data_raw[i_ss][i_se].ext_mat_data.npages(); i_za_sca++)
          {
            for (Index i_aa_sca = 0; i_aa_sca < scat_data_raw[i_ss][i_se].ext_mat_data.nrows(); i_aa_sca++)
            {
              for (Index i = 0; i < scat_data_raw[i_ss][i_se].ext_mat_data.ncols(); i++)
              {
                interp(scat_data[i_ss][i_se].ext_mat_data(joker, t_index,
                                                          i_za_sca, i_aa_sca, i),
                       itw,
                       scat_data_raw[i_ss][i_se].ext_mat_data(joker, t_index,
                                                              i_za_sca, i_aa_sca, i),
                       freq_gp);
              }
            }
          }
        }

        //Absorption vector data
        for (Index t_index = 0; t_index < scat_data_raw[i_ss][i_se].abs_vec_data.nbooks(); t_index ++)
        {
          for (Index i_za_sca = 0; i_za_sca < scat_data_raw[i_ss][i_se].abs_vec_data.npages(); i_za_sca++)
          {
            for (Index i_aa_sca = 0; i_aa_sca < scat_data_raw[i_ss][i_se].abs_vec_data.nrows(); i_aa_sca++)
            {
              for (Index i = 0; i < scat_data_raw[i_ss][i_se].abs_vec_data.ncols(); i++)
              {
                interp(scat_data[i_ss][i_se].abs_vec_data(joker, t_index,
                                                          i_za_sca, i_aa_sca, i),
                       itw,
                       scat_data_raw[i_ss][i_se].abs_vec_data(joker, t_index,
                                                              i_za_sca, i_aa_sca, i),
                       freq_gp);
              }
            }
          }
        }
      }
      else
      {
        assert( nf==1 );
        // we do only have one f_grid value in old and new data (and they have
        // been confirmed to be the same), hence only need to copy over/reassign
        // the data.
        scat_data[i_ss][i_se].pha_mat_data = scat_data_raw[i_ss][i_se].pha_mat_data;
        scat_data[i_ss][i_se].ext_mat_data = scat_data_raw[i_ss][i_se].ext_mat_data;
        scat_data[i_ss][i_se].abs_vec_data = scat_data_raw[i_ss][i_se].abs_vec_data;
      }
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void scat_dataReduceT(ArrayOfArrayOfSingleScatteringData& scat_data,
                const Index& i_ss,
                const Numeric& T,
                const Index& interp_order,
                const Index& phamat_only,
                const Numeric& threshold,
                const Verbosity&)
{
  // We are directly acting on the scat_data entries, modifying them
  // individually. That is, we don't need to resize these arrays. Only the
  // pha_mat and probably ext_mat and abs_vec Tensors (in the latter case also
  // T_grid!).

  // Check that species i_ss exists at all in scat_data
  const Index nss=scat_data.nelem();
  if( nss <= i_ss )
  {
    ostringstream os;
    os << "Can not T-reduce scattering species #" << i_ss << ".\n"
       << "*scat_data* contains only " << nss << " scattering species.";
    throw runtime_error( os.str() );
  }

  // Loop over the included scattering elements
  for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++)
  {
    // At very first check validity of the scatt elements ptype (so far we only
    // handle PTYPE_TOTAL_RND and PTYPE_AZIMUTH_RND).
    if( PART_TYPE != PTYPE_TOTAL_RND and PART_TYPE != PTYPE_AZIMUTH_RND )
    {
      ostringstream os;
      os << "Only ptypes " << PTYPE_TOTAL_RND << " and " << PTYPE_AZIMUTH_RND
         << " can be handled.\n"
         << "Scattering element #" << i_se << " has ptype " << PART_TYPE
         << ".";
      throw runtime_error( os.str() );
    }

    // If ssd.T_grid already has only a single point, we do nothing.
    // This is not necessarily expected behaviour. BUT, it is in line with
    // previous use (that if nT==1, then assume ssd constant in T).
    Index nT = T_DATAGRID.nelem();
    if( nT>1 )
    {
      // Check, that we not have data that has already been T-reduced (in
      // pha_mat only. complete ssd T-reduce should have been sorted away
      // already above).
      if( PHA_MAT_DATA.nvitrines()!=nT )
      {
        ostringstream os;
        os << "Single scattering data of scat element #" << i_se
           << " of scat species #" << i_ss << "\n"
           << "seems to have undergone some temperature grid manipulation in\n"
           << "*pha_mat_data* already. That can not be done twice!";
        throw runtime_error( os.str() );
      }

      // Check that ext_mat and abs_vec have the same temp dimensions as T_grid.
      // This should always be true, if not it's a bug not a user mistake, hence
      // use assert.
      assert( EXT_MAT_DATA.nbooks()==nT and ABS_VEC_DATA.nbooks()==nT );
      
      // Check that T_grid is consistent with requested interpolation order
      ostringstream ost;
      ost << "Scattering data temperature interpolation for\n"
         << "scat element #" << i_se << " of scat species #" << i_ss << ".";
      chk_interpolation_grids( ost.str(), T_DATAGRID, T, interp_order );

      // Gridpositions:
      GridPosPoly gp_T;
      gridpos_poly( gp_T,  T_DATAGRID, T, interp_order );
      Vector itw(interp_order+1);
      interpweights(itw, gp_T);

      //Sizing of temporary SSD data containers
      Tensor7 phamat_tmp(PHA_MAT_DATA.nlibraries(),
                         1,
                         PHA_MAT_DATA.nshelves(),
                         PHA_MAT_DATA.nbooks(),
                         PHA_MAT_DATA.npages(),
                         PHA_MAT_DATA.nrows(),
                         PHA_MAT_DATA.ncols(),
                         0.);
      Tensor5 extmat_tmp(EXT_MAT_DATA.nshelves(),
                         1,
                         EXT_MAT_DATA.npages(),
                         EXT_MAT_DATA.nrows(),
                         EXT_MAT_DATA.ncols(),
                         0.);
      Tensor5 absvec_tmp(ABS_VEC_DATA.nshelves(),
                         1,
                         ABS_VEC_DATA.npages(),
                         ABS_VEC_DATA.nrows(),
                         ABS_VEC_DATA.ncols(),
                         0.);

      // a1) temp interpol of pha mat
      //We have to apply the interpolation separately for each of the pha_mat
      //entries, i.e. loop over all remaining size dimensions
      //We don't apply any transformation here, but interpolate the actual
      //stored ssd (i.e. not the 4x4matrices, but the 7-16 elements separately).
      for( Index i_f=0; i_f<PHA_MAT_DATA.nlibraries(); i_f++ )
        for( Index i_za1=0; i_za1<PHA_MAT_DATA.nshelves(); i_za1++ )
          for( Index i_aa1=0; i_aa1<PHA_MAT_DATA.nbooks(); i_aa1++ )
            for( Index i_za2=0; i_za2<PHA_MAT_DATA.npages(); i_za2++ )
              for( Index i_aa2=0; i_aa2<PHA_MAT_DATA.nrows(); i_aa2++ )
                for( Index i_st=0; i_st<PHA_MAT_DATA.ncols(); i_st++ )
                  phamat_tmp(i_f,0,i_za1,i_aa1,i_za2,i_aa2,i_st) =
                    interp(itw,
                           PHA_MAT_DATA(i_f,joker,i_za1,i_aa1,i_za2,i_aa2,i_st),
                           gp_T);

      // a2) temp interpol of ext and abs.
      //We do that regardless of whether they should be reduced or not, because
      //we need them also for norm checking / renorming.
      for( Index i_f=0; i_f<EXT_MAT_DATA.nshelves(); i_f++ )
        for( Index i_za=0; i_za<EXT_MAT_DATA.npages(); i_za++ )
          for( Index i_aa=0; i_aa<EXT_MAT_DATA.nrows(); i_aa++ )
          {
            for( Index i_st=0; i_st<EXT_MAT_DATA.ncols(); i_st++ )
              extmat_tmp(i_f,0,i_za,i_aa,i_st) =
                interp(itw,
                       EXT_MAT_DATA(i_f,joker,i_za,i_aa,i_st),
                       gp_T);
            for( Index i_st=0; i_st<ABS_VEC_DATA.ncols(); i_st++ )
              absvec_tmp(i_f,0,i_za,i_aa,i_st) =
                interp(itw,
                       ABS_VEC_DATA(i_f,joker,i_za,i_aa,i_st),
                       gp_T);
          }

      // Norm & other consistency checks.
      // All done separately for PTYPE_TOTAL_RND and PTYPE_AZIMUTH_RND (the
      // latter needs to loop over scat_za_inc).
      //
      // b) calculate norm of T-reduced pha mat
      // c) check pha mat norm vs. sca xs from ext-abs at T_interpol
      // d) Ensure that T-reduced data is consistent/representative of all data.
      //    and throw error/disallow reduction if sca xs varying too much.
      // d1) in case of pha_mat only reduction, the scat xs (pha_mat norm) needs
      // to be consistent with the sca xs from over the ext/abs T_grid. This is
      // essentially an energy conservation issue. That is, we should be as
      // strict here as with pha_mat norm deviations in general (however, should
      // we allow the norm at T_grid to deviate by a threshold from the
      // T_interpol norm (that should be a little looser) or from the ext-abs
      // derived expected norm?).
      // d2) in case of all-ssp reduction, the data should still be
      // representative. b)&c) ensure data consistency in itself, making this
      // rather an error on the SSP as such. Hence, we might be a little more
      // loose here.
      // d) the resulting check for d1) and d2) is the same (ext-abs sca xs at
      // T_interpol vs ext-abs sca xs at T_grid), but we use different
      // thresholds.
      //
      // FIXME? 
      // Regarding b)&c) should we also calc norm of original-T pha mats? To get
      // a measure how strong they already deviate from expected norm (As we can
      // not be more exact here than what is already in the original data...).
      // On the other hand, a certain accuracy should be guaranteed from
      // scat_dataCheck already.
      // Hence, for now we skip that (but maybe added later when proves
      // necessary).
      //
      // FIXME?
      // Regarding d1), we could alternatively make sure here that norm at
      // T_interpol is good. And later on ignore any deviations between norm and
      // ext-abs sca xs and instead blindly renorm to expected norm (would that
      // be ok? correct norm here, doesn't imply correct norm at whatever scat
      // angle grids the user is applying. for that, we could in place also calc
      // the original-data norm. but that might be expensive (as we can't do
      // that from ext-abs sca xs, because we don't know to which T that refers.
      // that would go away if we'd actually store pha_mat normed to 1 or 4Pi.
      // but that's prob not going to happen. is it? Another option would be to
      // introduce an additional T_grid, eg T_grid_phamat.). which we actually
      // want to avoid :-/

      Numeric this_threshold;
      String errmsg;
      if( phamat_only )
      {
        this_threshold = threshold;
        errmsg = "T-reduced *pha_mat_data* norm (=sca xs) deviates too "
                 "much from non-reduced *ext_mat_data* and *abs_vec_data*:";
      }
      else
      {
        this_threshold = 2*threshold;
        errmsg = "T-reduced *scat_data* deviates too much from original "
                 "*scat_data*:";
      }

      // The norm-check code is copied and slightly adapted from scat_dataCheck.
      // Might be better to make a functon out of this and use in both places
      // for consistency.
      //
      // FIXME: no checks on higher Stokes elements are done. Should there?
      // Which?
      switch (PART_TYPE)
      {
        case PTYPE_TOTAL_RND:
        {
          for (Index f = 0; f < F_DATAGRID.nelem(); f++)
          {
            // b) calculate norm of T-reduced pha mat
            Numeric Csca = AngIntegrate_trapezoid(
                             phamat_tmp(f, 0, joker, 0, 0, 0, 0),
                             ZA_DATAGRID);
            Numeric Cext_data = extmat_tmp(f,0,0,0,0);
            //Numeric Cabs = Cext_data - Csca;
            Numeric Cabs_data = absvec_tmp(f,0,0,0,0);
            Numeric Csca_data = Cext_data - Cabs_data;

            /*
            cout << "  Coefficients in data: "
                 << "Cext: " << Cext_data << " Cabs: " << Cabs_data
                 << " Csca: " << Csca_data << "\n"
                 << "  Calculated coefficients: "
                 << "Cabs calc: " << Cabs
                 << " Csca calc: " << Csca << "\n"
                 << "  Deviations "
                 << "Cabs: " << 1e2*Cabs/Cabs_data-1e2
                 << "% Csca: " << 1e2*Csca/Csca_data-1e2
                 << "% Alb: " << (Csca-Csca_data)/Cext_data << "\n";
            */
                  
            // c) check pha mat norm vs. sca xs from ext-abs at T_interpol (as
            // albedo dev check)
            if (abs(Csca-Csca_data)/Cext_data > threshold)
            {
              ostringstream os;
              os << "  Deviations in T-reduced scat_data too large:\n"
                 << "  scat dev [%] " << 1e2*Csca/Csca_data-1e2
                 << " at nominal (actual) albedo of "
                 << Csca_data/Cext_data << " ("
                 << Csca/Cext_data << ").\n"
                 << "  Problem occurs for scattering element #" << i_se
                 << " at " << f << ".frequency!\n";
              throw runtime_error( os.str() );
            }
            Numeric norm_dev = (Csca-Csca)/Cext_data;

            // d) Ensure that T-reduced data is consistent/representative of all data.
            // below use theoretical (ext-abs derived) sca xs as reference.
            Csca = Csca_data;
            for (Index t = 0; t < T_DATAGRID.nelem(); t++)
            {
              Cext_data = EXT_MAT_DATA(f,t,0,0,0);
              Csca_data = Cext_data - ABS_VEC_DATA(f,t,0,0,0);
              Numeric xs_dev = (Csca-Csca_data)/Cext_data;
              if (abs(norm_dev+(Csca-Csca_data)/Cext_data) > this_threshold)
                cout << "Accumulated deviation (abs(" << norm_dev << "+" << xs_dev
                     << ")=" << abs(norm_dev+xs_dev) << " exceeding threshold ("
                     << this_threshold << ").\n";
              if (abs(Csca-Csca_data)/Cext_data > this_threshold)
              {
                ostringstream os;
                os << "  " << errmsg << "\n"
                   << "  scat dev [%] " << 1e2*Csca/Csca_data-1e2
                   << " at nominal (actual) albedo of "
                   << Csca_data/Cext_data << " ("
                   << Csca/Cext_data << ").\n"
                   << "  Problem occurs for scattering element #" << i_se
                   << " at " << f << ".frequency and " << t << ".temperature!\n";
                throw runtime_error( os.str() );
              }
            }
          }
          break;
        }
                    
        case PTYPE_AZIMUTH_RND:
        {
          for (Index f = 0; f < F_DATAGRID.nelem(); f++)
          {
            for (Index iza = 0; iza < ABS_VEC_DATA.npages(); iza++)
            {
              // b) calculate norm of T-reduced pha mat
              Numeric Csca = 2 * AngIntegrate_trapezoid(
                               phamat_tmp(f, 0, joker, joker, iza, 0, 0),
                               ZA_DATAGRID, AA_DATAGRID);
              Numeric Cext_data = extmat_tmp(f,0,iza,0,0);
              //Numeric Cabs = Cext_data - Csca;
              Numeric Cabs_data = absvec_tmp(f,0,iza,0,0);
              Numeric Csca_data = Cext_data - Cabs_data;

              /*
              cout << "  Coefficients in data: "
                   << "Cext: " << Cext_data << " Cabs: " << Cabs_data
                   << " Csca: " << Csca_data << "\n"
                   << "  Calculated coefficients: "
                   << "Cabs calc: " << Cabs
                   << " Csca calc: " << Csca << "\n"
                   << "  Deviations "
                   << "Cabs: " << 1e2*Cabs/Cabs_data-1e2
                   << "% Csca: " << 1e2*Csca/Csca_data-1e2
                   << "% Alb: " << (Csca-Csca_data)/Cext_data << "\n";
              */

              // c) check pha mat norm vs. sca xs from ext-abs at T_interpol (as
              // albedo dev check)
              if (abs(Csca-Csca_data)/Cext_data > threshold)
              {
                ostringstream os;
                os << "  Deviations in T-reduced scat_data too large:\n"
                   << "  scat dev [%] " << 1e2*Csca/Csca_data-1e2
                   << " at nominal (actual) albedo of "
                   << Csca_data/Cext_data << " ("
                   << Csca/Cext_data << ").\n"
                   << "  Problem occurs for scattering element #" << i_se
                   << " at " << f << ".frequency, and "
                   << iza << ". incident polar angle!\n";
                throw runtime_error( os.str() );
              }

              // d) Ensure that T-reduced data is consistent/representative of all data.
              // below use theoretical (ext-abs derived) sca xs as reference.
              Csca = Csca_data;
              for (Index t = 0; t < T_DATAGRID.nelem(); t++)
              {
                Cext_data = EXT_MAT_DATA(f,t,0,0,0);
                Csca_data = Cext_data - ABS_VEC_DATA(f,t,0,0,0);
                if (abs(Csca-Csca_data)/Cext_data > this_threshold)
                {
                  ostringstream os;
                  os << "  " << errmsg << "\n"
                     << "  scat dev [%] " << 1e2*Csca/Csca_data-1e2
                     << " at nominal (actual) albedo of "
                     << Csca_data/Cext_data << " ("
                     << Csca/Cext_data << ").\n"
                     << "  Problem occurs for scattering element #" << i_se
                     << " at " << f << ".frequency and " << t << ".temperature, and "
                     << iza << ". incident polar angle!\n";
                  throw runtime_error( os.str() );
                }
              }
            }
          }
          break;
        }

        default:
        {
          // other ptype cases already excluded above. i.e. we shouldn't end up
          // here. If we do, that's a bug.
          assert(0);
        }
      }

      PHA_MAT_DATA = phamat_tmp;
      //We don't need to reset the scat element's grids!
      //Except for T_grid in the case that we reduce ALL three ssd variables.
      if( !phamat_only )
      {
        T_DATAGRID.resize(1);
        T_DATAGRID = T;
        EXT_MAT_DATA = extmat_tmp;
        ABS_VEC_DATA = absvec_tmp;
      }
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void scat_data_monoCalc(ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                        const ArrayOfArrayOfSingleScatteringData& scat_data,
                        const Vector& f_grid,
                        const Index& f_index,
                        const Verbosity&)
{
  //Extrapolation factor:
  //const Numeric extpolfac = 0.5;
  
  // Check, whether single scattering data contains the right frequencies:
  // The check was changed to allow extrapolation at the boundaries of the 
  // frequency grid.
  for (Index h = 0; h<scat_data.nelem(); h++)
  {
    for (Index i = 0; i<scat_data[h].nelem(); i++)
    {
      // check with extrapolation
      chk_interpolation_grids("scat_data.f_grid to f_grid",
                              scat_data[h][i].f_grid,
                              f_grid[f_index]);
    }
  }

  //Initialise scat_data_mono
  scat_data_mono.resize(scat_data.nelem());

  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss<scat_data.nelem(); i_ss++)
  {
      const Index N_se = scat_data[i_ss].nelem();

      //Initialise scat_data_mono
      scat_data_mono[i_ss].resize(N_se);

      // Loop over the included scattering elements
      for (Index i_se = 0; i_se < N_se; i_se++)
      {
          // Gridpositions:
          GridPos freq_gp;
          gridpos(freq_gp, F_DATAGRID, f_grid[f_index]);

          // Interpolation weights:
          Vector itw(2);
          interpweights(itw, freq_gp);

          //Stuff that doesn't need interpolating
          scat_data_mono[i_ss][i_se].ptype=PART_TYPE;
          scat_data_mono[i_ss][i_se].f_grid.resize(1);
          scat_data_mono[i_ss][i_se].f_grid=f_grid[f_index];
          scat_data_mono[i_ss][i_se].T_grid=scat_data[i_ss][i_se].T_grid;
          scat_data_mono[i_ss][i_se].za_grid=ZA_DATAGRID;
          scat_data_mono[i_ss][i_se].aa_grid=AA_DATAGRID;

          //Phase matrix data
          scat_data_mono[i_ss][i_se].pha_mat_data.resize(1,
                                                   PHA_MAT_DATA.nvitrines(),
                                                   PHA_MAT_DATA.nshelves(),
                                                   PHA_MAT_DATA.nbooks(),
                                                   PHA_MAT_DATA.npages(),
                                                   PHA_MAT_DATA.nrows(),
                                                   PHA_MAT_DATA.ncols());

          for (Index t_index = 0; t_index < PHA_MAT_DATA.nvitrines(); t_index ++)
          {
              for (Index i_za_sca = 0; i_za_sca < PHA_MAT_DATA.nshelves();
                   i_za_sca++)
              {
                  for (Index i_aa_sca = 0; i_aa_sca < PHA_MAT_DATA.nbooks();
                       i_aa_sca++)
                  {
                      for (Index i_za_inc = 0; i_za_inc <
                           PHA_MAT_DATA.npages();
                           i_za_inc++)
                      {
                          for (Index i_aa_inc = 0;
                               i_aa_inc < PHA_MAT_DATA.nrows();
                               i_aa_inc++)
                          {
                              for (Index i = 0; i < PHA_MAT_DATA.ncols(); i++)
                              {
                                  scat_data_mono[i_ss][i_se].pha_mat_data(0, t_index,
                                                                    i_za_sca,
                                                                    i_aa_sca,
                                                                    i_za_inc,
                                                                    i_aa_inc, i) =
                                  interp(itw,
                                         PHA_MAT_DATA(joker, t_index,
                                                          i_za_sca,
                                                          i_aa_sca, i_za_inc,
                                                          i_aa_inc, i),
                                         freq_gp);
                              }
                          }
                      }
                  }
              }
              //Extinction matrix data
              scat_data_mono[i_ss][i_se].ext_mat_data.resize(1, T_DATAGRID.nelem(),
                                                       EXT_MAT_DATA.npages(),
                                                       EXT_MAT_DATA.nrows(),
                                                       EXT_MAT_DATA.ncols());
              for (Index i_za_sca = 0; i_za_sca < EXT_MAT_DATA.npages();
                   i_za_sca++)
              {
                  for(Index i_aa_sca = 0; i_aa_sca < EXT_MAT_DATA.nrows();
                      i_aa_sca++)
                  {
                      //
                      // Interpolation of extinction matrix:
                      //
                      for (Index i = 0; i < EXT_MAT_DATA.ncols(); i++)
                      {
                          scat_data_mono[i_ss][i_se].ext_mat_data(0, t_index,
                                                            i_za_sca, i_aa_sca, i)
                          = interp(itw, EXT_MAT_DATA(joker, t_index, i_za_sca,
                                                         i_aa_sca, i),
                                   freq_gp);
                      }
                  }
              }
              //Absorption vector data
              scat_data_mono[i_ss][i_se].abs_vec_data.resize(1, T_DATAGRID.nelem(),
                                                       ABS_VEC_DATA.npages(),
                                                       ABS_VEC_DATA.nrows(),
                                                       ABS_VEC_DATA.ncols());
              for (Index i_za_sca = 0; i_za_sca < ABS_VEC_DATA.npages() ;
                   i_za_sca++)
              {
                  for(Index i_aa_sca = 0; i_aa_sca < ABS_VEC_DATA.nrows();
                      i_aa_sca++)
                  {
                      //
                      // Interpolation of absorption vector:
                      //
                      for (Index i = 0; i < ABS_VEC_DATA.ncols(); i++)
                      {
                          scat_data_mono[i_ss][i_se].abs_vec_data(0, t_index, i_za_sca,
                                                            i_aa_sca, i) =
                          interp(itw, ABS_VEC_DATA(joker, t_index, i_za_sca,
                                                       i_aa_sca, i),
                                 freq_gp);
                      }
                  }
              }
          }
      }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void scat_data_monoExtract(ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                        const ArrayOfArrayOfSingleScatteringData& scat_data,
                        const Index& f_index,
                        const Verbosity&)
{
  //Initialise scat_data_mono
  scat_data_mono.resize(scat_data.nelem());

  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss<scat_data.nelem(); i_ss++)
    {
      const Index N_se = scat_data[i_ss].nelem();

      //Initialise scat_data_mono
      scat_data_mono[i_ss].resize(N_se);

      // Loop over the included scattering elements
      for (Index i_se = 0; i_se < N_se; i_se++)
      {

        Index nf = F_DATAGRID.nelem();
        if( nf == 1 )
        {
          scat_data_mono[i_ss][i_se] = scat_data[i_ss][i_se];
        }
        else
        {
          //Stuff that doesn't need interpolating
          scat_data_mono[i_ss][i_se].ptype=PART_TYPE;
          scat_data_mono[i_ss][i_se].T_grid=T_DATAGRID;
          scat_data_mono[i_ss][i_se].za_grid=ZA_DATAGRID;
          scat_data_mono[i_ss][i_se].aa_grid=AA_DATAGRID;

          scat_data_mono[i_ss][i_se].f_grid.resize(1);
          scat_data_mono[i_ss][i_se].f_grid = F_DATAGRID[f_index];

          Index this_f_index;

          //Phase matrix data
          /*scat_data_mono[i_ss][i_se].pha_mat_data.resize(1,
                                                 PHA_MAT_DATA.nvitrines(),
                                                 PHA_MAT_DATA.nshelves(),
                                                 PHA_MAT_DATA.nbooks(),
                                                 PHA_MAT_DATA.npages(),
                                                 PHA_MAT_DATA.nrows(),
                                                 PHA_MAT_DATA.ncols());*/
          this_f_index = (PHA_MAT_DATA.nlibraries()==1) ? 0 : f_index;
          scat_data_mono[i_ss][i_se].pha_mat_data =
              PHA_MAT_DATA(Range(this_f_index,1),joker,joker,joker,joker,joker,joker);

          //Extinction matrix data
          /*scat_data_mono[i_ss][i_se].ext_mat_data.resize(1, T_DATAGRID.nelem(),
                                                     EXT_MAT_DATA.npages(),
                                                     EXT_MAT_DATA.nrows(),
                                                     EXT_MAT_DATA.ncols());*/
          this_f_index = (EXT_MAT_DATA.nshelves()==1) ? 0 : f_index;
          scat_data_mono[i_ss][i_se].ext_mat_data =
              EXT_MAT_DATA(Range(this_f_index,1),joker,joker,joker,joker);

          //Absorption vector data
          /*scat_data_mono[i_ss][i_se].abs_vec_data.resize(1, T_DATAGRID.nelem(),
                                                     ABS_VEC_DATA.npages(),
                                                     ABS_VEC_DATA.nrows(),
                                                     ABS_VEC_DATA.ncols());*/
          this_f_index = (ABS_VEC_DATA.nshelves()==1) ? 0 : f_index;
          scat_data_mono[i_ss][i_se].abs_vec_data =
              ABS_VEC_DATA(Range(this_f_index,1),joker,joker,joker,joker);

        }
      }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void opt_prop_sptFromMonoData(// Output and Input:
                              ArrayOfPropagationMatrix& ext_mat_spt,
                              ArrayOfStokesVector& abs_vec_spt,
                              // Input:
                              const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                              const Vector& scat_za_grid,
                              const Vector& scat_aa_grid,
                              const Index& scat_za_index, // propagation directions
                              const Index& scat_aa_index,
                              const Numeric& rtp_temperature,
                              const Tensor4& pnd_field, 
                              const Index& scat_p_index,
                              const Index& scat_lat_index,
                              const Index& scat_lon_index,
                              const Verbosity& verbosity)
{
  DEBUG_ONLY(const Index N_se_total = TotalNumberOfElements(scat_data_mono);)
  const Index stokes_dim = ext_mat_spt[0].StokesDimensions();
  const Numeric za_sca = scat_za_grid[scat_za_index];
  const Numeric aa_sca = scat_aa_grid[scat_aa_index];
  
  if (stokes_dim > 4 or stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                         "must be 1,2,3 or 4");
  }
  
  assert( ext_mat_spt.nelem() == N_se_total );
  assert( abs_vec_spt.nelem() == N_se_total );

  // Check that we do indeed have scat_data_mono here. Only checking the first
  // scat element, assuming the other elements have been processed in the same
  // manner. That's save against having scat_data here if that originated from
  // scat_data_raw reading routines (ScatSpecies/Element*Add/Read), it's not safe
  // against data read by ReadXML directly or if scat_data(_raw) has been (partly)
  // produced from scat_data_singleTmatrix. That would be too costly here,
  // though.
  // Also, we can't check here whether data is at the correct frequency since we
  // don't know f_grid and f_index here (we could pass it in, though).
  if( scat_data_mono[0][0].f_grid.nelem() > 1 )
  {
      ostringstream os;
      os << "Scattering data seems to be *scat_data* (several freq points),\n"
         << "but *scat_data_mono* (1 freq point only) is expected here.";
      throw runtime_error( os.str() );
  }

  // Initialisation
  for(auto& pm : ext_mat_spt)
    pm.SetZero();
  for(auto& av : abs_vec_spt)
    av.SetZero();

  GridPos t_gp;
  
  Vector itw(2);

  Index i_se_flat = 0;
  // Loop over the included scattering elements
  for (Index i_ss = 0; i_ss < scat_data_mono.nelem(); i_ss++)
  {
      // Loop over the included scattering elements
      for (Index i_se = 0; i_se < scat_data_mono[i_ss].nelem(); i_se++)
      {
          // If the particle number density at a specific point in the
          // atmosphere for the i_se scattering element is zero, we don't need
          // to do the transformation!
          if (pnd_field(i_se_flat, scat_p_index, scat_lat_index, scat_lon_index)
              > PND_LIMIT)
          {

              // First we have to transform the data from the coordinate system
              // used in the database (depending on the kind of ptype) to the
              // laboratory coordinate system.

              //
              // Do the transformation into the laboratory coordinate system.
              //
              // Extinction matrix:
              //
              Index ext_npages = scat_data_mono[i_ss][i_se].ext_mat_data.npages();
              Index ext_nrows = scat_data_mono[i_ss][i_se].ext_mat_data.nrows();
              Index ext_ncols = scat_data_mono[i_ss][i_se].ext_mat_data.ncols();
              Index abs_npages = scat_data_mono[i_ss][i_se].abs_vec_data.npages();
              Index abs_nrows = scat_data_mono[i_ss][i_se].abs_vec_data.nrows();
              Index abs_ncols = scat_data_mono[i_ss][i_se].abs_vec_data.ncols();

              //Check that scattering data temperature range covers required temperature
              ConstVectorView t_grid = scat_data_mono[i_ss][i_se].T_grid;

              if (t_grid.nelem() > 1)
              {
                  ostringstream os;
                  os << "In opt_prop_sptFromMonoData.\n"
                     << "The temperature grid of the scattering data does not\n"
                     << "cover the atmospheric temperature at cloud location.\n"
                     << "The data should include the value T = "
                     << rtp_temperature << " K.";
                  chk_interpolation_grids( os.str(), t_grid, rtp_temperature );

                  //interpolate over temperature
                  Tensor3 ext_mat_data1temp(ext_npages,ext_nrows,ext_ncols);
                  gridpos(t_gp, t_grid, rtp_temperature);
                  interpweights(itw, t_gp);
                  for (Index i_p = 0; i_p < ext_npages ; i_p++)
                  {
                      for (Index i_r = 0; i_r < ext_nrows ; i_r++)
                      {
                          for (Index i_c = 0; i_c < ext_ncols ; i_c++)
                          {
                              ext_mat_data1temp(i_p,i_r,i_c) = 
                                interp(itw,
                                       scat_data_mono[i_ss][i_se].ext_mat_data( 
                                       0,joker,i_p,i_r,i_c),
                                       t_gp);
                          }
                      }
                  }
                  ext_matTransform(ext_mat_spt[i_se_flat],
                                   ext_mat_data1temp,
                                   scat_data_mono[i_ss][i_se].za_grid,
                                   scat_data_mono[i_ss][i_se].aa_grid,
                                   scat_data_mono[i_ss][i_se].ptype,
                                   za_sca, aa_sca,
                                   verbosity);
              }
              else
              {
                  ext_matTransform(ext_mat_spt[i_se_flat],
                                   scat_data_mono[i_ss][i_se].ext_mat_data(0, 0, joker, joker, joker),
                                   scat_data_mono[i_ss][i_se].za_grid,
                                   scat_data_mono[i_ss][i_se].aa_grid,
                                   scat_data_mono[i_ss][i_se].ptype,
                                   za_sca, aa_sca,
                                   verbosity);
              }
              //
              // Absorption vector:
              //

              if (t_grid.nelem() > 1)
              {
                  Tensor3 abs_vec_data1temp(abs_npages, abs_nrows, abs_ncols);
                  //interpolate over temperature
                  for (Index i_p = 0; i_p < abs_npages; i_p++)
                  {
                      for (Index i_r = 0; i_r < abs_nrows; i_r++)
                      {
                          for (Index i_c = 0; i_c < abs_ncols; i_c++)
                          {
                              abs_vec_data1temp(i_p, i_r, i_c) =
                                      interp(itw,
                                             scat_data_mono[i_ss][i_se].abs_vec_data(
                                                     0, joker, i_p, i_r, i_c),
                                             t_gp);
                          }
                      }
                  }
                  abs_vecTransform(abs_vec_spt[i_se_flat],
                                   abs_vec_data1temp,
                                   scat_data_mono[i_ss][i_se].za_grid,
                                   scat_data_mono[i_ss][i_se].aa_grid,
                                   scat_data_mono[i_ss][i_se].ptype,
                                   za_sca, aa_sca,
                                   verbosity);
              }
              else
              {
                  abs_vecTransform(abs_vec_spt[i_se_flat],
                                   scat_data_mono[i_ss][i_se].abs_vec_data(0, 0, joker, joker, joker),
                                   scat_data_mono[i_ss][i_se].za_grid,
                                   scat_data_mono[i_ss][i_se].aa_grid,
                                   scat_data_mono[i_ss][i_se].ptype,
                                   za_sca, aa_sca,
                                   verbosity);
              }
          }
          
          i_se_flat++;
      }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pha_mat_sptFromMonoData(// Output:
                             Tensor5& pha_mat_spt,
                             // Input:
                             const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                             const Index& doit_za_grid_size,
                             const Vector& scat_aa_grid,
                             const Index& scat_za_index, // propagation directions
                             const Index& scat_aa_index,
                             const Numeric& rtp_temperature,
                             const Tensor4& pnd_field, 
                             const Index& scat_p_index,
                             const Index& scat_lat_index,
                             const Index& scat_lon_index,
                             const Verbosity& verbosity)
{
  Vector za_grid;
  nlinspace(za_grid, 0, 180, doit_za_grid_size); 

  const Index N_se_total = TotalNumberOfElements(scat_data_mono);
  if( N_se_total != pnd_field.nbooks() )
    {
      ostringstream os;
      os << "Total number of scattering elements in *scat_data_mono* "
         << "inconsistent with size of pnd_field.";
      throw runtime_error(os.str());
    }
  // as pha_mat_spt is typically initialized from pnd_field, this theoretically
  // checks the same as the runtime_error above. Still, we keep it to be on the
  // save side.
  assert( pha_mat_spt.nshelves() == N_se_total );

  const Index stokes_dim = pha_mat_spt.ncols();
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                        "must be 1,2,3 or 4");
  }
  
  // Check that we do indeed have scat_data_mono here. Only checking the first
  // scat element, assuming the other elements have been processed in the same
  // manner. That's save against having scat_data here if that originated from
  // scat_data_raw reading routines (ScatSpecies/Element*Add/Read), it's not safe
  // against data read by ReadXML directly or if scat_data(_raw) has been (partly)
  // produced from scat_data_singleTmatrix. That would be too costly here,
  // though.
  // Also, we can't check here whether data is at the correct frequency since we
  // don't know f_grid and f_index here (we could pass it in, though).
  if( scat_data_mono[0][0].f_grid.nelem() > 1 )
  {
      ostringstream os;
      os << "Scattering data seems to be *scat_data* (several freq points),\n"
         << "but *scat_data_mono* (1 freq point only) is expected here.";
      throw runtime_error( os.str() );
  }

  GridPos T_gp, Tred_gp;
  Vector itw(2);

  // Initialisation
  pha_mat_spt = 0.;

  Index i_se_flat = 0;
  for (Index i_ss = 0; i_ss < scat_data_mono.nelem(); i_ss ++)
  {
      for (Index i_se = 0; i_se < scat_data_mono[i_ss].nelem(); i_se ++)
      {
          // If the particle number density at a specific point in the
          // atmosphere for scattering element i_se is zero, we don't need to
          // do the transformation!
          if (pnd_field(i_se_flat, scat_p_index, scat_lat_index, scat_lon_index)
              > PND_LIMIT)
          {
              // Temporary phase matrix which includes all temperatures.
              Index nT = scat_data_mono[i_ss][i_se].pha_mat_data.nvitrines();
              Tensor3 pha_mat_spt_tmp(nT,
                                      pha_mat_spt.nrows(), pha_mat_spt.ncols());

              pha_mat_spt_tmp = 0.;

              Index ti=-1;
              if( nT == 1 ) // just 1 T_grid element
              {
                  ti=0;
              }
              else if( rtp_temperature < 0. ) // coding for 'not interpolate, but
                                              // pick one temperature'
              {
                if( rtp_temperature > -10. )      // lowest T-point
                {
                  ti = 0;
                }
                else if( rtp_temperature > -20. ) // highest T-point
                {
                  ti = nT-1;
                }
                else                              // median T-point
                {
                  ti = nT/2;
                }
              }
              else
              {
                  ostringstream os;
                  os << "In pha_mat_sptFromMonoData.\n"
                     << "The temperature grid of the scattering data does not\n"
                     << "cover the atmospheric temperature at cloud location.\n"
                     << "The data should include the value T = "
                     << rtp_temperature << " K.";
                  chk_interpolation_grids( os.str(), 
                                           scat_data_mono[i_ss][i_se].T_grid, 
                                           rtp_temperature );

                  // Gridpositions:
                  gridpos( T_gp, scat_data_mono[i_ss][i_se].T_grid, 
                           rtp_temperature );
                  gridpos_copy( Tred_gp, T_gp );
                  Tred_gp.idx = 0;
                  // Interpolation weights:
                  interpweights(itw, Tred_gp);
              }

              // Do the transformation into the laboratory coordinate system.
              for (Index za_inc_idx = 0; za_inc_idx < doit_za_grid_size;
                   za_inc_idx ++)
              {
                  for (Index aa_inc_idx = 0; aa_inc_idx < scat_aa_grid.nelem();
                       aa_inc_idx ++)
                  {
                      if( ti<0 ) // Temperature interpolation
                      {
                          for (Index t_idx = 0; t_idx<2; t_idx ++)
                          {
                              pha_matTransform( pha_mat_spt_tmp(t_idx, joker, joker),
                                           scat_data_mono[i_ss][i_se].
                                           pha_mat_data
                                           (0,t_idx+T_gp.idx,joker,joker,joker,
                                            joker,joker),
                                           scat_data_mono[i_ss][i_se].za_grid,
                                           scat_data_mono[i_ss][i_se].aa_grid,
                                           scat_data_mono[i_ss][i_se].ptype,
                                           scat_za_index, scat_aa_index,
                                           za_inc_idx,
                                           aa_inc_idx, za_grid, scat_aa_grid,
                                           verbosity );
                          }

                          for (Index i = 0; i< stokes_dim; i++)
                          {
                              for (Index j = 0; j< stokes_dim; j++)
                              {
                                  pha_mat_spt(i_se_flat,
                                              za_inc_idx, aa_inc_idx, i, j)=
                                  interp(itw,
                                         pha_mat_spt_tmp(joker, i, j), Tred_gp);
                              }
                          }
                    }
                    else // no temperature interpolation required
                    {
                          pha_matTransform( pha_mat_spt(i_se_flat, za_inc_idx,
                                            aa_inc_idx, joker, joker),
                                           scat_data_mono[i_ss][i_se].
                                           pha_mat_data
                                           (0,ti,joker,joker,joker,
                                            joker,joker),
                                           scat_data_mono[i_ss][i_se].za_grid,
                                           scat_data_mono[i_ss][i_se].aa_grid,
                                           scat_data_mono[i_ss][i_se].ptype,
                                           scat_za_index, scat_aa_index,
                                           za_inc_idx,
                                           aa_inc_idx, za_grid, scat_aa_grid,
                                           verbosity );
                    }
                  }
              }
         }

          i_se_flat++;
      }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pha_mat_sptFromScat_data( // Output:
                         Tensor5& pha_mat_spt,
                         // Input:
                         const ArrayOfArrayOfSingleScatteringData& scat_data,
                         const Index& scat_data_checked,
                         const Vector& scat_za_grid,
                         const Vector& scat_aa_grid,
                         const Index& scat_za_index, // propagation directions
                         const Index& scat_aa_index,
                         const Index& f_index,
                         const Numeric& rtp_temperature,
                         const Tensor4& pnd_field, 
                         const Index& scat_p_index,
                         const Index& scat_lat_index,
                         const Index& scat_lon_index,
                         const Verbosity& verbosity
                         )
{
  if( scat_data_checked != 1 )
    throw runtime_error( "The scattering data must be flagged to have "
                         "passed a consistency check (scat_data_checked=1)." );

  const Index stokes_dim = pha_mat_spt.ncols();
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                        "must be 1,2,3 or 4");
  }

  // Determine total number of scattering elements
  const Index N_se_total = TotalNumberOfElements(scat_data);
  if( N_se_total != pnd_field.nbooks() )
    {
      ostringstream os;
      os << "Total number of scattering elements in scat_data "
         << "inconsistent with size of pnd_field.";
      throw runtime_error(os.str());
    }
  // as pha_mat_spt is typically initialized from pnd_field, this theoretically
  // checks the same as the runtime_error above. Still, we keep it to be on the
  // save side.
  assert( pha_mat_spt.nshelves() == N_se_total );
  
  const Index N_ss = scat_data.nelem();

  // Phase matrix in laboratory coordinate system. Dimensions:
  // [frequency, za_inc, aa_inc, stokes_dim, stokes_dim]
  Tensor5 pha_mat_data_int;

  Index this_f_index;

  Index i_se_flat = 0;
  // Loop over scattering species
  for (Index i_ss = 0; i_ss < N_ss; i_ss++)
  {
      const Index N_se = scat_data[i_ss].nelem();

      // Loop over the included scattering elements
      for (Index i_se = 0; i_se < N_se; i_se++)
      {
          // If the particle number density at a specific point in the
          // atmosphere for the i_se scattering element is zero, we don't need
          // to do the transfromation!
          if (abs(pnd_field(i_se_flat, scat_p_index,
                            scat_lat_index, scat_lon_index))
              > PND_LIMIT)
          {

              // First we have to transform the data from the coordinate system
              // used in the database (depending on the kind of ptype) to the
              // laboratory coordinate system.

              // Resize the variables for the interpolated data (1freq, 1T):
              pha_mat_data_int.resize(PHA_MAT_DATA.nshelves(),
                                      PHA_MAT_DATA.nbooks(),
                                      PHA_MAT_DATA.npages(),
                                      PHA_MAT_DATA.nrows(),
                                      PHA_MAT_DATA.ncols());

              // Frequency extraction and temperature interpolation

              // Gridpositions and interpolation weights;
              GridPos t_gp;
              Vector itw;
              Index this_T_index = -1;
              if ( PHA_MAT_DATA.nvitrines()==1 )
                {
                  this_T_index = 0;
                }
              else if( rtp_temperature < 0. ) // coding for 'not interpolate, but
                                              // pick one temperature'
                {
                  if( rtp_temperature > -10. )      // lowest T-point
                    {
                      this_T_index = 0;
                    }
                  else if( rtp_temperature > -20. ) // highest T-point
                    {
                      this_T_index = PHA_MAT_DATA.nvitrines()-1;
                    }
                  else                              // median T-point
                    {
                      this_T_index = PHA_MAT_DATA.nvitrines()/2;
                    }
                }
              else
                {
                  ostringstream os;
                  os << "In pha_mat_sptFromScat_data.\n"
                     << "The temperature grid of the scattering data does not\n"
                     << "cover the atmospheric temperature at cloud location.\n"
                     << "The data should include the value T = "
                     << rtp_temperature << " K.";
                  chk_interpolation_grids( os.str(), T_DATAGRID, 
                                           rtp_temperature );

                  gridpos(t_gp, T_DATAGRID, rtp_temperature);

                  // Interpolation weights:
                  itw.resize(2);
                  interpweights(itw, t_gp);
                }

              if( PHA_MAT_DATA.nlibraries()==1 )
                this_f_index = 0;
              else
                this_f_index = f_index;

              if ( this_T_index < 0 )
                {
                  // Interpolation of scattering matrix:
                  for (Index i_za_sca = 0;
                       i_za_sca < PHA_MAT_DATA.nshelves(); i_za_sca++)
                    for (Index i_aa_sca = 0;
                         i_aa_sca < PHA_MAT_DATA.nbooks(); i_aa_sca++)
                      for (Index i_za_inc = 0;
                           i_za_inc < PHA_MAT_DATA.npages(); i_za_inc++)
                        for (Index i_aa_inc = 0;
                             i_aa_inc < PHA_MAT_DATA.nrows(); i_aa_inc++)
                          for (Index i = 0; i < PHA_MAT_DATA.ncols(); i++)
                            pha_mat_data_int(i_za_sca, i_aa_sca,
                                             i_za_inc, i_aa_inc, i) =
                              interp(itw,
                                     PHA_MAT_DATA(this_f_index, joker,
                                                      i_za_sca, i_aa_sca,
                                                      i_za_inc, i_aa_inc, i),
                                     t_gp);
                }
              else
                {
                  pha_mat_data_int = PHA_MAT_DATA(this_f_index, this_T_index,
                                                      joker, joker,
                                                      joker, joker, joker);
                  /*
                  for (Index i_za_sca = 0;
                       i_za_sca < PHA_MAT_DATA.nshelves(); i_za_sca++)
                    for (Index i_aa_sca = 0;
                         i_aa_sca < PHA_MAT_DATA.nbooks(); i_aa_sca++)
                      for (Index i_za_inc = 0;
                           i_za_inc < PHA_MAT_DATA.npages(); i_za_inc++)
                        for (Index i_aa_inc = 0;
                             i_aa_inc < PHA_MAT_DATA.nrows(); i_aa_inc++)
                          for (Index i = 0; i < PHA_MAT_DATA.ncols(); i++)
                            // Interpolation of phase matrix:
                            pha_mat_data_int(i_za_sca, i_aa_sca,
                                             i_za_inc, i_aa_inc, i) =
                                PHA_MAT_DATA(this_f_index, this_T_index,
                                                 i_za_sca, i_aa_sca,
                  */
                }


              // Do the transformation into the laboratory coordinate system.
              for (Index za_inc_idx = 0; za_inc_idx < scat_za_grid.nelem();
                   za_inc_idx ++)
              {
                  for (Index aa_inc_idx = 0; aa_inc_idx < scat_aa_grid.nelem();
                       aa_inc_idx ++)
                  {
                      pha_matTransform(pha_mat_spt(i_se_flat,
                                                   za_inc_idx, aa_inc_idx,
                                                   joker, joker),
                                       pha_mat_data_int,
                                       ZA_DATAGRID, AA_DATAGRID,
                                       PART_TYPE, scat_za_index, scat_aa_index,
                                       za_inc_idx, 
                                       aa_inc_idx, scat_za_grid, scat_aa_grid,
                                       verbosity);
                  }
              }
          }
          i_se_flat++;
      }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesMerge(//WS Output:
                      Tensor4& pnd_field,
                      ArrayOfArrayOfSingleScatteringData& scat_data,
                      ArrayOfArrayOfScatteringMetaData& scat_meta,
                      ArrayOfString& scat_species,
                      Index& cloudbox_checked,
                      //WS Input:
                      const Index& atmosphere_dim,
                      const Index& cloudbox_on,
                      const ArrayOfIndex& cloudbox_limits,
                      const Tensor3& t_field,
                      const Tensor3& z_field,
                      const Matrix& z_surface,
                      const Verbosity& /*verbosity*/)
{
    // FIXME:
    // so far, this works for both scat_data and scat_data_raw. Needs to be
    // adapted, though, once we have WSM that can create Z/K/a with different
    // f/T dimensions than scat_data_single.f/T_grid.

    // cloudbox variable state should be ok before entering here
    if (!cloudbox_checked)
        throw std::runtime_error(
              "You must call *cloudbox_checkedCalc* before this method.");
    //however, we modify cloudbox variables. hence force re-checking the new
    //variables by resetting cloudbox_checked to False.
    cloudbox_checked = 0;
    
    if (atmosphere_dim != 1)
        throw std::runtime_error(
              "Merging scattering elements only works with a 1D atmoshere");

    // Cloudbox on/off?
    if ( !cloudbox_on )
    {
        /* Must initialise pnd_field anyway; but empty */
        pnd_field.resize(0, 0, 0, 0);
        return;
    }
    
    // ------- setup new pnd_field and scat_data -------------------
    ArrayOfIndex limits(2);
    //pressure
    limits[0] = cloudbox_limits[0];
    limits[1] = cloudbox_limits[1] + 1;

    Tensor4 pnd_field_merged(limits[1] - limits[0],
                          limits[1] - limits[0],
                          1,
                          1,
                          0.);

    ArrayOfArrayOfSingleScatteringData scat_data_merged;
    scat_data_merged.resize(1);
    scat_data_merged[0].resize(pnd_field_merged.nbooks());
    ArrayOfArrayOfScatteringMetaData scat_meta_merged;
    scat_meta_merged.resize(1);
    scat_meta_merged[0].resize(pnd_field_merged.nbooks());
    ArrayOfString scat_species_merged;
    scat_species_merged.resize(1);
    scat_species_merged[0] = "mergedfield-mergedpsd";
    for (Index sp = 0; sp < scat_data_merged[0].nelem(); sp++)
    {
        SingleScatteringData &this_part = scat_data_merged[0][sp];
        this_part.ptype = scat_data[0][0].ptype;
        this_part.description = "Merged scattering elements";
        this_part.f_grid = scat_data[0][0].f_grid;
        this_part.za_grid = scat_data[0][0].za_grid;
        this_part.aa_grid = scat_data[0][0].aa_grid;
        this_part.pha_mat_data.resize(scat_data[0][0].pha_mat_data.nlibraries(),
                                      1,
                                      scat_data[0][0].pha_mat_data.nshelves(),
                                      scat_data[0][0].pha_mat_data.nbooks(),
                                      scat_data[0][0].pha_mat_data.npages(),
                                      scat_data[0][0].pha_mat_data.nrows(),
                                      scat_data[0][0].pha_mat_data.ncols());
        this_part.ext_mat_data.resize(scat_data[0][0].ext_mat_data.nshelves(),
                                      1,
                                      scat_data[0][0].ext_mat_data.npages(),
                                      scat_data[0][0].ext_mat_data.nrows(),
                                      scat_data[0][0].ext_mat_data.ncols());
        this_part.abs_vec_data.resize(scat_data[0][0].abs_vec_data.nshelves(),
                                      1,
                                      scat_data[0][0].abs_vec_data.npages(),
                                      scat_data[0][0].abs_vec_data.nrows(),
                                      scat_data[0][0].abs_vec_data.ncols());
        this_part.pha_mat_data = 0.;
        this_part.ext_mat_data = 0.;
        this_part.abs_vec_data = 0.;
        this_part.T_grid.resize(1);
        this_part.T_grid[0] = t_field(sp, 0, 0);

        ScatteringMetaData &this_meta = scat_meta_merged[0][sp];
        ostringstream os;
        os << "Merged scattering element of cloudbox-level #" << sp;
        this_meta.description = os.str();
        this_meta.source = "ARTS internal";
        this_meta.refr_index = "Unknown";
        this_meta.mass = -1.;
        this_meta.diameter_max = -1.;
        this_meta.diameter_volume_equ = -1.;
        this_meta.diameter_area_equ_aerodynamical = -1.;
    }

    // Check that all scattering elements have same ptype and data dimensions
    SingleScatteringData &first_part = scat_data[0][0];
    for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++)
    {
        for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++)
        {
            SingleScatteringData &orig_part = scat_data[i_ss][i_se];

            if (orig_part.ptype != first_part.ptype)
                throw std::runtime_error(
                  "All scattering elements must have the same type");

            if (orig_part.f_grid.nelem() != first_part.f_grid.nelem())
                throw std::runtime_error(
                  "All scattering elements must have the same f_grid");

            if (!is_size(orig_part.pha_mat_data(joker, 0, joker, joker,
                                                joker, joker, joker),
                         first_part.pha_mat_data.nlibraries(),
                         first_part.pha_mat_data.nshelves(),
                         first_part.pha_mat_data.nbooks(),
                         first_part.pha_mat_data.npages(),
                         first_part.pha_mat_data.nrows(),
                         first_part.pha_mat_data.ncols()
                         ))
                throw std::runtime_error(
                  "All scattering elements must have the same pha_mat_data size"
                  " (except for temperature).");

            if (!is_size(orig_part.ext_mat_data(joker, 0, joker, joker, joker),
                         first_part.ext_mat_data.nshelves(),
                         first_part.ext_mat_data.npages(),
                         first_part.ext_mat_data.nrows(),
                         first_part.ext_mat_data.ncols()
                         ))
                throw std::runtime_error(
                  "All scattering elements must have the same ext_mat_data size"
                  " (except for temperature).");

            if (!is_size(orig_part.abs_vec_data(joker, 0, joker, joker, joker),
                         first_part.abs_vec_data.nshelves(),
                         first_part.abs_vec_data.npages(),
                         first_part.abs_vec_data.nrows(),
                         first_part.abs_vec_data.ncols()
                         ))
                throw std::runtime_error(
                  "All scattering elements must have the same abs_vec_data size"
                  " (except for temperature).");
        }
    }
    
    //----- Start pnd_field_merged and scat_data_array_merged calculations -----
    
    GridPos T_gp;
    Vector itw(2);
    
    Index nlevels = pnd_field_merged.nbooks();
    // loop over pressure levels in cloudbox
    for (Index i_lv = 0; i_lv < nlevels-1; i_lv++)
    {
        pnd_field_merged(i_lv,i_lv,0,0) = 1.;

        SingleScatteringData &this_part = scat_data_merged[0][i_lv];
        for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++)
        {
            for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++)
            {
                SingleScatteringData &orig_part = scat_data[i_ss][i_se];
                const Index pnd_index = FlattenedIndex(scat_data, i_ss, i_se);

                // If the particle number density at a specific point in the
                // atmosphere for the i_se scattering element is zero, we don't
                // need to do the transformation!
                if (pnd_field(pnd_index, i_lv, 0, 0) > PND_LIMIT) //TRS
                {
                    Numeric temperature = this_part.T_grid[0];
                    if( orig_part.T_grid.nelem() > 1)
                    {
                        ostringstream os;
                        os << "The temperature grid of the scattering data "
                           << "does not cover the\n"
                           << "atmospheric temperature at cloud location. "
                           << "The data should\n"
                           << "include the value T = "<< temperature << " K.\n"
                           << "Offending particle is scat_data[" << i_ss
                           << "][" << i_se << "]:\n"
                           << "Description: " << orig_part.description << "\n";
                        chk_interpolation_grids(os.str(), orig_part.T_grid,
                                                temperature);

                        // Gridpositions:
                        gridpos( T_gp, orig_part.T_grid, temperature );
                        // Interpolation weights:
                        interpweights(itw, T_gp);
                    }

                    ////////// Extinction matrix and absorption vector

                    // Loop over frequencies
                    for (Index i_f = 0;
                         i_f < orig_part.pha_mat_data.nlibraries(); i_f++)
                    {

                        // Loop over zenith angles
                        for (Index i_za = 0;
                             i_za < orig_part.ext_mat_data.npages(); i_za++)
                        {
                            // Loop over azimuth angles
                            for (Index i_aa = 0;
                                 i_aa < orig_part.ext_mat_data.nrows(); i_aa++)
                            {
                                // Weighted sum of ext_mat_data and abs_vec_data
                                if( orig_part.T_grid.nelem() == 1)
                                {
                                    Vector v = orig_part.ext_mat_data(i_f, 0, i_za, i_aa, joker);
                                    v *= pnd_field(pnd_index, i_lv, 0, 0);
                                    this_part.ext_mat_data(i_f, 0, i_za, 0, joker) += v;

                                    v = orig_part.abs_vec_data(i_f, 0, i_za, i_aa, joker);
                                    v *= pnd_field(pnd_index, i_lv, 0, 0);
                                    this_part.abs_vec_data(i_f, 0, i_za, i_aa, joker) += v;
                                }
                                else
                                {
                                    for (Index i = 0;
                                         i < orig_part.ext_mat_data.ncols(); i++)
                                    {
                                        // Temperature interpolation
                                        this_part.ext_mat_data(i_f, 0, i_za, i_aa, i) +=
                                        pnd_field(pnd_index, i_lv, 0, 0)
                                        * interp(itw,
                                                 orig_part.ext_mat_data(i_f, joker, i_za, i_aa, i),
                                                 T_gp);
                                    }
                                    for (Index i = 0;
                                         i < orig_part.abs_vec_data.ncols(); i++)
                                    {
                                        // Temperature interpolation
                                        this_part.abs_vec_data(i_f, 0, i_za, i_aa, i) +=
                                        pnd_field(pnd_index, i_lv, 0, 0)
                                        * interp(itw,
                                                 orig_part.abs_vec_data(i_f, joker, i_za, i_aa, i),
                                                 T_gp);
                                    }
                                }
                            }
                        }

                        ////////// Phase matrix

                        // Loop over outgoing zenith angles
                        for (Index i_za_out = 0;
                             i_za_out < orig_part.pha_mat_data.nshelves(); i_za_out++)
                        {
                            // Loop over outgoing azimuth angles
                            for (Index i_aa_out = 0;
                                 i_aa_out < orig_part.pha_mat_data.nbooks(); i_aa_out++)
                            {
                                // Loop over incoming zenith angles
                                for (Index i_za_inc = 0;
                                     i_za_inc < orig_part.pha_mat_data.npages(); i_za_inc++)
                                {
                                    // Loop over incoming azimuth angles
                                    for (Index i_aa_inc = 0;
                                         i_aa_inc < orig_part.pha_mat_data.nrows(); i_aa_inc++)
                                    {
                                        // Weighted sum of pha_mat_data
                                        if( orig_part.T_grid.nelem() == 1)
                                        {
                                            Vector v = orig_part.pha_mat_data(i_f, 0,
                                                                              i_za_out, i_aa_out,
                                                                              i_za_inc, i_aa_inc,
                                                                              joker);
                                            v *= pnd_field(pnd_index, i_lv, 0, 0);
                                            this_part.pha_mat_data(i_f, 0,
                                                                   i_za_out, i_aa_out,
                                                                   i_za_inc, i_aa_inc, joker) = v;
                                        }
                                        else
                                        {
                                            // Temperature interpolation
                                            for (Index i = 0;
                                                 i < orig_part.pha_mat_data.ncols(); i++)
                                            {
                                                this_part.pha_mat_data(i_f, 0, i_za_out,
                                                                       i_aa_out, i_za_inc, i_aa_inc, i) +=
                                                pnd_field(pnd_index, i_lv, 0, 0)
                                                * interp(itw,
                                                         orig_part.pha_mat_data(i_f, joker,
                                                                                i_za_out, i_aa_out,
                                                                                i_za_inc, i_aa_inc, i),
                                                         T_gp);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Set new pnd_field at lowest altitude to 0 if the cloudbox doesn't touch
    // the ground.
    // The consistency for the original pnd_field has already been ensured by
    // cloudbox_checkedCalc
    if (z_field(cloudbox_limits[0], 0, 0) > z_surface(0, 0))
        pnd_field_merged(0, 0, 0, 0) = 0.;

    pnd_field = pnd_field_merged;
    scat_data = scat_data_merged;
    scat_meta = scat_meta_merged;
    scat_species = scat_species_merged;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void ExtractFromMetaSingleScatSpecies(
                     //WS Output:
                     Vector& meta_param,
                     //WS Input:
                     const ArrayOfArrayOfScatteringMetaData& scat_meta,
                     const String& meta_name,
                     const Index& scat_species_index,
                     const Verbosity& /*verbosity*/)
{
    if ( scat_species_index<0 )
    {
      ostringstream os;
      os << "scat_species_index can't be <0!";
      throw runtime_error( os.str() );
    }

    const Index nss = scat_meta.nelem();

    // check that scat_meta actually has at least scat_species_index elements
    if ( !(nss>scat_species_index) )
    {
      ostringstream os;
      os << "Can not extract data for scattering species #"
         << scat_species_index << "\n"
         << "because scat_meta has only " << nss << " elements.";
      throw runtime_error( os.str() );
    }

    const Index nse = scat_meta[scat_species_index].nelem();
    meta_param.resize(nse);

    for ( Index i=0; i<nse; i++ )
      {
        if ( meta_name=="mass" )
          meta_param[i] =
            scat_meta[scat_species_index][i].mass;
        else if ( meta_name=="diameter_max" )
          meta_param[i] =
            scat_meta[scat_species_index][i].diameter_max;
        else if ( meta_name=="diameter_volume_equ" )
          meta_param[i] =
            scat_meta[scat_species_index][i].diameter_volume_equ;
        else if ( meta_name=="diameter_area_equ_aerodynamical" )
          meta_param[i] =
            scat_meta[scat_species_index][i].diameter_area_equ_aerodynamical;
        else
          {
            ostringstream os;
            os << "Meta parameter \"" << meta_name << "\"is unknown.";
            throw runtime_error( os.str() );
          }
      }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void TestScatDataInterp(
//                        Workspace& ws,
                        const ArrayOfArrayOfSingleScatteringData&   scat_data,
                        const Index&       stokes_dim,
                        const Vector&      f_grid,
                        const Vector&      rtp_los,
                        const Numeric&     rtp_temperature,
//                        const ArrayOfIndex& cloudbox_limits,
//                        const Agenda&      pha_mat_spt_agenda,
                        const Index&       scat_elem_index,
                        const ArrayOfIndex&  compare,
                        const Numeric&     tolerance,
                        const Index&       za_printinfo_index,
                        const Index&       aa_printinfo_index,
                        const Index&       mirror,
                        const Verbosity&   verbosity )
{

  cout << endl;
  cout << "========== START ==========" << endl;
  cout << endl;
  cout << "LOS direction: " << rtp_los << endl;
  cout << endl;

  if ( (rtp_los[0] > 180. ) || (rtp_los[0] < 0.) )
    throw runtime_error( "LOS zenith angle must be between 0 and 180deg." );
  if ( (rtp_los[1] > 180. ) || (rtp_los[1] < -180.) )
    throw runtime_error( "LOS azimuth angle must be between -180 and 180deg." );
  
  // Hard-coded grids for phase matrix incidence directions
  const Vector pha_mat_za( 0.0, 37, 5.0 );
  const Vector pha_mat_aa( -180.0, 37, 10.0 );

  // Some sizes
  Index N_se = TotalNumberOfElements(scat_data);
  const Index n_za        = pha_mat_za.nelem();
  const Index n_aa        = pha_mat_aa.nelem();
  const Index n_f = f_grid.nelem();

  // Input checks
  if( f_grid.nelem() != 1 )
    throw runtime_error( "Only length 1 *f_grid* allowed." );
  if( rtp_los.nelem() != 2 )
    throw runtime_error( "Only length 2 *rtp_los* allowed." );
  if( scat_elem_index < 0  ||  scat_elem_index >= N_se )
    throw runtime_error( "Invalid choice for *scat_element_index*" );

  Index printinfo=1;
  if ( za_printinfo_index<0 || aa_printinfo_index<0 )
    printinfo=0;
  else if ( za_printinfo_index>=n_za || aa_printinfo_index>=n_aa )
    {
      ostringstream os;
      os << "Printout indices too large (requirement: za<" << n_za << ", aa<"
         << n_aa << ").\n";
      throw runtime_error( os.str() );
    }

  // Common variables
  const Index f_index = 0;
  //
  Vector pnd_vec(N_se,0); pnd_vec[scat_elem_index] = 1; // 1-pt pnd as needed for MC

  Tensor4 pnd(N_se, 1, 1, 1);                           // 3D field (containing 1pt only) as used in DO methods)
  pnd(joker,0,0,0) = pnd_vec;


  
  ArrayOfArrayOfSingleScatteringData scat_data_mono;
  //
  scat_data_monoCalc( scat_data_mono, scat_data,
                      f_grid, f_index, verbosity );
  //

  ////// unified system //////
  //

  // making containers
  ArrayOfArrayOfTensor6 pha_mat_Nse;
  ArrayOfArrayOfTensor5 ext_mat_Nse;
  ArrayOfArrayOfTensor4 abs_vec_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;
  ArrayOfTensor6 pha_mat_ssbulk;
  ArrayOfTensor5 ext_mat_ssbulk;
  ArrayOfTensor4 abs_vec_ssbulk;
  ArrayOfIndex ptype_ssbulk;
  Tensor6 pha_mat_bulk;
  Tensor5 ext_mat_bulk;
  Tensor4 abs_vec_bulk;
  Index ptype_bulk;

  // preparing input in format needed
  Vector T_array(1,rtp_temperature);
  Matrix dir_array(1,2);
  dir_array(0,joker) = rtp_los;

  opt_prop_NScatElems( ext_mat_Nse, abs_vec_Nse, ptypes_Nse, t_ok,
                       scat_data_mono, stokes_dim, T_array, dir_array, f_index );
  opt_prop_ScatSpecBulk( ext_mat_ssbulk, abs_vec_ssbulk, ptype_ssbulk,
                         ext_mat_Nse, abs_vec_Nse, ptypes_Nse,
                         pnd(joker,joker,0,0), t_ok );
  opt_prop_Bulk( ext_mat_bulk, abs_vec_bulk, ptype_bulk,
                 ext_mat_ssbulk, abs_vec_ssbulk, ptype_ssbulk );

  Matrix idir_array(pha_mat_za.nelem()*pha_mat_aa.nelem(),2);
  ArrayOfArrayOfIndex printinfo_index;
  printinfo_index.resize(pha_mat_za.nelem());
  Index j=0;
  for( Index iz=0; iz<n_za; iz++ )
  {
    printinfo_index[iz].resize(pha_mat_za.nelem());
    for( Index ia=0; ia<n_aa; ia++ )
    {
      idir_array(j,0) = pha_mat_za[iz];
      idir_array(j,1) = pha_mat_aa[ia];
      printinfo_index[iz][ia] = j;
      j++;
    }
  }

  pha_mat_NScatElems( pha_mat_Nse, ptypes_Nse, t_ok,
                      scat_data_mono, stokes_dim, T_array, dir_array,
                      idir_array, f_index );
  pha_mat_ScatSpecBulk( pha_mat_ssbulk, ptype_ssbulk,
                        pha_mat_Nse, ptypes_Nse,
                        pnd(joker,joker,0,0), t_ok );
  pha_mat_Bulk( pha_mat_bulk, ptype_bulk,
                pha_mat_ssbulk, ptype_ssbulk );

  if (printinfo)
    {
      cout << "----- unified system -----" << endl;
      cout << "absorption vector:\n" << abs_vec_bulk(f_index,0,0,joker) << endl;
      cout << "extinction matrix:\n" << ext_mat_bulk(f_index,0,0,joker,joker) << endl;
      cout << "phase matrix (" << za_printinfo_index
           << "," << aa_printinfo_index << "):\n"
           << pha_mat_bulk(f_index,0,
                           0,printinfo_index[za_printinfo_index][aa_printinfo_index],
                           joker,joker)
           << endl;
      cout << endl;
    }



  ////// Monte Carlo //////
  //
  Matrix ext_mat_mc( stokes_dim, stokes_dim, 0.0 );
  Vector abs_vec_mc( stokes_dim, 0.0 );
  Tensor4 pha_mat_mc( n_za, n_aa, stokes_dim, stokes_dim, 0.0 );
  //
  Vector out;
  if (mirror)
    mirror_los( out, rtp_los, 3 );
  else
    out=rtp_los;

  //  
  opt_propCalc( ext_mat_mc, abs_vec_mc, out[0], out[1], scat_data_mono,
                stokes_dim, pnd_vec, rtp_temperature, verbosity );
  //
  for( Index iz=0; iz<n_za; iz++ )
    {
      Vector in_los(2);
      in_los[0] = pha_mat_za[iz];
      
      for( Index ia=0; ia<n_aa; ia++ )
        {
          Vector inc;
          in_los[1] = pha_mat_aa[ia];
          if (mirror)
            mirror_los( inc, in_los, 3 );
          else
            inc=in_los;
          pha_mat_singleCalc( pha_mat_mc(iz,ia,joker,joker),
                              out[0], out[1], inc[0], inc[1],
                              scat_data_mono, stokes_dim, pnd_vec,
                              rtp_temperature, verbosity );
        }
    }

  if (printinfo)
    {
      Vector in_los(2);
      in_los[0] = pha_mat_za[za_printinfo_index];
      in_los[1] = pha_mat_aa[aa_printinfo_index];
      cout << "Incident LOS direction: " << in_los << endl;
      cout << endl;

      Vector inc;
      if (mirror)
        mirror_los( inc, in_los, 3 );
      else
        inc=in_los;
      cout << "----- MC -----" << endl;
      cout << "photon propagation direction: " << out << endl;
      cout << "incident scattered photon propagation direction: " << inc << endl;
      cout << endl;

      cout << "absorption vector:\n" << abs_vec_mc << endl;
      cout << "extinction matrix:\n" << ext_mat_mc << endl;
      cout << "phase matrix (" << za_printinfo_index
           << "," << aa_printinfo_index << "):\n"
           << pha_mat_mc(za_printinfo_index,aa_printinfo_index,joker,joker)
           << endl;
      cout << endl;
    }


  ////// DO variable prep //////
  //

  // for abs_vec and ext_mat, the same implementation is used in RT4 and DOIT.
  ArrayOfStokesVector abs_vec_spt(N_se);
  for(auto& av : abs_vec_spt)
  {
    av = StokesVector(1, stokes_dim);
    av.SetZero();
  }
  
  ArrayOfPropagationMatrix ext_mat_spt(N_se);
  for(auto& pm : ext_mat_spt)
  {
    pm = PropagationMatrix(1, stokes_dim);
    pm.SetZero();
  }
  
/*
  ////// DOIT //////
  //

  Tensor5 pha_mat_spt( N_se, n_za, n_aa, stokes_dim, stokes_dim, 0.0 );
  
  StokesVector abs_vec_doit(n_f, stokes_dim);
  abs_vec_doit.SetZero();
  
  PropagationMatrix ext_mat_doit(n_f, stokes_dim);
  ext_mat_doit.SetZero();
  
  Tensor4 pha_mat_doit( n_za, n_aa, stokes_dim, stokes_dim, 0.0 );

  Tensor4 pnd3D(N_se, 1, 2, 2, 0.);
  pnd3D(scat_elem_index,0,joker,joker) = 1.;
  Tensor3 t3D(1, 2, 2, rtp_temperature);
  
  //doit_mono_agenda
  ArrayOfArrayOfSingleScatteringData scat_data_mono_optdoit;
  ArrayOfTensor7 pha_mat_sptDOITOpt;
  Tensor7 dummy;  // This is to catch output from DoitScatterinDataPrepare
                  // that is not needed here

  // JM180119:
  // tried to update this - remove the additional TestScatDataInterp interface
  // parameters introduced 170516 as they are inconsistent with the aim of this
  // WSM - to test the different opt prop extraction methods on identical input
  // I failed to go through with it - we'd also need some assumptions on
  // cloudbox_limits. But moreover, i couldn't see a way to get rid of
  // pha_mat_spt_agenda, or to replace that properly (before it was mimicked
  // below. but we can't pass that into DoitScatteringDataPrepare that way. That
  // is, we'd have to find another replacement. Or we indeed use
  // pha_mat_spt_agenda from the user. but than have to adapt the code below
  // accordingly (to also use this user-set agenda).
  DoitScatteringDataPrepare( ws, pha_mat_sptDOITOpt, scat_data_mono_optdoit,
                             dummy, n_za, pha_mat_aa,
                             scat_data, f_grid, f_index, 3, stokes_dim,
                             t3D, cloudbox_limits, pnd3D,
                             pha_mat_spt_agenda,
                             verbosity);

  //pha_mat_spt_agenda
  Numeric ANG_TOL=1e-3;
  Index scat_za_index=0;
  Index not_found=1;
  while( not_found && scat_za_index<n_za )
    {
      if( abs(pha_mat_za[scat_za_index]-rtp_los[0]) < ANG_TOL )
        not_found=0;
      else
        ++scat_za_index;
     }
  if( not_found )
    {
      ostringstream os;
      os << "For DOIT-like SSP extraction (using pha_mat_sptDOITOpt),\n"
         << "rtp_los[0] needs to be included in pha_mat_za\n"
         << "(within " << ANG_TOL << "deg tolerance), but is not.";
    throw runtime_error( os.str() );
    }
  Index scat_aa_index=0;
  not_found=1;
  while( not_found && scat_aa_index<n_aa )
    {
      if( abs(pha_mat_aa[scat_aa_index]-rtp_los[1]) < ANG_TOL )
        not_found=0;
      else
        ++scat_aa_index;
     }
  if( not_found )
    {
      ostringstream os;
      os << "For DOIT-like SSP extraction (using pha_mat_sptDOITOpt),\n"
         << "rtp_los[1] needs to be included in pha_mat_aa\n"
         << "(within " << ANG_TOL << "deg tolerance), but is not.";
    throw runtime_error( os.str() );
    }

  pha_mat_sptFromDataDOITOpt(pha_mat_spt,
                             pha_mat_sptDOITOpt, scat_data_mono_optdoit,
                             n_za, pha_mat_aa, scat_za_index, scat_aa_index,
                             rtp_temperature, pnd3D, 0, 0, 0, verbosity);

  //doit_scat_fieldCalc
  pha_matCalc(pha_mat_doit, pha_mat_spt, pnd, 3, 0, 0, 0, verbosity );

  //spt_calc_agenda
  opt_prop_sptFromMonoData( ext_mat_spt, abs_vec_spt,
                            scat_data_mono_optdoit, pha_mat_za, pha_mat_aa,
                            scat_za_index, scat_aa_index, rtp_temperature,
                            pnd, 0, 0, 0, verbosity );

  opt_prop_bulkCalc(ext_mat_doit, abs_vec_doit,
                    ext_mat_spt, abs_vec_spt,
                    pnd, 0, 0, 0, verbosity );

  if (printinfo)
    {
      Matrix tmp1(stokes_dim, stokes_dim);
      ext_mat_doit.MatrixAtPosition(tmp1, 0);
      Vector tmp2(stokes_dim);
      abs_vec_doit.VectorAtPosition(tmp2, 0);
      cout << "----- DOIT -----" << endl;
      cout << "absorption vector:\n" << tmp2 << endl;
      cout << "extinction matrix:\n" << tmp1 << endl;
      cout << "phase matrix (" << za_printinfo_index
           << "," << aa_printinfo_index << "):\n"
           << pha_mat_doit(za_printinfo_index,aa_printinfo_index,joker,joker)
           << endl;
      cout << endl;
    }
*/
    
  ////// RT4 //////
  //
  // FIXME: shouldn't we use (somehow; smartly setup) par_optpropCalc and
  // sca_optpropCalc?
  StokesVector abs_vec_rt4(n_f, stokes_dim);
  abs_vec_rt4.SetZero();
  
  PropagationMatrix ext_mat_rt4(n_f, stokes_dim);
  ext_mat_rt4.SetZero();
  
  Tensor4 pha_mat_rt4( n_za, n_aa, stokes_dim, stokes_dim, 0.0 );
  Vector sza_grid(1);
  sza_grid[0] = rtp_los[0];
  Vector saa_grid(1);
  saa_grid[0] = rtp_los[1];
  Vector siza_grid(n_za+1);
  Vector siaa_grid(n_aa+1);
  siza_grid[0] = sza_grid[0];
  siza_grid[Range(1,n_za)] = pha_mat_za;
  siaa_grid[0] = saa_grid[0];
  siaa_grid[Range(1,n_aa)] = pha_mat_aa;


  //spt_calc_agenda
  opt_prop_sptFromMonoData( ext_mat_spt, abs_vec_spt,
                            scat_data_mono, sza_grid, saa_grid,
                            0, 0, rtp_temperature,
                            pnd, 0, 0, 0, verbosity );

  opt_prop_bulkCalc(ext_mat_rt4, abs_vec_rt4, 
                    ext_mat_spt, abs_vec_spt,
                    pnd, 0, 0, 0, verbosity );

  // RT4-like Z-extraction (note: not exactly the same, as RT4 does not
  // require the explicit aa_grid value of Z(za_sca,aa_sca=delta_aa,za_inc), but
  // just the azimuthally averaged value of z(za_sca,za_inc), aka the 0th
  // component of the Fourier series of the azimuthal dependence of
  // Z(za_sca,za_inc).
  Index i_se_flat=0;
  for (Index i_ss = 0; i_ss < scat_data_mono.nelem(); i_ss++)
    for (Index i_se = 0; i_se < scat_data_mono[i_ss].nelem(); i_se++)
    {
      if (i_se_flat == scat_elem_index)
      {
        SingleScatteringData ssd=scat_data_mono[i_ss][i_se];
        Index i_pfct = ssd.T_grid.nelem()/2;

        if (ssd.ptype == PTYPE_TOTAL_RND)
          for (Index iza=0; iza<n_za; iza++)
            for (Index iaa=0; iaa<n_aa; iaa++)
              pha_matTransform( pha_mat_rt4(iza,iaa,joker,joker),
                                ssd.pha_mat_data(0,i_pfct,joker,
                                                 joker,joker,joker,joker),
                                ssd.za_grid, ssd.aa_grid, ssd.ptype,
                                0, 0, iza+1, iaa+1,
                                siza_grid, siaa_grid,
                                verbosity );

        else if (ssd.ptype == PTYPE_AZIMUTH_RND)
          {
            Index nza_se = ssd.za_grid.nelem();
            assert(nza_se==ssd.pha_mat_data.npages());
            ConstVectorView za_datagrid = ssd.za_grid;

            // in the actual SSP prep for RT4, the phase matrix is first
            // extracted and the azimuthal average (0th Fourier mode) is derived
            // at and based on the scatt elements' own polar angle grid.
            // That's not appropriate here for comparison to MC extraction,
            // hence here we interpolate to the requested azimuthal angles at
            // the scatt elements' own polar angles.
            // That is, we keep the original scat angle grids (inc and sca!)
            // meaning we loop over them and derive a pha_mat for each of them
            // (and each requested incident azimuth for fixed scattered azimuth)
            for (Index iaa=0; iaa<n_aa; iaa++)
              {
                Tensor4 pha_mat_interp(nza_se,nza_se,stokes_dim,stokes_dim,0.);

                GridPos daa_gp;
                Vector itw_aa(2);
                Numeric daa = abs(saa_grid[0]-pha_mat_aa[iaa]);
                if (daa>180.)
                  daa = 360.-daa;

                gridpos(daa_gp,ssd.aa_grid,daa);
                interpweights(itw_aa,daa_gp);

                for (Index iza=0; iza<nza_se; iza++)
                  for (Index sza=0; sza<nza_se; sza++)
                    for (Index ist1=0; ist1<stokes_dim; ist1++)
                      for (Index ist2=0; ist2<stokes_dim; ist2++)
                        pha_mat_interp(sza,iza,ist1,ist2) = interp(itw_aa,
                          ssd.pha_mat_data(0,i_pfct,sza,Range(joker),
                                           iza,0,ist1*4+ist2),
                          daa_gp);

            // in the actual SSP prep for RT4, the extracted azimuthal mode is
            // then interpolated to the RT4 solver polar angles.
            // Here, we're interpolating to the fixed scattered polar angle
            // and the hardcoded incidence angle grid. Need to do this for each
            // (delta) azimuth angle.
                for (Index iza=0; iza<n_za; iza++)
                {
                  GridPos za_sca_gp;
                  GridPos za_inc_gp;
                  Vector itw_za(4);
                  Numeric za_sca = sza_grid[0]; 
                  Numeric za_inc = pha_mat_za[iza]; 

                  gridpos(za_inc_gp,za_datagrid,za_inc);
                  gridpos(za_sca_gp,za_datagrid,za_sca);

                  interpweights(itw_za,za_sca_gp,za_inc_gp);
      
                  for (Index ist1=0; ist1<stokes_dim; ist1++)
                    for (Index ist2=0; ist2<stokes_dim; ist2++)
                      pha_mat_rt4(iza,iaa,ist1,ist2) = interp(itw_za,
                        pha_mat_interp(Range(joker),Range(joker),ist1,ist2),
                                       za_sca_gp,za_inc_gp);
                }
              }
          }

        else
          {
            ostringstream os;
            os << "Unsuitable particle type encountered.";
            throw runtime_error( os.str() );
          }
      }
      i_se_flat++;
    }

  if (printinfo)
    {
      Matrix tmp1(stokes_dim, stokes_dim);
      ext_mat_rt4.MatrixAtPosition(tmp1, 0);
      Vector tmp2(stokes_dim);
      abs_vec_rt4.VectorAtPosition(tmp2, 0);
      cout << "----- RT4 -----" << endl;
      cout << "absorption vector:\n" << tmp2 << endl;
      cout << "extinction matrix:\n" << tmp1 << endl;
      cout << "phase matrix (" << za_printinfo_index
           << "," << aa_printinfo_index << "):\n"
           << pha_mat_rt4(za_printinfo_index,aa_printinfo_index,joker,joker)
           << endl;
      cout << endl;
    }

  if (printinfo)
    {
      Vector dabs_vec(stokes_dim, -999.);
      Matrix dext_mat(stokes_dim,stokes_dim, -999.);
      Matrix dpha_mat(stokes_dim,stokes_dim, -999.);

      Matrix tmp1(stokes_dim, stokes_dim);
      Vector tmp2(stokes_dim);

      for (Index ist=0; ist<stokes_dim; ist++)
        dabs_vec[ist] = abs_vec_mc[ist]-abs_vec_bulk(f_index,0,0,ist);
      for (Index ist1=0; ist1<stokes_dim; ist1++)
        for (Index ist2=0; ist2<stokes_dim; ist2++)
          dext_mat(ist1,ist2) =
            ext_mat_mc(ist1,ist2)-ext_mat_bulk(f_index,0,0,ist1,ist2);

      for (Index ist1=0; ist1<stokes_dim; ist1++)
        for (Index ist2=0; ist2<stokes_dim; ist2++)
          dpha_mat(ist1,ist2) = 
              pha_mat_mc(za_printinfo_index,aa_printinfo_index,ist1,ist2)
            - pha_mat_bulk(f_index,0,
                           0,printinfo_index[za_printinfo_index][aa_printinfo_index],
                           ist1,ist2);

      cout << "----- differences (MC-US) -----" << endl;
      cout << "absorption vector:\n" << dabs_vec << endl;
      cout << "extinction matrix:\n" << dext_mat << endl;
      cout << "phase matrix (" << za_printinfo_index
           << "," << aa_printinfo_index << "):\n"
           << dpha_mat << endl;
      cout << endl;


      ext_mat_rt4.MatrixAtPosition(tmp1, 0);
      abs_vec_rt4.VectorAtPosition(tmp2, 0);
      for (Index ist=0; ist<stokes_dim; ist++)
        dabs_vec[ist] = tmp2[ist]-abs_vec_bulk(f_index,0,0,ist);
      for (Index ist1=0; ist1<stokes_dim; ist1++)
        for (Index ist2=0; ist2<stokes_dim; ist2++)
          dext_mat(ist1,ist2) = tmp1(ist1,ist2)-ext_mat_bulk(f_index,0,0,ist1,ist2);
      for (Index ist1=0; ist1<stokes_dim; ist1++)
        for (Index ist2=0; ist2<stokes_dim; ist2++)
          dpha_mat(ist1,ist2) = 
              pha_mat_rt4(za_printinfo_index,aa_printinfo_index,ist1,ist2)
            - pha_mat_bulk(f_index,0,
                           0,printinfo_index[za_printinfo_index][aa_printinfo_index],
                           ist1,ist2);

      cout << "----- differences (RT4-US) -----" << endl;
      cout << "absorption vector:\n" << dabs_vec << endl;
      cout << "extinction matrix:\n" << dext_mat << endl;
      cout << "phase matrix (" << za_printinfo_index
           << "," << aa_printinfo_index << "):\n"
           << dpha_mat << endl;
      cout << endl;

/*
      ext_mat_doit.MatrixAtPosition(tmp1, 0);
      abs_vec_doit.VectorAtPosition(tmp2, 0);
      for (Index ist=0; ist<stokes_dim; ist++)
        dabs_vec[ist] = tmp2[ist]-abs_vec_mc[ist];
      for (Index ist1=0; ist1<stokes_dim; ist1++)
        for (Index ist2=0; ist2<stokes_dim; ist2++)
          dext_mat(ist1,ist2) = tmp1(ist1,ist2)-ext_mat_mc(ist1,ist2);
      for (Index ist1=0; ist1<stokes_dim; ist1++)
        for (Index ist2=0; ist2<stokes_dim; ist2++)
          dpha_mat(ist1,ist2) = 
              pha_mat_doit(za_printinfo_index,aa_printinfo_index,ist1,ist2)
            - pha_mat_mc(za_printinfo_index,aa_printinfo_index,ist1,ist2);

      cout << "----- differences (DOIT-MC) -----" << endl;
      cout << "absorption vector:\n" << dabs_vec << endl;
      cout << "extinction matrix:\n" << dext_mat << endl;
      cout << "phase matrix (" << za_printinfo_index
           << "," << aa_printinfo_index << "):\n"
           << dpha_mat << endl;
      cout << endl;
*/
      
/*
      ext_mat_rt4.MatrixAtPosition(tmp1, 0);
      abs_vec_rt4.VectorAtPosition(tmp2, 0);
      for (Index ist=0; ist<stokes_dim; ist++)
        dabs_vec[ist] = tmp2[ist]-abs_vec_mc[ist];
      for (Index ist1=0; ist1<stokes_dim; ist1++)
        for (Index ist2=0; ist2<stokes_dim; ist2++)
          dext_mat(ist1,ist2) = tmp1(ist1,ist2)-ext_mat_mc(ist1,ist2);
      for (Index ist1=0; ist1<stokes_dim; ist1++)
        for (Index ist2=0; ist2<stokes_dim; ist2++)
          dpha_mat(ist1,ist2) = 
              pha_mat_rt4(za_printinfo_index,aa_printinfo_index,ist1,ist2)
            - pha_mat_mc(za_printinfo_index,aa_printinfo_index,ist1,ist2);

      cout << "----- differences (RT4-MC) -----" << endl;
      cout << "absorption vector:\n" << dabs_vec << endl;
      cout << "extinction matrix:\n" << dext_mat << endl;
      cout << "phase matrix (" << za_printinfo_index
           << "," << aa_printinfo_index << "):\n"
           << dpha_mat << endl;
      cout << endl;
*/
    }
  cout << "========== END ==========" << endl;
  cout << endl;


  ////// Compare //////
  //
  for( Index ic=0; ic<compare.nelem(); ic++ )
  {
    if (compare[ic]==1)
    {
      Numeric dmax;
      dmax = 0.5*tolerance*(abs_vec_bulk(f_index,0,0,0)+abs_vec_mc[0]);
      Compare(abs_vec_bulk(f_index,0,0,joker), abs_vec_mc, dmax, "Deviation in abs_vec",
              "US", "MC", "", "", verbosity);
      dmax = 0.5*tolerance*(ext_mat_bulk(f_index,0,0,0,0)+ext_mat_mc(0,0));
      Compare(ext_mat_bulk(f_index,0,0,joker,joker), ext_mat_mc, dmax, "Deviation in ext_mat",
              "US", "MC", "", "", verbosity);
      for (Index iza=0; iza<n_za; iza++)
        for (Index iaa=0; iaa<n_aa; iaa++)
        {
          ostringstream os;
          os << "Deviation in pha_mat at za[" << iza << "]=" << pha_mat_za[iza]
             << "deg and aa[" << iaa << "]=" << pha_mat_aa[iaa] << "deg.";
          dmax = 0.5*tolerance*(pha_mat_bulk(f_index,0,0,printinfo_index[iza][iaa],0,0)+
                         pha_mat_mc(iza,iaa,0,0));
          Compare(pha_mat_bulk(f_index,0,0,printinfo_index[iza][iaa],joker,joker),
                  pha_mat_mc(iza,iaa,joker,joker),
                  dmax, os.str(), "US", "MC", "", "", verbosity);
        }
    }

    if (compare[ic]==2)
    {
      Matrix tmp1(stokes_dim, stokes_dim);
      ext_mat_rt4.MatrixAtPosition(tmp1, 0);
      Vector tmp2(stokes_dim);
      abs_vec_rt4.VectorAtPosition(tmp2, 0);

      Numeric dmax;
      dmax = 0.5*tolerance*(abs_vec_bulk(f_index,0,0,0)+tmp2[0]);
      Compare(abs_vec_bulk(f_index,0,0,joker), tmp2, dmax, "Deviation in abs_vec",
              "US", "RT4", "", "", verbosity);
      dmax = 0.5*tolerance*(ext_mat_bulk(f_index,0,0,0,0)+tmp1(0,0));
      Compare(ext_mat_bulk(f_index,0,0,joker,joker), tmp1, dmax, "Deviation in ext_mat",
              "US", "RT4", "", "", verbosity);
      for (Index iza=0; iza<n_za; iza++)
        for (Index iaa=0; iaa<n_aa; iaa++)
        {
          ostringstream os;
          os << "Deviation in pha_mat at za[" << iza << "]=" << pha_mat_za[iza]
             << "deg and aa[" << iaa << "]=" << pha_mat_aa[iaa] << "deg.";
          dmax = 0.5*tolerance*(pha_mat_bulk(f_index,0,0,printinfo_index[iza][iaa],0,0)+
                         pha_mat_rt4(iza,iaa,0,0));
          Compare(pha_mat_bulk(f_index,0,0,printinfo_index[iza][iaa],joker,joker),
                  pha_mat_rt4(iza,iaa,joker,joker),
                  dmax, os.str(), "US", "RT4", "", "", verbosity);
        }
    }

    if (compare[ic]==3)
    {
      Matrix tmp1(stokes_dim, stokes_dim);
      ext_mat_rt4.MatrixAtPosition(tmp1, 0);
      Vector tmp2(stokes_dim);
      abs_vec_rt4.VectorAtPosition(tmp2, 0);
    
      Numeric dmax;
      dmax = 0.5*tolerance*(tmp2[0]+abs_vec_mc[0]);
      Compare(tmp2, abs_vec_mc, dmax, "Deviation in abs_vec",
              "RT4", "MC", "", "", verbosity);
      dmax = 0.5*tolerance*(tmp1(0,0)+ext_mat_mc(0,0));
      Compare(tmp1, ext_mat_mc, dmax, "Deviation in ext_mat",
              "RT4", "MC", "", "", verbosity);
      for (Index iza=0; iza<n_za; iza++)
        for (Index iaa=0; iaa<n_aa; iaa++)
        {
          ostringstream os;
          os << "Deviation in pha_mat at za[" << iza << "]=" << pha_mat_za[iza]
             << "deg and aa[" << iaa << "]=" << pha_mat_aa[iaa] << "deg.";
          dmax = 0.5*tolerance*(pha_mat_rt4(iza,iaa,0,0)+pha_mat_mc(iza,iaa,0,0));
          Compare(pha_mat_rt4(iza,iaa,joker,joker), pha_mat_mc(iza,iaa,joker,joker),
                  dmax, os.str(), "RT4", "MC", "", "", verbosity);
        }
    }

  /*
  if (compare[ic]==...)
    {
      Matrix tmp1(stokes_dim, stokes_dim, 0.0);
      ext_mat_doit.MatrixAtPosition(tmp1, 0);
      Vector tmp2(stokes_dim);
      abs_vec_doit.VectorAtPosition(tmp2, 0);

      Numeric dmax;
      dmax = 0.5*tolerance*(tmp2[0]+abs_vec_mc[0]);
      Compare(tmp2, abs_vec_mc, dmax, "Deviation in abs_vec",
              "DOIT", "MC", "", "", verbosity);
      dmax = 0.5*tolerance*(tmp1(0,0)+ext_mat_mc(0,0));
      Compare(tmp1, ext_mat_mc, dmax, "Deviation in ext_mat",
              "DOIT", "MC", "", "", verbosity);
      for (Index iza=0; iza<n_za; iza++)
        for (Index iaa=0; iaa<n_aa; iaa++)
        {
          ostringstream os;
          os << "Deviation in pha_mat at za[" << iza << "]=" << pha_mat_za[iza]
             << "deg and aa[" << iaa << "]=" << pha_mat_aa[iaa] << "deg.";
          dmax = 0.5*tolerance*(pha_mat_doit(iza,iaa,0,0)+pha_mat_mc(iza,iaa,0,0));
          Compare(pha_mat_doit(iza,iaa,joker,joker), pha_mat_mc(iza,iaa,joker,joker),
                  dmax, os.str(), "DOIT", "MC", "", "", verbosity);
        }
    }
  */
  }
}

