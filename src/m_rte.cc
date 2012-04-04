/* Copyright (C) 2002-2010
   Patrick Eriksson <patrick.eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
                            
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */



/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_rte.cc
  \author Patrick Eriksson <patrick.eriksson@chalmers.se>
  \date   2002-05-11 

  \brief  Workspace functions for solving clear sky radiative transfer.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "check_input.h"
#include "geodetic.h"
#include "jacobian.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "montecarlo.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"

extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT;
extern const String ABSSPECIES_MAINTAG;
extern const String TEMPERATURE_MAINTAG;





/*===========================================================================
  === Help sub-functions to handle analytical jacobians (in alphabetical order)
  ===========================================================================*/


//! diy_from_path_to_rgrids
/*!
    Maps jacobian data for points along the propagation path, to
    jacobian retrieval grid data.

    \param   diy_dx              Out: Jacobians for selected retrieval grids.
    \param   jacobian_quantity   As the WSV.
    \param   diy_dpath           Jacobians along the propagation path.
    \param   atmosphere_dim      As the WSV.
    \param   ppath               As the WSV.
    \param   ppath_p             The pressure at each ppath point.

    \author Patrick Eriksson 
    \date   2009-10-08
*/
// A small help function, to make the code below cleaner
void from_dpath_to_dx(
         MatrixView   diy_dx,
    ConstMatrixView   diy_dq,
    const Numeric&    w )
{
  for( Index irow=0; irow<diy_dx.nrows(); irow++ )
    { 
      for( Index icol=0; icol<diy_dx.ncols(); icol++ )
        { diy_dx(irow,icol) += w * diy_dq(irow,icol); }
    }
}
// Kept as a local function as long as just used here.
// We trust here gridpos, and "extrapolation points" are identified simply
// by a fd outside [0,1].
void add_extrap( ArrayOfGridPos&   gp )
{
  for( Index i=0; i<gp.nelem(); i++ )
    { 
      if( gp[i].fd[0] < 0 ) 
        {  
          gp[i].fd[0] = 0; 
          gp[i].fd[1] = 1; 
        }
      else if( gp[i].fd[0] > 1 ) 
        {  
          gp[i].fd[0] = 1; 
          gp[i].fd[1] = 0; 
        }
    }
}
//
void diy_from_path_to_rgrids(
         Tensor3View          diy_dx,
   const RetrievalQuantity&   jacobian_quantity,
   ConstTensor3View           diy_dpath,
   const Index&               atmosphere_dim,
   const Ppath&               ppath,
   ConstVectorView            ppath_p )
{
  // We want here an extrapolation to infinity -> 
  //                                        extremly high extrapolation factor
  const Numeric   extpolfac = 1.0e99;

  // Retrieval grid of interest
  Vector r_grid;

  if( ppath.np > 1 )  // Otherwise nothing to do here
    {
      // Pressure
      r_grid = jacobian_quantity.Grids()[0];
      Index            nr1 = r_grid.nelem();
      ArrayOfGridPos   gp_p(ppath.np);
      p2gridpos( gp_p, r_grid, ppath_p, extpolfac );
      add_extrap( gp_p );

      // Latitude
      Index            nr2 = 1;
      ArrayOfGridPos   gp_lat;
      if( atmosphere_dim > 1 )
        {
          gp_lat.resize(ppath.np);
          r_grid = jacobian_quantity.Grids()[1];
          nr2    = r_grid.nelem();
          gridpos( gp_lat, r_grid, ppath.pos(joker,1), extpolfac );
          add_extrap( gp_lat );
        }

      // Longitude
      ArrayOfGridPos   gp_lon;
      if( atmosphere_dim > 2 )
        {
          gp_lon.resize(ppath.np);
          r_grid = jacobian_quantity.Grids()[2];
          gridpos( gp_lon, r_grid, ppath.pos(joker,2), extpolfac );
          add_extrap( gp_lon );
        }
      
      //- 1D
      if( atmosphere_dim == 1 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              if( gp_p[ip].fd[0] < 1 )
                {
                  from_dpath_to_dx( diy_dx(gp_p[ip].idx,joker,joker),
                                    diy_dpath(ip,joker,joker), gp_p[ip].fd[1] );
                }
              if( gp_p[ip].fd[0] > 0 )
                {
                  from_dpath_to_dx( diy_dx(gp_p[ip].idx+1,joker,joker),
                                    diy_dpath(ip,joker,joker), gp_p[ip].fd[0] );
                }
            }
        }

      //- 2D
      else if( atmosphere_dim == 2 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              Index   ix = nr1*gp_lat[ip].idx + gp_p[ip].idx;
              // Low lat, low p
              from_dpath_to_dx( diy_dx(ix,joker,joker),
                                diy_dpath(ip,joker,joker), 
                                gp_lat[ip].fd[1]*gp_p[ip].fd[1] );
              // Low lat, high p
              from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                                gp_lat[ip].fd[1]*gp_p[ip].fd[0] );
              // High lat, low p
              from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                                gp_lat[ip].fd[0]*gp_p[ip].fd[1] );
              // High lat, high p
              from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                                gp_lat[ip].fd[0]*gp_p[ip].fd[0] );
            }
        }

      //- 3D
      else if( atmosphere_dim == 3 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              Index   ix = nr2*nr1*gp_lon[ip].idx +
                           nr1*gp_lat[ip].idx + gp_p[ip].idx;
              // Low lon, low lat, low p
              from_dpath_to_dx( diy_dx(ix,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[1]*gp_p[ip].fd[1]);
              // Low lon, low lat, high p
              from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[1]*gp_p[ip].fd[0]);
              // Low lon, high lat, low p
              from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[0]*gp_p[ip].fd[1]);
              // Low lon, high lat, high p
              from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[0]*gp_p[ip].fd[0]);

              // Increase *ix* (to be valid for high lon level)
              ix += nr2*nr1;

              // High lon, low lat, low p
              from_dpath_to_dx( diy_dx(ix,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[1]*gp_p[ip].fd[1]);
              // High lon, low lat, high p
              from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[1]*gp_p[ip].fd[0]);
              // High lon, high lat, low p
              from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[0]*gp_p[ip].fd[1]);
              // High lon, high lat, high p
              from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[0]*gp_p[ip].fd[0]);
            }
        }
    }
}



//! Help function for analytical jacobian calculations
/*!
    The function determines which terms in jacobian_quantities that are 
    analytical absorption species and temperature jacobians. 

    *abs_species_i* and *is_t* shall be sized to have the same length
    as *jacobian_quantities*. For analytical absorption species
    jacobians, *abs_species_i* is set to the matching index in
    *abs_species*. Otherwise, to -1. For analytical temperature
    jacobians, *is_t* is set to 1. Otherwise to 0.

    \param   abs_species_i         Out: Matching index in abs_species 
    \param   is_t                  Out: Flag for: Is a temperature jacobian?
    \param   jacobian_quantities   As the WSV.
    \param   abs_species           As the WSV.


    \author Patrick Eriksson 
    \date   2009-10-07
*/
void get_pointers_for_analytical_jacobians( 
        ArrayOfIndex&               abs_species_i, 
        ArrayOfIndex&               is_t,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfSpeciesTag&   abs_species )
{

  FOR_ANALYTICAL_JACOBIANS_DO( 
    //
    if( jacobian_quantities[iq].MainTag() == TEMPERATURE_MAINTAG )
      { is_t[iq] = 1; }
    else
      { is_t[iq] = 0; }
    //
    if( jacobian_quantities[iq].MainTag() == ABSSPECIES_MAINTAG )
      {
        ArrayOfSpeciesTag  atag;
        array_species_tag_from_string( atag, jacobian_quantities[iq].Subtag() );
        abs_species_i[iq] = chk_contains( "abs_species", abs_species, atag );
      }
    else
      { abs_species_i[iq] = -1; }
  )
}






/*===========================================================================
  === Workspace methods 
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void iyBeerLambertStandardClearsky(
        Workspace&            ws,
        Matrix&               iy,
        Matrix&               iy_error,
        Index&                iy_error_type,
        Matrix&               iy_aux,
        ArrayOfTensor3&       diy_dx,
  const Index&                iy_agenda_call1,
  const Tensor3&              iy_transmission,
  const Vector&               rte_pos,      
  const Vector&               rte_los,      
  const Index&                jacobian_do,
  const Index&                atmosphere_dim,
  const Vector&               p_grid,
  const Vector&               lat_grid,
  const Vector&               lon_grid,
  const Tensor3&              z_field,
  const Tensor3&              t_field,
  const Tensor4&              vmr_field,
  const Tensor3&              wind_u_field,
  const Tensor3&              wind_v_field,
  const Tensor3&              wind_w_field,
  const Tensor3&              edensity_field,
  const Vector&               refellipsoid,
  const Matrix&               z_surface,
  const Index&                cloudbox_on,
  const Index&                stokes_dim,
  const Vector&               f_grid,
  const ArrayOfArrayOfSpeciesTag&   abs_species,
  const Index&                mblock_index,
  const Agenda&               ppath_agenda,
  const Agenda&               abs_scalar_gas_agenda,
  const Agenda&               iy_clearsky_agenda,
  const Agenda&               iy_space_agenda,
  const Agenda&               surface_prop_agenda,
  const Agenda&               iy_cloudbox_agenda,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&  jacobian_indices,
  const Verbosity&            verbosity )
{
  // See initial comments of iyEmissionStandardClearsky

  // Sizes
  //
  const Index nf = f_grid.nelem();
  const Index nq = jacobian_quantities.nelem();

  // Determine if there are any analytical jacobians to handle, and if primary
  // call, resize diy_dx.
  //
  Index j_analytical_do = 0;
  //
  if( jacobian_do ) { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }
  //
  if( iy_agenda_call1 && j_analytical_do )
    {
      diy_dx.resize( jacobian_indices.nelem() ); 
      //
      FOR_ANALYTICAL_JACOBIANS_DO( 
        diy_dx[iq].resize( jacobian_indices[iq][1] - jacobian_indices[iq][0] + 
                           1, nf, stokes_dim ); 
        diy_dx[iq] = 0.0;
      ) 
    }

  //- Determine propagation path
  //
  Ppath  ppath;
  //
  ppath_agendaExecute( ws, ppath, rte_pos, rte_los, cloudbox_on, 0,
                       mblock_index, t_field, z_field, vmr_field, 
                       edensity_field, -1, ppath_agenda );

  // Get atmospheric and RT quantities for each ppath point/step
  //
  // If np = 1, we only need to determine the radiative background
  //
  // "atmvars"
  Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w;
  Matrix    ppath_vmr;
  // "rtvars"
  Vector    total_tau;
  Matrix    emission_dummy, ppath_tau;
  Tensor3   ppath_abs_scalar, iy_trans_new;
  Agenda    agenda_dummy;
  //
  const Index np  = ppath.np;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      get_ppath_atmvars( ppath_p, ppath_t, ppath_vmr, 
                         ppath_wind_u, ppath_wind_v, ppath_wind_w,
                         ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                         wind_u_field, wind_v_field, wind_w_field );

      // Absorption and optical thickness for each step
      get_ppath_rtvars( ws, ppath_abs_scalar, ppath_tau, total_tau, 
                   emission_dummy,  abs_scalar_gas_agenda, agenda_dummy,
                   ppath, ppath_p, ppath_t, ppath_vmr, ppath_wind_u, 
                   ppath_wind_v, ppath_wind_w, f_grid, -1, atmosphere_dim, 0 );
    }
  else // To handle cases inside the cloudbox, or totally outside the atmosphere
    {
      total_tau.resize( nf );
      total_tau = 0;
    }

  // iy_transmission
  //
  if( iy_agenda_call1 )
    { iy_transmission_for_scalar_tau( iy_trans_new, stokes_dim, total_tau ); }
  else
    { iy_transmission_mult_scalar_tau( iy_trans_new, iy_transmission, 
                                                                total_tau ); }

  // Radiative background
  //
  get_iy_of_background( ws, iy, iy_error, iy_error_type, iy_aux, diy_dx, 
                        iy_trans_new, jacobian_do, ppath, atmosphere_dim, 
                        lat_grid, lon_grid, t_field, z_field, 
                        vmr_field, refellipsoid, z_surface, cloudbox_on, 
                        stokes_dim, f_grid, iy_clearsky_agenda, iy_space_agenda,
                        surface_prop_agenda, iy_cloudbox_agenda, verbosity );
  
  // Do RT calculations
  //
  if( np > 1 )
    {
      // Create container for the derivatives with respect to changes
      // at each ppath point. And find matching abs_species-index and 
      // "temperature flag" (for analytical jacobians).
      //
      ArrayOfTensor3  diy_dpath; 
      ArrayOfIndex    abs_species_i, is_t; 
      //
      const Numeric   dt = 0.1;
            Tensor3   ppath_as2;
      //
      if( j_analytical_do )
        { 
          diy_dpath.resize( nq ); 
          abs_species_i.resize( nq ); 
          is_t.resize( nq ); 
          //
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_dpath[iq].resize( np, nf, stokes_dim ); 
            diy_dpath[iq] = 0.0;
          )
          get_pointers_for_analytical_jacobians( abs_species_i, is_t, 
                                            jacobian_quantities, abs_species );
          //
          // Determine if temperature is among the analytical jac. quantities.
          // If yes, get emission and absorption for disturbed temperature
          Index do_t=0;
          for( Index iq=0; iq<is_t.nelem() && do_t==0; iq++ )
            { if( is_t[iq] ) { do_t = 1; } }
          if( do_t )
            {
              Matrix tau_dummy; 
              Vector total_tau_dummy, t2 = ppath_t;
              t2 += dt;
              get_ppath_rtvars( ws, ppath_as2, tau_dummy, total_tau_dummy,
                                emission_dummy, abs_scalar_gas_agenda, 
                                agenda_dummy, ppath, ppath_p, t2, ppath_vmr, 
                                ppath_wind_u, ppath_wind_v, ppath_wind_w, 
                                f_grid, -1, atmosphere_dim, 0 );
            }
        }

      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {

          // Jacobian quantities
          if( j_analytical_do )
            {
              // Common terms introduced for efficiency and clarity
              Vector X(nf);   // See AUG
              //
              for( Index iv=0; iv<nf; iv++ )
                {
                  X[iv] = 0.5 * ppath.lstep[ip] * exp( -total_tau[iv] );
                }

              // Loop quantities
              // Loop quantities
              for( Index iq=0; iq<nq; iq++ ) 
                {
                  if( jacobian_quantities[iq].Analytical() )
                    {
                      // Absorbing species
                      const Index isp = abs_species_i[iq]; 
                      if( isp >= 0 )
                        {
                          // Scaling factors to handle retrieval unit
                          // (gives the cross-section to apply)
                          Numeric unitscf1, unitscf2;
                          vmrunitscf( unitscf1, 
                                      jacobian_quantities[iq].Mode(), 
                                      ppath_vmr(isp,ip), ppath_p[ip], 
                                      ppath_t[ip] );
                          vmrunitscf( unitscf2, 
                                      jacobian_quantities[iq].Mode(), 
                                      ppath_vmr(isp,ip+1), ppath_p[ip+1], 
                                      ppath_t[ip+1] );
                          //
                          for( Index iv=0; iv<nf; iv++ )
                            {
                              // All Stokes components equal
                              for( Index is=0; is<stokes_dim; is++ )
                                { 
                                  const Numeric Z = -X[iv] * iy(iv,is);
                                  diy_dpath[iq](ip  ,iv,is) += Z * unitscf1 *
                                                 ppath_abs_scalar(iv,isp,ip);
                                  diy_dpath[iq](ip+1,iv,is) += Z * unitscf2 *
                                                 ppath_abs_scalar(iv,isp,ip+1);
                                }
                            }
                        }

                      // Temperature
                      else if( is_t[iq] )
                        {
                          for( Index iv=0; iv<nf; iv++ )
                            {
                              // The terms associated with Dtau/Dk:
                              const Numeric k1 = 
                                         ppath_abs_scalar(iv,joker,ip  ).sum();
                              const Numeric k2 = 
                                         ppath_abs_scalar(iv,joker,ip+1).sum();
                              const Numeric dkdt1 =
                                  ( ppath_as2(iv,joker,ip  ).sum() - k1 ) / dt;
                              const Numeric dkdt2 =
                                  ( ppath_as2(iv,joker,ip+1).sum() - k2 ) / dt;
                              for( Index is=0; is<stokes_dim; is++ )
                                { 
                                  const Numeric Z = -X[iv] * iy(iv,is);
                                  diy_dpath[iq](ip  ,iv,is) += Z * dkdt1;
                                  diy_dpath[iq](ip+1,iv,is) += Z * dkdt2;
                                }
                              //
                              // The terms associated with Delta-s:
                              if( jacobian_quantities[iq].Subtag() == "HSE on" )
                                {
                                  const Numeric kbar = 0.5 * ( k1 + k2 );
                                  for( Index is=0; is<stokes_dim; is++ )
                                    { 
                                      const Numeric Z = -X[iv] * iy(iv,is);
                                      diy_dpath[iq](ip  ,iv,is) += Z * kbar /
                                                                  ppath_t[ip];
                                      diy_dpath[iq](ip+1,iv,is) += Z * kbar /
                                                                 ppath_t[ip+1];
                                    }
                                } //hse
                            }  // frequency
                        }  // if is_t
                    } // if analytical
                } // for iq
              
              // Update total_tau (not used for spectrum calculation)
              for( Index iv=0; iv<nf; iv++ )
                { total_tau[iv] -= ppath_tau(iv,ip); }
            }

          // Spectrum
          //
          for( Index iv=0; iv<nf; iv++ )
            {
              const Numeric step_tr = exp( -ppath_tau(iv,ip) );
              //
              for( Index is=0; is<stokes_dim; is++ )
                { iy(iv,is) *= step_tr; }
            }
        } 

      // Map jacobians from ppath to retrieval grids
      // (this operation corresponds to the term Dx_i/Dx)
      if( j_analytical_do )
        { 
          // Weight with iy_transmission
          if( !iy_agenda_call1 )
            {
              Matrix X, Y(stokes_dim,diy_dpath[0].npages()); 
              //
              FOR_ANALYTICAL_JACOBIANS_DO( 
                for( Index iv=0; iv<nf; iv++ )
                  { 
                    X = transpose( diy_dpath[iq](joker,iv,joker) );
                    mult( Y, iy_transmission(iv,joker,joker), X );
                    diy_dpath[iq](joker,iv,joker) = transpose( Y );
                  }
               )
            }

          // Map to retrieval grids
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_from_path_to_rgrids( diy_dx[iq], jacobian_quantities[iq], 
                               diy_dpath[iq], atmosphere_dim, ppath, ppath_p );
          )
        }
    }
}






/* Workspace method: Doxygen documentation will be auto-generated */
void iyBeerLambertStandardCloudbox(
        Workspace&                     ws,
        Matrix&                        iy,
        Matrix&                        iy_error,
        Index&                         iy_error_type,
        Matrix&                        iy_aux,
        ArrayOfTensor3&                diy_dx,
  const Tensor3&                       iy_transmission,
  const Vector&                        rte_pos,
  const Vector&                        rte_los,
  const Index&                         jacobian_do,
  const Index&                         atmosphere_dim,
  const Vector&                        p_grid,
  const Tensor3&                       z_field,
  const Tensor3&                       t_field,
  const Tensor4&                       vmr_field,
  const Tensor3&                       edensity_field,
  const Index&                         cloudbox_on,
  const ArrayOfIndex&                  cloudbox_limits,
  const Index&                         stokes_dim,
  const Vector&                        f_grid,
  const Agenda&                        ppath_agenda,
  const Agenda&                        abs_scalar_gas_agenda,
  const Agenda&                        iy_clearsky_agenda,
  const Tensor4&                       pnd_field,
  const Index&                         use_mean_scat_data,
  const ArrayOfSingleScatteringData&   scat_data_raw,
  const Agenda&                        opt_prop_gas_agenda,
  const Verbosity&                     verbosity)
{
  // Input checks
  if( !cloudbox_on )
    throw runtime_error( "The cloudbox must be defined to use this method." );
  if( jacobian_do )
    throw runtime_error( 
       "This method does not provide any jacobians (jacobian_do must be 0)." );

  // Sizes
  //
  const Index nf = f_grid.nelem();

  // Determine ppath through the cloudbox
  //
  Ppath  ppath;
  //
  ppath_agendaExecute( ws, ppath, rte_pos, rte_los, cloudbox_on, 1, -1,
                       t_field, z_field, vmr_field, edensity_field, -1,
                       ppath_agenda );

  // Check radiative background
  const Index bkgr = ppath_what_background( ppath );
  if( bkgr == 2 )
    throw runtime_error( "Observations where (unscattered) propagation path "
                         "hits the surface inside the cloudbox are not yet "
                         "handled by this method." );
  assert( bkgr == 3 );

  // Get atmospheric and RT quantities for each ppath point/step (inside box)
  // 
  // If np = 1, we only need to determine the radiative background
  //
  // "atmvars"
  Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w;
  Matrix    ppath_vmr, ppath_pnd;
  // "rtvars"
  Matrix    emission_dummy, ppath_tau;
  Tensor3   wind_field_dummy(0,0,0), iy_trans_new;
  Tensor3   ppath_asp_abs_vec, ppath_pnd_abs_vec, total_transmission;
  Tensor4   ppath_asp_ext_mat, ppath_pnd_ext_mat, ppath_transmission;
  Agenda    agenda_dummy;
  Array<ArrayOfSingleScatteringData>  scat_data;
  //
  const Index np  = ppath.np;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      get_ppath_atmvars( ppath_p, ppath_t, ppath_vmr, 
                         ppath_wind_u, ppath_wind_v, ppath_wind_w,
                         ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                         wind_field_dummy, wind_field_dummy, wind_field_dummy );

      // Particle number densities
      get_ppath_pnd( ppath_pnd, ppath, atmosphere_dim, cloudbox_limits, 
                     pnd_field );

      // Absorption and optical thickness for each step
      get_ppath_cloudrtvars( ws, ppath_asp_abs_vec, ppath_asp_ext_mat,  
                             ppath_pnd_abs_vec, ppath_pnd_ext_mat, 
                             ppath_transmission, total_transmission, 
                             emission_dummy, scat_data, abs_scalar_gas_agenda, 
                             agenda_dummy, opt_prop_gas_agenda,
                             ppath, ppath_p, ppath_t, ppath_vmr, ppath_wind_u,  
                             ppath_wind_v, ppath_wind_w, ppath_pnd, 
                             use_mean_scat_data, scat_data_raw, stokes_dim, 
                             f_grid, atmosphere_dim, 0, verbosity );

      // *scat_data* not used. Free the memory.
      scat_data.resize( 0 );
    }
  else // Just in case, should not happen
    { assert( 0 ); }

  // iy_transmission
  //
  iy_transmission_mult( iy_trans_new, iy_transmission, total_transmission );

  // Get iy for unscattered direction
  //
  // Note that the Ppath positions (ppath.pos) for 1D have one column more
  // than expected by most functions. Only the first atmosphere_dim values
  // shall be copied.
  //
  {
    Vector  rte_pos2, rte_los2;
    rte_pos2 = ppath.pos(ppath.np-1,Range(0,atmosphere_dim));
    rte_los2 = ppath.los(ppath.np-1,joker);
    //
    iy_clearsky_agendaExecute( ws, iy, iy_error, iy_error_type,
                               iy_aux, diy_dx, 0, iy_trans_new,
                               rte_pos2, rte_los2, 0, jacobian_do, t_field, 
                               z_field, vmr_field, -1, iy_clearsky_agenda );
  }

  // Without jacobians the complete RT calculations are just multiplication
  // with *total_transmission*
  for( Index iv=0; iv<nf; iv++ )
    {
      const Vector tmp = iy(iv,joker);
      mult( iy(iv,joker), total_transmission(iv,joker,joker), tmp );
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void iyEmissionStandardClearsky(
        Workspace&                  ws,
        Matrix&                     iy,
        Matrix&                     iy_error,
        Index&                      iy_error_type,
        Matrix&                     iy_aux,
        ArrayOfTensor3&             diy_dx,
  const Index&                      iy_agenda_call1,
  const Tensor3&                    iy_transmission,
  const Vector&                     rte_pos,      
  const Vector&                     rte_los,      
  const Index&                      jacobian_do,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Tensor3&                    z_field,
  const Tensor3&                    t_field,
  const Tensor4&                    vmr_field,
  const Tensor3&                    wind_u_field,
  const Tensor3&                    wind_v_field,
  const Tensor3&                    wind_w_field,
  const Tensor3&                    edensity_field,
  const Vector&                     refellipsoid,
  const Matrix&                     z_surface,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const ArrayOfArrayOfSpeciesTag&   abs_species,
  const Index&                      mblock_index,
  const Agenda&                     ppath_agenda,
  const Agenda&                     emission_agenda,
  const Agenda&                     abs_scalar_gas_agenda,
  const Agenda&                     iy_clearsky_agenda,
  const Agenda&                     iy_space_agenda,
  const Agenda&                     surface_prop_agenda,
  const Agenda&                     iy_cloudbox_agenda,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices,
  const Verbosity&                  verbosity )
{
  // The method can in principle be used "stand-alone", but for efficiency
  // reasons we skip all checks needed to handle such usage. Those checks are
  // done by yCalc.

  // iy_error, iy_error_type, iy_aux and diy_dx are initiated by iyb_calc
  // to be empty/zero. Accordingly, if any radiative background wants
  // to flag an error, iy_error must be resized. No subsequent scaling of the
  // error is made. Thus the iy_transmission should be considered when adding
  // a term to the error. A setting of iy_aux must also include a resizing. 

  // The WSV *iy_transmission* holds the transmission between the sensor
  // and the end of the propagation path.

  // Sizes
  //
  const Index nf = f_grid.nelem();
  const Index nq = jacobian_quantities.nelem();

  // Determine if there are any analytical jacobians to handle, and if primary
  // call, resize diy_dx.
  //
  Index j_analytical_do = 0;
  //
  if( jacobian_do ) { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }
  //
  if( iy_agenda_call1 && j_analytical_do )
    {
      diy_dx.resize( nq ); 
      //
      FOR_ANALYTICAL_JACOBIANS_DO( 
        diy_dx[iq].resize( jacobian_indices[iq][1] - jacobian_indices[iq][0] +
                           1, nf, stokes_dim ); 
        diy_dx[iq] = 0.0;
      )
    }

  //- Determine propagation path
  //
  Ppath  ppath;
  //
  ppath_agendaExecute( ws, ppath, rte_pos, rte_los, cloudbox_on, 0, 
                       mblock_index, t_field, z_field, vmr_field, 
                       edensity_field, -1, ppath_agenda );

  // Get atmospheric and RT quantities for each ppath point/step
  //
  // If np = 1, we only need to determine the radiative background
  //
  // "atmvars"
  Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w;
  Matrix    ppath_vmr;
  // "rtvars"
  Vector    total_tau;
  Matrix    ppath_emission, ppath_tau;
  Tensor3   ppath_abs_scalar, iy_trans_new;
  //
  const Index np  = ppath.np;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      get_ppath_atmvars( ppath_p, ppath_t, ppath_vmr, 
                         ppath_wind_u, ppath_wind_v, ppath_wind_w,
                         ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                         wind_u_field, wind_v_field, wind_w_field );

      // Get emission, absorption and optical thickness for each step
      get_ppath_rtvars( ws, ppath_abs_scalar, ppath_tau, total_tau, 
                   ppath_emission, abs_scalar_gas_agenda, emission_agenda, 
                   ppath, ppath_p, ppath_t, ppath_vmr, ppath_wind_u, 
                   ppath_wind_v, ppath_wind_w, f_grid, -1, atmosphere_dim, 1 );
    }
  else // To handle cases inside the cloudbox, or totally outside the atmosphere
    { 
      total_tau.resize( nf );
      total_tau = 0;
    }

  // iy_transmission
  //
  if( iy_agenda_call1 )
    { iy_transmission_for_scalar_tau( iy_trans_new, stokes_dim, total_tau ); }
  else
    { iy_transmission_mult_scalar_tau( iy_trans_new, iy_transmission, 
                                                                total_tau ); }

  // Radiative background
  //
  get_iy_of_background( ws, iy, iy_error, iy_error_type, iy_aux, diy_dx, 
                        iy_trans_new, jacobian_do, ppath, atmosphere_dim, 
                        lat_grid, lon_grid, t_field, z_field, 
                        vmr_field, refellipsoid, z_surface, cloudbox_on, 
                        stokes_dim, f_grid, iy_clearsky_agenda, iy_space_agenda,
                        surface_prop_agenda, iy_cloudbox_agenda, verbosity);
  
  // Do RT calculations
  //
  if( np > 1 )
    {
      // Create container for the derivatives with respect to changes
      // at each ppath point. And find matching abs_species-index and 
      // "temperature flag" (for analytical jacobians).
      //
      ArrayOfTensor3  diy_dpath; 
      ArrayOfIndex    abs_species_i, is_t; 
      //
      const Numeric   dt = 0.1;
            Matrix    ppath_e2;
            Tensor3   ppath_as2;
      //
      if( j_analytical_do )
        { 
          diy_dpath.resize( nq ); 
          abs_species_i.resize( nq ); 
          is_t.resize( nq ); 
          //
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_dpath[iq].resize( np, nf, stokes_dim ); 
            diy_dpath[iq] = 0.0;
          )
          get_pointers_for_analytical_jacobians( abs_species_i, is_t, 
                                            jacobian_quantities, abs_species );
          //
          // Determine if temperature is among the analytical jac. quantities.
          // If yes, get emission and absorption for disturbed temperature
          Index do_t = 0;
          for( Index iq=0; iq<is_t.nelem() && do_t==0; iq++ )
            { if( is_t[iq] ) { do_t = 1; } }
          if( do_t )
            {
              Matrix tau_dummy; 
              Vector total_tau_dummy, t2 = ppath_t;
              t2 += dt;
              get_ppath_rtvars( ws, ppath_as2, tau_dummy, total_tau_dummy, 
                                ppath_e2, abs_scalar_gas_agenda, 
                                emission_agenda, ppath, ppath_p, t2, ppath_vmr, 
                                ppath_wind_u, ppath_wind_v, ppath_wind_w, 
                                f_grid, -1, atmosphere_dim, 1 );
            }
        }

      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {
          Vector esource(nf);            // B-bar
          Matrix iy_new(nf,stokes_dim);  // Intensity at end of step

          // Spectrum at end of ppath step 
          // (iy updated to iy_new at end of ip-loop)
          for( Index iv=0; iv<nf; iv++ )  
            {
              const Numeric step_tr = exp( -ppath_tau(iv,ip) );
              //
              // First Stokes element:
              esource[iv] = 0.5 * ( ppath_emission(iv,ip) + 
                                    ppath_emission(iv,ip+1) );    
              iy_new(iv,0) = iy(iv,0) * step_tr + esource[iv] * (1-step_tr); 
              //
              // Higher Stokes:
              for( Index is=1; is<stokes_dim; is++ )
                { iy_new(iv,is) = step_tr * iy(iv,is); }
            }

          // Jacobian quantities
          if( j_analytical_do )
            {
              // Common terms introduced for efficiency and clarity
              Vector X(nf), Y(nf);   // See AUG
              //
              for( Index iv=0; iv<nf; iv++ )
                {
                  X[iv] = 0.5 * ppath.lstep[ip] * exp( -total_tau[iv] );
                  Y[iv] = X[iv] * ( esource[iv] - iy(iv,0) );
                }

              // Loop quantities
              for( Index iq=0; iq<nq; iq++ ) 
                {
                  if( jacobian_quantities[iq].Analytical() )
                    {
                      // Absorbing species
                      const Index isp = abs_species_i[iq]; 
                      if( isp >= 0 )
                        {
                          // Scaling factors to handle retrieval unit
                          // (gives the cross-section to apply)
                          Numeric unitscf1, unitscf2;
                          vmrunitscf( unitscf1, 
                                      jacobian_quantities[iq].Mode(), 
                                      ppath_vmr(isp,ip), ppath_p[ip], 
                                      ppath_t[ip] );
                          vmrunitscf( unitscf2, 
                                      jacobian_quantities[iq].Mode(), 
                                      ppath_vmr(isp,ip+1), ppath_p[ip+1], 
                                      ppath_t[ip+1] );
                          //
                          for( Index iv=0; iv<nf; iv++ )
                            {
                              // Stokes component 1
                              diy_dpath[iq](ip  ,iv,0) += Y[iv] * unitscf1 *
                                                 ppath_abs_scalar(iv,isp,ip);
                              diy_dpath[iq](ip+1,iv,0) += Y[iv] * unitscf2 * 
                                                 ppath_abs_scalar(iv,isp,ip+1);
                              // Higher stokes components
                              for( Index is=1; is<stokes_dim; is++ )
                                { 
                                  const Numeric Z = -X[iv] * iy(iv,is); 
                                  diy_dpath[iq](ip  ,iv,is) += Z * unitscf1 *
                                                 ppath_abs_scalar(iv,isp,ip);
                                  diy_dpath[iq](ip+1,iv,is) += Z * unitscf2 *
                                                 ppath_abs_scalar(iv,isp,ip+1);
                                }
                            }
                        }

                      // Temperature
                      else if( is_t[iq] )
                        {
                          for( Index iv=0; iv<nf; iv++ )
                            {
                              // The terms associated with Dtau/Dk:
                              const Numeric k1 = 
                                         ppath_abs_scalar(iv,joker,ip  ).sum();
                              const Numeric k2 = 
                                         ppath_abs_scalar(iv,joker,ip+1).sum();
                              const Numeric dkdt1 =
                                  ( ppath_as2(iv,joker,ip  ).sum() - k1 ) / dt;
                              const Numeric dkdt2 =
                                  ( ppath_as2(iv,joker,ip+1).sum() - k2 ) / dt;
                              // Stokes 1:
                              diy_dpath[iq](ip  ,iv,0) += Y[iv] * dkdt1;
                              diy_dpath[iq](ip+1,iv,0) += Y[iv] * dkdt2;
                              // Higher Stokes
                              for( Index is=1; is<stokes_dim; is++ )
                                { 
                                  const Numeric Z = -X[iv] * iy(iv,is);
                                  diy_dpath[iq](ip  ,iv,is) += Z * dkdt1;
                                  diy_dpath[iq](ip+1,iv,is) += Z * dkdt2;
                                }
                              //
                              // The terms associated with B-bar:
                              const Numeric V = 0.5 * 
                                      exp(-(total_tau[iv]-ppath_tau(iv,ip))) *
                                              ( 1.0 - exp(-ppath_tau(iv,ip)) );
                              diy_dpath[iq](ip  ,iv,0) += V *
                                              ( ppath_e2(iv,ip) -
                                                ppath_emission(iv,ip) ) / dt;
                              diy_dpath[iq](ip+1,iv,0) += V * 
                                              ( ppath_e2(iv,ip+1) -
                                                ppath_emission(iv,ip+1) ) / dt;
                              // Zero for higher Stokes
                              //
                              // The terms associated with Delta-s:
                              if( jacobian_quantities[iq].Subtag() == "HSE on" )
                                {
                                  // Stokes 1:
                                  const Numeric kbar = 0.5 * ( k1 + k2 );
                                  diy_dpath[iq](ip  ,iv,0) += Y[iv] * kbar /
                                                                  ppath_t[ip];
                                  diy_dpath[iq](ip+1,iv,0) += Y[iv] * kbar /
                                                                 ppath_t[ip+1];
                                  // Higher Stokes
                                  for( Index is=1; is<stokes_dim; is++ )
                                    { 
                                      const Numeric Z = -X[iv] * iy(iv,is);
                                      diy_dpath[iq](ip  ,iv,is) += Z * kbar /
                                                                  ppath_t[ip];
                                      diy_dpath[iq](ip+1,iv,is) += Z * kbar /
                                                                 ppath_t[ip+1];
                                    }
                                } //hse
                            }  // frequency
                        }  // if is_t
                    } // if analytical
                } // for iq
              
              // Update total_tau (not used for spectrum calculation)
              for( Index iv=0; iv<nf; iv++ )
                { total_tau[iv] -= ppath_tau(iv,ip); }
            }

          // Update iy
          iy = iy_new;
        } 

      // Map jacobians from ppath to retrieval grids
      // (this operation corresponds to the term Dx_i/Dx)
      if( j_analytical_do )
        { 
          // Weight with iy_transmission
          if( !iy_agenda_call1 )
            {
              Matrix X, Y(stokes_dim,diy_dpath[0].npages()); 
              //
              FOR_ANALYTICAL_JACOBIANS_DO( 
                for( Index iv=0; iv<nf; iv++ )
                  { 
                    X = transpose( diy_dpath[iq](joker,iv,joker) );
                    mult( Y, iy_transmission(iv,joker,joker), X );
                    diy_dpath[iq](joker,iv,joker) = transpose( Y );
                  }
               )
            }

          // Map to retrieval grids
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_from_path_to_rgrids( diy_dx[iq], jacobian_quantities[iq], 
                               diy_dpath[iq], atmosphere_dim, ppath, ppath_p );
          )
        }
    } // if np>1

  // Set iy_aux, if background is TOA or surface
  //
  if( ppath_what_background( ppath ) <= 2 )
    {
      iy_aux.resize( nf, stokes_dim ); 
      //
      for( Index iv=0; iv<nf; iv++ )
        { 
          for( Index is=0; is<stokes_dim; is++ )
            { 
              iy_aux( iv, is ) = iy_trans_new( iv, is, is ); 
            }
        }
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void iyEmissionStandardClearskyBasic(
        Workspace&     ws,
        Matrix&        iy,
  const Vector&        rte_pos,      
  const Vector&        rte_los,      
  const Index&         jacobian_do,
  const Index&         atmosphere_dim,
  const Vector&        p_grid,
  const Vector&        lat_grid,
  const Vector&        lon_grid,
  const Tensor3&       z_field,
  const Tensor3&       t_field,
  const Tensor4&       vmr_field,
  const Tensor3&       wind_u_field,
  const Tensor3&       wind_v_field,
  const Tensor3&       wind_w_field,
  const Tensor3&       edensity_field,
  const Vector&        refellipsoid,
  const Matrix&        z_surface,
  const Index&         cloudbox_on,
  const Index&         stokes_dim,
  const Vector&        f_grid,
  const Agenda&        ppath_agenda,
  const Agenda&        emission_agenda,
  const Agenda&        abs_scalar_gas_agenda,
  const Agenda&        iy_clearsky_basic_agenda,
  const Agenda&        iy_space_agenda,
  const Agenda&        surface_prop_agenda,
  const Agenda&        iy_cloudbox_agenda,
  const Verbosity&     verbosity )
{
  // See initial comments of iyEmissionStandardClearsky

  if( jacobian_do )
    throw runtime_error( 
        "This method does not provide any jacobians (jacobian_do must be 0)" );

  // Sizes
  //
  const Index nf = f_grid.nelem();

  //- Determine propagation path
  //
  Ppath  ppath;
  //
  ppath_agendaExecute( ws, ppath, rte_pos, rte_los, cloudbox_on, 0, -1,
                       t_field, z_field, vmr_field, edensity_field, -1,
                       ppath_agenda );

  // Get quantities required for each ppath point
  //
  // If np = 1, we only need to determine the radiative background
  // "atmvars"
  Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w;
  Matrix    ppath_vmr;
  // "rtvars"
  Vector    total_tau;
  Matrix    ppath_emission, ppath_tau;
  Tensor3   ppath_abs_scalar, iy_trans_new;
  //
  const Index np  = ppath.np;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      get_ppath_atmvars( ppath_p, ppath_t, ppath_vmr, 
                         ppath_wind_u, ppath_wind_v, ppath_wind_w,
                         ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                         wind_u_field, wind_v_field, wind_w_field );

      // Get emission, absorption and optical thickness for each step
      get_ppath_rtvars( ws, ppath_abs_scalar, ppath_tau, total_tau, 
                   ppath_emission, abs_scalar_gas_agenda, emission_agenda, 
                   ppath, ppath_p, ppath_t, ppath_vmr, ppath_wind_u,
                   ppath_wind_v, ppath_wind_w, f_grid, -1, atmosphere_dim, 1 );
    }

  // Radiative background
  //
  Matrix          dummy1, dummy3;
  Index           dummy2;
  ArrayOfTensor3  dummy4;
  Tensor3         dummy5;
  //  
  get_iy_of_background( ws, iy, dummy1, dummy2, dummy3, dummy4, dummy5,  
                        0, ppath, atmosphere_dim, lat_grid, lon_grid,
                        t_field, z_field, vmr_field, refellipsoid, z_surface, 
                        cloudbox_on, stokes_dim, f_grid, 
                        iy_clearsky_basic_agenda, iy_space_agenda, 
                        surface_prop_agenda, iy_cloudbox_agenda, verbosity );

  // Do RT calculations
  //
  if( np > 1 )
    {
      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {
          // Loop frequencies
          for( Index iv=0; iv<nf; iv++ )
            { 
              const Numeric step_tr = exp( -ppath_tau(iv,ip) );
              //
              iy(iv,0) = iy(iv,0) * step_tr + 0.5 * (1-step_tr) * 
                           ( ppath_emission(iv,ip) + ppath_emission(iv,ip+1) ); 
              //
              for( Index is=1; is<stokes_dim; is++ )
                { iy(iv,is) *= step_tr; }
            }
        } 
    }
}


 


/* Workspace method: Doxygen documentation will be auto-generated */
void iyMC(
        Workspace&                  ws,
        Matrix&                     iy,
        Matrix&                     iy_error,
        Index&                      iy_error_type,
        Matrix&                     iy_aux,
        ArrayOfTensor3&             diy_dx,
  const Index&                      iy_agenda_call1,
  const Tensor3&                    iy_transmission,
  const Vector&                     rte_pos,      
  const Vector&                     rte_los,      
  const Index&                      jacobian_do,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Tensor3&                    z_field,
  const Tensor3&                    t_field,
  const Tensor4&                    vmr_field,
  const Vector&                     refellipsoid,
  const Matrix&                     z_surface,
  const Index&                      cloudbox_on,
  const ArrayOfIndex&               cloudbox_limits,
  const Index&                      cloudbox_checked,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const ArrayOfSingleScatteringData&   scat_data_raw,
  const Agenda&                     iy_space_agenda,
  const Agenda&                     surface_prop_agenda,
  const Agenda&                     abs_scalar_gas_agenda, 
  const Agenda&                     opt_prop_gas_agenda,
  const Tensor4&                    pnd_field,
  const String&                     y_unit,
  const Numeric&                    mc_std_err,
  const Index&                      mc_max_time,
  const Index&                      mc_max_iter,
  const Verbosity&                  verbosity)
{
  // Throw error if unsupported features are requested
  if( atmosphere_dim != 3 )
    throw runtime_error( 
                "Only 3D atmospheres are allowed (atmosphere_dim must be 3)" );
  if( !cloudbox_on )
    throw runtime_error( 
                    "The cloudbox must be activated (cloudbox_on must be 1)" );
  if( jacobian_do )
    throw runtime_error( 
        "This method does not provide any jacobians (jacobian_do must be 0)" );
  if( !iy_agenda_call1 )
    throw runtime_error( 
                  "Recursive usage not possible (iy_agenda_call1 must be 1)" );
  if( iy_transmission.ncols() )
    throw runtime_error( "*iy_transmission* must be empty" );

  // See initial comments of iyEmissionStandardClearsky

  // Size output variables
  //
  const Index   nf  = f_grid.nelem();
  //
  iy.resize( nf, stokes_dim );
  iy_error.resize( nf, stokes_dim );
  iy_error_type = 1;

  // To avoid warning about unused variables
  iy_aux.resize(0,0);
  diy_dx.resize(0);

  // Some MC variables are only local here
  Tensor3  mc_points;
  Index    mc_iteration_count, mc_seed;
  //
  MCAntenna mc_antenna;
  mc_antenna.set_pencil_beam();

  // Pos and los must be matrices 
  Matrix pos(1,3), los(1,2);
  //
  pos(0,joker) = rte_pos;
  los(0,joker) = rte_los;

  for( Index f_index=0; f_index<nf; f_index++ )
    {
      ArrayOfSingleScatteringData   scat_data_mono;
      
      scat_data_monoCalc( scat_data_mono, scat_data_raw, 
                                          f_grid, f_index, verbosity );

      // Seed reset for each loop. If not done, the errors 
      // appear to be highly correlated.
      MCSetSeedFromTime( mc_seed, verbosity );

      Vector y, mc_error;
                  
      MCGeneral( ws, y, mc_iteration_count, mc_error, mc_points, mc_antenna, 
                 f_grid, f_index, pos, los, stokes_dim, atmosphere_dim,
                 iy_space_agenda, surface_prop_agenda, opt_prop_gas_agenda,
                 abs_scalar_gas_agenda, p_grid, lat_grid, lon_grid, z_field, 
                 refellipsoid, z_surface, t_field, vmr_field, 
                 cloudbox_on, cloudbox_limits, 
                 pnd_field, scat_data_mono, 1, cloudbox_checked,
                 mc_seed, y_unit, mc_std_err, mc_max_time, mc_max_iter,
                 verbosity); 
                 // GH 2011-06-17, mc_z_field_is_1D);

      assert( y.nelem() == stokes_dim );

      // Data returned can not be in Tb
      if ( y_unit == "RJBT" )
        { 
          const Numeric scfac = invrayjean( 1, f_grid[f_index] );
          y /= scfac;
          mc_error /= scfac;
        }
     
      iy(f_index,joker) = y;
      iy_error(f_index,joker) = mc_error;
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void iyRadioLink(
        Workspace&            ws,
        Matrix&               iy,
        Matrix&               iy_error _U_,
        Index&                iy_error_type _U_,
        Matrix&               iy_aux,
        ArrayOfTensor3&       diy_dx _U_,
  const Index&                iy_agenda_call1,
  const Tensor3&              iy_transmission _U_,
  const Vector&               rte_pos,      
  const Index&                jacobian_do,
  const Index&                atmosphere_dim,
  const Vector&               p_grid,
  const Vector&               lat_grid,
  const Vector&               lon_grid,
  const Tensor3&              z_field,
  const Tensor3&              t_field,
  const Tensor4&              vmr_field,
  const Tensor3&              wind_u_field,
  const Tensor3&              wind_v_field,
  const Tensor3&              wind_w_field,
  const Tensor3&              edensity_field,
  const Vector&               refellipsoid,
  const Index&                cloudbox_on,
  const Index&                stokes_dim,
  const Vector&               f_grid,
  const Index&                dispersion_do,
  const Index&                mblock_index,
  const Agenda&               ppath_agenda,
  const Agenda&               abs_scalar_gas_agenda,
  const Agenda&               iy_space_agenda,
  const String&               aux_var,
  const Verbosity&            verbosity )
{
  // See initial comments of iyEmissionStandardClearsky

  if( !iy_agenda_call1 )
    throw runtime_error( 
                  "Recursive usage not possible (iy_agenda_call1 must be 1)" );

  // Jacobians not yet handled. Set vars to be empty.
  if( jacobian_do )
    throw runtime_error( "This method does not yet provide any jacobians and "
                         "*jacobian_do* must be 0." );
  CREATE_OUT2
  if( dispersion_do )
    { out2 << "  Dispersion mode selected !!!\n"; }

  // Sizes
  //
  Index nf = f_grid.nelem();
  //
  iy_aux.resize( nf, stokes_dim ); 
  iy_aux = 0;  // Remove when iy_aux handled

  // Dispersion?
  Index nloops=1; 
  if( dispersion_do )
    { 
      iy.resize( nf, stokes_dim );  // iy_space_handles this if no dispersion
      nloops = nf; 
      nf     = 1;
    }  
  
  // Different ppath variables
  Ppath  ppath;
  Numeric lbg=0;  // Bent geometrical length of ray path
  Numeric lba=0;  // Bent apparent length of ray path

  for( Index i=0; i<nloops; i++ )
    {
      Index f_index = -1; if( dispersion_do ){ f_index = i; }

      //- Determine propagation path
      Vector rte_los(0);  // Dummy value
      //
      ppath_agendaExecute( ws, ppath, rte_pos, rte_los, cloudbox_on, 0,
                           mblock_index, t_field, z_field, vmr_field,
                           edensity_field, f_index, ppath_agenda );
      if( ppath_what_background(ppath) > 2 )
        { throw runtime_error( "Radiative background not set to \"space\" by "
                      "*ppath_agenda*. Is correct WSM used in the agenda?" ); }

      // Get atmospheric and RT quantities for each ppath point/step
      //
      // "atmvars":
      Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w;
      Matrix    ppath_vmr;
      // "rtvars":
      Vector    total_tau;
      Matrix    emission_dummy, ppath_tau;
      Tensor3   ppath_abs_scalar, iy_trans_new;
      Agenda    agenda_dummy;
      //
      const Index np  = ppath.np;
      //
      if( np > 1 )
        {
          // Get pressure, temperature and VMRs
          get_ppath_atmvars( ppath_p, ppath_t, ppath_vmr, 
                             ppath_wind_u, ppath_wind_v, ppath_wind_w,
                             ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                             wind_u_field, wind_v_field, wind_w_field );

          // Absorption and optical thickness for each step
          get_ppath_rtvars( ws, ppath_abs_scalar, ppath_tau, total_tau, 
                            emission_dummy,  abs_scalar_gas_agenda, 
                            agenda_dummy, ppath, ppath_p, ppath_t, ppath_vmr, 
                            ppath_wind_u, ppath_wind_v, ppath_wind_w, f_grid, 
                            f_index, atmosphere_dim, 0 );
        }

      // Radiative background
      //
      if( dispersion_do )
        { 
          Matrix iy1;  // Single frequency version of iy
          iy_space_agendaExecute( ws, iy1, rte_pos, rte_los, 
                                Vector(1,f_grid[f_index]), iy_space_agenda ); 
          if( iy1.ncols() != stokes_dim  ||  iy1.nrows() != 1 )
            { throw runtime_error( "The size of *iy* returned from "
                                   "*iy_space_agenda* is not correct." ); }
          iy(f_index,joker) = iy1(0,joker);
        }
      else
        { 
          iy_space_agendaExecute( ws, iy, rte_pos, rte_los, f_grid,
                                                           iy_space_agenda ); 
          if( iy.ncols() != stokes_dim  ||  iy.nrows() != nf )
            { throw runtime_error( "The size of *iy* returned from "
                                   "*iy_space_agenda* is not correct." ); }
        }
      //      

      // Do RT calculations
      //
      lbg = ppath.lspace;
      lba = ppath.lspace;
      //
      if( np > 1 )
        {
          // Loop ppath steps
          for( Index ip=np-2; ip>=0; ip-- )
            {
              lbg += ppath.lstep[ip];
              lba += ppath.lstep[ip] * (ppath.nreal[ip]+ppath.nreal[ip+1])/2.0;

              if( dispersion_do )
                {
                  const Numeric step_tr = exp( -ppath_tau(0,ip) );
                  for( Index is=0; is<stokes_dim; is++ )
                    { iy(f_index,is) *= step_tr; }
                }
              else
                {
                  // Loop frequencies
                  for( Index iv=0; iv<nf; iv++ )
                    {
                      const Numeric step_tr = exp( -ppath_tau(iv,ip) );
                      for( Index is=0; is<stokes_dim; is++ )
                        { iy(iv,is) *= step_tr; }
                    }
                }
            } 
        }
      
      // Geomtrical distance between start and end point
      Numeric lgd = 0;
      {
        // Radius of rte_pos and rte_pos2
        const Numeric r1 = pos2refell_r( atmosphere_dim, refellipsoid, lat_grid,
                                  lon_grid, ppath.end_pos ) + ppath.end_pos[0];
        const Numeric r2 = pos2refell_r( atmosphere_dim, refellipsoid, lat_grid,
                              lon_grid, ppath.start_pos ) + ppath.start_pos[0];
        if( atmosphere_dim <= 2 )
          { distance2D( lgd, r1, ppath.end_pos[1], r2, ppath.start_pos[1] ); }
        else 
          { distance3D( lgd, r1, ppath.end_pos[1],   ppath.end_pos[2],
                         r2, ppath.start_pos[1], ppath.start_pos[2] ); }
      }

      // Geometric LOS from rte_pos to rte_pos2
      Vector rte_los_geom;
      rte_losGeometricFromRtePosToRtePos2( rte_los_geom, atmosphere_dim, 
                                         lat_grid, lon_grid, refellipsoid, 
                                         rte_pos, ppath.start_pos, verbosity );

      // Free space loss
      const Numeric fspl = 1 / ( 4 * PI * lbg*lbg );

      // Extra path delay
      const Numeric epd = ( lba - lgd ) / SPEED_OF_LIGHT;

      // Zenith bending angle
      const Numeric zba = rte_los_geom[0]-ppath.end_los[0];

      // Azimuth bending angle
      Numeric aba = 0;
      if( atmosphere_dim == 3 )
        { aba = rte_los_geom[1]-ppath.end_los[1]; }

      out2 << "  Atmospheric attenuation : " << iy(joker,0) << "\n";
      out2 << "          Free space loss : " << fspl << "\n";
      out2 << "         Extra path delay : " << 1e9*epd << " ns\n";
      out2 << "     Zenith bending angle : " << zba << " deg.\n";
      if( atmosphere_dim == 3 )
        { out2 << "    Azimuth bending angle : " << aba << " deg.\n"; }

      // Apply loss terms and fill iy_aux
      if( dispersion_do )
        {
          if( aux_var == "AtmosAtten" )
            { iy_aux(f_index,joker) = iy(f_index,joker); }
          else if( aux_var == "FreeSpaceLoss" )
            { iy_aux(f_index,0) = fspl; }
          else if( aux_var == "DefocusingLoss" )
            { iy_aux(f_index,0) = 0; }
          else if( aux_var == "ExtraPathDelay" )
            { iy_aux(f_index,0) = epd; }
          else if( aux_var == "ZenithBendingAngle" )
            { iy_aux(f_index,0) = zba; }
          else if( aux_var == "AzimuthBendingAngle" )
            { iy_aux(f_index,0) = aba; }
          //
          for( Index is=0; is<stokes_dim; is++ )
            { iy(f_index,is) *= fspl; }
        }
      else
        {
          if( aux_var == "AtmosAtten" )
            { iy_aux = iy; }
          else if( aux_var == "FreeSpaceLoss" )
            { iy_aux(joker,0) = fspl; }
          else if( aux_var == "DefocusingLoss" )
            { iy_aux(joker,0) = 0; }
          else if( aux_var == "ExtraPathDelay" )
            { iy_aux(joker,0) = epd; }
          else if( aux_var == "ZenithBendingAngle" )
            { iy_aux(joker,0) = zba; }
          else if( aux_var == "AzimuthBendingAngle" )
            { iy_aux(joker,0) = aba; }
          //
          iy *= fspl;
        }
    }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void iyCalc(
        Workspace&   ws,
        Matrix&      iy,
        Matrix&      iy_aux,
        Matrix&      iy_error,
        Index&       iy_error_type,
  ArrayOfTensor3&    diy_dx,
  const Index&       basics_checked,
  const Tensor3&     t_field,
  const Tensor3&     z_field,
  const Tensor4&     vmr_field,
  const Index&       cloudbox_on,
  const Index&       cloudbox_checked,
  const Vector&      rte_pos,
  const Vector&      rte_los,
  const Index&       jacobian_do,
  const Index&       mblock_index,
  const Agenda&      iy_clearsky_agenda,
  const Verbosity& )
{
  // Basics and cloudbox OK?
  //
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );

  // Make sure that output is empty if not set by the agenda
  iy.resize(0,0);
  iy_aux.resize(0,0);
  iy_error.resize(0,0);
  iy_error_type = 0;
  diy_dx.resize(0);

  // iy_transmission is just input and can be left empty for first call
  Tensor3   iy_transmission(0,0,0);

  iy_clearsky_agendaExecute( ws, iy, iy_error, iy_error_type, iy_aux, diy_dx, 
                             1, iy_transmission, rte_pos, rte_los, cloudbox_on,
                             jacobian_do, t_field, z_field, vmr_field, 
                             mblock_index, iy_clearsky_agenda );
}





/* Workspace method: Doxygen documentation will be auto-generated */
void yCalc(
        Workspace&                  ws,
        Vector&                     y,
        Vector&                     y_f,
        ArrayOfIndex&               y_pol,
        Matrix&                     y_pos,
        Matrix&                     y_los,
        Vector&                     y_error,
        Vector&                     y_aux,
        Matrix&                     jacobian,
  const Index&                      basics_checked,
  const Index&                      atmosphere_dim,
  const Tensor3&                    t_field,
  const Tensor3&                    z_field,
  const Tensor4&                    vmr_field,
  const Index&                      cloudbox_on,
  const Index&                      cloudbox_checked,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const Matrix&                     sensor_pos,
  const Matrix&                     sensor_los,
  const Vector&                     mblock_za_grid,
  const Vector&                     mblock_aa_grid,
  const Index&                      antenna_dim,
  const Sparse&                     sensor_response,
  const Vector&                     sensor_response_f,
  const ArrayOfIndex&               sensor_response_pol,
  const Vector&                     sensor_response_za,
  const Vector&                     sensor_response_aa,
  const Agenda&                     iy_clearsky_agenda,
  const String&                     y_unit,
  const Agenda&                     jacobian_agenda,
  const Index&                      jacobian_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices,
  const Verbosity&                  verbosity )
{
  // Some sizes
  const Index   nf      = f_grid.nelem();
  const Index   nza     = mblock_za_grid.nelem();
        Index   naa     = mblock_aa_grid.nelem();   
  if( antenna_dim == 1 )  
    { naa = 1; }
  const Index   n1y     = sensor_response.nrows();
  const Index   nmblock = sensor_pos.nrows();
  const Index   niyb    = nf * nza * naa * stokes_dim;


  //---------------------------------------------------------------------------
  // Input checks
  //---------------------------------------------------------------------------

  // Basics and cloudbox OK?
  //
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );

  // Sensor position and LOS.
  //
  if( sensor_pos.ncols() != atmosphere_dim )
    throw runtime_error( "The number of columns of sensor_pos must be "
                         "equal to the atmospheric dimensionality." );
  if( atmosphere_dim <= 2  &&  sensor_los.ncols() != 1 )
    throw runtime_error( "For 1D and 2D, sensor_los shall have one column." );
  if( atmosphere_dim == 3  &&  sensor_los.ncols() != 2 )
    throw runtime_error( "For 3D, sensor_los shall have two columns." );
  if( sensor_los.nrows() != nmblock )
    {
      ostringstream os;
      os << "The number of rows of sensor_pos and sensor_los must be "
         << "identical, but sensor_pos has " << nmblock << " rows,\n"
         << "while sensor_los has " << sensor_los.nrows() << " rows.";
      throw runtime_error( os.str() );
    }
  if( max( sensor_los(joker,0) ) > 180 )
    throw runtime_error( 
     "First column of *sensor_los* is not allowed to have values above 180." );
  if( atmosphere_dim == 2 )
    {
      if( min( sensor_los(joker,0) ) < -180 )
          throw runtime_error( "For atmosphere_dim = 2, first column of "
                    "*sensor_los* is not allowed to have values below -180." );
    }     
  else
    {
      if( min( sensor_los(joker,0)  ) < 0 )
          throw runtime_error( "For atmosphere_dim != 2, first column of "
                       "*sensor_los* is not allowed to have values below 0." );
    }    
  if( atmosphere_dim == 3  &&  max( sensor_los(joker,1) ) > 180 )
    throw runtime_error( 
    "Second column of *sensor_los* is not allowed to have values above 180." );

  // Antenna
  //
  chk_if_in_range( "antenna_dim", antenna_dim, 1, 2 );
  //
  if( nza == 0 )
    throw runtime_error( "The measurement block zenith angle grid is empty." );
  chk_if_increasing( "mblock_za_grid", mblock_za_grid );
  //
  if( antenna_dim == 1 )
    {
      if( mblock_aa_grid.nelem() != 0 )
        throw runtime_error( 
          "For antenna_dim = 1, the azimuthal angle grid must be empty." );
    }
  else
    {
      if( atmosphere_dim < 3 )
        throw runtime_error( "2D antennas (antenna_dim=2) can only be "
                                                 "used with 3D atmospheres." );
      if( mblock_aa_grid.nelem() == 0 )
        throw runtime_error(
                      "The measurement block azimuthal angle grid is empty." );
      chk_if_increasing( "mblock_aa_grid", mblock_aa_grid );
    }

  // Sensor
  //
  if( sensor_response.ncols() != niyb ) 
    {
      ostringstream os;
      os << "The *sensor_response* matrix does not have the right size,\n"
         << "either the method *sensor_responseInit* has not been run or some\n"
         << "of the other sensor response methods has not been correctly\n"
         << "configured.";
      throw runtime_error( os.str() );
    }

  // Sensor aux variables
  //
  if( n1y != sensor_response_f.nelem()  || n1y != sensor_response_pol.nelem() ||
      n1y != sensor_response_za.nelem() || n1y != sensor_response_za.nelem() )
    {
      ostringstream os;
      os << "Sensor auxiliary variables do not have the correct size.\n"
         << "The following variables should all have same size:\n"
         << "length(y) for one block      : " << n1y << "\n"
         << "sensor_response_f.nelem()    : " << sensor_response_f.nelem()
         << "\nsensor_response_pol.nelem(): " << sensor_response_pol.nelem()
         << "\nsensor_response_za.nelem() : " << sensor_response_za.nelem() 
         << "\nsensor_response_aa.nelem() : " << sensor_response_za.nelem() 
         << "\n";
      throw runtime_error( os.str() );
    }


  //---------------------------------------------------------------------------
  // Allocations and resizing
  //---------------------------------------------------------------------------

  // Resize *y* and *y_XXX*
  //
  y.resize( nmblock*n1y );
  y_f.resize( nmblock*n1y );
  y_pol.resize( nmblock*n1y );
  y_pos.resize( nmblock*n1y, sensor_pos.ncols() );
  y_los.resize( nmblock*n1y, sensor_los.ncols() );
  y_error.resize( nmblock*n1y );
  y_aux.resize( nmblock*n1y );
  //
  y_error = 0;
  y_aux = 999;

  // Jacobian variables
  //
  Index  j_analytical_do = 0;
  //
  if( jacobian_do )
    {
      jacobian.resize( nmblock*n1y, 
                       jacobian_indices[jacobian_indices.nelem()-1][1]+1 );
      jacobian = 0;
      //
      FOR_ANALYTICAL_JACOBIANS_DO(
        j_analytical_do  = 1; 
      )
    }
  else
    { jacobian.resize(0, 0); }


  //---------------------------------------------------------------------------
  // The calculations
  //---------------------------------------------------------------------------

  // We have to make a local copy of the Workspace and the agendas because
  // only non-reference types can be declared firstprivate in OpenMP
  Workspace l_ws (ws);
  Agenda l_jacobian_agenda (jacobian_agenda);
  Agenda l_iy_clearsky_agenda (iy_clearsky_agenda);

/*#pragma omp parallel for                                           \
  if(!arts_omp_in_parallel() && nmblock>1 && nmblock>=nza)        \
    default(none)                                                   \
    firstprivate(l_ws, l_jacobian_agenda, l_iy_clearsky_agenda)     \
    shared(j_analytical_do, sensor_los, mblock_za_grid, mblock_aa_grid, \
           vmr_field, t_field, lon_grid, lat_grid, p_grid, f_grid,      \
           sensor_pos, joker, naa)*/
#pragma omp parallel for                                          \
  if(!arts_omp_in_parallel() && nmblock>1 && nmblock>=nza)        \
  firstprivate(l_ws, l_jacobian_agenda, l_iy_clearsky_agenda)
  for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
    {
      // Calculate monochromatic pencil beam data for 1 measurement block
      //
      Vector          iyb, iyb_aux, iyb_error, yb(n1y);
      Index           iy_error_type;
      ArrayOfMatrix   diyb_dx;      
      //
      iyb_calc( l_ws, iyb, iyb_error, iy_error_type, iyb_aux, diyb_dx, 
                mblock_index, atmosphere_dim, t_field, 
                z_field, vmr_field, cloudbox_on, stokes_dim, f_grid, sensor_pos,
                sensor_los, mblock_za_grid, mblock_aa_grid, antenna_dim, 
                l_iy_clearsky_agenda, y_unit, j_analytical_do, 
                jacobian_quantities, jacobian_indices, verbosity );


      // Apply sensor response matrix on iyb, and put into y
      //
      const Range rowind = get_rowindex_for_mblock( sensor_response, 
                                                                mblock_index );
      const Index row0   = rowind.get_start();
      //
      mult( yb, sensor_response, iyb );
      //
      y[rowind] = yb;  // *yb* also used below, as input to jacobian_agenda

      // Apply sensor response matrix on iyb_error, and put into y_error
      // (nothing to do if type is 0).
      if( iy_error_type == 1 )
        { 
          for( Index i=0; i<n1y; i++ )
            {
              const Index ithis = row0+i;
              y_error[ithis] = 0;
              for( Index icol=0; icol<niyb; icol++ )
                { 
                  y_error[ithis] += pow( 
                       sensor_response(i,icol)*iyb_error[icol], (Numeric)2.0 ); 
                }
              y_error[ithis] = sqrt( y_error[ithis] );              
            }
        }
      else if( iy_error_type == 2 )
        { mult( y_error[rowind], sensor_response, iyb_error ); }

      // Apply sensor response matrix on iyb_aux, and put into y_aux
      // (if filled)
      if( iyb_aux.nelem() )
        { mult( y_aux[rowind], sensor_response, iyb_aux ); }

      // Fill information variables
      //
      for( Index i=0; i<n1y; i++ )
        { 
          y_f[row0+i]         = sensor_response_f[i];
          y_pol[row0+i]       = sensor_response_pol[i]; 
          y_pos(row0+i,joker) = sensor_pos(mblock_index,joker);
          y_los(row0+i,0)     = sensor_los(mblock_index,0) + 
                                sensor_response_za[i];
          if( sensor_response_aa.nelem() )
            { 
              y_los(row0+i,1) = sensor_los(mblock_index,0) + 
                                sensor_response_aa[i]; 
            }
        }

      // Apply sensor response matrix on diyb_dx, and put into jacobian
      // (that is, analytical jacobian part)
      //
      if( j_analytical_do )
        {
          FOR_ANALYTICAL_JACOBIANS_DO(
            mult( jacobian(rowind, Range(jacobian_indices[iq][0],
                          jacobian_indices[iq][1]-jacobian_indices[iq][0]+1)),
                                                sensor_response, diyb_dx[iq] );
          )
        }

      // Rest of *jacobian*
      //
      if( jacobian_do )
        { 
          jacobian_agendaExecute( l_ws, jacobian, mblock_index, iyb, yb, 
                                                            l_jacobian_agenda );
        }
    }  // End mblock loop
}





/* Workspace method: Doxygen documentation will be auto-generated */
void yCalc2(
        Workspace&                  ws,
        Vector&                     y,
        Vector&                     y_f,
        ArrayOfIndex&               y_pol,
        Matrix&                     y_pos,
        Matrix&                     y_los,
        Vector&                     y_error,
        Vector&                     y_aux,
        Matrix&                     jacobian,
  const Index&                      basics_checked,
  const Index&                      atmosphere_dim,
  const Tensor3&                    t_field,
  const Tensor3&                    z_field,
  const Tensor4&                    vmr_field,
  const Index&                      cloudbox_on,
  const Index&                      cloudbox_checked,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const Matrix&                     sensor_pos,
  const Matrix&                     sensor_los,
  const Vector&                     mblock_za_grid,
  const Vector&                     mblock_aa_grid,
  const Index&                      antenna_dim,
  const Agenda&                     sensor_response_agenda,
  const Agenda&                     iy_clearsky_agenda,
  const String&                     y_unit,
  const Agenda&                     jacobian_agenda,
  const Index&                      jacobian_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices,
  const Verbosity&                  verbosity )
{
  // Some sizes
  const Index   nf      = f_grid.nelem();
  const Index   nza     = mblock_za_grid.nelem();
        Index   naa     = mblock_aa_grid.nelem();   
  if( antenna_dim == 1 )  
    { naa = 1; }
  const Index   nmblock = sensor_pos.nrows();
  const Index   niyb    = nf * nza * naa * stokes_dim;


  //---------------------------------------------------------------------------
  // Input checks
  //---------------------------------------------------------------------------

  // Basics and cloudbox OK?
  //
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control varaibles must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );

  // Sensor position and LOS.
  //
  if( sensor_pos.ncols() != atmosphere_dim )
    throw runtime_error( "The number of columns of sensor_pos must be "
                         "equal to the atmospheric dimensionality." );
  if( atmosphere_dim <= 2  &&  sensor_los.ncols() != 1 )
    throw runtime_error( "For 1D and 2D, sensor_los shall have one column." );
  if( atmosphere_dim == 3  &&  sensor_los.ncols() != 2 )
    throw runtime_error( "For 3D, sensor_los shall have two columns." );
  if( sensor_los.nrows() != nmblock )
    {
      ostringstream os;
      os << "The number of rows of sensor_pos and sensor_los must be "
         << "identical, but sensor_pos has " << nmblock << " rows,\n"
         << "while sensor_los has " << sensor_los.nrows() << " rows.";
      throw runtime_error( os.str() );
    }
  if( max( sensor_los(joker,0) ) > 180 )
    throw runtime_error( 
     "First column of *sensor_los* is not allowed to have values above 180." );
  if( atmosphere_dim == 2 )
    {
      if( min( sensor_los(joker,0) ) < -180 )
          throw runtime_error( "For atmosphere_dim = 2, first column of "
                    "*sensor_los* is not allowed to have values below -180." );
    }     
  else
    {
      if( min( sensor_los(joker,0)  ) < 0 )
          throw runtime_error( "For atmosphere_dim != 2, first column of "
                       "*sensor_los* is not allowed to have values below 0." );
    }    
  if( atmosphere_dim == 3  &&  max( sensor_los(joker,1) ) > 180 )
    throw runtime_error( 
    "Second column of *sensor_los* is not allowed to have values above 180." );

  // Antenna
  //
  chk_if_in_range( "antenna_dim", antenna_dim, 1, 2 );
  //
  if( nza == 0 )
    throw runtime_error( "The measurement block zenith angle grid is empty." );
  chk_if_increasing( "mblock_za_grid", mblock_za_grid );
  //
  if( antenna_dim == 1 )
    {
      if( mblock_aa_grid.nelem() != 0 )
        throw runtime_error( 
          "For antenna_dim = 1, the azimuthal angle grid must be empty." );
    }
  else
    {
      if( atmosphere_dim < 3 )
        throw runtime_error( "2D antennas (antenna_dim=2) can only be "
                                                 "used with 3D atmospheres." );
      if( mblock_aa_grid.nelem() == 0 )
        throw runtime_error(
                      "The measurement block azimuthal angle grid is empty." );
      chk_if_increasing( "mblock_aa_grid", mblock_aa_grid );
    }


  //---------------------------------------------------------------------------
  // Allocations 
  //---------------------------------------------------------------------------

  // Temporary storage for *y* etc.
  //
  ArrayOfVector ya(nmblock), yf(nmblock), yerror(nmblock), yaux(nmblock);
  ArrayOfMatrix ypos(nmblock), ylos(nmblock), ja(nmblock);
  ArrayOfArrayOfIndex ypol(nmblock);

  // Jacobian variables
  //
  Index  j_analytical_do = 0;
  //
  if( jacobian_do )
    {
      FOR_ANALYTICAL_JACOBIANS_DO(
        j_analytical_do  = 1; 
      )
    }


  //---------------------------------------------------------------------------
  // The calculations
  //---------------------------------------------------------------------------

  // We have to make a local copy of the Workspace and the agendas because
  // only non-reference types can be declared firstprivate in OpenMP
  Workspace l_ws (ws);
  Agenda l_jacobian_agenda (jacobian_agenda);
  Agenda l_iy_clearsky_agenda (iy_clearsky_agenda);
  Agenda l_sensor_response_agenda (sensor_response_agenda);

/*#pragma omp parallel for                                              \
  if(!arts_omp_in_parallel() && nmblock>1 && nmblock>=nza)              \
    default(none)                                                       \
    firstprivate(l_ws, l_jacobian_agenda, l_iy_clearsky_agenda,         \
                 l_sensor_response_agenda)                              \
    shared(j_analytical_do, sensor_los, mblock_za_grid, mblock_aa_grid, \
           vmr_field, t_field, lon_grid, lat_grid, p_grid, f_grid,      \
           sensor_pos, joker, naa)*/
#pragma omp parallel for                                                \
  if(!arts_omp_in_parallel() && nmblock>1 && nmblock>=nza)              \
  firstprivate(l_ws, l_jacobian_agenda, l_iy_clearsky_agenda)
  for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
    {
      // Calculate monochromatic pencil beam data for 1 measurement block
      //
      Vector          iyb, iyb_aux, iyb_error;
      Index           iy_error_type;
      ArrayOfMatrix   diyb_dx;      
      //
      iyb_calc( l_ws, iyb, iyb_error, iy_error_type, iyb_aux, diyb_dx, 
                mblock_index, atmosphere_dim, t_field, 
                z_field, vmr_field, cloudbox_on, stokes_dim, f_grid, sensor_pos,
                sensor_los, mblock_za_grid, mblock_aa_grid, antenna_dim, 
                l_iy_clearsky_agenda, y_unit, j_analytical_do, 
                jacobian_quantities, jacobian_indices, verbosity );

      // Obtain sensor response data
      //
      Sparse         sensor_response;
      Vector         sensor_response_f;
      ArrayOfIndex   sensor_response_pol;
      Vector         sensor_response_za;
      Vector         sensor_response_aa;
      //
      sensor_response_agendaExecute( l_ws, sensor_response, sensor_response_f, 
                                     sensor_response_pol, sensor_response_za,
                                     sensor_response_aa, 
                                     mblock_index, l_sensor_response_agenda );

      if( sensor_response.ncols() != niyb ) 
        {
          ostringstream os;
          os << "The number of columns in *sensor_response* and the length\n"
             << "of *iyb* do not match.";
          throw runtime_error( os.str() );
        }
      //
      const Index yl = sensor_response.nrows(); // Length of this part of y
      //
      if( yl!=sensor_response_f.nelem() || yl!=sensor_response_pol.nelem() ||
          yl!= sensor_response_za.nelem() || yl!=sensor_response_za.nelem() )
        {
          ostringstream os;
          os << "Sensor auxiliary variables do not have the correct size.\n"
             << "The following variables should all have same size:\n"
             << "length(y) for one block      : " << yl << "\n"
             << "sensor_response_f.nelem()    : " << sensor_response_f.nelem()
             << "\nsensor_response_pol.nelem(): " << sensor_response_pol.nelem()
             << "\nsensor_response_za.nelem() : " << sensor_response_za.nelem() 
             << "\nsensor_response_aa.nelem() : " << sensor_response_za.nelem() 
             << "\n";
          throw runtime_error( os.str() );
        }

      // Set sizes of ya and associated variables
      //
      ya[mblock_index].resize( yl );
      yf[mblock_index].resize( yl );
      ypol[mblock_index].resize( yl );
      ypos[mblock_index].resize( yl, sensor_pos.ncols() );
      ylos[mblock_index].resize( yl, sensor_los.ncols() );
      yerror[mblock_index].resize( yl );
      yaux[mblock_index].resize( yl );
      if( jacobian_do )
        { ja[mblock_index].resize( yl, 
                         jacobian_indices[jacobian_indices.nelem()-1][1]+1 ); }

      // Apply sensor response matrix on iyb, and put into ya
      //
      mult( ya[mblock_index], sensor_response, iyb );

      // Apply sensor response matrix on iyb_error, and put into yerror
      //
      if( iy_error_type == 0 )
        { yerror[mblock_index] = 0; }
      else if( iy_error_type == 1 )
        { 
          for( Index i=0; i<yl; i++ )
            {
              yerror[mblock_index][i] = 0;
              for( Index icol=0; icol<niyb; icol++ )
                { 
                  yerror[mblock_index][i] += pow( 
                       sensor_response(i,icol)*iyb_error[icol], (Numeric)2.0 ); 
                }
              yerror[mblock_index][i] = sqrt( yerror[mblock_index][i] );              
            }
        }
      else if( iy_error_type == 2 )
        { mult( yerror[mblock_index], sensor_response, iyb_error ); }

      // Apply sensor response matrix on iyb_aux, and put into y_aux
      // (if filled)
      if( iyb_aux.nelem() )
        { mult( yaux[mblock_index], sensor_response, iyb_aux ); }
      else
        { yaux[mblock_index] = 999; }

      // Fill information variables
      //
      for( Index i=0; i<yl; i++ )
        { 
          yf[mblock_index][i]         = sensor_response_f[i];
          ypol[mblock_index][i]       = sensor_response_pol[i]; 
          ypos[mblock_index](i,joker) = sensor_pos(mblock_index,joker);
          ylos[mblock_index](i,0)     = sensor_los(mblock_index,0) + 
                                        sensor_response_za[i];
          if( sensor_response_aa.nelem() )
            { 
              ylos[mblock_index](i,1) = sensor_los(mblock_index,0) + 
                                        sensor_response_aa[i];
            }
        }

      // Apply sensor response matrix on diyb_dx, and put into jacobian
      // (that is, analytical jacobian part)
      //
      if( j_analytical_do )
        {
          FOR_ANALYTICAL_JACOBIANS_DO(
            mult( ja[mblock_index](joker, Range(jacobian_indices[iq][0],
                          jacobian_indices[iq][1]-jacobian_indices[iq][0]+1)),
                                                sensor_response, diyb_dx[iq] );
          )
        }

      // Rest of *jacobian*
      //
      if( jacobian_do )
        { 
          jacobian_agendaExecute( l_ws, ja[mblock_index], mblock_index, iyb, 
                                  ya[mblock_index], l_jacobian_agenda );
        }
    }  // End mblock loop


  // Combine data from the mblocks to form output variables
  //
  Index i0 =0 , n = 0;
  //
  for( Index m=0; m<nmblock; m++ )
    { n += ya[m].nelem(); }
  // 
  y.resize( n );
  y_f.resize( n );
  y_pol.resize( n );
  y_pos.resize( n, sensor_pos.ncols() );
  y_los.resize( n, sensor_los.ncols() );
  y_error.resize( n );
  y_aux.resize( n ); 
  if( jacobian_do )
    { jacobian.resize( n, jacobian_indices[jacobian_indices.nelem()-1][1]+1 );}
  //
  for( Index m=0; m<nmblock; m++ )
    {
      for( Index i=0; i<ya[m].nelem(); i++ )
        {
          const Index r = i0 + i;
          y[r]              = ya[m][i];
          y_f[r]            = yf[m][i];
          y_pol[r]          = ypol[m][i];
          y_pos(r,joker)    = ypos[m](i,joker);
          y_los(r,joker)    = ylos[m](i,joker);
          y_error[r]        = yerror[m][i];
          y_aux[r]          = yaux[m][i];
          if( jacobian_do )
            { jacobian(r,joker) = ja[m](i,joker); }
        }
      i0 += ya[m].nelem();
    }  
}





/* Workspace method: Doxygen documentation will be auto-generated */
void yFromIy(
        Vector&         y,
        Vector&         y_f,
        ArrayOfIndex&   y_pol,
        Matrix&         y_pos,
        Matrix&         y_los,
        Vector&         y_error,
        Vector&         y_aux,
        Matrix&         jacobian,
  const Index&          stokes_dim,
  const Vector&         f_grid,
  const Index&                 jacobian_do,
  const ArrayOfArrayOfIndex&   jacobian_indices,
  const Vector&         rte_pos,
  const Vector&         rte_los,
  const Matrix&         iy,
  const Matrix&         iy_aux,
  const Matrix&         iy_error,
  const Index&          iy_error_type,
  const ArrayOfTensor3& diy_dx,
  const Verbosity& )
{
  const Index nf = f_grid.nelem();

  if( iy.nrows() != nf  ||  iy.ncols() != stokes_dim )
    throw runtime_error( "The size of *iy* does not match length of *f_grid* "
                         "and/or *stokes_dim*." );

  // Size and init y-vars
  const Index n = nf*stokes_dim;
  y.resize( n );
  y_f.resize( n );
  y_aux.resize( n );     y_aux = 999;
  y_error.resize( n );   y_error = 0;
  y_pol.resize( n );
  y_pos.resize( n, rte_pos.nelem() );
  y_los.resize( n, rte_los.nelem() );

  // Size jacobian
  Index nq = 0;
  if( jacobian_do )
    {
      nq = jacobian_indices.nelem();
      jacobian.resize( n, jacobian_indices[nq-1][1]+1 );
      jacobian = 0;
    }
  else
    { jacobian.resize( 0, 0); }

  for( Index i=0; i <nf; i++ )
    {
      const Index i0 = i*stokes_dim;

      for( Index j=0; j<stokes_dim; j++ )
        {
          Index ii = i0 + j;
          
          y[ii]       = iy(i,j);
          y_f[ii]     = f_grid[i];
          y_pol[ii]   = j+1;
          if( iy_error_type )  // Nothing to do if 0
            { y_error[ii] = iy_error(i,j); }
          if( iy_aux.nrows() )
            { y_aux[ii]   = iy_aux(i,j); }

          y_pos(ii,joker) = rte_pos;
          y_los(ii,joker) = rte_los;

          if( jacobian_do )
            {
              for( Index q=0; q<nq; q++ )
                {
                  if( diy_dx[q].npages() )
                    {
                      for( Index r=0; r<diy_dx[q].npages(); r++ )
                        {
                          jacobian(ii,jacobian_indices[q][0]+r) = 
                                                              diy_dx[q](r,i,j);
                        }
                    }
                }
            }
        }
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void y_unitApply(
        Vector&         y,
        Vector&         y_error,
        Matrix&         jacobian,
  const Vector&         y_f,
  const ArrayOfIndex&   y_pol,
  const String&         y_unit,
  const Verbosity&)
{
  if( y_unit == "1" )
    { throw runtime_error(
        "No need to use this method with *y_unit* = \"1\"." ); }

  if( max(y) > 1e-3 )
    {
      ostringstream os;
      os << "The spectrum vector *y* is required to have original radiance\n"
         << "unit, but this seems not to be the case. This as a value above\n"
         << "1e-3 is found in *y*.";
      throw runtime_error( os.str() );      
    }

  // Are jacobian and y_error set?
  //
  const Index ny = y.nelem();
  //
  const bool do_j = jacobian.nrows() == ny;
  const bool do_e = y_error.nelem() == ny;

  // Some jacobian quantities can not be handled
  if( do_j  &&  max(jacobian) > 1e-3 )
    {
      ostringstream os;
      os << "The method can not be used with jacobian quantities that are not\n"
         << "obtained through radiative transfer calculations. One example on\n"
         << "quantity that can not be handled is *jacobianAddPolyfit*.\n"
         << "The maximum value of *jacobian* indicates that one or several\n"
         << "such jacobian quantities are included.";
      throw runtime_error( os.str() );      
    }

  // All conversions could be handled as Planck-BT, but before everything is
  // tested properly and it is clear which version that is most efficient, both
  // versions are kept

  // Planck-Tb 
  //-------------------------------------------------------------------------- 
  if( y_unit == "PlanckBT" )
    {
      // Hard to use telescoping here as the data are sorted differently in y
      // and jacobian, than what is expected apply_y_unit. Copy to temporary
      // variables instead.

      // Handle the elements in "frequency chunks"

      Index i0 = 0;           // Index of first element for present chunk
      //
      while( i0 < ny )
        { 
          // Find number values for this chunk
          Index n = 1;
          //
          while( i0+n < ny &&  y_f[i0] == y_f[i0+n] ) 
            { n++; }                              

          Matrix yv(1,n);  
          ArrayOfIndex i_pol(n);
          bool any_quv = false;
          //
          for( Index i=0; i<n; i++ )
            { 
              const Index ix=i0+i;   
              yv(0,i) = y[ix];   
              i_pol[i] = y_pol[ix]; 
              if( i_pol[i] > 1  &&  i_pol[i] < 5 )
                { any_quv = true; }
            }

          // Index of elements to convert
          Range ii( i0, n );

          if( do_e || do_j )
            {
              if( any_quv  &&  i_pol[0] != 1 )
                {
                  ostringstream os;
                  os << "The conversion to PlanckBT, of the Jacobian and "
                     << "errors for Q, U and V, requires that I (first Stokes "
                     << "element) is at hand and that the data are sorted in "
                     << "such way that I comes first for each frequency.";
                  throw runtime_error( os.str() );      
                }

              // y_error
              if( do_e )
                {
                  Tensor3 e(1,1,n);
                  e(0,0,joker) = y_error[ii];
                  apply_y_unit2( e, yv, y_unit, y_f[i0], i_pol ); 
                  y_error[ii] = e(0,0,joker);
                }

              // Jacobian
              if( do_j )
                {
                  Tensor3 J(jacobian.ncols(),1,n);
                  J(joker,0,joker) = transpose( jacobian(ii,joker) );
                  apply_y_unit2( J, yv, y_unit, y_f[i0], i_pol ); 
                  jacobian(ii,joker) = transpose( J(joker,0,joker) );
                }
            }

          // y (must be done last)
          apply_y_unit( yv, y_unit, y_f[i0], i_pol );
          y[ii] = yv(0,joker);

          i0 += n;
        }
    }


  // Other conversions
  //-------------------------------------------------------------------------- 
  else
    {
      // Here we take each element of y separately. 

      Matrix yv(1,1);
      ArrayOfIndex i_pol(1);

      for( Index i=0; i<ny; i++ )
        {
          yv(0,0)  = y[i];      // To avoid repeated telescoping
          i_pol[0] = y_pol[i];

          // y_error
          if( do_e )
            {
              apply_y_unit2( MatrixView(VectorView( y_error[i] )), yv, 
                             y_unit, y_f[i], i_pol ); 
            }

          // Jacobian
          if( do_j )
            {
              apply_y_unit2( MatrixView( jacobian(i,joker) ), yv, 
                             y_unit, y_f[i], i_pol ); 
            }

          // y (must be done last)
          apply_y_unit( yv, y_unit, y_f[i], i_pol );
          y[i] = yv(0,0);
        }
    }
}






