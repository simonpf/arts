/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>

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
  ===  File description 
  ===========================================================================*/

/*!
  \file   ppath.h
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-02
  
  \brief  Propagation path structure and functions.
  
   This file contains the definition of the Ppath structure and the
   functions in ppath.cc that are of interest elsewhere.
*/


#ifndef ppath_h
#define ppath_h


#include "agenda_class.h"
#include "array.h"
#include "arts.h"
#include "interpolation.h"
#include "matpackI.h"
#include "mystring.h"



/*===========================================================================
  === The Ppath structure
  ===========================================================================*/

//! The structure to describe a propagation path and releated quantities.
/*! 
   The fields of the structure are described in the ARTS user guide (AUG).
   It is listed as a sub-entry to "data structures".  
*/

struct Ppath {
  Index             dim;
  Index             np;
  Index             refraction;
  String            method;
  Numeric           constant;
  Matrix            pos;
  Vector            z;
  Vector            l_step;
  ArrayOfGridPos    gp_p;
  ArrayOfGridPos    gp_lat;
  ArrayOfGridPos    gp_lon;
  Matrix            los;
  String            background;
  Vector            tan_pos;
  Vector            geom_tan_pos;
};



/*===========================================================================
  === Functions from ppath.cc
  ===========================================================================*/

double geometrical_ppc( const double& r, const double& za );

double geompath_za_at_r(
        const double&   ppc,
        const double&   a_za,
        const double&   r );

double geompath_lat_at_za(
        const double&   za0,
        const double&   lat0,
        const double&   za );

bool is_los_downwards( 
        const double&   za,
        const double&   tilt );

double psurface_slope_2d(
        ConstVectorView   lat_grid,           
        ConstVectorView   r_geoid,
        ConstVectorView   z_surf,
        const GridPos&    gp,
        const double&     za );

double psurface_slope_3d(
        const double&   lat1,
        const double&   lat3,
        const double&   lon5,
        const double&   lon6,
        const double&   r15,
        const double&   r35,
        const double&   r36,
        const double&   r16,
        const double&   lat,
        const double&   lon,
        const double&   aa );

double psurface_angletilt(
        const double&   r,
        const double&   c );

void ppath_init_structure( 
              Ppath&      ppath,
        const Index&      atmosphere_dim,
        const Index&      np );

void ppath_set_background( 
              Ppath&      ppath,
        const Index&      case_nr );

Index ppath_what_background( const Ppath&   ppath );

void ppath_step_geom_1d(
              Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   z_field,
        const double&     r_geoid,
        const double&     z_ground,
        const double&     lmax );

void ppath_step_geom_2d(
              Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstVectorView   r_geoid,
        ConstVectorView   z_ground,
        const double&     lmax );

void ppath_step_geom_3d(
              Ppath&       ppath,
        ConstVectorView    p_grid,
        ConstVectorView    lat_grid,
        ConstVectorView    lon_grid,
        ConstTensor3View   z_field,
        ConstMatrixView    r_geoid,
        ConstMatrixView    z_ground,
	const double&      lmax );

void ppath_step_refr_1d(
              Ppath&      ppath,
              Numeric&    a_pressure,
              Numeric&    a_temperature,
              Vector&     a_vmr_list,
              Numeric&    refr_index,
        const Agenda&     refr_index_agenda,
        ConstVectorView   p_grid,
        ConstVectorView   z_field,
        ConstVectorView   t_field,
        ConstMatrixView   vmr_field,
        const double&     r_geoid,
        const double&     z_ground,
        const String&     rtrace_method,
        const double&     lraytrace,
        const double&     lmax );

void ppath_step_refr_2d(
              Ppath&      ppath,
              Numeric&    a_pressure,
              Numeric&    a_temperature,
              Vector&     a_vmr_list,
              Numeric&    refr_index,
        const Agenda&     refr_index_agenda,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstMatrixView   t_field,
        ConstTensor3View  vmr_field,
        ConstVectorView   r_geoid,
        ConstVectorView   z_ground,
        const String&     rtrace_method,
        const double&     lraytrace,
        const double&     lmax );

void ppath_step_refr_3d(
              Ppath&      ppath,
              Numeric&    a_pressure,
              Numeric&    a_temperature,
              Vector&     a_vmr_list,
              Numeric&    refr_index,
        const Agenda&     refr_index_agenda,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid,
        ConstTensor3View  z_field,
        ConstTensor3View  t_field,
        ConstTensor4View  vmr_field,
        ConstMatrixView   r_geoid,
        ConstMatrixView   z_ground,
        const String&     rtrace_method,
        const double&     lraytrace,
        const double&     lmax );

void ppath_calc(
              Ppath&          ppath,
              Ppath&          ppath_step,
        const Agenda&         ppath_step_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Matrix&         r_geoid,
        const Matrix&         z_ground,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         a_pos,
        const Vector&         a_los,
        const bool&           outside_cloudbox,
        const Index&          agenda_verb );

#endif  // ppath_h
