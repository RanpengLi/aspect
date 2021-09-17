/*
  Copyright (C) 2019 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_material_model_rheology_friction_options_h
#define _aspect_material_model_rheology_friction_options_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parsed_function.h>
#include <aspect/utilities.h>

#include <aspect/material_model/rheology/drucker_prager.h>
#include <deal.II/fe/component_mask.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;
    namespace Rheology
    {
      /**
       * Enumeration for selecting which type of friction dependence to use.
       *
       * For the type 'independent', the user-supplied internal angle of friction is used.
       *
       * For the type 'dynamic friction', the friction angle is rate dependent following
       * Equation 13 from \\cite{van_dinther_seismic_2013}.
       */
      enum FrictionDependenceMechanism
      {
        independent,
        dynamic_friction,
      };

      template <int dim>
      class FrictionOptions : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm);

          /**
           * A function that computes the new friction angle when rate and/or state
           * dependence is taken into account. Given a volume fraction with
           * the index j and a vector of all compositional fields, it returns
           * the newly calculated friction angle.
           */
          double
          compute_dependent_friction_angle(const double current_edot_ii,
                                           const unsigned int j,
                                           double current_friction) const;

          /**
           * A function that returns the selected type of friction dependence.
           */
          FrictionDependenceMechanism
          get_friction_dependence_mechanism () const;
          
        private:
          /**
           * Select the mechanism to be used for the friction dependence.
           * Possible options: none | dynamic friction | rate and state dependent friction
           * | rate and state dependent friction plus linear slip weakening | slip rate
           * dependent rate and state dependent friction | regularized rate and state
           * dependent friction
           */
          FrictionDependenceMechanism friction_dependence_mechanism;

          /**
           * Dynamic friction input parameters
           */

          /**
           * Dynamic angles of internal friction that are used at high strain rates.
           */
          std::vector<double> dynamic_angles_of_internal_friction;

          /**
           * The characteristic strain rate value at which the angle of friction is taken as
           * the mean of the dynamic and the static angle of friction. When the effective
           * strain rate in a cell is very high the dynamic angle of friction is taken, when
           * it is very low the static angle of internal friction is chosen.
           */
          double dynamic_characteristic_strain_rate;

          /**
           * An exponential factor in the equation for the calculation of the friction angle
           * to make the transition between static and dynamic friction angle more smooth or
           * more step-like.
           */
          double dynamic_friction_smoothness_exponent;
      };
    }
  }
}
#endif
