/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_compositing_h
#define _aspect_material_model_compositing_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

#include <map>
#include <vector>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Property
    {
      /**
       * An enum to define what material properties are available.
       */
      enum MaterialProperty
      {
        viscosity,
        density,
        thermal_expansion_coefficient,
        specific_heat,
        thermal_conductivity,
        compressibility,
        entropy_derivative_pressure,
        entropy_derivative_temperature,
        reaction_terms
      };
    }

    /**
     * A material model that selects properties from other
     * (non-compositing) material models.
     *
     * @ingroup MaterialModels
     */

    template <int dim>
    class Compositing : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * @copydoc MaterialModel::Interface::initialize()
         */
        void
        initialize () override;

        /**
         * @copydoc MaterialModel::Interface::evaluate()
         */
        void
        evaluate (const typename Interface<dim>::MaterialModelInputs &in,
                  typename Interface<dim>::MaterialModelOutputs &out) const override;

        /**
         * If this material model can produce additional named outputs
         * that are derived from NamedAdditionalOutputs, create them in here.
         */
        void
        create_additional_named_outputs (typename Interface<dim>::MaterialModelOutputs &outputs) const override;

        /**
         * @copydoc MaterialModel::Interface::declare_parameters()
         */
        static void
        declare_parameters (ParameterHandler &prm);

        /**
         * @copydoc MaterialModel::Interface::parse_parameters()
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * @copydoc MaterialModel::Interface::is_compressible()
         *
         * Returns value from material model providing compressibility.
         */
        bool is_compressible () const override;



      private:
        /**
         * Copy desired properties from material model outputs
         * produced by another material model,
         *
         * @param model_index Internal index for pointer to the material model evaluated
         * @param base_output Properties generated by the material model specified
         * @param out MaterialModelOutputs to be used.
         */
        void
        copy_required_properties(const unsigned int model_index,
                                 const typename Interface<dim>::MaterialModelOutputs &base_output,
                                 typename Interface<dim>::MaterialModelOutputs &out) const;

        /**
         * Map specifying which material model is responsible for a
         * property.
         */
        std::map<Property::MaterialProperty, unsigned int> model_property_map;

        /**
         * Names of and pointers to the material models used for
         * compositing.
         */
        std::vector<std::string>                             model_names;
        std::vector<std::unique_ptr<Interface<dim>>> models;
    };
  }
}

#endif
