/*
  Copyright (C) 2016 - 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/entropy_model.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>
#include <aspect/material_model/rheology/visco_plastic.h>
#include <aspect/material_model/steinberger.h>
#include <aspect/material_model/equation_of_state/interface.h>

namespace aspect //TEST
{
  namespace MaterialModel
  {
    namespace
    {
      template <int dim>
      bool solver_scheme_is_supported(const Parameters<dim> &parameters)
      {
        // If we solve advection equations, we need to iterate them, because this material
        // models splits temperature diffusion from entropy advection.
        switch (parameters.nonlinear_solver)
          {
            case Parameters<dim>::NonlinearSolver::Kind::iterated_Advection_and_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::iterated_Advection_and_defect_correction_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::iterated_Advection_and_Newton_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::no_Advection_no_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::no_Advection_iterated_defect_correction_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::no_Advection_iterated_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::no_Advection_single_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::first_timestep_only_single_Stokes:
              return true;

            case Parameters<dim>::NonlinearSolver::Kind::single_Advection_single_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::single_Advection_iterated_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::single_Advection_iterated_defect_correction_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::single_Advection_iterated_Newton_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::single_Advection_no_Stokes:
              return false;
          }
        Assert(false, ExcNotImplemented());
        return false;
      }
    }



    template <int dim>
    void
    EntropyModel<dim>::initialize()
    {
      AssertThrow (this->get_parameters().formulation_mass_conservation ==
                   Parameters<dim>::Formulation::MassConservation::projected_density_field,
                   ExcMessage("The 'entropy model' material model was only tested with the "
                              "'projected density field' approximation "
                              "for the mass conservation equation, which is not selected."));

      AssertThrow (this->introspection().composition_type_exists(CompositionalFieldDescription::Type::entropy),
                   ExcMessage("The 'entropy model' material model requires the existence of a compositional field "
                              "named 'entropy'. This field does not exist."));

      AssertThrow(solver_scheme_is_supported(this->get_parameters()) == true,
                  ExcMessage("The 'entropy model' material model requires the use of a solver scheme that "
                             "iterates over the advection equations but a non iterating solver scheme was selected. "
                             "Please check the consistency of your solver scheme."));

// Uncomment this line to enable the use of the entropy averaging for multiple compositions,
// which has NOT been tested yet.
      /*      AssertThrow(material_file_names.size() == 1 || SimulatorAccess<dim>::get_end_time () == 0,
                       ExcMessage("The 'entropy model' material model can only handle one composition, "
                                  "and can therefore only read one material lookup table."));
      */


      for (unsigned int i = 0; i < material_file_names.size(); ++i)
        {
          entropy_reader.push_back(std::make_unique<MaterialUtilities::Lookup::EntropyReader>());
          entropy_reader[i]->initialize(this->get_mpi_communicator(), data_directory, material_file_names[i]);
        }

      lateral_viscosity_prefactor_lookup = std::make_unique<internal::LateralViscosityLookup>(data_directory+lateral_viscosity_file_name,
                                           this->get_mpi_communicator());
    }



    template <int dim>
    bool
    EntropyModel<dim>::
    is_compressible () const
    {
      return true;
    }



    template <int dim>
    double
    EntropyModel<dim>::
    thermal_conductivity (const double temperature,
                          const double pressure,
                          const Point<dim> &position) const
    {
      if (conductivity_formulation == constant)
        return thermal_conductivity_value;

      else if (conductivity_formulation == p_T_dependent)
        {
          // Find the conductivity layer that corresponds to the depth of the evaluation point.
          const double depth = this->get_geometry_model().depth(position);
          unsigned int layer_index = std::distance(conductivity_transition_depths.begin(),
                                                   std::lower_bound(conductivity_transition_depths.begin(),conductivity_transition_depths.end(), depth));

          const double p_dependence = reference_thermal_conductivities[layer_index] + conductivity_pressure_dependencies[layer_index] * pressure;

          // Make reasonably sure we will not compute any invalid values due to the temperature-dependence.
          // Since both the temperature-dependence and the saturation time scale with (Tref/T), we have to
          // make sure we can compute the square of this number. If the temperature is small enough to
          // be close to yielding NaN values, the conductivity will be set to the maximum value anyway.
          const double T = std::max(temperature, std::sqrt(std::numeric_limits<double>::min()) * conductivity_reference_temperatures[layer_index]);
          const double T_dependence = std::pow(conductivity_reference_temperatures[layer_index] / T, conductivity_exponents[layer_index]);

          // Function based on the theory of Roufosse and Klemens (1974) that accounts for saturation.
          // For the Tosi formulation, the scaling should be zero so that this term is 1.
          double saturation_function = 1.0;
          if (1./T_dependence > 1.)
            saturation_function = (1. - saturation_scaling[layer_index])
                                  + saturation_scaling[layer_index] * (2./3. * std::sqrt(T_dependence) + 1./3. * 1./T_dependence);

          return std::min(p_dependence * saturation_function * T_dependence, maximum_conductivity);
        }
      else
        {
          AssertThrow(false, ExcNotImplemented());
          return numbers::signaling_nan<double>();
        }
    }



    template <int dim>
    double
    EntropyModel<dim>::
    equilibrate_temperature (std::vector<double> &composition_equalibrated_S,
                             const std::vector<double> &temperature,
                             const std::vector<double> &mass_fractions,
                             const std::vector<double> &entropy,
                             const std::vector<double> &Cp,
                             const double pressure) const
    {
      AssertThrow(material_file_names.size() == temperature.size() && temperature.size() == mass_fractions.size() &&  temperature.size() == entropy.size() && temperature.size() == Cp.size(),
                  ExcMessage("The temperature, chemical composition, entropy and specific heat capacity vectors"
                             " must all have the same size as the number of look-up tables."));

      std::vector<double> composition_initial_S  = entropy;
      std::vector<double> composition_initial_T  = temperature;
      std::vector<double> composition_initial_Cp = Cp;
      std::vector<double> composition_lookup_T(temperature.size());
      std::vector<double> composition_lookup_Cp(temperature.size());

      bool equalibration = false;
      unsigned int iteration = 0;
      double ln_equalibrated_T = 0;

      // Step1

      // TODO: set the iteration number as a parameter
      while (equalibration == false || iteration == 500)
        {
          double T_numerator = 0;
          double T_denominator = 0;

          iteration += 1;

          for (unsigned int i = 0; i < material_file_names.size(); ++i)
            {
              T_numerator += mass_fractions[i] * composition_initial_Cp[i] * std::log(composition_initial_T[i]);
              T_denominator += mass_fractions[i] * composition_initial_Cp[i];
            }

          ln_equalibrated_T = T_numerator/T_denominator;

          // step2
          for (unsigned int i = 0; i < material_file_names.size(); ++i)
            {
              composition_equalibrated_S[i] = composition_initial_S[i] + composition_initial_Cp[i] * (ln_equalibrated_T - std::log (composition_initial_T[i]));
              // step3
              composition_lookup_T[i] = entropy_reader[i]->temperature(composition_equalibrated_S[i], pressure);

              // composition_lookup_Cp[i] = entropy_reader[i]->specific_heat(composition_equalibrated_S[i], pressure);
            }
          // step4
          // update the T0 and S0 to prepare for another iteration
          composition_initial_T = composition_lookup_T;
          composition_initial_S = composition_equalibrated_S;
          //   composition_initial_Cp = composition_lookup_Cp;

          equalibration = true;
          for (unsigned int i = 0; i < material_file_names.size(); ++i)
            {
              // TODO: set the small value (currently 1e-5) as a parameter
              if (std::abs (composition_lookup_T[i] - std::exp(ln_equalibrated_T)) >= 1e-8)
                {
                  equalibration = false;
                  break;
                }
            }
        }
//(0920)      std::cout << "S for component = " << composition_equalibrated_S[0] <<" "<<composition_equalibrated_S[1] <<std::endl;
      //  entropy = composition_equalibrated_S;
      return exp(ln_equalibrated_T); // vector composition_equalibrated_S could be modified while reading in reference

    }

    template <int dim>
    void
    EntropyModel<dim>::evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                                MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      const unsigned int projected_density_index = this->introspection().compositional_index_for_name("density_field");
      //TODO : need to make it work for more than one field
      const std::vector<unsigned int> &entropy_indices = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::entropy);
      const std::vector<unsigned int> &composition_indices = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::chemical_composition);

      AssertThrow(composition_indices.size() == material_file_names.size() - 1,
                  ExcMessage("The 'entropy model' material model assumes that there exists a background field in addition to the compositional fields, "
                             "and therefore it requires one more lookup table than there are chemical compositional fields."));

      EquationOfStateOutputs<dim> eos_outputs (material_file_names.size());
      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();

      std::vector<double> volume_fractions (material_file_names.size());
      std::vector<double> mass_fractions (material_file_names.size());

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          // Use the adiabatic pressure instead of the real one,
          // to stabilize against pressure oscillations in phase transitions.
          // This is a requirement of the projected density approximation for
          // the Stokes equation and not related to the entropy formulation.
          // Also convert pressure from Pa to bar, bar is used in the table.
          //const double entropy = in.composition[i][entropy_index];
          std::vector<double> component_entropy (material_file_names.size());
          std::vector<double> composition_temperature_lookup (material_file_names.size()); // NEED TO CHANGE
          const double pressure = this->get_adiabatic_conditions().pressure(in.position[i]) / 1.e5;

          // Loop over all material files, and store the looked-up values for all compositions.
          for (unsigned int j=0; j<material_file_names.size(); ++j)
            {
              component_entropy[j] = in.composition[i][entropy_indices[j]];
              composition_temperature_lookup[j] = entropy_reader[j]->temperature(component_entropy[j], pressure);
              // std::cout << "component_entropy = " <<component_entropy[j]<<" " << std::endl;
              eos_outputs.densities[j] = entropy_reader[j]->density(component_entropy[j], pressure);
              // std::cout << "densities = " << eos_outputs.densities[j]<<" " << std::endl;
              eos_outputs.thermal_expansion_coefficients[j] = entropy_reader[j]->thermal_expansivity(component_entropy[j],pressure);
              eos_outputs.specific_heat_capacities[j] = entropy_reader[j]->specific_heat(component_entropy[j],pressure);

              const Tensor<1, 2> pressure_unit_vector({0.0, 1.0});
              eos_outputs.compressibilities[j] = ((entropy_reader[j]->density_gradient(component_entropy[j],pressure)) * pressure_unit_vector) / eos_outputs.densities[j];
            }



          // Calculate volume fractions from mass fractions
          // If there is only one lookup table, set the mass and volume fractions to 1
          if (material_file_names.size() == 1)
            mass_fractions [0] = 1.0;

          else
            {
              // We only want to compute mass/volume fractions for fields that are chemical compositions.
              mass_fractions = MaterialUtilities::compute_only_composition_fractions(in.composition[i], this->introspection().chemical_composition_field_indices());
            }

          volume_fractions = MaterialUtilities::compute_volumes_from_masses(mass_fractions,
                                                                            eos_outputs.densities,
                                                                            true);

          out.densities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);
          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.thermal_expansion_coefficients, MaterialUtilities::arithmetic);

          out.specific_heat[i] = MaterialUtilities::average_value (mass_fractions, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic);

          out.compressibilities[i] = MaterialUtilities::average_value (mass_fractions, eos_outputs.compressibilities, MaterialUtilities::arithmetic);


          // Thermal conductivity can be pressure temperature dependent
          std::vector<double> composition_equalibrated_S(material_file_names.size());



          const double equilibrated_T = equilibrate_temperature (composition_equalibrated_S, composition_temperature_lookup, mass_fractions, component_entropy, eos_outputs.specific_heat_capacities, pressure);
          const double temperature_lookup = equilibrated_T;

          /*
          std::cout << "equilibrated_T = " << equilibrated_T<<" " << std::endl;
          std::cout << "equilibrated_S = " << composition_equalibrated_S[0]<<" " <<composition_equalibrated_S[1]<<" " << std::endl;
          std::cout << "equilibrated_T_lookup = " << composition_temperature_lookup[0]<<" " << composition_temperature_lookup[1]<<" " <<std::endl;
          */
          //entropy_reader[0]->temperature(entropy,pressure);
          //   const std::vector<double> temp_temperature_lookup (material_file_names.size(), temperature_lookup); // NEED TO CHANGE


          out.thermal_conductivities[i] = thermal_conductivity(temperature_lookup, in.pressure[i], in.position[i]);

          out.entropy_derivative_pressure[i]    = 0.;
          out.entropy_derivative_temperature[i] = 0.;

          // Calculate the reaction terms

          if (material_file_names.size()==1)
            {
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                {
                  out.reaction_terms[i][c] = 0.;
                }
            }

          else
            {
              // Calculate the reaction rates for the operator splitting
              for (unsigned int c = 0; c < in.composition[i].size(); ++c)
                {
                  if (this->get_parameters().use_operator_splitting)
                    {
                      if (reaction_rate_out != nullptr)
                        {
                          //AssertThrow(this->get_parameters().use_operator_splitting == 1,
                          //ExcMessage("The 'entropy model' material model requires the use of operator splitting for multiple chemical composition."));



                          reaction_rate_out->reaction_rates[i][c] = 0.0;



                          for (unsigned int c = 0; c < in.composition[i].size(); ++c)
                            {
                              bool c_is_entropy_field = false;
                              unsigned int c_is_nth_entropy_field = 0;

                              unsigned int nth_entropy_index = 0;
                              for (unsigned int entropy_index : entropy_indices)
                                {
                                  if (c == entropy_index)
                                    {
                                      c_is_entropy_field = true;

                                      c_is_nth_entropy_field = nth_entropy_index;

                                    }
                                  ++nth_entropy_index;
                                }

//                          out.reaction_terms[i][c] = (composition_equalibrated_S[c_is_nth_entropy_field] - in.composition[i][entropy_indices[c_is_nth_entropy_field]]); //


                              const unsigned int timestep_number = this->simulator_is_past_initialization()
                                                                   ?
                                                                   this->get_timestep_number()
                                                                   :
                                                                   0;

                              if (c_is_entropy_field == true && timestep_number > 0)
                                //      const unsigned int dif = (composition_equalibrated_S[c_is_nth_entropy_field] - in.composition[i][entropy_indices[c_is_nth_entropy_field]])/ this->get_timestep();
                                reaction_rate_out->reaction_rates[i][c] = (composition_equalibrated_S[c_is_nth_entropy_field] - in.composition[i][entropy_indices[c_is_nth_entropy_field]]) / this->get_timestep();

                              //     std::cout << "reaction_rate_out = " << composition_equalibrated_S[c_is_nth_entropy_field] <<" " << std::endl;
                            }
                        }

                      out.reaction_terms[i][c] = 0.0;


                    }
                }
            }




          // set up variable to interpolate prescribed field outputs onto compositional fields
          if (PrescribedFieldOutputs<dim> *prescribed_field_out = out.template get_additional_output<PrescribedFieldOutputs<dim>>())
            {
              prescribed_field_out->prescribed_field_outputs[i][projected_density_index] = out.densities[i];
            }

          // set up variable to interpolate prescribed field outputs onto temperature field
          if (PrescribedTemperatureOutputs<dim> *prescribed_temperature_out = out.template get_additional_output<PrescribedTemperatureOutputs<dim>>())
            {
              prescribed_temperature_out->prescribed_temperature_outputs[i] = temperature_lookup;
            }

          // Calculate Viscosity
          if (in.requests_property(MaterialProperties::viscosity))
            {
              // read in the viscosity profile
              const double depth = this->get_geometry_model().depth(in.position[i]);
              const double viscosity_profile = depth_dependent_rheology->compute_viscosity(depth);

              // lateral viscosity variations
              const double reference_temperature = this->get_adiabatic_conditions().is_initialized()
                                                   ?
                                                   this->get_adiabatic_conditions().temperature(in.position[i])
                                                   :
                                                   this->get_parameters().adiabatic_surface_temperature;

              const double delta_temperature = temperature_lookup-reference_temperature;

              // Steinberger & Calderwood viscosity
              if (temperature_lookup*reference_temperature == 0)
                out.viscosities[i] = max_eta;
              else
                {
                  double vis_lateral = std::exp(-lateral_viscosity_prefactor_lookup->lateral_viscosity(depth)*delta_temperature/(temperature_lookup*reference_temperature));
                  // lateral vis variation
                  vis_lateral = std::max(std::min((vis_lateral),max_lateral_eta_variation),1/max_lateral_eta_variation);

                  if (std::isnan(vis_lateral))
                    vis_lateral = 1.0;

                  double effective_viscosity = vis_lateral * viscosity_profile;

                  const double strain_rate_effective = std::fabs(second_invariant(deviator(in.strain_rate[i])));

                  if (std::sqrt(strain_rate_effective) >= std::numeric_limits<double>::min())
                    {
                      const double pressure =  this->get_adiabatic_conditions().pressure(in.position[i]);
                      const double eta_plastic = drucker_prager_plasticity.compute_viscosity(cohesion,
                                                                                             angle_of_internal_friction,
                                                                                             pressure,
                                                                                             std::sqrt(strain_rate_effective),
                                                                                             std::numeric_limits<double>::infinity());

                      effective_viscosity = 1.0 / ( ( 1.0 /  eta_plastic  ) + ( 1.0 / (vis_lateral * viscosity_profile) ) );

                      PlasticAdditionalOutputs<dim> *plastic_out = out.template get_additional_output<PlasticAdditionalOutputs<dim>>();
                      if (plastic_out != nullptr)
                        {
                          plastic_out->cohesions[i] = cohesion;
                          plastic_out->friction_angles[i] = angle_of_internal_friction;
                          plastic_out->yielding[i] = eta_plastic < (vis_lateral * viscosity_profile) ? 1 : 0;
                        }
                    }
                  out.viscosities[i] = std::max(std::min(effective_viscosity,max_eta),min_eta);
                }
            }

          // fill seismic velocities outputs if they exist
          if (SeismicAdditionalOutputs<dim> *seismic_out = out.template get_additional_output<SeismicAdditionalOutputs<dim>>())
            {

              std::vector<double> vp (material_file_names.size());
              std::vector<double> vs (material_file_names.size());
              for (unsigned int j=0; j<material_file_names.size(); ++j)
                {
                  vp[j] = entropy_reader[j]->seismic_vp(component_entropy[j],pressure);
                  vs[j] = entropy_reader[j]->seismic_vs(component_entropy[j],pressure);
                }
              seismic_out->vp[i] = MaterialUtilities::average_value (volume_fractions, vp, MaterialUtilities::arithmetic);
              seismic_out->vs[i] = MaterialUtilities::average_value (volume_fractions, vs, MaterialUtilities::arithmetic);
            }
        }
    }




    template <int dim>
    void
    EntropyModel<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Entropy model");
        {
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/entropy-table/opxtable/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the `data/' subdirectory of ASPECT.");
          prm.declare_entry ("Material file name", "material_table.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file name of the material data. The first material data file is intended for the background composition. ");
          prm.declare_entry ("Reference viscosity", "1e22",
                             Patterns::Double(0),
                             "The viscosity that is used in this model. "
                             "\n\n"
                             "Units: \\si{\\pascal\\second}");
          prm.declare_entry ("Lateral viscosity file name", "temp-viscosity-prefactor.txt",
                             Patterns::Anything (),
                             "The file name of the lateral viscosity prefactor.");
          prm.declare_entry ("Minimum viscosity", "1e19",
                             Patterns::Double (0.),
                             "The minimum viscosity that is allowed in the viscosity "
                             "calculation. Smaller values will be cut off.");
          prm.declare_entry ("Maximum viscosity", "1e23",
                             Patterns::Double (0.),
                             "The maximum viscosity that is allowed in the viscosity "
                             "calculation. Larger values will be cut off.");
          prm.declare_entry ("Maximum lateral viscosity variation", "1e2",
                             Patterns::Double (0.),
                             "The relative cutoff value for lateral viscosity variations "
                             "caused by temperature deviations. The viscosity may vary "
                             "laterally by this factor squared.");
          prm.declare_entry ("Angle of internal friction", "0.",
                             Patterns::Double (0.),
                             "The value of the angle of internal friction, $\\phi$."
                             "For a value of zero, in 2D the von Mises criterion is retrieved. "
                             "Angles higher than 30 degrees are harder to solve numerically."
                             "Units: degrees.");
          prm.declare_entry ("Cohesion", "1e20",
                             Patterns::Double (0.),
                             "The value of the cohesion, $C$. The extremely large default"
                             "cohesion value (1e20 Pa) prevents the viscous stress from "
                             "exceeding the yield stress. Units: \\si{\\pascal}.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Thermal conductivity formulation", "constant",
                             Patterns::Selection("constant|p-T-dependent"),
                             "Which law should be used to compute the thermal conductivity. "
                             "The 'constant' law uses a constant value for the thermal "
                             "conductivity. The 'p-T-dependent' formulation uses equations "
                             "from Stackhouse et al. (2015): First-principles calculations "
                             "of the lattice thermal conductivity of the lower mantle "
                             "(https://doi.org/10.1016/j.epsl.2015.06.050), and Tosi et al. "
                             "(2013): Mantle dynamics with pressure- and temperature-dependent "
                             "thermal expansivity and conductivity "
                             "(https://doi.org/10.1016/j.pepi.2013.02.004) to compute the "
                             "thermal conductivity in dependence of temperature and pressure. "
                             "The thermal conductivity parameter sets can be chosen in such a "
                             "way that either the Stackhouse or the Tosi relations are used. "
                             "The conductivity description can consist of several layers with "
                             "different sets of parameters. Note that the Stackhouse "
                             "parametrization is only valid for the lower mantle (bridgmanite).");
          prm.declare_entry ("Thermal conductivity transition depths", "410000, 520000, 660000",
                             Patterns::List(Patterns::Double (0.)),
                             "A list of depth values that indicate where the transitions between "
                             "the different conductivity parameter sets should occur in the "
                             "'p-T-dependent' Thermal conductivity formulation (in most cases, "
                             "this will be the depths of major mantle phase transitions). "
                             "Units: \\si{\\meter}.");
          prm.declare_entry ("Reference thermal conductivities", "2.47, 3.81, 3.52, 4.9",
                             Patterns::List(Patterns::Double (0.)),
                             "A list of base values of the thermal conductivity for each of the "
                             "horizontal layers in the 'p-T-dependent' thermal conductivity "
                             "formulation. Pressure- and temperature-dependence will be applied"
                             "on top of this base value, according to the parameters 'Pressure "
                             "dependencies of thermal conductivity' and 'Reference temperatures "
                             "for thermal conductivity'. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}");
          prm.declare_entry ("Pressure dependencies of thermal conductivity", "3.3e-10, 3.4e-10, 3.6e-10, 1.05e-10",
                             Patterns::List(Patterns::Double ()),
                             "A list of values that determine the linear scaling of the "
                             "thermal conductivity with the pressure in the 'p-T-dependent' "
                             "thermal conductivity formulation. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin\\per\\pascal}.");
          prm.declare_entry ("Reference temperatures for thermal conductivity", "300, 300, 300, 1200",
                             Patterns::List(Patterns::Double (0.)),
                             "A list of values of reference temperatures used to determine "
                             "the temperature-dependence of the thermal conductivity in the "
                             "'p-T-dependent' thermal conductivity formulation. "
                             "Units: \\si{\\kelvin}.");
          prm.declare_entry ("Thermal conductivity exponents", "0.48, 0.56, 0.61, 1.0",
                             Patterns::List(Patterns::Double (0.)),
                             "A list of exponents in the temperature-dependent term of the "
                             "'p-T-dependent' thermal conductivity formulation. Note that this "
                             "exponent is not used (and should have a value of 1) in the "
                             "formulation of Stackhouse et al. (2015). "
                             "Units: none.");
          prm.declare_entry ("Saturation prefactors", "0, 0, 0, 1",
                             Patterns::List(Patterns::Double (0., 1.)),
                             "A list of values that indicate how a given layer in the "
                             "conductivity formulation should take into account the effects "
                             "of saturation on the temperature-dependence of the thermal "
                             "conducitivity. This factor is multiplied with a saturation function "
                             "based on the theory of Roufosse and Klemens, 1974. A value of 1 "
                             "reproduces the formulation of Stackhouse et al. (2015), a value of "
                             "0 reproduces the formulation of Tosi et al., (2013). "
                             "Units: none.");
          prm.declare_entry ("Maximum thermal conductivity", "1000",
                             Patterns::Double (0.),
                             "The maximum thermal conductivity that is allowed in the "
                             "model. Larger values will be cut off.");
          prm.leave_subsection();
        }

        // Depth-dependent parameters from the rheology plugin
        Rheology::AsciiDepthProfile<dim>::declare_parameters(prm,
                                                             "Depth dependent viscosity");

        prm.leave_subsection();
      }
    }



    template <int dim>
    void
    EntropyModel<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Entropy model");
        {
          data_directory              = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          material_file_names          = Utilities::split_string_list(prm.get ("Material file name"));
          lateral_viscosity_file_name  = prm.get ("Lateral viscosity file name");
          min_eta                     = prm.get_double ("Minimum viscosity");
          max_eta                     = prm.get_double ("Maximum viscosity");
          max_lateral_eta_variation    = prm.get_double ("Maximum lateral viscosity variation");
          thermal_conductivity_value  = prm.get_double ("Thermal conductivity");

          if (prm.get ("Thermal conductivity formulation") == "constant")
            conductivity_formulation = constant;
          else if (prm.get ("Thermal conductivity formulation") == "p-T-dependent")
            conductivity_formulation = p_T_dependent;
          else
            AssertThrow(false, ExcMessage("Not a valid thermal conductivity formulation"));

          conductivity_transition_depths = Utilities::string_to_double
                                           (Utilities::split_string_list(prm.get ("Thermal conductivity transition depths")));
          const unsigned int n_conductivity_layers = conductivity_transition_depths.size() + 1;
          AssertThrow (std::is_sorted(conductivity_transition_depths.begin(), conductivity_transition_depths.end()),
                       ExcMessage("The list of 'Thermal conductivity transition depths' must "
                                  "be sorted such that the values increase monotonically."));

          reference_thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference thermal conductivities"))),
                                                                                     n_conductivity_layers,
                                                                                     "Reference thermal conductivities");
          conductivity_pressure_dependencies = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Pressure dependencies of thermal conductivity"))),
                                                                                       n_conductivity_layers,
                                                                                       "Pressure dependencies of thermal conductivity");
          conductivity_reference_temperatures = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference temperatures for thermal conductivity"))),
                                                                                        n_conductivity_layers,
                                                                                        "Reference temperatures for thermal conductivity");
          conductivity_exponents = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivity exponents"))),
                                                                           n_conductivity_layers,
                                                                           "Thermal conductivity exponents");
          saturation_scaling = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Saturation prefactors"))),
                                                                       n_conductivity_layers,
                                                                       "Saturation prefactors");
          maximum_conductivity = prm.get_double ("Maximum thermal conductivity");

          angle_of_internal_friction = prm.get_double ("Angle of internal friction") * constants::degree_to_radians;
          cohesion = prm.get_double("Cohesion");

          prm.leave_subsection();
        }

        depth_dependent_rheology = std::make_unique<Rheology::AsciiDepthProfile<dim>>();
        depth_dependent_rheology->initialize_simulator (this->get_simulator());
        depth_dependent_rheology->parse_parameters(prm, "Depth dependent viscosity");
        depth_dependent_rheology->initialize();

        prm.leave_subsection();

        // Declare dependencies on solution variables
        this->model_dependence.viscosity = NonlinearDependence::temperature;
        this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.compressibility = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.specific_heat = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      }
    }



    template <int dim>
    void
    EntropyModel<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<SeismicAdditionalOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>> (n_points));
        }

      if (out.template get_additional_output<PrescribedFieldOutputs<dim>>() == NULL)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>>
            (n_points, this->n_compositional_fields()));
        }

      if (out.template get_additional_output<PrescribedTemperatureOutputs<dim>>() == NULL)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedTemperatureOutputs<dim>>
            (n_points));
        }

      if (out.template get_additional_output<PlasticAdditionalOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<PlasticAdditionalOutputs<dim>> (n_points));
        }

      if (this->get_parameters().use_operator_splitting
          && out.template get_additional_output<ReactionRateOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::ReactionRateOutputs<dim>> (n_points, this->n_compositional_fields()));
        }

    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(EntropyModel,
                                   "entropy model",
                                   "A material model that is designed to use pressure and entropy (rather "
                                   "than pressure and temperature) as independent variables. "
                                   "It requires a thermodynamic data table that contains "
                                   "all relevant properties in a specific format as illustrated in "
                                   "the data/material-model/entropy-table/opxtable example folder. "
                                   "The material model requires the use of the projected density "
                                   "approximation for compressibility, and the existence of a "
                                   "compositional field called 'entropy'.")
  }
}
