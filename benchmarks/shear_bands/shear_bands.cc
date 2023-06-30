/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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
#include <aspect/initial_composition/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/postprocess/interface.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/melt.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>
#include <array>


namespace aspect
{
  /**
   * This is the "Shear bands" benchmark defined in the following paper:
   * @code
   *  @Article{DMGT11,
   *    author =       {R. F. Katz and M. Spiegelman and B. Holtzman},
   *    title =        {The dynamics of melt and shear localization in
   *                    partially molten aggregates},
   *    journal =      {Nature},
   *    year =         2006,
   *    volume =       442,
   *    pages =        {676-679}
   * @endcode
   *
   * Magmatic shear bands are generated by a combination of porosity- and
   * strain-rate dependent viscosity in a simple shear regime; and the
   * angle of the shear bands is measured.
   *
   */
  namespace ShearBands
  {
    using namespace dealii;


    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class ShearBandsMaterial : public MaterialModel::MeltInterface<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:
        bool is_compressible () const override
        {
          return false;
        }

        double reference_viscosity () const
        {
          return eta_0;
        }

        double reference_darcy_coefficient () const override
        {
          return reference_permeability * pow(0.01, permeability_exponent) / eta_f;
        }


        double
        get_background_porosity () const;

        double
        get_reference_compaction_viscosity () const;

        double
        get_porosity_exponent () const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                      typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const override
        {
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
          const double strain_rate_dependence = (1.0 - dislocation_creep_exponent) / dislocation_creep_exponent;

          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
              double porosity = std::max(in.composition[i][porosity_idx],1e-4);
              out.viscosities[i] = eta_0 * std::exp(alpha*(porosity - background_porosity));
              if (in.requests_property(MaterialModel::MaterialProperties::viscosity))
                {
                  const SymmetricTensor<2,dim> shear_strain_rate = in.strain_rate[i]
                                                                   - 1./dim * trace(in.strain_rate[i]) * unit_symmetric_tensor<dim>();
                  const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

                  if (std::abs(second_strain_rate_invariant) > 1e-30)
                    out.viscosities[i] *= std::pow(second_strain_rate_invariant,strain_rate_dependence);
                }


              out.densities[i] = reference_rho_s;
              out.thermal_expansion_coefficients[i] = 0.0;
              out.specific_heat[i] = 1.0;
              out.thermal_conductivities[i] = 0.0;
              out.compressibilities[i] = 0.0;
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                out.reaction_terms[i][c] = 0.0;
            }

          // fill melt outputs if they exist
          aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim>>();

          if (melt_out != NULL)
            for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
              {
                double porosity = std::max(in.composition[i][porosity_idx],1e-4);

                melt_out->compaction_viscosities[i] = xi_0 * pow(porosity/background_porosity,-compaction_viscosity_exponent);
                melt_out->fluid_viscosities[i]= eta_f;
                melt_out->permeabilities[i]= reference_permeability * std::pow(porosity,permeability_exponent);
                melt_out->fluid_densities[i]= reference_rho_f;
                melt_out->fluid_density_gradients[i] = 0.0;
              }
        }

      private:
        double reference_rho_s;
        double reference_rho_f;
        double eta_0;
        double xi_0;
        double eta_f;
        double reference_permeability;
        double dislocation_creep_exponent;
        double alpha;
        double background_porosity;
        double compaction_viscosity_exponent;
        double permeability_exponent;
    };

    template <int dim>
    double
    ShearBandsMaterial<dim>::get_background_porosity () const
    {
      return background_porosity;
    }

    template <int dim>
    double
    ShearBandsMaterial<dim>::get_reference_compaction_viscosity () const
    {
      return xi_0;
    }

    template <int dim>
    double
    ShearBandsMaterial<dim>::get_porosity_exponent () const
    {
      return alpha;
    }


    template <int dim>
    void
    ShearBandsMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Shear bands material");
        {
          prm.declare_entry ("Reference solid density", "3000",
                             Patterns::Double (0),
                             "Reference density of the solid $\\rho_{s,0}$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Reference melt density", "3000",
                             Patterns::Double (0),
                             "Reference density of the melt/fluid$\\rho_{f,0}$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Reference shear viscosity", "1.41176e7",
                             Patterns::Double (0),
                             "The value of the constant viscosity $\\eta_0$ of the solid matrix. "
                             "Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Reference compaction viscosity", "1.41176e8",
                             Patterns::Double (0),
                             "The value of the constant volumetric viscosity $\\xi_0$ of the solid matrix. "
                             "Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Reference melt viscosity", "100.0",
                             Patterns::Double (0),
                             "The value of the constant melt viscosity $\\eta_f$. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Reference permeability", "5e-9",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: \\si{\\meter\\squared}.");
          prm.declare_entry ("Dislocation creep exponent", "6.0",
                             Patterns::Double(0),
                             "Power-law exponent $n_{dis}$ for dislocation creep. "
                             "Units: none.");
          prm.declare_entry ("Porosity exponent", "-27.0",
                             Patterns::Double(),
                             "Exponent $\alpha$ for the exponential porosity-dependence "
                             "of the viscosity. "
                             "Units: none.");
          prm.declare_entry ("Background porosity", "0.05",
                             Patterns::Double (0),
                             "Background porosity used in the viscosity law. Units: none.");
          prm.declare_entry ("Compaction viscosity exponent", "0.0",
                             Patterns::Double(0),
                             "Power-law exponent $m$ for the compaction viscosity. "
                             "Units: none.");
          prm.declare_entry ("Permeability exponent", "3.0",
                             Patterns::Double(0),
                             "Power-law exponent $n$ for the porosity-dependence of permeability. "
                             "Units: none.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ShearBandsMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Shear bands material");
        {
          reference_rho_s            = prm.get_double ("Reference solid density");
          reference_rho_f            = prm.get_double ("Reference melt density");
          eta_0                      = prm.get_double ("Reference shear viscosity");
          xi_0                       = prm.get_double ("Reference compaction viscosity");
          eta_f                      = prm.get_double ("Reference melt viscosity");
          reference_permeability     = prm.get_double ("Reference permeability");
          dislocation_creep_exponent = prm.get_double ("Dislocation creep exponent");
          alpha                      = prm.get_double ("Porosity exponent");
          background_porosity        = prm.get_double ("Background porosity");
          compaction_viscosity_exponent = prm.get_double ("Compaction viscosity exponent");
          permeability_exponent      = prm.get_double ("Permeability exponent");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    /**
     * An initial conditions model for the shear bands benchmark
     * that implements a background porosity plus white noise with
     * a certain amplitude.
     */
    template <int dim>
    class ShearBandsInitialCondition : public InitialComposition::Interface<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Return the initial porosity as a function of position.
         */
        double
        initial_composition (const Point<dim> &position, const unsigned int n_comp) const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Initialization function.
         */
        void
        initialize () override;

      private:
        double noise_amplitude;
        double background_porosity;
        std::array<unsigned int,dim> grid_intervals;
        Functions::InterpolatedUniformGridData<dim> *interpolate_noise;
    };



    template <int dim>
    void
    ShearBandsInitialCondition<dim>::initialize ()
    {
      AssertThrow(Plugins::plugin_type_matches<const ShearBandsMaterial<dim>>(this->get_material_model()),
                  ExcMessage("Initial condition shear bands only works with the material model shear bands."));

      const ShearBandsMaterial<dim> &
      material_model
        = Plugins::get_plugin_as_type<const ShearBandsMaterial<dim>>(this->get_material_model());

      background_porosity = material_model.get_background_porosity();

      AssertThrow(noise_amplitude < background_porosity,
                  ExcMessage("Amplitude of the white noise must be smaller "
                             "than the background porosity."));

      Point<dim> extents;
      TableIndices<dim> size_idx;
      for (unsigned int d=0; d<dim; ++d)
        size_idx[d] = grid_intervals[d]+1;

      Table<dim,double> white_noise;
      white_noise.TableBase<dim,double>::reinit(size_idx);
      std::array<std::pair<double,double>,dim> grid_extents;

      AssertThrow(Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("Initial condition shear bands only works with the box geometry model."));

      const GeometryModel::Box<dim> &
      geometry_model
        = Plugins::get_plugin_as_type<const GeometryModel::Box<dim>>(this->get_geometry_model());

      extents = geometry_model.get_extents();

      for (unsigned int d=0; d<dim; ++d)
        {
          grid_extents[d].first=0;
          grid_extents[d].second=extents[d];
        }

      // use a fixed number as seed for random generator
      // this is important if we run the code on more than 1 processor
      std::srand(0);

      TableIndices<dim> idx;

      for (unsigned int i=0; i<white_noise.size()[0]; ++i)
        {
          idx[0] = i;
          for (unsigned int j=0; j<white_noise.size()[1]; ++j)
            {
              idx[1] = j;
              white_noise(idx) = noise_amplitude * ((std::rand() % 10000) / 5000.0 - 1.0);
            }
        }

      interpolate_noise = new Functions::InterpolatedUniformGridData<dim> (grid_extents,
                                                                           grid_intervals,
                                                                           white_noise);
    }


    template <int dim>
    double
    ShearBandsInitialCondition<dim>::
    initial_composition (const Point<dim> &position, const unsigned int /*n_comp*/) const
    {
      return background_porosity + interpolate_noise->value(position);
    }


    template <int dim>
    void
    ShearBandsInitialCondition<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Shear bands initial condition");
        {
          prm.declare_entry ("Noise amplitude", "0.0005",
                             Patterns::Double (0),
                             "Amplitude of the white noise added to the initial "
                             "porosity. Units: none.");
          prm.declare_entry ("Grid intervals for noise X", "100",
                             Patterns::Integer (0),
                             "Grid intervals in X directions for the white noise added to "
                             "the initial background porosity that will then be interpolated "
                             "to the model grid. "
                             "Units: none.");
          prm.declare_entry ("Grid intervals for noise Y", "25",
                             Patterns::Integer (0),
                             "Grid intervals in Y directions for the white noise added to "
                             "the initial background porosity that will then be interpolated "
                             "to the model grid. "
                             "Units: none.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ShearBandsInitialCondition<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Shear bands initial condition");
        {
          noise_amplitude      = prm.get_double ("Noise amplitude");
          grid_intervals[0]    = prm.get_integer ("Grid intervals for noise X");
          grid_intervals[1]    = prm.get_integer ("Grid intervals for noise Y");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    /**
     * An initial conditions model for the plane wave melt bands benchmark
     * that implements a background porosity plus a plave wave in the porosity
     * field with a certain amplitude.
     */
    template <int dim>
    class PlaneWaveMeltBandsInitialCondition : public InitialComposition::Interface<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Return the initial porosity as a function of position.
         */
        virtual
        double
        initial_composition (const Point<dim> &position, const unsigned int n_comp) const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * Initialization function.
         */
        void
        initialize ();

        double
        get_wave_amplitude () const;

        double
        get_wave_number () const;

        double
        get_initial_band_angle () const;

      private:
        double amplitude;
        double background_porosity;
        double wave_number;
        double initial_band_angle;
    };



    template <int dim>
    void
    PlaneWaveMeltBandsInitialCondition<dim>::initialize ()
    {
      AssertThrow(Plugins::plugin_type_matches<ShearBandsMaterial<dim>>(this->get_material_model()),
                  ExcMessage("Initial condition shear bands only works with the material model shear bands."));

      const ShearBandsMaterial<dim> &
      material_model
        = Plugins::get_plugin_as_type<const ShearBandsMaterial<dim>>(this->get_material_model());

      background_porosity = material_model.get_background_porosity();

      AssertThrow(amplitude < 1.0,
                  ExcMessage("Amplitude of the melt bands must be smaller "
                             "than the background porosity."));
    }


    template <int dim>
    double
    PlaneWaveMeltBandsInitialCondition<dim>::get_wave_amplitude () const
    {
      return amplitude;
    }

    template <int dim>
    double
    PlaneWaveMeltBandsInitialCondition<dim>::get_wave_number () const
    {
      return wave_number;
    }

    template <int dim>
    double
    PlaneWaveMeltBandsInitialCondition<dim>::get_initial_band_angle () const
    {
      return initial_band_angle;
    }


    template <int dim>
    double
    PlaneWaveMeltBandsInitialCondition<dim>::
    initial_composition (const Point<dim> &position, const unsigned int /*n_comp*/) const
    {
      return background_porosity * (1.0 + amplitude * cos(wave_number*position[0]*sin(initial_band_angle)
                                                          + wave_number*position[1]*cos(initial_band_angle)));
    }


    template <int dim>
    void
    PlaneWaveMeltBandsInitialCondition<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Plane wave melt bands initial condition");
        {
          prm.declare_entry ("Wave amplitude", "1e-4",
                             Patterns::Double (0),
                             "Amplitude of the plane wave added to the initial "
                             "porosity. Units: none.");
          prm.declare_entry ("Wave number", "2000",
                             Patterns::Double (0),
                             "Wave number of the plane wave added to the initial "
                             "porosity. Is multiplied by 2 pi internally. "
                             "Units: 1/m.");
          prm.declare_entry ("Initial band angle", "30",
                             Patterns::Double (0),
                             "Initial angle of the plane wave added to the initial "
                             "porosity. Units: degrees.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    PlaneWaveMeltBandsInitialCondition<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Plane wave melt bands initial condition");
        {
          amplitude          = prm.get_double ("Wave amplitude");
          wave_number        = 2.0 * numbers::PI * prm.get_double ("Wave number");
          initial_band_angle = 2.0 * numbers::PI / 360.0 * prm.get_double ("Initial band angle");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }




    /**
      * A postprocessor that evaluates the angle of the shear bands.
      */
    template <int dim>
    class ShearBandsPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };


    template <int dim>
    std::pair<std::string,std::string>
    ShearBandsPostprocessor<dim>::execute (TableHandler &/*statistics*/)
    {
      // write output that can be used to calculate the angle of the shear bands
      const unsigned int max_lvl = this->get_triangulation().n_global_levels();

      std::stringstream output;
      output.precision (16);
      output << std::scientific;

      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          // write header
          output << "x                      y                      porosity" << std::endl;
        }

      // we want to have equidistant points in the output
      const QMidpoint<1> mp_rule;
      const QIterated<dim> quadrature_formula (mp_rule, 2);
      const unsigned int n_q_points =  quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(), this->get_fe(),  quadrature_formula,
                               update_JxW_values | update_values | update_quadrature_points);

      std::vector<double>         porosity_values (quadrature_formula.size());
      const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell != endc; ++cell)
        {
          if (!cell->is_locally_owned())
            continue;

          fe_values.reinit (cell);
          fe_values[this->introspection().extractors.compositional_fields[porosity_idx]].get_function_values (this->get_solution(), porosity_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              output << fe_values.quadrature_point (q) (0) << " "
                     << fe_values.quadrature_point (q) (1) << " "
                     << porosity_values[q] << std::endl;
            }
        }

      std::string filename = this->get_output_directory() + "shear_bands_" +
                             Utilities::int_to_string(max_lvl) +
                             ".csv";

      Utilities::collect_and_write_file_content(filename, output.str(), this->get_mpi_communicator());

      return std::make_pair("writing:", filename);
    }



    /**
      * A postprocessor that evaluates the groeth rate of the shear bands.
      *
      * The implementation of error evaluators correspond to the
      * benchmarks defined in the paper:
      * Laura Alisic, John F. Rudge, Richard F. Katz, Garth N. Wells, Sander Rhebergen:
      * "Compaction around a rigid, circular inclusion in partially molten rock."
      * Journal of Geophysical Research: Solid Earth 119.7 (2014): 5903-5920.
      */
    template <int dim>
    class ShearBandsGrowthRate : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

        /**
         * Initialization function. Take references to the material model and
         * initial conditions model to get the parameters necessary for computing
         * the analytical solution for the growth rate and store them.
         */
        void
        initialize ();

      private:
        double amplitude;
        double background_porosity;
        double initial_band_angle;
        double eta_0;
        double xi_0;
        double alpha;
    };


    template <int dim>
    void
    ShearBandsGrowthRate<dim>::initialize ()
    {
      AssertThrow(Plugins::plugin_type_matches<const ShearBandsMaterial<dim>>(this->get_material_model()),
                  ExcMessage("Postprocessor shear bands growth rate only works with the material model shear bands."));

      const ShearBandsMaterial<dim> &
      material_model
        = Plugins::get_plugin_as_type<const ShearBandsMaterial<dim>>(this->get_material_model());

      background_porosity = material_model.get_background_porosity();
      eta_0               = material_model.reference_viscosity();
      xi_0                = material_model.get_reference_compaction_viscosity();
      alpha               = material_model.get_porosity_exponent();

      const PlaneWaveMeltBandsInitialCondition<dim> &initial_composition
        = this->get_initial_composition_manager().template
          get_matching_initial_composition_model<PlaneWaveMeltBandsInitialCondition<dim>> ();

      amplitude           = initial_composition.get_wave_amplitude();
      initial_band_angle  = initial_composition.get_initial_band_angle();
    }


    template <int dim>
    std::pair<std::string,std::string>
    ShearBandsGrowthRate<dim>::execute (TableHandler &/*statistics*/)
    {
      // compute analytical melt band growth rate
      const double time = this->get_time();
      const Point<dim> upper_boundary_point = this->get_geometry_model().representative_point(0.0);
      const Point<dim> lower_boundary_point = this->get_geometry_model().representative_point(this->get_geometry_model().maximal_depth());

      types::boundary_id upper_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
      types::boundary_id lower_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id("bottom");

      const BoundaryVelocity::Manager<dim> &bm = this->get_boundary_velocity_manager();
      const double max_velocity = bm.boundary_velocity(upper_boundary, upper_boundary_point).norm();
      const double min_velocity = bm.boundary_velocity(lower_boundary, lower_boundary_point).norm();

      const double strain_rate = 0.5 * (max_velocity + min_velocity) / this->get_geometry_model().maximal_depth();
      const double theta = std::atan(std::sin(initial_band_angle) / (std::cos(initial_band_angle) - time * strain_rate/sqrt(2.0) * std::sin(initial_band_angle)));
      const double analytical_growth_rate = - eta_0 / (xi_0 + 4.0/3.0 * eta_0) * alpha * (1.0 - background_porosity)
                                            * 2.0 * strain_rate * std::sin(2.0 * theta);

      // integrate velocity divergence
      const QGauss<dim> quadrature_formula (this->get_fe()
                                            .base_element(this->introspection().base_elements.velocities).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_gradients   |
                               update_quadrature_points |
                               update_JxW_values);
      std::vector<double> velocity_divergences(n_q_points);
      std::vector<Point<dim>> position(n_q_points);

      double local_velocity_divergence_max = 0.0;
      double local_velocity_divergence_min = 0.0;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.velocities].get_function_divergences (this->get_solution(),
                velocity_divergences);
            position = fe_values.get_quadrature_points();

            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                const double relative_depth = this->get_geometry_model().depth(position[q]) / this->get_geometry_model().maximal_depth();
                const double relative_x = 0.25 * position[q](0) / this->get_geometry_model().maximal_depth();
                if (relative_depth < 0.55 && relative_depth > 0.45 && relative_x < 0.75 && relative_x > 0.25)
                  {
                    local_velocity_divergence_max = std::max (velocity_divergences[q], local_velocity_divergence_max);
                    local_velocity_divergence_min = std::min (velocity_divergences[q], local_velocity_divergence_min);
                  }
              }
          }

      const double global_velocity_divergence_max
        = Utilities::MPI::max (local_velocity_divergence_max, this->get_mpi_communicator());
      const double global_velocity_divergence_min
        = Utilities::MPI::min (local_velocity_divergence_min, this->get_mpi_communicator());

      // compute modelled melt band growth rate
      const double numerical_growth_rate = (0 <= initial_band_angle && initial_band_angle < numbers::PI /2.
                                            ?
                                            (1.0 - background_porosity) / (amplitude * background_porosity) * global_velocity_divergence_max
                                            :
                                            (1.0 - background_porosity) / (amplitude * background_porosity) * global_velocity_divergence_min);

      const double relative_error = (analytical_growth_rate != 0.0) ?
                                    std::abs(numerical_growth_rate - analytical_growth_rate)/analytical_growth_rate
                                    :
                                    1.0;

      // compute and output error
      std::ostringstream os;
      os << std::scientific << initial_band_angle
         << ", " << analytical_growth_rate
         << ", " << numerical_growth_rate
         << ", " << relative_error;

      return std::make_pair("Initial angle, Analytical growth rate, Modelled growth rate, Error:", os.str());
    }


  }
}



// explicit instantiations
namespace aspect
{
  namespace ShearBands
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ShearBandsMaterial,
                                   "shear bands material",
                                   "A material model that corresponds to the setup to"
                                   "generate magmatic shear bands described in Katz et al., "
                                   "Nature, 2006.")

    ASPECT_REGISTER_POSTPROCESSOR(ShearBandsPostprocessor,
                                  "shear bands statistics",
                                  "A postprocessor that writes output for the porosity field, "
                                  "which can be used to calculate the angle of the shear "
                                  "bands in the model.")

    ASPECT_REGISTER_POSTPROCESSOR(ShearBandsGrowthRate,
                                  "shear bands growth rate",
                                  "A postprocessor computes the growth rate of melt bands as "
                                  "presicted from linear stability analysis (Spiegelman, 2003), "
                                  "and compares modelled and analytical solution.")

    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(ShearBandsInitialCondition,
                                              "shear bands initial condition",
                                              "Composition is set to background porosity plus "
                                              "white noise.")

    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(PlaneWaveMeltBandsInitialCondition,
                                              "plane wave melt bands initial condition",
                                              "Composition is set to background porosity plus "
                                              "plane wave.")
  }
}
