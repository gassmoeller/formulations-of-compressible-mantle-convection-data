/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/simulator/assemblers/interface.h>
#include <aspect/simulator/assemblers/stokes.h>
#include <aspect/material_model/simple_compressible.h>
#include <aspect/simulator_access.h>
#include <aspect/simulator_signals.h>
#include <aspect/adiabatic_conditions/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that is identical to the simple compressible model,
     * except that the density is tracked in a compositional field using
     * the reactions.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class ProjectedDensity : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        virtual
        void initialize();

        virtual
        void update();

        virtual
        bool is_compressible () const;

        virtual
        double reference_viscosity () const;

        virtual
        void
        evaluate (const MaterialModelInputs<dim> &in,
                  MaterialModelOutputs<dim> &out) const;

        virtual
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;

        static void
        declare_parameters (ParameterHandler &prm);

        virtual void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Pointer to the material model used as the base model
         */
        std::shared_ptr<MaterialModel::Interface<dim> > base_model;

        /**
         * The reference density
         */
        double reference_rho;

        /**
         * The constant thermal expansivity
         */
        double thermal_alpha;
        bool use_exponential_alpha;

        /**
         * The constant compressibility.
         */
        double reference_compressibility;

        /**
         * Whether to use drucker-prager plasticity to
         * limit the viscosity.
         */
        bool use_drucker_prager_viscosity;

        /**
         * The angle of internal friction
         */
        double angle_of_internal_friction;

        /**
         * The cohesion
         */
        double cohesion;

        /**
         * The applied viscosity bounds
         */
        double minimum_viscosity;
        double maximum_viscosity;

        /**
         * The reference strain rate used as a first estimate
         */
        double reference_strain_rate;

        /**
         * Reduce viscosity according to composition.
         */
        bool reduce_viscosity_by_composition;

        /**
         * Use the adiabatic instead of the full pressure for the
         * density calculation.
         */
        bool use_adiabatic_pressure_for_density;
    };



    template <int dim>
    void
    ProjectedDensity<dim>::initialize()
    {
      base_model->initialize();
    }



    template <int dim>
    void
    ProjectedDensity<dim>::update()
    {
      base_model->update();
    }



    template <int dim>
    double
    ProjectedDensity<dim>::
    reference_viscosity () const
    {
      return base_model->reference_viscosity();
    }



    template <int dim>
    bool
    ProjectedDensity<dim>::
    is_compressible () const
    {
      return this->get_parameters().formulation_mass_conservation != Parameters<dim>::Formulation::MassConservation::incompressible;
    }



    template <int dim>
    void
    ProjectedDensity<dim>::evaluate(const MaterialModelInputs<dim> &in,
             MaterialModelOutputs<dim> &out) const
    {
      base_model->evaluate(in,out);

      const unsigned int projected_density_index = this->introspection().compositional_index_for_name("projected_density");

      unsigned int fault_index = numbers::invalid_unsigned_int;
      unsigned int crust_index = numbers::invalid_unsigned_int;

      if (this->introspection().compositional_name_exists("fault"))
        fault_index = this->introspection().compositional_index_for_name("fault");

      if (this->introspection().compositional_name_exists("crust"))
        crust_index = this->introspection().compositional_index_for_name("crust");

      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const double transition_pressure = 4.e9;
          const double transition_temperature = 870.;
          const double pressure_width = 5.e7; //2.5e8;
          const double slope = 0;
          const double density_jump = 600.0;

          double pressure_for_density = in.pressure[i];

          // If using the PDA avoid using the real pressure for the density calculation to prevent
          // pressure waves from forming
          if (use_adiabatic_pressure_for_density)
        	  pressure_for_density = this->get_adiabatic_conditions().pressure(in.position[i]);

          // then calculate the deviation from the transition point (both in temperature
          // and in pressure)
          const double pressure_deviation = pressure_for_density - transition_pressure
              - slope * (in.temperature[i] - transition_temperature);

          // last, calculate the percentage of material that has undergone the transition
          // (also in dependence of the phase transition width - this is an input parameter)
          const double phase_dependence =  (this->introspection().compositional_name_exists("crust"))
        		      ?
                                  0.5*(-1.0 + std::tanh(pressure_deviation / pressure_width)) * density_jump * in.composition[i][crust_index]
        		                                                                                                         :
        		          0.0;

          double thermal_expansivity = thermal_alpha;
          if (use_exponential_alpha)
            {
        	  thermal_expansivity *= std::exp(-1.117979323e-11*pressure_for_density);
        	  out.thermal_expansion_coefficients[i] = thermal_expansivity;
            }

          // Use a thermodynamically consistent thermal expansivity
          // (if alpha is constant, rho needs to depend exponentially on it, not linearly)
          out.densities[i] = (reference_rho + phase_dependence) * std::exp(reference_compressibility * (pressure_for_density - this->get_surface_pressure()) -
        		  thermal_expansivity * (in.temperature[i] - this->get_adiabatic_surface_temperature()));

          // For the ICA, we have to provide the isentropic rather than isothermal compressibility
          // for the mass conservation equation.
          if (this->get_parameters().formulation_mass_conservation ==
        		  Parameters<dim>::Formulation::MassConservation::isentropic_compression)
            {
        	  out.compressibilities[i] = reference_compressibility - std::pow(thermal_expansivity, 2) * in.temperature[i]
															         / (out.densities[i] * out.specific_heat[i]);
            }

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            if (c == projected_density_index)
              out.reaction_terms[i][c] = out.densities[i] - in.composition[i][c];
            else
              out.reaction_terms[i][c] = 0.0;

          if (reduce_viscosity_by_composition && this->introspection().compositional_name_exists("fault"))
          out.viscosities[i] = out.viscosities[i] * (1.0 - in.composition[i][fault_index]) +
              out.viscosities[i] * 1e-2 * in.composition[i][fault_index];

          // calculate effective viscosity
          if (in.strain_rate.size() > 0 && use_drucker_prager_viscosity)
            {
              const SymmetricTensor<2,dim> strain_rate_deviator = deviator(in.strain_rate[i]);

              // For the very first time this function is called
              // (the first iteration of the first timestep), this function is called
              // with a zero input strain rate. We provide a representative reference
              // strain rate for this case, which avoids division by zero and produces
              // a representative first guess of the viscosities.
              //
              // In later iterations and timesteps we calculate the second moment
              // invariant of the deviatoric strain rate tensor.
              // This is equal to the negative of the second principle
              // invariant of the deviatoric strain rate (calculated with the function second_invariant),
              // as shown in Appendix A of Zienkiewicz and Taylor (Solid Mechanics, 2000).
              //
              // The negative of the second principle invariant is equal to 0.5 e_dot_dev_ij e_dot_dev_ji,
              // where e_dot_dev is the deviatoric strain rate tensor. The square root of this quantity
              // gives the common definition of effective strain rate.
              const double edot_ii_strict = ((this->get_timestep_number() == 0)
                  &&
                  (in.strain_rate[i].norm() <= std::numeric_limits<double>::min())
                  ?
                      reference_strain_rate * reference_strain_rate
                      :
                      std::fabs(second_invariant(strain_rate_deviator)));

                      const double strain_rate_effective = edot_ii_strict;

                      if (std::sqrt(strain_rate_effective) > std::numeric_limits<double>::min())
                        {
                          // plasticity
                          const MaterialUtilities::DruckerPragerInputs plastic_in(cohesion, angle_of_internal_friction, in.pressure[i], std::sqrt(strain_rate_effective));
                          MaterialUtilities::DruckerPragerOutputs plastic_out;
                          MaterialUtilities::compute_drucker_prager_yielding<dim> (plastic_in, plastic_out);

                          const double eta_plastic = plastic_out.plastic_viscosity;

                          // Cut off the viscosity between a minimum and maximum value to avoid
                          // a numerically unfavourable large viscosity range.
                          const double effective_viscosity = std::min(std::max(eta_plastic, minimum_viscosity), maximum_viscosity);

                          Assert(dealii::numbers::is_finite(effective_viscosity), ExcMessage ("Error: Viscosity is not finite."));

                          out.viscosities[i] = std::min(out.viscosities[i],effective_viscosity);

                          Assert(dealii::numbers::is_finite(out.viscosities[i]),
                                 ExcMessage ("Error: Averaged viscosity is not finite."));
                        }
            }
        }

      // set up variable to interpolate prescribed field outputs onto compositional fields
      PrescribedFieldOutputs<dim> *prescribed_field_out = out.template get_additional_output<PrescribedFieldOutputs<dim> >();

      if (prescribed_field_out != NULL)
        for (unsigned int i=0; i < in.position.size(); ++i)
          prescribed_field_out->prescribed_field_outputs[i][projected_density_index] = out.densities[i];
    }



    template <int dim>
    void
    ProjectedDensity<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<PrescribedFieldOutputs<dim> >() == NULL)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std::unique_ptr<MaterialModel::AdditionalMaterialOutputs<dim> >
            (new MaterialModel::PrescribedFieldOutputs<dim> (n_points, this->n_compositional_fields())));
        }
    }



    template <int dim>
    void
    ProjectedDensity<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Projected density model");
        {
          prm.declare_entry ("Minimum viscosity", "1e19",
              Patterns::Double (0),
              "The value of the minimum viscosity cutoff $\\eta_min$. Units: $Pa\\;s$.");
          prm.declare_entry ("Maximum viscosity", "1e24",
              Patterns::Double (0),
              "The value of the maximum viscosity cutoff $\\eta_max$. Units: $Pa\\;s$.");
          prm.declare_entry ("Reference strain rate", "1e-15",
              Patterns::Double (0),
              "The value of the initial strain rate prescribed during the "
              "first nonlinear iteration $\\dot{\\epsilon}_ref$. Units: $1/s$.");
          prm.declare_entry ("Angle of internal friction", "0",
              Patterns::Double (0),
              "The value of the angle of internal friction $\\phi$. "
              "For a value of zero, in 2D the von Mises "
              "criterion is retrieved. Angles higher than 30 degrees are "
              "harder to solve numerically. Units: degrees.");
          prm.declare_entry ("Cohesion", "2e7",
              Patterns::Double (0),
              "The value of the cohesion $C$. Units: $Pa$.");
          prm.declare_entry("Use drucker prager viscosity","false",
              Patterns::Bool(),
              "Whether to use drucker prager plasticity for the viscosity.");
          prm.declare_entry ("Reduce viscosity for composition 1", "false",
              Patterns::Bool(),
              "");
          prm.declare_entry ("Use adiabatic pressure for density", "false",
              Patterns::Bool(),
              "");
          prm.declare_entry ("Use exponential thermal expansivity", "false",
              Patterns::Bool(),
              "");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      SimpleCompressible<dim>::declare_parameters (prm);
    }



    template <int dim>
    void
    ProjectedDensity<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Projected density model");
        {
          minimum_viscosity          = prm.get_double ("Minimum viscosity");
          maximum_viscosity          = prm.get_double ("Maximum viscosity");
          reference_strain_rate      = prm.get_double ("Reference strain rate");
          // Convert degrees to radians
          angle_of_internal_friction = prm.get_double ("Angle of internal friction") * numbers::PI/180.0;
          cohesion                   = prm.get_double ("Cohesion");

          use_drucker_prager_viscosity = prm.get_bool("Use drucker prager viscosity");

          reduce_viscosity_by_composition = prm.get_bool("Reduce viscosity for composition 1");

          use_adiabatic_pressure_for_density = prm.get_bool("Use adiabatic pressure for density");
          use_exponential_alpha = prm.get_bool("Use exponential thermal expansivity");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();


      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple compressible model");
        {
          reference_rho              = prm.get_double ("Reference density");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          reference_compressibility  = prm.get_double ("Reference compressibility");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // create the base model and initialize its SimulatorAccess base
      // class; it will get a chance to read its parameters below after we
      // leave the current section
      base_model.reset(create_material_model<dim>("simple compressible"));
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
        sim->initialize_simulator (this->get_simulator());
      base_model->parse_parameters(prm);
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ProjectedDensity,
                                   "projected density",
                                   "")
  }
}

