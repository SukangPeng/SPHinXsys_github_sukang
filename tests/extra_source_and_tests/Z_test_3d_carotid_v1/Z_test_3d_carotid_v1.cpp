/**
 * @file     Z_test_3d_v8_vessel_solid_test1.cpp
 * @brief 	 test 正常的血管流动
 * @details  test
 *
 * @author 	Sukang Peng
 */

#include "Z_test_3d_carotid_v1.h"
#include "sphinxsys.h"
using namespace SPH;
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, dp_0);
    //sph_system.setRunParticleRelaxation(true); // Tag for run particle relaxation for body-fitted distribution
    //sph_system.setReloadParticles(false);      // Tag for computation with save particles distribution
    sph_system.setRunParticleRelaxation(false); // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(true);        // Tag for computation with save particles distribution
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.cd
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape()->cleanLevelSet();
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticlesWithReserve<BaseParticles, Reload>(in_outlet_particle_buffer, water_block.getName())
        : water_block.generateParticles<BaseParticles, Lattice>();

     BodyStatesRecordingToVtp write_water_block_to_vtp(water_block);
     write_water_block_to_vtp.writeToFile(0);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
   // wall_boundary.defineAdaptationRatios(1.15, 2.0);
    wall_boundary.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    wall_boundary.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

     BodyStatesRecordingToVtp write_vessel_wall_to_vtp(wall_boundary);
     write_vessel_wall_to_vtp.writeToFile(0);
        
      RealBody test_body_in(
          sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_cut_translation)), inlet_half_cut, "TestBodyCutIn"));
      test_body_in.generateParticles<BaseParticles, Lattice>();
      BodyStatesRecordingToVtp write_body_in_to_vtp(test_body_in);
      write_body_in_to_vtp.writeToFile(0);

      RealBody test_body_out_up(
          sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_up_emitter_rotation), Vec3d(outlet_up_cut_translation)), outlet_up_half_cut, "TestBodyCutOutUp"));
      test_body_out_up.generateParticles<BaseParticles, Lattice>();
      BodyStatesRecordingToVtp write_body_out_up_to_vtp(test_body_out_up);
      write_body_out_up_to_vtp.writeToFile(0);

      RealBody test_body_out_down(
          sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_down_emitter_rotation), Vec3d(outlet_down_cut_translation)), outlet_down_half_cut, "TestBodyCutOutDown"));
      test_body_out_down.generateParticles<BaseParticles, Lattice>();
      BodyStatesRecordingToVtp write_body_out_down_to_vtp(test_body_out_down);
      write_body_out_down_to_vtp.writeToFile(0);


      RealBody test_body_in2(
          sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_buffer_translation)), inlet_half, "TestBodyCutIn2"));
      test_body_in2.generateParticles<BaseParticles, Lattice>();
      BodyStatesRecordingToVtp write_body_in2_to_vtp(test_body_in2);
      write_body_in2_to_vtp.writeToFile(0);

      RealBody test_body_out_up2(
          sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_up_emitter_rotation), Vec3d(outlet_up_buffer_translation)), outlet_up_half, "TestBodyCutOutUp2"));
      test_body_out_up2.generateParticles<BaseParticles, Lattice>();
      BodyStatesRecordingToVtp write_body_out_up2_to_vtp(test_body_out_up2);
      write_body_out_up2_to_vtp.writeToFile(0);

      RealBody test_body_out_down2(
          sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_down_emitter_rotation), Vec3d(outlet_down_buffer_translation)), outlet_down_half, "TestBodyCutOutDown2"));
      test_body_out_down2.generateParticles<BaseParticles, Lattice>();
      BodyStatesRecordingToVtp write_body_out_down2_to_vtp(test_body_out_down2);
      write_body_out_down2_to_vtp.writeToFile(0);
      //----------------------------------------------------------------------
      //	SPH Particle relaxation section
      //----------------------------------------------------------------------
      /** check whether run particle relaxation for body fitted particle distribution. */
      if (sph_system.RunParticleRelaxation())
      {
          InnerRelation wall_inner(wall_boundary);
          InnerRelation blood_inner(water_block);
          using namespace relax_dynamics;
          SimpleDynamics<RandomizeParticlePosition> random_particles(wall_boundary);
          SimpleDynamics<RandomizeParticlePosition> random_blood_particles(water_block);
          RelaxationStepLevelSetCorrectionInner relaxation_step_inner(wall_inner);
          RelaxationStepLevelSetCorrectionInner relaxation_step_inner_blood(blood_inner);
          /** Write the body state to Vtp file. */
          BodyStatesRecordingToVtp write_state_to_vtp(sph_system);
          /** Write the particle reload files. */
          ReloadParticleIO write_particle_reload_files({&wall_boundary, &water_block});
          //----------------------------------------------------------------------
          //	Physics relaxation starts here.
          //----------------------------------------------------------------------
          random_particles.exec(0.25);
          random_blood_particles.exec(0.25);
          relaxation_step_inner.SurfaceBounding().exec();
          relaxation_step_inner_blood.SurfaceBounding().exec();
          write_state_to_vtp.writeToFile(0.0);
          //----------------------------------------------------------------------
          // From here the time stepping begins.
          //----------------------------------------------------------------------
          int ite = 0;
          int relax_step = 1000;
          while (ite < relax_step)
          {
              relaxation_step_inner.exec();
              relaxation_step_inner_blood.exec();
              ite++;
              if (ite % 200 == 0)
              {
                  std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                  write_state_to_vtp.writeToFile(ite);
              }
          }

          std::cout << "The physics relaxation process of wall particles finish !" << std::endl;
          write_state_to_vtp.writeToFile(ite);
          write_particle_reload_files.writeToFile(0);

          return 0;
      }
      //----------------------------------------------------------------------
      //	Define body relation map.
      //	The contact map gives the topological connections between the bodies.
      //	Basically the the range of bodies to build neighbor particle lists.
      //  Generally, we first define all the inner relations, then the contact relations.
      //  At last, we define the complex relaxations by combining previous defined
      //  inner and contact relations.
      //----------------------------------------------------------------------
      InnerRelation water_block_inner(water_block);
      ContactRelation water_wall_contact(water_block, {&wall_boundary});
      //----------------------------------------------------------------------
      // Combined relations built from basic relations
      // which is only used for update configuration.
      //----------------------------------------------------------------------
      ComplexRelation water_block_complex(water_block_inner, water_wall_contact);
      //----------------------------------------------------------------------
      //	Define the main numerical methods used in the simulation.
      //	Note that there may be data dependence on the constructors of these methods.
      //----------------------------------------------------------------------
      SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
      InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_wall_contact);
      InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_wall_contact);

      Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
      Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_wall_contact);
      InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_wall_contact);
      InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
      ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
      ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
      //----------------------------------------------------------------------
      // Bottom buffer
      //----------------------------------------------------------------------
      BodyAlignedBoxByCell left_emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_buffer_translation)), inlet_half));
      fluid_dynamics::BidirectionalBuffer<SPH::fluid_dynamics::NonPrescribedPressure> left_emitter_inflow_injection(left_emitter, in_outlet_particle_buffer);
      BodyAlignedBoxByCell right_up_emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_up_emitter_rotation), Vec3d(outlet_up_buffer_translation)), outlet_up_half));
      fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_up_emitter_inflow_injection(right_up_emitter, in_outlet_particle_buffer);
      BodyAlignedBoxByCell right_down_emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_down_emitter_rotation), Vec3d(outlet_down_buffer_translation)), outlet_down_half));
      fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_down_emitter_inflow_injection(right_down_emitter, in_outlet_particle_buffer);

      BodyAlignedBoxByCell left_disposer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(inlet_disposer_rotation), Vec3d(inlet_buffer_translation)), inlet_half));
      SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> left_disposer_outflow_deletion(left_disposer);
      BodyAlignedBoxByCell right_up_disposer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_up_disposer_rotation), Vec3d(outlet_up_buffer_translation)), outlet_up_half));
      SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> right_up_disposer_outflow_deletion(right_up_disposer);
      BodyAlignedBoxByCell right_down_disposer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_down_disposer_rotation), Vec3d(outlet_down_buffer_translation)), outlet_down_half));
      SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> right_down_disposer_outflow_deletion(right_down_disposer);

      InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_wall_contact);
      SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
      SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_up_inflow_pressure_condition(right_up_emitter);
      SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_down_inflow_pressure_condition(right_down_emitter);
      SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
      //----------------------------------------------------------------------
      //	Define the configuration related particles dynamics.
      //----------------------------------------------------------------------
      ParticleSorting particle_sorting(water_block);
      //----------------------------------------------------------------------
      //	Define the methods for I/O operations, observations
      //	and regression tests of the simulation.
      //----------------------------------------------------------------------
      BodyStatesRecordingToVtp body_states_recording(sph_system);
      body_states_recording.addToWrite<Real>(water_block, "Pressure");
      body_states_recording.addToWrite<int>(water_block, "Indicator");
      body_states_recording.addToWrite<Real>(water_block, "Density");
      body_states_recording.addToWrite<int>(water_block, "BufferParticleIndicator");
      body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
      //----------------------------------------------------------------------
      //	Prepare the simulation with cell linked list, configuration
      //	and case specified initial condition if necessary.
      //----------------------------------------------------------------------
      sph_system.initializeSystemCellLinkedLists();
      sph_system.initializeSystemConfigurations();
      water_block_complex.updateConfiguration();
      boundary_indicator.exec();
      left_emitter_inflow_injection.tag_buffer_particles.exec();
      right_up_emitter_inflow_injection.tag_buffer_particles.exec();
      right_down_emitter_inflow_injection.tag_buffer_particles.exec();
      wall_boundary_normal_direction.exec();
      //----------------------------------------------------------------------
      //	Setup for time-stepping control
      //----------------------------------------------------------------------
      Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
      size_t number_of_iterations = 0.0;
      int screen_output_interval = 100;
      int observation_sample_interval = screen_output_interval * 2;
      Real end_time = 2.5;               /**< End time. */
      Real Output_Time = end_time / 250; /**< Time stamps for output of body states. */
      Real dt = 0.0;                     /**< Default acoustic time step sizes. */
      Real dt_s = 0.0;                   /**< Default acoustic time step sizes for solid. */
      //----------------------------------------------------------------------
      //	Statistics for CPU time
      //----------------------------------------------------------------------
      TickCount t1 = TickCount::now();
      TimeInterval interval;
      TimeInterval interval_computing_time_step;
      TimeInterval interval_computing_pressure_relaxation;
      TimeInterval interval_updating_configuration;
      TickCount time_instance;
      //----------------------------------------------------------------------
      //	First output before the main loop.
      //----------------------------------------------------------------------
      body_states_recording.writeToFile();
      //----------------------------------------------------------------------
      //	Main loop starts here.
      //----------------------------------------------------------------------
      while (physical_time < end_time)
      {
          Real integration_time = 0.0;
          /** Integrate time (loop) until the next output time. */
          while (integration_time < Output_Time)
          {
              time_instance = TickCount::now();
              Real Dt = get_fluid_advection_time_step_size.exec();
              update_fluid_density.exec();
              viscous_acceleration.exec();
              transport_velocity_correction.exec();
              interval_computing_time_step += TickCount::now() - time_instance;
              time_instance = TickCount::now();
              Real relaxation_time = 0.0;
              while (relaxation_time < Dt)
              {
                  dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                  pressure_relaxation.exec(dt);
                  kernel_summation.exec();
                  left_inflow_pressure_condition.exec(dt);
                  right_up_inflow_pressure_condition.exec(dt);
                  right_down_inflow_pressure_condition.exec(dt);
                  inflow_velocity_condition.exec();
                  density_relaxation.exec(dt);
                  relaxation_time += dt;
                  integration_time += dt;
                  physical_time += dt;
              }
              interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
              {
                  std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                            << physical_time
                            << "	Dt = " << Dt << "	dt = " << dt << "\n";
                  body_states_recording.writeToFile();
              }
              number_of_iterations++;

              time_instance = TickCount::now();

              left_emitter_inflow_injection.injection.exec();
              right_up_emitter_inflow_injection.injection.exec();
              right_down_emitter_inflow_injection.injection.exec();

              left_disposer_outflow_deletion.exec();
              right_up_disposer_outflow_deletion.exec();
              right_down_disposer_outflow_deletion.exec();

              if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
              {
                  particle_sorting.exec();
              }
              water_block.updateCellLinkedList();
              water_block_complex.updateConfiguration();

              interval_updating_configuration += TickCount::now() - time_instance;
              boundary_indicator.exec();
              left_emitter_inflow_injection.tag_buffer_particles.exec();
              right_up_emitter_inflow_injection.tag_buffer_particles.exec();
              right_down_emitter_inflow_injection.tag_buffer_particles.exec();
          }
          TickCount t2 = TickCount::now();
          body_states_recording.writeToFile();
          TickCount t3 = TickCount::now();
          interval += t3 - t2;
      }
      TickCount t4 = TickCount::now();

      TimeInterval tt;
      tt = t4 - t1 - interval;
      std::cout << "Total wall time for computation: " << tt.seconds()
                << " seconds." << std::endl;
      std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
                << interval_computing_time_step.seconds() << "\n";
      std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
                << interval_computing_pressure_relaxation.seconds() << "\n";
      std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
                << interval_updating_configuration.seconds() << "\n";

    return 0;
}
