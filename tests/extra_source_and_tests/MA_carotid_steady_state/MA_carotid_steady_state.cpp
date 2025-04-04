/**
 * @file     carotid_steady_state_simulation.cpp
 * @brief    Simulation of steady-state blood flow in the carotid artery
 * @details  This simulation models blood flow within a carotid artery using SPH methods.
 * @author   Sukang Peng
 */
#include "MA_carotid_steady_state.h"
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
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(true); // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(false);      // Tag for computation with save particles distribution
    //sph_system.setRunParticleRelaxation(false); // Tag for run particle relaxation for body-fitted distribution
    //sph_system.setReloadParticles(true);        // Tag for computation with save particles distribution
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.cd
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet();
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticlesWithReserve<BaseParticles, Reload>(in_outlet_particle_buffer, water_block.getName())
        : water_block.generateParticles<BaseParticles, Lattice>();

     BodyStatesRecordingToVtp write_water_block_to_vtp(water_block);
     write_water_block_to_vtp.writeToFile(0);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
   // wall_boundary.defineAdaptationRatios(1.15, 2.0);
    wall_boundary.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet();
    wall_boundary.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

     BodyStatesRecordingToVtp write_vessel_wall_to_vtp(wall_boundary);
     write_vessel_wall_to_vtp.writeToFile(0);

     Vec3d observer_inlet_location = Vec3d(1.5611, 5.8559, -30.8885) * length_scale + inlet_normal * 5.0 * resolution_ref;
     Vec3d observer_outlet_large_location = Vec3d(-2.6975, -0.4330, 21.7855) * length_scale - outlet_large_normal * 5.0 * resolution_ref;
     Vec3d observer_outlet_small_location = Vec3d(9.0220, 0.9750, 18.6389) * length_scale - outlet_small_normal * 5.0 * resolution_ref;
     Vecd oberver_point1 = Vecd(2.51680152 * length_scale, 5.20993528 * length_scale, -20.95315931 * length_scale);
     StdVec<Vecd> observer_location = {observer_inlet_location, observer_outlet_large_location, observer_outlet_small_location, oberver_point1}; /**< Displacement observation point. */
     ObserverBody velocity_observer(sph_system, "VelocityObserver");
     velocity_observer.generateParticles<ObserverParticles>(observer_location);

     BodyStatesRecordingToVtp write_velocity_observer_to_vtp(velocity_observer);
     write_velocity_observer_to_vtp.writeToFile(0);

    // Cut Position Test
     RealBody test_body_inlet_cut(
     sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_cut_translation)), inlet_cut_half, "TestBodyInletCut"));
     test_body_inlet_cut.generateParticles<BaseParticles, Lattice>();
     BodyStatesRecordingToVtp write_body_inlet_cut_to_vtp(test_body_inlet_cut);
     write_body_inlet_cut_to_vtp.writeToFile(0);

     RealBody test_body_outlet_large_cut(
     sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_large_emitter_rotation), Vec3d(outlet_large_cut_translation)), outlet_large_cut_half, "TestBodyOutletLargeCut"));
     test_body_outlet_large_cut.generateParticles<BaseParticles, Lattice>();
     BodyStatesRecordingToVtp write_body_outlet_large_cut_to_vtp(test_body_outlet_large_cut);
     write_body_outlet_large_cut_to_vtp.writeToFile(0);

     RealBody test_body_outlet_small_cut(
     sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_small_emitter_rotation), Vec3d(outlet_small_cut_translation)), outlet_small_cut_half, "TestBodyOutletSmallCut"));
     test_body_outlet_small_cut.generateParticles<BaseParticles, Lattice>();
     BodyStatesRecordingToVtp write_body_outlet_small_cut_to_vtp(test_body_outlet_small_cut);
     write_body_outlet_small_cut_to_vtp.writeToFile(0);

     // Buffer Position Test
     RealBody test_body_inlet_buffer(
     sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_buffer_translation)), inlet_half, "TestBodyInletBuffer"));
     test_body_inlet_buffer.generateParticles<BaseParticles, Lattice>();
     BodyStatesRecordingToVtp write_body_inlet_buffer_to_vtp(test_body_inlet_buffer);
     write_body_inlet_buffer_to_vtp.writeToFile(0);

     RealBody test_body_outlet_large_buffer(
     sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_large_emitter_rotation), Vec3d(outlet_large_buffer_translation)), outlet_large_half, "TestBodyOutletLargeBuffer"));
     test_body_outlet_large_buffer.generateParticles<BaseParticles, Lattice>();
     BodyStatesRecordingToVtp write_body_outlet_large_buffer_to_vtp(test_body_outlet_large_buffer);
     write_body_outlet_large_buffer_to_vtp.writeToFile(0);

     RealBody test_body_outlet_small_buffer(
     sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_small_emitter_rotation), Vec3d(outlet_small_buffer_translation)), outlet_small_half, "TestBodyOutletSmallBuffer"));
     test_body_outlet_small_buffer.generateParticles<BaseParticles, Lattice>();
     BodyStatesRecordingToVtp write_body_outlet_small_buffer_to_vtp(test_body_outlet_small_buffer);
     write_body_outlet_small_buffer_to_vtp.writeToFile(0);

     // Delete Position Test
     BodyAlignedBoxByCell inlet_detection_box(wall_boundary, AlignedBox(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_delete_translation)), inlet_delete_half));
     BodyAlignedBoxByCell outlet_large_detection_box(wall_boundary, AlignedBox(xAxis, Transform(Rotation3d(outlet_large_emitter_rotation), Vec3d(outlet_large_delete_translation)), outlet_large_delete_half));
     BodyAlignedBoxByCell outlet_small_detection_box(wall_boundary, AlignedBox(xAxis, Transform(Rotation3d(outlet_small_emitter_rotation), Vec3d(outlet_small_delete_translation)), outlet_small_delete_half));

     RealBody test_body_inlet_delete(
     sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_delete_translation)), inlet_delete_half, "TestBodyInletDelete"));
     test_body_inlet_delete.generateParticles<BaseParticles, Lattice>();
     BodyStatesRecordingToVtp write_body_inlet_delete_to_vtp(test_body_inlet_delete);
     write_body_inlet_delete_to_vtp.writeToFile(0);

     RealBody test_body_outlet_large_delete(
     sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_large_emitter_rotation), Vec3d(outlet_large_delete_translation)), outlet_large_delete_half, "TestBodyOutletLargeDelete"));
     test_body_outlet_large_delete.generateParticles<BaseParticles, Lattice>();
     BodyStatesRecordingToVtp write_body_outlet_large_delete_to_vtp(test_body_outlet_large_delete);
     write_body_outlet_large_delete_to_vtp.writeToFile(0);

     RealBody test_body_outlet_small_delete(
     sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_small_emitter_rotation), Vec3d(outlet_small_delete_translation)), outlet_small_delete_half, "TestBodyOutletSmallDelete"));
     test_body_outlet_small_delete.generateParticles<BaseParticles, Lattice>();
     BodyStatesRecordingToVtp write_body_outlet_small_delete_to_vtp(test_body_outlet_small_delete);
     write_body_outlet_large_delete_to_vtp.writeToFile(0);

     Vec3d delete_half(3.0 * resolution_ref, 3.0 * resolution_ref, 3.0 * resolution_ref);
     BodyAlignedBoxByCell outlet_detection_box(wall_boundary, AlignedBox(xAxis, Transform(Vec3d(0.007, -0.00065, -0.0315)), delete_half));

     RealBody test_body_delete1(
     sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Vec3d(0.006, -0.00065, -0.0315)), delete_half, "TestBodyDelete1"));
     test_body_delete1.generateParticles<BaseParticles, Lattice>();
     BodyStatesRecordingToVtp write_body_delete1_to_vtp(test_body_delete1);
     write_body_delete1_to_vtp.writeToFile(0);

     SimpleDynamics<DeleteParticlesInAlignedBoxByCell> inlet_particles_detection(inlet_detection_box);
     SimpleDynamics<DeleteParticlesInAlignedBoxByCell> outlet_large_particles_detection(outlet_large_detection_box);
     SimpleDynamics<DeleteParticlesInAlignedBoxByCell> outlet_small_particles_detection(outlet_small_detection_box);
     SimpleDynamics<DeleteParticlesInAlignedBoxByCell> particles_detection1(outlet_detection_box);

     ParticleSorting particle_sorting_wall(wall_boundary);
     //----------------------------------------------------------------------
     //	SPH Particle relaxation section
     //----------------------------------------------------------------------
     /** check whether run particle relaxation for body fitted particle distribution. */
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

         inlet_particles_detection.exec();
         particle_sorting_wall.exec();
         wall_boundary.updateCellLinkedList();

         outlet_large_particles_detection.exec();
         particle_sorting_wall.exec();
         wall_boundary.updateCellLinkedList();

         outlet_small_particles_detection.exec();
         particle_sorting_wall.exec();
         wall_boundary.updateCellLinkedList();

         particles_detection1.exec();
         particle_sorting_wall.exec();
         wall_boundary.updateCellLinkedList();

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
     ContactRelation velocity_observer_contact(velocity_observer, {&water_block});
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
     InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);

     Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
     Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_wall_contact);
     InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_wall_contact);
     InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
     ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
     ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
     //----------------------------------------------------------------------
     // Definition Buffer
     //----------------------------------------------------------------------
     BodyAlignedBoxByCell inlet_emitter(water_block, AlignedBox(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_buffer_translation)), inlet_half));
     fluid_dynamics::BidirectionalBuffer<SPH::fluid_dynamics::NonPrescribedPressure> inlet_buffer(inlet_emitter, in_outlet_particle_buffer);
     BodyAlignedBoxByCell outlet_large_emitter(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_large_emitter_rotation), Vec3d(outlet_large_buffer_translation)), outlet_large_half));
     fluid_dynamics::BidirectionalBuffer<OutletInflowPressure> outlet_large_buffer(outlet_large_emitter, in_outlet_particle_buffer);
     BodyAlignedBoxByCell outlet_small_emitter(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_small_emitter_rotation), Vec3d(outlet_small_buffer_translation)), outlet_small_half));
     fluid_dynamics::BidirectionalBuffer<OutletInflowPressure> outlet_small_buffer(outlet_small_emitter, in_outlet_particle_buffer);

     InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_wall_contact);
     SimpleDynamics<fluid_dynamics::PressureCondition<InflowPressure>> inlet_inflow_pressure_condition(inlet_emitter);
     SimpleDynamics<fluid_dynamics::PressureCondition<OutletInflowPressure>> outlet_large_inflow_pressure_condition(outlet_large_emitter);
     SimpleDynamics<fluid_dynamics::PressureCondition<OutletInflowPressure>> outlet_small_inflow_pressure_condition(outlet_small_emitter);
     SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(inlet_emitter);
     //----------------------------------------------------------------------
     // FlowRateCalculator
     //----------------------------------------------------------------------
     Vec3d inlet_flow_rate_translation = Vec3d(1.5611, 5.8559, -30.8885) * length_scale + inlet_normal * 1.355 * resolution_ref;
     Rotation3d inlet_flow_rate_rotation = inlet_emitter_rotation;
     Vec3d inlet_flow_rate_half(1.355 * resolution_ref, 3.5 * length_scale, 3.5 * length_scale);
     BodyAlignedBoxByCell inlet_flow_rate_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(inlet_flow_rate_rotation), Vec3d(inlet_flow_rate_translation)), inlet_flow_rate_half));

     Vec3d outlet_large_flow_rate_translation = Vec3d(-2.6975, -0.4330, 21.7855) * length_scale - outlet_large_normal * 1.0 * resolution_ref;
     Rotation3d outlet_large_flow_rate_rotation = outlet_large_emitter_rotation;
     Vec3d outlet_large_flow_rate_half(1.0 * resolution_ref, 0.6 * length_scale, 0.6 * length_scale);
     BodyAlignedBoxByCell outlet_large_flow_rate_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_large_flow_rate_rotation), Vec3d(outlet_large_flow_rate_translation)), outlet_large_flow_rate_half));

     Vec3d outlet_small_flow_rate_translation = Vec3d(9.0220, 0.9750, 18.6389) * length_scale - outlet_small_normal * 1.0 * resolution_ref;
     Rotation3d outlet_small_flow_rate_rotation = outlet_small_disposer_rotation;
     Vec3d outlet_small_flow_rate_half(0.5 * resolution_ref, 0.5 * length_scale, 0.5 * length_scale);
     //Vec3d outlet_small_flow_rate_half(0.5 * dp_0, 2.0 * length_scale, 2.0 * length_scale);
     BodyAlignedBoxByCell outlet_small_flow_rate_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_small_flow_rate_rotation), Vec3d(outlet_small_flow_rate_translation)), outlet_small_flow_rate_half));

     RealBody test_body_inlet_flow(
     sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(inlet_flow_rate_rotation), Vec3d(inlet_flow_rate_translation)), inlet_flow_rate_half, "TestBodyFlowInlet"));
     test_body_inlet_flow.generateParticles<BaseParticles, Lattice>();
     BodyStatesRecordingToVtp write_body_inlet_flow_to_vtp(test_body_inlet_flow);
     write_body_inlet_flow_to_vtp.writeToFile(0);

     RealBody test_body_outlet_large_flow(
     sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_large_flow_rate_rotation), Vec3d(outlet_large_flow_rate_translation)), outlet_large_flow_rate_half, "TestBodyFlowOutletLarge"));
     test_body_outlet_large_flow.generateParticles<BaseParticles, Lattice>();
     BodyStatesRecordingToVtp write_body_outlet_large_flow_to_vtp(test_body_outlet_large_flow);
     write_body_outlet_large_flow_to_vtp.writeToFile(0);

     RealBody test_body_outlet_small_flow(
     sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_small_flow_rate_rotation), Vec3d(outlet_small_flow_rate_translation)), outlet_small_flow_rate_half, "TestBodyFlowOutletSmall"));
     test_body_outlet_small_flow.generateParticles<BaseParticles, Lattice>();
     BodyStatesRecordingToVtp write_body_outlet__small_flow_to_vtp(test_body_outlet_small_flow);
     write_body_outlet__small_flow_to_vtp.writeToFile(0);

     ReduceDynamics<FlowRateCalculator> compute_inlet_flow_rate(inlet_flow_rate_box, inlet_normal_flow, A_inlet, "inlet");
     ReduceDynamics<FlowRateCalculator> compute_outlet_large_flow_rate(outlet_large_flow_rate_box, outlet_large_normal_flow, A_outlet_large, "outlet_large");
     ReduceDynamics<FlowRateCalculator> compute_outlet_small_flow_rate(outlet_small_flow_rate_box, outlet_small_normal_flow, A_outlet_small, "outlet_small");
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
     body_states_recording.addToWrite<Vecd>(water_block, "Force");
     body_states_recording.addToWrite<int>(water_block, "BufferParticleIndicator");
     body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
     RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_point_velocity("Velocity", velocity_observer_contact);
     // RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_point_displacement("Position", velocity_observer_contact);
     //----------------------------------------------------------------------
     //	Prepare the simulation with cell linked list, configuration
     //	and case specified initial condition if necessary.
     //----------------------------------------------------------------------
     sph_system.initializeSystemCellLinkedLists();
     sph_system.initializeSystemConfigurations();
     water_block_complex.updateConfiguration();
     inlet_outlet_surface_particle_indicator.exec();
     inlet_buffer.tag_buffer_particles.exec();
     outlet_large_buffer.tag_buffer_particles.exec();
     outlet_small_buffer.tag_buffer_particles.exec();
     wall_boundary_normal_direction.exec();
     //----------------------------------------------------------------------
     //	Setup for time-stepping control
     //----------------------------------------------------------------------
     Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
     size_t number_of_iterations = 0.0;
     int screen_output_interval = 100;
     int observation_sample_interval = screen_output_interval * 2;
     Real end_time = 2.0;               /**< End time. */
     Real Output_Time = end_time / 200; /**< Time stamps for output of body states. */
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
     write_point_velocity.writeToFile(number_of_iterations);
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
                 inlet_inflow_pressure_condition.exec(dt);
                 outlet_large_inflow_pressure_condition.exec(dt);
                 outlet_small_inflow_pressure_condition.exec(dt);
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
                 //body_states_recording.writeToFile();

                 //if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                 //{
                 //    write_point_velocity.writeToFile(number_of_iterations);
                 //}
             }
             number_of_iterations++;

             time_instance = TickCount::now();

            inlet_buffer.injection.exec();
             outlet_large_buffer.injection.exec();
             outlet_small_buffer.injection.exec();

             inlet_buffer.deletion.exec();
             outlet_large_buffer.deletion.exec();
             outlet_small_buffer.deletion.exec();

             if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
             {
                 particle_sorting.exec();
             }
             water_block.updateCellLinkedList();
             water_block_complex.updateConfiguration();

             interval_updating_configuration += TickCount::now() - time_instance;

             inlet_outlet_surface_particle_indicator.exec();
             inlet_buffer.tag_buffer_particles.exec();
             outlet_large_buffer.tag_buffer_particles.exec();
             outlet_small_buffer.tag_buffer_particles.exec();
         }
         TickCount t2 = TickCount::now();
         compute_inlet_flow_rate.exec(dt);
         compute_outlet_large_flow_rate.exec(dt);
         compute_outlet_small_flow_rate.exec(dt);
         velocity_observer_contact.updateConfiguration();
         write_point_velocity.writeToFile(number_of_iterations);
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
