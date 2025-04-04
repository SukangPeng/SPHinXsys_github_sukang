/**
 * @file Z_test_3d_v8_vessel_solid_post_stent.cpp
 * @brief SPH-based simulation of post-stent hemodynamics in a stenotic coronary artery.
 * @details This file implements a high-fidelity numerical simulation of coronary blood flow
 *          after stent implantation using the Smoothed Particle Hydrodynamics (SPH) method
 *          within the SPHinXsys framework. 
 * @author Sukang Peng
 * @date March 2025
 */
#include "MA_coronary_with_stenosis_pre_stent.h"
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
    //wall_boundary.defineAdaptationRatios(1.15, 2.0);
    wall_boundary.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet();
    wall_boundary.defineMaterial<NeoHookeanSolid>(rho0_s_vessel, Youngs_modulus_vessel, poisson_vessel);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

    BodyStatesRecordingToVtp write_vessel_wall_to_vtp(wall_boundary);
    write_vessel_wall_to_vtp.writeToFile(0);

    ObserverBody stenosis_observer(sph_system, "StenosisObserver");
    Vecd oberver_point_stenosis_before = Vecd(15.60849846, 0.25177237, -4.89181985) * length_scale;
    Vecd oberver_point_stenosis_after = Vecd(20.36823748, 1.49065448, -14.73102995) * length_scale;
    StdVec<Vecd> observer_location = {oberver_point_stenosis_before, oberver_point_stenosis_after}; /**< Displacement observation point. */
    stenosis_observer.generateParticles<ObserverParticles>(observer_location);

    BodyStatesRecordingToVtp write_stenosis_observer_to_vtp(stenosis_observer);
    write_stenosis_observer_to_vtp.writeToFile(0);

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

    RealBody test_body_outlet_middle_buffer(
    sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_middle_emitter_rotation), Vec3d(outlet_middle_buffer_translation)), outlet_middle_half, "TestBodyOutletMiddleBuffer"));
    test_body_outlet_middle_buffer.generateParticles<BaseParticles, Lattice>();
    BodyStatesRecordingToVtp write_body_outlet_middle_buffer_to_vtp(test_body_outlet_middle_buffer);
    write_body_outlet_middle_buffer_to_vtp.writeToFile(0);

    RealBody test_body_outlet_small_buffer(
        sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_small_emitter_rotation), Vec3d(outlet_small_buffer_translation)), outlet_small_half, "TestBodyOutletSmallBuffer"));
    test_body_outlet_small_buffer.generateParticles<BaseParticles, Lattice>();
    BodyStatesRecordingToVtp write_body_outlet_small_buffer_to_vtp(test_body_outlet_small_buffer);
    write_body_outlet_small_buffer_to_vtp.writeToFile(0);
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
    InnerRelation wall_boundary_inner(wall_boundary);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelation wall_water_contact(wall_boundary, {&water_block});
    ContactRelation stenosis_observer_contact(stenosis_observer, {&water_block});
    //----------------------------------------------------------------------
    //  Combined relations built from basic relations
    //  which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);
    // solid dynamics
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> wall_corrected_configuration(wall_boundary_inner);

    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> wall_stress_relaxation_first_half(wall_boundary_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> wall_stress_relaxation_second_half(wall_boundary_inner);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> wall_computing_time_step_size(wall_boundary);

    // Constrain
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_mass_center(wall_boundary);

    ConstrainGeometryAlignedBox inlet_boundary(wall_boundary, "InletConstrain",
    AlignedBox(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_buffer_translation)), inlet_fix_half));

    ConstrainGeometryAlignedBox outlet_large_boundary(wall_boundary, "OutletLargeConstrain",
    AlignedBox(xAxis, Transform(Rotation3d(outlet_large_emitter_rotation), Vec3d(outlet_large_buffer_translation)), outlet_large_fix_half));

    ConstrainGeometryAlignedBox outlet_middle_boundary(wall_boundary, "OutletMiddleConstrain",
    AlignedBox(xAxis, Transform(Rotation3d(outlet_middle_emitter_rotation), Vec3d(outlet_middle_buffer_translation)), outlet_middle_fix_half));

    ConstrainGeometryAlignedBox outlet_small_boundary(wall_boundary, "OutletSmallConstrain",
    AlignedBox(xAxis, Transform(Rotation3d(outlet_small_emitter_rotation), Vec3d(outlet_small_buffer_translation)), outlet_small_fix_half));
    // Apply constraints
    SimpleDynamics<FixBodyPartConstraint> constrain_inlet(inlet_boundary);
    SimpleDynamics<FixBodyPartConstraint> constrain_outlet_large(outlet_large_boundary);
    SimpleDynamics<FixBodyPartConstraint> constrain_outlet_middle(outlet_middle_boundary);
    SimpleDynamics<FixBodyPartConstraint> constrain_outlet_small(outlet_small_boundary);

    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        wall_velocity_damping(0.2, wall_boundary_inner, "Velocity", physical_viscosity_vessel);
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> kernel_correction_complex(water_block_inner, water_wall_contact);
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_wall_contact);
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
    BodyAlignedBoxByCell outlet_middle_emitter(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_middle_emitter_rotation), Vec3d(outlet_middle_buffer_translation)), outlet_middle_half));
    fluid_dynamics::BidirectionalBuffer<OutletInflowPressure> outlet_middle_buffer(outlet_middle_emitter, in_outlet_particle_buffer);
    BodyAlignedBoxByCell outlet_small_emitter(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_small_emitter_rotation), Vec3d(outlet_small_buffer_translation)), outlet_small_half));
    fluid_dynamics::BidirectionalBuffer<OutletInflowPressure> outlet_small_buffer(outlet_small_emitter, in_outlet_particle_buffer);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_wall_contact);
    SimpleDynamics<fluid_dynamics::PressureCondition<InflowPressure>> inlet_inflow_pressure_condition(inlet_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<OutletInflowPressure>> outlet_large_inflow_pressure_condition(outlet_large_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<OutletInflowPressure>> outlet_middle_inflow_pressure_condition(outlet_middle_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<OutletInflowPressure>> outlet_small_inflow_pressure_condition(outlet_small_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(inlet_emitter);
    //----------------------------------------------------------------------
    // FlowRateCalculator
    //----------------------------------------------------------------------
    Vec3d inlet_flow_rate_translation = Vec3d(0, 0, 0) * length_scale + inlet_normal * 2.0 * resolution_ref;
    Rotation3d inlet_flow_rate_rotation = inlet_emitter_rotation;
    Vec3d inlet_flow_rate_half(2.0 * resolution_ref, 2.15 * length_scale, 2.15 * length_scale);
    BodyAlignedBoxByCell inlet_flow_rate_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(inlet_flow_rate_rotation), Vec3d(inlet_flow_rate_translation)), inlet_flow_rate_half));

    Vec3d outlet_large_flow_rate_translation = Vec3d(24.26695907, -9.65881011, -3.31277237) * length_scale - outlet_large_normal * 1.0 * resolution_ref;
    Rotation3d outlet_large_flow_rate_rotation = outlet_large_emitter_rotation;
    Vec3d outlet_large_flow_rate_half(1.0 * resolution_ref, 2.0 * length_scale, 2.0 * length_scale);
    BodyAlignedBoxByCell outlet_large_flow_rate_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_large_flow_rate_rotation), Vec3d(outlet_large_flow_rate_translation)), outlet_large_flow_rate_half));

    Vec3d outlet_middle_flow_rate_translation = Vec3d(27.00875603, 12.17813813, -28.45272815) * length_scale - outlet_middle_normal * 1.0 * resolution_ref;
    Rotation3d outlet_middle_flow_rate_rotation = outlet_small_disposer_rotation;
    Vec3d outlet_middle_flow_rate_half(1.0 * resolution_ref, 2.0 * length_scale, 2.0 * length_scale);
    BodyAlignedBoxByCell outlet_middle_flow_rate_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_middle_flow_rate_rotation), Vec3d(outlet_middle_flow_rate_translation)), outlet_middle_flow_rate_half));

    Vec3d outlet_small_flow_rate_translation = Vec3d(27.00875603, 12.17813813, -28.45272815) * length_scale - outlet_small_normal * 1.0 * resolution_ref;
    Rotation3d outlet_small_flow_rate_rotation = outlet_small_disposer_rotation;
    Vec3d outlet_small_flow_rate_half(1.0 * resolution_ref, 2.0 * length_scale, 2.0 * length_scale);
    BodyAlignedBoxByCell outlet_small_flow_rate_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_small_flow_rate_rotation), Vec3d(outlet_small_flow_rate_translation)), outlet_small_flow_rate_half));

    Vec3d stent_start_point = Vec3d(15.60849846, 0.25177237, -4.89181985) * length_scale;
    Vec3d stent_end_point = Vec3d(20.36823748, 1.49065448, -14.73102995) * length_scale;
    Vecd stent_normal = computeNormalDirectional(stent_start_point, stent_end_point);
    Real A_stenosis = M_PI * pow(3.3 * length_scale / 2.0, 2.0);

    Vec3d stenosis_start_flow_rate_translation = Vec3d(15.60849846, 0.25177237, -4.89181985) * length_scale;
    RotationCalculator stenosis_start_rotation_calculator(stent_normal, x_Axis);
    Rotation3d stenosis_start_rotation(stenosis_start_rotation_calculator.getRotationAngle(), stenosis_start_rotation_calculator.getRotationAxis());
    Vec3d stenosis_start_flow_rate_half(1.0 * resolution_ref, 2.0 * length_scale, 2.0 * length_scale);
    BodyAlignedBoxByCell stenosis_start_flow_rate_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(stenosis_start_rotation), Vec3d(stenosis_start_flow_rate_translation)), stenosis_start_flow_rate_half));

    Vec3d stenosis_end_flow_rate_translation = Vec3d(20.36823748, 1.49065448, -14.73102995) * length_scale;
    Rotation3d stenosis_end_rotation = stenosis_start_rotation;
    Vec3d stenosis_end_flow_rate_half(1.0 * resolution_ref, 2.0 * length_scale, 2.0 * length_scale);
    BodyAlignedBoxByCell stenosis_end_flow_rate_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(stenosis_end_rotation), Vec3d(stenosis_end_flow_rate_translation)), stenosis_end_flow_rate_half));

    ReduceDynamics<FlowRateCalculator> compute_inlet_flow_rate(inlet_flow_rate_box, inlet_normal, A_inlet, "Inlet");
    ReduceDynamics<FlowRateCalculator> compute_outlet_large_flow_rate(outlet_large_flow_rate_box, outlet_large_normal, A_outlet_large, "Outlet_large");
    ReduceDynamics<FlowRateCalculator> compute_outlet_middle_flow_rate(outlet_middle_flow_rate_box, outlet_middle_normal, A_outlet_middle, "Outlet_middle");
    ReduceDynamics<FlowRateCalculator> compute_outlet_small_flow_rate(outlet_small_flow_rate_box, outlet_small_normal, A_outlet_small, "Outlet_small");
    ReduceDynamics<FlowRateCalculator> compute_stenosis_start_flow_rate(stenosis_start_flow_rate_box, stent_normal, A_stenosis, "Stenosis_start");
    ReduceDynamics<FlowRateCalculator> compute_stenosis_end_flow_rate(stenosis_end_flow_rate_box, stent_normal, A_stenosis, "Stenosis_end");

    ReduceDynamics<PressureCalculator> compute_inlet_pressure(inlet_flow_rate_box, "Inlet_pre");
    ReduceDynamics<PressureCalculator> compute_stenosis_start_pressure(stenosis_start_flow_rate_box, "Stenosis_start");
    ReduceDynamics<PressureCalculator> compute_stenosis_end_pressure(stenosis_end_flow_rate_box, "Stenosis_end");
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(wall_boundary);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> wall_update_normal(wall_boundary);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_on_wall(wall_water_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(wall_water_contact);
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
    SimpleDynamics<VonMisesStress> vessel_stress(wall_boundary);
    body_states_recording.addToWrite<Real>(wall_boundary, "VonMisesStress");
    MaxStressCalculator max_stress_calculator(wall_boundary);
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_states_recording.addToWrite<Vecd>(wall_boundary, "PressureForceFromFluid");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_stenosis_velocity("Velocity", stenosis_observer_contact);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>> write_stenosis_pressure("Pressure", stenosis_observer_contact);
    RegressionTestTimeAverage<ReducedQuantityRecording<QuantitySummation<Vecd>>> write_total_viscous_force_from_fluid(wall_boundary, "ViscousForceFromFluid");
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
    wall_corrected_configuration.exec();
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
    write_stenosis_velocity.writeToFile(number_of_iterations);
    write_stenosis_pressure.writeToFile(number_of_iterations);
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

            /** FSI for viscous force. */
            viscous_force_on_wall.exec();
            /** Update normal direction on elastic body.*/
            // wall_update_normal.exec();

            interval_computing_time_step += TickCount::now() - time_instance;
            time_instance = TickCount::now();
            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;
            Real relaxation_time = 0.0;

            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
                pressure_force_from_fluid.exec();

                kernel_summation.exec();
                inlet_inflow_pressure_condition.exec(dt);
                outlet_large_inflow_pressure_condition.exec(dt);
                outlet_middle_inflow_pressure_condition.exec(dt);
                outlet_small_inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                /** Solid dynamics. */
                inner_ite_dt_s = 0.0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(wall_computing_time_step_size.exec(), dt - dt_s_sum);

                    wall_stress_relaxation_first_half.exec(dt_s);
                    constrain_inlet.exec();
                    constrain_outlet_large.exec();
                    constrain_outlet_middle.exec();
                    constrain_outlet_small.exec();
                    constrain_mass_center.exec();
                    wall_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    inner_ite_dt_s++;
                    // body_states_recording.writeToFile();
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                inner_ite_dt++;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "\tTime = "
                          << physical_time
                          << "\tDt = " << Dt << "\tdt = " << dt << "\tDt / dt = " << inner_ite_dt << "\tdt / dt_s = " << inner_ite_dt_s << "\n";

                // body_states_recording.writeToFile();

                // if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                //{
                //     write_point_velocity.writeToFile(number_of_iterations);
                // }
            }
            number_of_iterations++;

            time_instance = TickCount::now();

            inlet_buffer.injection.exec();
            outlet_large_buffer.injection.exec();
            outlet_middle_buffer.injection.exec();
            outlet_small_buffer.injection.exec();

            inlet_buffer.deletion.exec();
            outlet_large_buffer.deletion.exec();
            outlet_middle_buffer.deletion.exec();
            outlet_small_buffer.deletion.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }

            water_block.updateCellLinkedList();
            wall_update_normal.exec();
            wall_boundary.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            wall_water_contact.updateConfiguration();
            stenosis_observer_contact.updateConfiguration();

            interval_updating_configuration += TickCount::now() - time_instance;

            inlet_outlet_surface_particle_indicator.exec();
            inlet_buffer.tag_buffer_particles.exec();
            outlet_large_buffer.tag_buffer_particles.exec();
            outlet_middle_buffer.tag_buffer_particles.exec();
            outlet_small_buffer.tag_buffer_particles.exec();
        }
        TickCount t2 = TickCount::now();
        compute_inlet_flow_rate.exec(dt);
        compute_outlet_large_flow_rate.exec(dt);
        compute_outlet_middle_flow_rate.exec(dt);
        compute_outlet_small_flow_rate.exec(dt);
        compute_stenosis_start_flow_rate.exec(dt);
        compute_stenosis_end_flow_rate.exec(dt);
        compute_inlet_pressure.exec(dt);
        compute_stenosis_start_pressure.exec(dt);
        compute_stenosis_end_pressure.exec(dt);
        vessel_stress.exec(dt);
        max_stress_calculator.exec(dt);
        write_stenosis_pressure.writeToFile(number_of_iterations);
        write_stenosis_velocity.writeToFile(number_of_iterations);
        write_total_viscous_force_from_fluid.writeToFile(number_of_iterations);
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
