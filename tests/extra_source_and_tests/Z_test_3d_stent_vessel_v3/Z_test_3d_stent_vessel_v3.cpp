/**
 * @file     Z_test_3d_stent_v2.cpp
 * @brief 	 test 支架
 * @details  test
 *
 * @author 	Sukang Peng
 */

#include "Z_test_3d_stent_vessel_v3.h"
#include "sphinxsys.h"
using namespace SPH; 
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);

#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody vessel_wall(sph_system, makeShared<VesselWall>("VesselWall"));
    vessel_wall.defineAdaptationRatios(1.15, 1.5);
    vessel_wall.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    //vessel_wall.defineMaterial<NeoHookeanSolid>(rho0_s_vessel, Youngs_modulus_vessel, poisson_vessel);
    vessel_wall.defineMaterial<Solid>();
    vessel_wall.generateParticles<BaseParticles, Reload>(vessel_wall.getName());
    // vessel_wall.generateParticles<BaseParticles, Lattice>();

    BodyStatesRecordingToVtp write_vessel_wall_to_vtp(vessel_wall);
    write_vessel_wall_to_vtp.writeToFile(0);

    //sph_system.setRunParticleRelaxation(true); // Tag for run particle relaxation for body-fitted distribution
    //sph_system.setReloadParticles(false);      // Tag for computation with save particles distribution
    sph_system.setRunParticleRelaxation(false); // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(true);      // Tag for computation with save particles distribution

    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(1.0);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticlesWithReserve<BaseParticles, Reload>(in_outlet_particle_buffer, water_block.getName())
        : water_block.generateParticles<BaseParticles, Lattice>();

    //BodyStatesRecordingToVtp write_water_block_to_vtp(water_block);
    //write_water_block_to_vtp.writeToFile(0);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation blood_inner(water_block);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_blood_particles(water_block);
        RelaxationStepInner relaxation_step_blood_inner(blood_inner);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_blood_state_to_vtp(water_block);
        /** Write the particle reload files. */
        ReloadParticleIO write_blood_particle_reload_files(water_block);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_blood_particles.exec(0.25);
        relaxation_step_blood_inner.SurfaceBounding().exec();
        write_blood_state_to_vtp.writeToFile(0.0);
        //----------------------------------------------------------------------
        //	Relax particles of the vessel wall.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_blood_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the vessel wall N = " << ite_p << "\n";
                write_blood_state_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of vessel wall finish !" << std::endl;
        /** Output results. */
        write_blood_particle_reload_files.writeToFile(0);
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
    ContactRelation water_wall_contact(water_block, {&vessel_wall});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(vessel_wall);
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);

    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_wall_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    //----------------------------------------------------------------------
    // Buffer
    //----------------------------------------------------------------------
    /*left buffer*/
    Vecd left_buffer_halfsize = Vecd(0.5 * BW, diameter / 2, diameter / 2);
    Vecd left_buffer_translation = Vecd(0.5 * BW - 0.5 * resolution_ref, 0.0, 0.0);
    AlignedBoxShape left_shape(xAxis, Transform(Vecd(left_buffer_translation)), left_buffer_halfsize);
    BodyAlignedBoxByCell left_emitter(water_block, left_shape);
    fluid_dynamics::BidirectionalBuffer<InflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);

    /*right buffer*/
    Vecd right_buffer_halfsize = Vecd(0.5 * BW, diameter / 2, diameter / 2);
    Vecd right_buffer_translation = Vecd(0.025 - 0.5 * BW, 0.0, 0.0);
    Real angle = M_PI;        // 180 度，弧度表示
    Vecd axis(0.0, 1.0, 0.0); // y 轴
    AlignedBoxShape right_shape(xAxis, Transform(Rotation3d(angle, axis), Vecd(right_buffer_translation)), right_buffer_halfsize);
    BodyAlignedBoxByCell right_emitter(water_block, right_shape);
    fluid_dynamics::BidirectionalBuffer<OutflowPressure> right_bidirection_buffer(right_emitter, in_outlet_particle_buffer);

    /*shape_test*/
    SolidBody Left_(sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Vecd(left_buffer_translation)), left_buffer_halfsize), "Left");
    Left_.defineMaterial<Solid>();
    Left_.generateParticles<BaseParticles, Lattice>();
    BodyStatesRecordingToVtp write_left_to_vtp(Left_);
    write_left_to_vtp.writeToFile(0);

    SolidBody Right_(sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Vecd(right_buffer_translation)), right_buffer_halfsize), "Right");
    Right_.defineMaterial<Solid>();
    Right_.generateParticles<BaseParticles, Lattice>();
    BodyStatesRecordingToVtp write_right_to_vtp(Right_);
    write_right_to_vtp.writeToFile(0);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_wall_contact);
    SimpleDynamics<fluid_dynamics::PressureCondition<InflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<OutflowPressure>> right_inflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_body_states(sph_system);
    write_body_states.addToWrite<Real>(water_block, "Pressure"); // output for debug
    write_body_states.addToWrite<int>(water_block, "Indicator"); // output for debug
    write_body_states.addToWrite<Real>(water_block, "Density");
    write_body_states.addToWrite<int>(water_block, "BufferParticleIndicator");
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_water_kinetic_energy(water_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    inlet_outlet_surface_particle_indicator.exec();
    left_bidirection_buffer.tag_buffer_particles.exec();
    right_bidirection_buffer.tag_buffer_particles.exec();
    wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0.0;
    int screen_output_interval = 500;
    Real end_time = 20.0;
    Real Output_Time = end_time / 200.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;                       /**< Default acoustic time step sizes. */
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
    write_body_states.writeToFile();
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();
            viscous_force.exec();
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
                right_inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();
                density_relaxation.exec(dt);
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                //write_body_states.writeToFile();

            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            time_instance = TickCount::now();

            left_bidirection_buffer.injection.exec();
            right_bidirection_buffer.injection.exec();

            left_bidirection_buffer.deletion.exec();
            right_bidirection_buffer.deletion.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();

            interval_updating_configuration += TickCount::now() - time_instance;
            inlet_outlet_surface_particle_indicator.exec();
            left_bidirection_buffer.tag_buffer_particles.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();
        }
        TickCount t2 = TickCount::now();
        write_body_states.writeToFile();
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