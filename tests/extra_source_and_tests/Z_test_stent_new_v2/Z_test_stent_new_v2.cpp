/**
 * @file     Z_test_3d_stent_v2.cpp
 * @brief 	 test 支架
 * @details  test
 *
 * @author 	Sukang Peng
 */

#include "Z_test_stent_new_v2.h"
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
    //sph_system.setRunParticleRelaxation(true); // Tag for run particle relaxation for body-fitted distribution
    //sph_system.setReloadParticles(false);      // Tag for computation with save particles distribution
     sph_system.setRunParticleRelaxation(false); // Tag for run particle relaxation for body-fitted distribution
     sph_system.setReloadParticles(true);        // Tag for computation with save particles distribution
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody stent_body(sph_system, makeShared<Stent>("Stent"));
    stent_body.defineAdaptationRatios(1.15, 7.0);
    stent_body.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    stent_body.defineMaterial<NeoHookeanSolid>(rho0_s_stent, youngs_modulus_stent, poisson_stent);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? stent_body.generateParticles<BaseParticles, Reload>(stent_body.getName())
        : stent_body.generateParticles<BaseParticles, Lattice>();

    BodyStatesRecordingToVtp write_stent_to_vtp(stent_body);
    write_stent_to_vtp.writeToFile(0);

    SolidBody vessel_wall(sph_system, makeShared<VesselWall>("VesselWall"));
    vessel_wall.defineAdaptationRatios(1.15, 1.5);
    vessel_wall.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    // vessel_wall.defineMaterial<Solid>();
    vessel_wall.defineMaterial<NeoHookeanSolid>(rho0_s_vessel, Youngs_modulus_vessel, poisson_vessel);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? vessel_wall.generateParticles<BaseParticles, Reload>(vessel_wall.getName())
        : vessel_wall.generateParticles<BaseParticles, Lattice>();

    BodyStatesRecordingToVtp write_vessel_wall_to_vtp(vessel_wall);
    write_vessel_wall_to_vtp.writeToFile(0);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation wall_relax_inner(vessel_wall);
        InnerRelation stent_relax_inner(stent_body);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_vessel_wall_particles(vessel_wall);
        SimpleDynamics<RandomizeParticlePosition> random_stent_particles(stent_body);
        RelaxationStepLevelSetCorrectionInner relaxation_step_wall_inner(wall_relax_inner);
        RelaxationStepLevelSetCorrectionInner relaxation_step_stent_inner(stent_relax_inner);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_wall_state_to_vtp(vessel_wall);
        BodyStatesRecordingToVtp write_stent_state_to_vtp(stent_body);
        /** Write the particle reload files. */
        ReloadParticleIO write_wall_particle_reload_files(vessel_wall);
        ReloadParticleIO write_stent_particle_reload_files(stent_body);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_vessel_wall_particles.exec(0.25);
        random_stent_particles.exec(0.25);
        relaxation_step_wall_inner.SurfaceBounding().exec();
        relaxation_step_stent_inner.SurfaceBounding().exec();
        write_wall_state_to_vtp.writeToFile(0.0);
        write_stent_state_to_vtp.writeToFile(0.0);
        //----------------------------------------------------------------------
        //	Relax particles of the vessel wall.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_wall_inner.exec();
            relaxation_step_stent_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the vessel wall N = " << ite_p << "\n";
                write_wall_state_to_vtp.writeToFile(ite_p);
                write_stent_state_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of vessel wall finish !" << std::endl;
        /** Output results. */
        write_wall_particle_reload_files.writeToFile(0);
        write_stent_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation stent_inner(stent_body);
    InnerRelation vessel_inner(vessel_wall);
    SurfaceContactRelation stent_vessel_contact(stent_body, {&vessel_wall});
    SurfaceContactRelation vessel_stent_contact(vessel_wall, {&stent_body});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    // RadialForce radial_force(23000.0, xAxis);
    StartupRadialForce radial_force(2000000.0, xAxis, 0.005);
    // 使用 SimpleDynamics 创建径向力应用对象
    // SimpleDynamics<RadialForceApplication<RadialForce>> apply_radial_force(stent_body, radial_force);
    SimpleDynamics<RadialForceApplication<StartupRadialForce>> apply_radial_force(stent_body, radial_force);

    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_stent(stent_inner);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_vessel(vessel_inner);
    /** active and passive stress relaxation. */
    Dynamics1Level<solid_dynamics::DecomposedIntegration1stHalf> stress_relaxation_first_half_stent(stent_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half_stent(stent_inner);
    Dynamics1Level<solid_dynamics::DecomposedIntegration1stHalf> stress_relaxation_first_half_vessel(vessel_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half_vessel(vessel_inner);
    /** Algorithms for stent-vessel contact. */
    InteractionDynamics<solid_dynamics::ContactFactorSummation> stent_update_contact_density(stent_vessel_contact);
    InteractionDynamics<solid_dynamics::ContactFactorSummation> vessel_update_contact_density(vessel_stent_contact);
    InteractionWithUpdate<solid_dynamics::ContactForce> stent_compute_solid_contact_forces(stent_vessel_contact);
    InteractionWithUpdate<solid_dynamics::ContactForce> vessel_compute_solid_contact_forces(vessel_stent_contact);

    /**Constrain  */
    BoundaryGeometry boundary_geometry(vessel_wall, "BoundaryGeometry", resolution_ref * 2.0);
    SimpleDynamics<FixBodyPartConstraint> constrain_holder(boundary_geometry);
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_mass_center_stent(stent_body);
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_mass_center_vessel(vessel_wall);

    /** Damping with the solid body*/
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> stent_damping(1.0, stent_inner, "Velocity", physical_viscosity_stent);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> vessel_damping(1.0, vessel_inner, "Velocity", physical_viscosity_vessel);

    ReduceDynamics<QuantitySummation<Real, SolidBody>> compute_total_mass_stent(stent_body, "Mass");
    ReduceDynamics<QuantityMassPosition<SolidBody>> compute_mass_position_stent(stent_body);
    Vecd mass_center_stent = compute_mass_position_stent.exec() / compute_total_mass_stent.exec();
    Matd moment_of_inertia_stent = Matd::Zero();

    // 计算惯性矩
    for (int i = 0; i != Dimensions; ++i)
    {
        for (int j = 0; j != Dimensions; ++j)
        {
            ReduceDynamics<QuantityMomentOfInertia<SolidBody>> compute_moment_of_inertia_stent(stent_body, mass_center_stent, i, j);
            moment_of_inertia_stent(i, j) = compute_moment_of_inertia_stent.exec();
        }
    }

    SimpleDynamics<Constrain3DSolidBodyRotation> constrain_rotation_stent(stent_body, mass_center_stent, moment_of_inertia_stent);

    ReduceDynamics<QuantitySummation<Real, SolidBody>> compute_total_mass_vessel(vessel_wall, "Mass");
    ReduceDynamics<QuantityMassPosition<SolidBody>> compute_mass_position_vessel(vessel_wall);
    Vecd mass_center_vessel = compute_mass_position_vessel.exec() / compute_total_mass_vessel.exec();
    Matd moment_of_inertia_vessel = Matd::Zero();

    // 计算惯性矩
    for (int i = 0; i != Dimensions; ++i)
    {
        for (int j = 0; j != Dimensions; ++j)
        {
            ReduceDynamics<QuantityMomentOfInertia<SolidBody>> compute_moment_of_inertia_vessel(vessel_wall, mass_center_vessel, i, j);
            moment_of_inertia_vessel(i, j) = compute_moment_of_inertia_vessel.exec();
        }
    }

    SimpleDynamics<Constrain3DSolidBodyRotation> constrain_rotation_vessel(vessel_wall, mass_center_vessel, moment_of_inertia_vessel);

    //    try
    //{
    //    // 计算总质量
    //    Real total_mass = compute_total_mass_.exec();
    //    std::cout << "Total mass of the stent: " << std::fixed << std::setprecision(12) << total_mass << " kg" << std::endl;

    //    // 计算质量加权的位置
    //    Vecd mass_weighted_position = compute_mass_position_.exec();
    //    std::cout << "Mass-weighted position of the stent: ("
    //              << mass_weighted_position[0] << ", "
    //              << mass_weighted_position[1] << ", "
    //              << mass_weighted_position[2] << ")" << std::endl;

    //    // 计算质心位置
    //    Vecd mass_center = mass_weighted_position / total_mass;
    //    std::cout << "Mass center of the stent: ("
    //              << mass_center[0] << ", "
    //              << mass_center[1] << ", "
    //              << mass_center[2] << ")" << std::endl;

    //    // 计算平均位置
    //    Vecd average_position = computeAveragePosition(stent_body);

    //    // 打印结果
    //    std::cout << "Average Position of Stent: ("
    //              << average_position[0] << ", "
    //              << average_position[1] << ", "
    //              << average_position[2] << ")" << std::endl;
    //}
    // catch (const std::exception &e)
    //{
    //    std::cerr << "Error during computation: " << e.what() << std::endl;
    //}

    ////----------------------------------------------------------------------
    ////  Prevent the console window from closing immediately
    ////----------------------------------------------------------------------
    // std::cout << "\nPress Enter to exit the program..." << std::endl;
    // std::cin.get(); // 等待用户按下回车键

    /** Update normal direction. */
    // SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> stent_update_normal_direction(stent_body);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    SimpleDynamics<VonMisesStress> vessel_stress(vessel_wall);
    write_states.addToWrite<Real>(vessel_wall, "VonMisesStress");
    SimpleDynamics<VonMisesStress> stent_stress(stent_body);
    write_states.addToWrite<Real>(stent_body, "VonMisesStress");
    ReloadParticleRecordingToXml write_particle_state(vessel_wall);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_stent_kinetic_energy(stent_body);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    corrected_configuration_stent.exec();
    corrected_configuration_vessel.exec();
    // apply_radial_force.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0.0;
    Real end_time = 0.1;
    Real output_period = end_time / 50.0;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    bool stop_simulation = false;
    while (physical_time < end_time && !stop_simulation)
    {
        Real integration_time = 0.0;

        while (integration_time < output_period && !stop_simulation)
        {
            if (ite % 50 == 0)
            {
                std::cout << "N=" << ite << "   Time: "
                          << physical_time << "\tdt: "
                          << dt << "\n";
                write_stent_kinetic_energy.writeToFile(ite);

                // 新增：输出当前径向力
                Vecd current_radial_force = radial_force.GetCurrentForce(physical_time);
                std::cout << "Current Radial Force: ("
                          << current_radial_force[0] << ", "
                          << current_radial_force[1] << ", "
                          << current_radial_force[2] << ")" << std::endl;

                BoundingBox current_bbox = getRealTimeBoundingBox(stent_body.getBaseParticles());
                BoundingBox aligned_bbox = getAlignedBoundingBox(stent_body.getBaseParticles(), Mat3d::Identity(), translation_stent);
                printBoundingBoxAndDelta(current_bbox);
                printBoundingBoxAndDelta(aligned_bbox);
                std::cout << "\n";

                vessel_stress.exec(dt);
                stent_stress.exec(dt);
                write_states.writeToFile();
            }

            apply_radial_force.exec(dt);

            // 接触模型计算
            stent_update_contact_density.exec();
            stent_compute_solid_contact_forces.exec();
            vessel_update_contact_density.exec();
            vessel_compute_solid_contact_forces.exec();

            /** Stress relaxation and damping. */
            stress_relaxation_first_half_stent.exec(dt);
            // constrain_rotation_stent.exec(dt);
            constrain_mass_center_stent.exec(dt);
            // stent_damping.exec(dt);
            // constrain_rotation_stent.exec(dt);
            // constrain_mass_center_stent.exec(dt);
            stress_relaxation_second_half_stent.exec(dt);

            stress_relaxation_first_half_vessel.exec(dt);
            // constrain_holder.exec(dt);
            // constrain_rotation_vessel.exec(dt);
            // constrain_mass_center_vessel.exec(dt);
            // vessel_damping.exec(dt);
            // constrain_rotation_vessel.exec(dt);
            // constrain_mass_center_vessel.exec(dt);
            // constrain_holder.exec(dt);
            stress_relaxation_second_half_vessel.exec(dt);

            ite++;
            dt = sph_system.getSmallestTimeStepAmongSolidBodies();
            integration_time += dt;
            physical_time += dt;

            stent_body.updateCellLinkedList();
            vessel_wall.updateCellLinkedList();

            stent_vessel_contact.updateConfiguration();
            vessel_stent_contact.updateConfiguration();

            // vessel_stress.exec(dt);
            // stent_stress.exec(dt);
            // write_states.writeToFile();

            // 计算当前的边界框
            BoundingBox current_bbox1 = getRealTimeBoundingBox(stent_body.getBaseParticles());
            BoundingBox aligned_bbox1 = getAlignedBoundingBox(stent_body.getBaseParticles(), Mat3d::Identity(), translation_stent);
            Vec3d delta = aligned_bbox1.second_ - aligned_bbox1.first_;

            // 检查 delta 是否已经达到阈值
            if (delta[1] >= 0.0036 || delta[2] >= 0.0036)
            {
                std::cout << "Delta reached 4 mm in Z or Y direction, stopping force application.\n";
                printBoundingBoxAndDelta(aligned_bbox1);
                // 设置标志位以提前结束模拟
                stop_simulation = true;
                // write_states.writeToFile();

                break;
            }
        }
        TickCount t2 = TickCount::now();
        vessel_stress.exec(dt);
        stent_stress.exec(dt);
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    write_particle_state.writeToFile(ite);
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    // 回归测试数据
    if (sph_system.GenerateRegressionData())
    {
        write_stent_kinetic_energy.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_stent_kinetic_energy.testResult();
    }

    return 0;
}