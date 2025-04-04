/**
 * @file     Z_test_3d_stent_v2.cpp
 * @brief 	 test 支架
 * @details  test
 *
 * @author 	Sukang Peng
 */

#include "Z_test_stent_v2.h"
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
    stent_body.defineAdaptationRatios(1.15, 3.5);
    stent_body.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet();
    //stent_body.defineMaterial<NeoHookeanSolid>(rho0_s_stent, youngs_modulus_stent, poisson_stent);
    stent_body.defineMaterial<NonLinearHardeningPlasticSolid>(
        rho0_s_stent,
        youngs_modulus_stent,
        poisson_stent,
        yield_stress_stent,           // 屈服应力
        hardening_modulus_stent,      // 线性硬化模量
        saturation_flow_stress_stent, // 饱和流动应力
        saturation_exponent_stent     // 硬化指数
    );
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? stent_body.generateParticles<BaseParticles, Reload>(stent_body.getName())
        : stent_body.generateParticles<BaseParticles, Lattice>();

    BodyStatesRecordingToVtp write_stent_to_vtp(stent_body);
    write_stent_to_vtp.writeToFile(0);

    SolidBody vessel_wall(sph_system, makeShared<VesselWall>("VesselWall"));
    //vessel_wall.defineAdaptationRatios(1.15, 2.0);
    vessel_wall.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet();
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
    Vec3d stent_start_point = Vec3d(15.60849846, 0.25177237, -4.89181985) * length_scale;
    Vec3d stent_end_point = Vec3d(20.36823748, 1.49065448, -14.73102995) * length_scale;
    // Vec3d stent_start_point = Vec3d(1,0,0) * length_scale;
    // Vec3d stent_end_point = Vec3d(-1,0,0) * length_scale;

    // 计算支架中心点
    Vec3d stent_center_point = (stent_start_point + stent_end_point) * 0.5;

    // 1. 初始化径向力
    Real radial_force_magnitude = 4.0e5; // N/mass
    int radial_force_axis = 0;         // x 轴
    Vec3d initial_direction(1.0, 0.0, 0.0);
    RotationCalculator rotation_calculator(stent_start_point, stent_end_point, initial_direction);
    Mat3d rotation_matrix = rotation_calculator.getRotationMatrix();
    RadialForce radial_force(radial_force_magnitude, stent_start_point, stent_end_point);

    // 2. 计算旋转矩阵并应用
    SimpleDynamics<RadialForceApplication<RadialForce>> apply_radial_force(stent_body, radial_force);

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
    // Constrain Stent
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_mass_center_stent(stent_body);
    /** Damping with the solid body*/
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> stent_damping(1.0, stent_inner, "Velocity", physical_viscosity_stent);
    // Rotation constrain stent
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

    // Constrain Vessel
    BoundaryGeometry boundary_geometry(vessel_wall, "BoundaryGeometry", resolution_ref * 2.0);
    SimpleDynamics<FixBodyPartConstraint> constrain_holder(boundary_geometry);
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_mass_center_vessel(vessel_wall);
    // Damping with the solid body
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> vessel_damping(0.5, vessel_inner, "Velocity", physical_viscosity_vessel);
    // Rotation constrain vessel
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
    Real output_period = end_time / 100.0;
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
                std::cout << "N=" << ite << " Time: " << physical_time << " dt: " << dt << "\n";
                write_stent_kinetic_energy.writeToFile(ite);

                // 计算 **全局 Bounding Box** 并 **打印**
                BoundingBoxCalculator::getGlobalBoundingBox(stent_body.getBaseParticles(), true);

                // 计算 **支架自身坐标系 Bounding Box** 并 **打印**
                BoundingBoxCalculator::getLocalBoundingBox(stent_body.getBaseParticles(), rotation_matrix, stent_center_point, true);

                // 打印施加在支架上的径向力（支架坐标系下）
                apply_radial_force.printAppliedForce();
                
                vessel_stress.exec(dt);
                stent_stress.exec(dt);
                write_states.writeToFile();
            }
            apply_radial_force.exec(dt);

            // 接触模型计算
            stent_update_contact_density.exec();
            //stent_update_contact_density.scaleOffsetW(1.5);
            stent_compute_solid_contact_forces.exec();
            vessel_update_contact_density.exec();
            //vessel_update_contact_density.scaleOffsetW(1.5);
            vessel_compute_solid_contact_forces.exec();

            /** Stress relaxation and damping. */
            stress_relaxation_first_half_stent.exec(dt);
            //constrain_rotation_stent.exec(dt);
            constrain_mass_center_stent.exec(dt);
            //stent_damping.exec(dt);
            //constrain_rotation_stent.exec(dt);
            //constrain_mass_center_stent.exec(dt);
            stress_relaxation_second_half_stent.exec(dt);

            stress_relaxation_first_half_vessel.exec(dt);
            // constrain_holder.exec(dt);
            // constrain_rotation_vessel.exec(dt);
            constrain_mass_center_vessel.exec(dt);
            //vessel_damping.exec(dt);
            // constrain_rotation_vessel.exec(dt);
            //constrain_mass_center_vessel.exec(dt);
            // constrain_holder.exec(dt);
            stress_relaxation_second_half_vessel.exec(dt);

            ite++;
            //dt = sph_system.getSmallestTimeStepAmongSolidBodies();
            dt = 0.5 * sph_system.getSmallestTimeStepAmongSolidBodies();

            integration_time += dt;
            physical_time += dt;

            stent_body.updateCellLinkedList();
            vessel_wall.updateCellLinkedList();

            stent_vessel_contact.updateConfiguration();
            vessel_stent_contact.updateConfiguration();

            // 计算全局 Bounding Box（只计算，不打印）
            BoundingBox global_bbox_test = BoundingBoxCalculator::getGlobalBoundingBox(stent_body.getBaseParticles(), false);
            // 计算支架自身坐标系下的 Bounding Box（只计算，不打印）
            BoundingBox local_bbox_test = BoundingBoxCalculator::getLocalBoundingBox(stent_body.getBaseParticles(), rotation_matrix, stent_center_point, false);

            // 获取 Δy 和 Δz
            Real delta_y = local_bbox_test.second_[1] - local_bbox_test.first_[1];
            Real delta_z = local_bbox_test.second_[2] - local_bbox_test.first_[2];
            // 检测是否达到了支架扩张条件
            if (delta_y >= 3.5 * length_scale || delta_z >= 3.5 * length_scale)
            {
                std::cout << "========================================\n";
                std::cout << "         Stent Expansion Reached!\n";
                std::cout << "----------------------------------------\n";
                std::cout << "  Δy = " << delta_y << " (Threshold: 3.5)\n";
                std::cout << "  Δz = " << delta_z << " (Threshold: 3.5)\n";
                std::cout << "========================================\n\n";

                stent_stress.exec(dt);
                write_states.writeToFile();
                stop_simulation = true;
                break; // 退出当前 while 循环
            }
        }
        TickCount t2 = TickCount::now();
        vessel_stress.exec(dt);
        stent_stress.exec(dt);
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
        // 直接跳出外层循环
        // goto END_SIMULATION;
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
    // END_SIMULATION:
    //  std::cout << "BoundingBox computation complete, exiting simulation.\n\n";
     return 0;
}