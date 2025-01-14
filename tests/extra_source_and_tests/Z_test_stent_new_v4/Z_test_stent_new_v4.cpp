/**
 * @file     Z_test_3d_stent_v2.cpp
 * @brief 	 test 支架
 * @details  test
 *
 * @author 	Sukang Peng
 */

#include "Z_test_stent_new_v4.h"
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
    stent_body.defineAdaptationRatios(1.15, 5.0);
    stent_body.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    stent_body.defineMaterial<NeoHookeanSolid>(rho0_s_stent, youngs_modulus_stent, poisson_stent);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? stent_body.generateParticles<BaseParticles, Reload>(stent_body.getName())
        : stent_body.generateParticles<BaseParticles, Lattice>();

    BodyStatesRecordingToVtp write_stent_to_vtp(stent_body);
    write_stent_to_vtp.writeToFile(0);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation stent_relax_inner(stent_body);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_stent_particles(stent_body);
        RelaxationStepLevelSetCorrectionInner relaxation_step_stent_inner(stent_relax_inner);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_stent_state_to_vtp(stent_body);
        /** Write the particle reload files. */
        ReloadParticleIO write_stent_particle_reload_files(stent_body);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_stent_particles.exec(0.25);
        relaxation_step_stent_inner.SurfaceBounding().exec();
        write_stent_state_to_vtp.writeToFile(0.0);
        //----------------------------------------------------------------------
        //	Relax particles of the vessel wall.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_stent_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the stent N = " << ite_p << "\n";
                write_stent_state_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of stent finish !" << std::endl;
        /** Output results. */
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
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    // RadialForce radial_force(23000.0, xAxis);
    StartupRadialForce radial_force(5000000.0, xAxis, 0.001);
    // 使用 SimpleDynamics 创建径向力应用对象
    // SimpleDynamics<RadialForceApplication<RadialForce>> apply_radial_force(stent_body, radial_force);
    SimpleDynamics<RadialForceApplication<StartupRadialForce>> apply_radial_force(stent_body, radial_force);

    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_stent(stent_inner);
    /** active and passive stress relaxation. */
    Dynamics1Level<solid_dynamics::DecomposedIntegration1stHalf> stress_relaxation_first_half_stent(stent_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half_stent(stent_inner);

    /** Damping with the solid body*/
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> stent_damping(0.2, stent_inner, "Velocity", physical_viscosity_stent);

    /**Constrain  */
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_mass_center_stent(stent_body);

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
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    SimpleDynamics<VonMisesStress> stent_stress(stent_body);
    write_states.addToWrite<Real>(stent_body, "VonMisesStress");
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_stent_kinetic_energy(stent_body);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    corrected_configuration_stent.exec();
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
                stent_stress.exec(dt);
                write_states.writeToFile();
            }

            apply_radial_force.exec(dt);

            /** Stress relaxation and damping. */
            stress_relaxation_first_half_stent.exec(dt);
            // constrain_rotation_stent.exec(dt);
            //constrain_mass_center_stent.exec(dt);
            //stent_damping.exec(dt);
            //constrain_rotation_stent.exec(dt);
            //constrain_mass_center_stent.exec(dt);
            stress_relaxation_second_half_stent.exec(dt);

            ite++;
            dt = sph_system.getSmallestTimeStepAmongSolidBodies();
            integration_time += dt;
            physical_time += dt;

            stent_body.updateCellLinkedList();

            // 计算当前的边界框
            BoundingBox current_bbox = getRealTimeBoundingBox(stent_body.getBaseParticles());
            BoundingBox aligned_bbox = getAlignedBoundingBox(stent_body.getBaseParticles(), Mat3d::Identity(), translation_stent);
            Vec3d delta = aligned_bbox.second_ - aligned_bbox.first_;

            // 检查 delta 是否已经达到阈值
            if (delta[1] >= 0.004 && delta[2] >= 0.004)
            {
                std::cout << "Delta reached 4 mm in Z or Y direction, stopping force application.\n";
                printBoundingBoxAndDelta(aligned_bbox);
                // 设置标志位以提前结束模拟
                stop_simulation = true;
                write_states.writeToFile();

                break;
            }
        }
        TickCount t2 = TickCount::now();
        stent_stress.exec(dt);
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
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