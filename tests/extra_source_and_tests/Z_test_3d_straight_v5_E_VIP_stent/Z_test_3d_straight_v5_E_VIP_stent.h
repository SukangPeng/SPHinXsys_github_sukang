/**
 * @file Z_test_3d_straight_v5_E_VIP_stent
 * @brief
 * @details
 * @author Sukang
 */

#pragma once
#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"
using namespace SPH;
constexpr Real PI = 3.14159265358979323846;
//----------------------------------------------------------------------
//	Set the file path to the data file
//----------------------------------------------------------------------
std::string vessel_wall = "./input/vessel_wall.stl";
std::string vessel_fluid = "./input/vessel_fluid.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Vec3d translation_water_block(0.0, 0.0, 0.0);                             /**< 水块的初始平移，单位：m (米) */
Vec3d translation_wall_boundary(0.0, 0.0, 0.0);                           /**< 血管壁的初始平移，单位：m (米) */
Real length_scale = 1.0;                                                  /**< 长度比例因子，无量纲 (无单位) */
Real resolution_ref = 0.00035;                                            /**< 初始参考粒子间距，单位：m (米) */
Real BW = resolution_ref * 4.0;                                           /**< 发射器的参考大小，单位：m (米) */
Real diameter = 0.004;                                                    /**< 血管外径，单位：m (米) */
Vec3d domain_lower_bound(-0.001, -0.003, -0.003);                         /**< 系统域的下边界，单位：m (米) */
Vec3d domain_upper_bound(0.026, 0.003, 0.003);                            /**< 系统域的上边界，单位：m (米) */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound); /**< 系统域的边界框，单位：m (米) */
Real full_length = 0.025;                                                 /**< 血管的总长度，单位：m (米) */
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1060.0; /**< 流体的参考密度，单位：kg/m³ (千克每立方米) */
Real Outlet_pressure = 0.0;
Real Re = 100.0; /**< 雷诺数，无量纲 (无单位) */
Real U_f = 0.1;  /**< Characteristic velocity. */
Real mu_f = rho0_f * U_f * diameter / Re;
Real c_f = 10.0 * U_f;
// Real mu_f = 3.6e-3;   /**< 动态粘度，单位：Pa·s (帕·秒) */
// const Real U_f = Re * mu_f / rho0_f / diameter;
// const Real U_max = 2.0 * U_f;  // parabolic inflow, Thus U_max = 2*U_f
// const Real c_f = 10.0 * U_max; /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties (血管壁参数)
//----------------------------------------------------------------------
Real rho0_s = 1265;              /**< 血管壁的密度，单位：kg/m³ (千克每立方米) */
Real poisson = 0.45;             /**< 血管壁的泊松比，无量纲 (无单位) */
Real Youngs_modulus = 50000.0;   /**< 血管壁的杨氏模量，单位：Pa (帕) */
Real physical_viscosity = 500.0; /**< 血管壁的物理粘度，单位：Pa·s (帕·秒) */
//----------------------------------------------------------------------
// Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(2.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = Vecd::Zero(); // 初始为零，确保边界条件不会受到上一步速度影响
        Real run_time = current_time;

        // 改进的速度平滑函数
        Real u_ave = run_time < t_ref_ ? u_ref_ * (1.0 - exp(-run_time / t_ref_)) : u_ref_;
        target_velocity[0] = u_ave; // 假设 x 方向是主要流入方向

        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Define SPH bodies.
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(vessel_fluid, translation_water_block, length_scale);
    }
};
class VesselWall : public ComplexShape
{
  public:
    explicit VesselWall(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(vessel_wall, translation_wall_boundary, length_scale);
    }
};
//----------------------------------------------------------------------
//	BoundaryGeometry.
//----------------------------------------------------------------------
class BoundaryGeometry : public BodyPartByParticle
{
  public:
    BoundaryGeometry(SPHBody &body, const std::string &body_part_name, Real constrain_len)
        : BodyPartByParticle(body, body_part_name), constrain_len_(constrain_len)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry(){};

  private:
    Real constrain_len_;

    void tagManually(size_t index_i)
    {
        if (base_particles_.ParticlePositions()[index_i][0] < constrain_len_ || base_particles_.ParticlePositions()[index_i][0] > full_length - constrain_len_)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};