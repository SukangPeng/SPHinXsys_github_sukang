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
std::string stent = "./input/ml_stent_modified_0.0034.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Vec3d translation_water_block(0.0, 0.0, 0.0);                             /**< 水块的初始平移，单位：m (米) */
Vec3d translation_wall_boundary(0.0, 0.0, 0.0);                           /**< 血管壁的初始平移，单位：m (米) */
Vec3d translation_stent(0.0, 0.0, 0.0);                                   /**< 支架的初始平移，单位：m (米) */
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
//	Global parameters for vessel (血管壁参数)
//----------------------------------------------------------------------
/*Vessel Material*/
Real rho0_s = 1265;              /**< 血管壁的密度，单位：kg/m³ (千克每立方米) */
Real poisson = 0.45;             /**< 血管壁的泊松比，无量纲 (无单位) */
Real Youngs_modulus = 50000.0;   /**< 血管壁的杨氏模量，单位：Pa (帕) */
Real physical_viscosity = 500.0; /**< 血管壁的物理粘度，单位：Pa·s (帕·秒) */
//----------------------------------------------------------------------
//	Global parameters for stent (血管壁参数)
//----------------------------------------------------------------------
/* Stent Material */
Real rho0_s_stent = 6450.0;            /**< 支架的密度，单位：kg/m³ (千克每立方米) */
Real poisson_stent = 0.33;             /**< 支架的泊松比，无量纲 (无单位) */
Real youngs_modulus_stent = 1e6;       /**< 支架的杨氏模量，单位：Pa (帕) */
Real physical_viscosity_stent = 100.0; /**< 支架的物理粘度，单位：Pa·s (帕·秒) */
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
class Stent : public ComplexShape
{
  public:
    explicit Stent(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(stent, translation_stent, length_scale);
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

//----------------------------------------------------------------------
// RadialForce.
//----------------------------------------------------------------------
/**
 * @class RadialForce
 * @brief This class defines a radial force that can be used to simulate the expansion of a stent.
 */
class RadialForce : public Gravity
{
  protected:
    Vecd translation_;      // 支架的平移向量
    Real force_magnitude_;  // 力的大小
    Mat3d rotation_matrix_; // 旋转矩阵，将力从全局坐标转换到局部坐标

  public:
    RadialForce(const Vecd &translation, Real force_magnitude, const Mat3d &rotation_matrix)
        : Gravity(Vecd::Zero()), translation_(translation), force_magnitude_(force_magnitude), rotation_matrix_(rotation_matrix) {}

    virtual ~RadialForce() {}

    /** This function computes the induced radial acceleration along支架坐标系中的径向力 */
    virtual Vecd InducedAcceleration(const Vecd &position, Real physical_time = 0.0) const
    {
        // 1. 将全局坐标转换为支架的局部坐标系
        Vecd global_position = position - translation_;                       // 应用平移
        Vecd local_position = rotation_matrix_.transpose() * global_position; // 转换到局部坐标系

        // 2. 计算局部坐标系中的径向力
        Vecd radial_direction = local_position; // 获取径向分量
        radial_direction[2] = 0;                // 将Z轴分量置为0，确保径向力仅在XY平面

        // 3. 归一化径向力方向
        if (radial_direction.norm() > Eps)
        {
            radial_direction.normalize();
        }
        else
        {
            radial_direction = Vecd::Zero();
        }

        // 4. 施加径向力并将其转换回全局坐标系
        Vecd global_radial_direction = rotation_matrix_ * radial_direction;

        return force_magnitude_ * global_radial_direction;
    }

    /** 更新力的大小 */
    void setForceMagnitude(Real new_force_magnitude)
    {
        force_magnitude_ = new_force_magnitude;
    }

    /** 获取当前力的大小 */
    Real getForceMagnitude() const
    {
        return force_magnitude_;
    }
};

//----------------------------------------------------------------------
// RadialExpansionForce.
//----------------------------------------------------------------------
template <class GravityType>
class RadialExpansionForce : public GravityForce<GravityType>
{
  protected:
    RadialForce &radial_force_; // 径向力
    StdLargeVec<Vecd> &pos_;    // 粒子的位置
    StdLargeVec<Real> &mass_;   // 粒子的质量
    Real *physical_time_;       // 物理时间

  public:
    explicit RadialExpansionForce(SPHBody &sph_body, RadialForce &radial_force)
        : GravityForce<GravityType>(sph_body, radial_force), radial_force_(radial_force),
          pos_(*this->particles_->getVariableDataByName<Vecd>("Position")),
          mass_(*this->particles_->registerStateVariable<Real>("Mass")),
          physical_time_(this->sph_system_.getSystemVariableDataByName<Real>("PhysicalTime")) {}
    virtual ~RadialExpansionForce() {}

    void update(size_t index_i, Real dt = 0.0) override
    {
        // 计算径向力
        this->current_force_[index_i] = mass_[index_i] * radial_force_.InducedAcceleration(pos_[index_i], *physical_time_);
        // 调用 ForcePrior 的 update 方法
        GravityForce<GravityType>::update(index_i, dt);
    }
};




