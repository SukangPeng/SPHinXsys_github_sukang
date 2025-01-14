/**
 * @file Z_test_3d_straight_v4_E_VI
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
std::string vessel_wall = "./input/carotid_wall_geo.stl";
std::string vessel_fluid = "./input/carotid_fluid_geo.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Vec3d translation_water_block(0.0, 0.0, 0.0);                             /**< 水块的初始平移，单位：m (米) */
Vec3d translation_wall_boundary(0.0, 0.0, 0.0);                           /**< 血管壁的初始平移，单位：m (米) */
Real length_scale_vessel = 1.0;                                   /**< 长度比例因子，无量纲 (无单位) */
Real length_scale = pow(10, -3);                                                  /**< 长度比例因子，无量纲 (无单位) */
Real resolution_ref = 0.2 * length_scale;                                /**< 初始参考粒子间距，单位：m (米) */
Real BW = resolution_ref * 4.0;                                           /**< 发射器的参考大小，单位：m (米) */
Real diameter = 0.004;                                                    /**< 血管外径，单位：m (米) */
Vec3d domain_lower_bound(-7.0 * length_scale, -4.0 * length_scale, -35.0 * length_scale);
Vec3d domain_upper_bound(20.0 * length_scale, 12.0 * length_scale, 30.0 * length_scale); /**< 系统域的上边界，单位：m (米) */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound); /**< 系统域的边界框，单位：m (米) */
Real full_length = 0.025;                                                 /**< 血管的总长度，单位：m (米) */
//----------------------------------------------------------------------
//	Buffer location.
//----------------------------------------------------------------------
struct RotationResult
{
    Vec3d axis;
    Real angle;
};

RotationResult RotationCalculator(Vecd target_normal, Vecd standard_direction)
{
    target_normal.normalize();

    Vec3d axis = standard_direction.cross(target_normal);
    Real angle = std::acos(standard_direction.dot(target_normal));

    if (axis.norm() < 1e-6)
    {
        if (standard_direction.dot(target_normal) < 0)
        {
            axis = Vec3d(1, 0, 0);
            angle = M_PI;
        }
        else
        {
            axis = Vec3d(0, 0, 1);
            angle = 0;
        }
    }
    else
    {
        axis.normalize();
    }

    return {axis, angle};
}

Vecd standard_direction(1, 0, 0);

// inlet R=2.9293, (1.5611, 5.8559, -30.8885), (-0.1034, 0.0458, -0.9935)
Real DW_in = 2.9293 * 2 * length_scale;
Vecd inlet_buffer_half = Vecd(2.0 * resolution_ref, 3.2 * length_scale, 3.2 * length_scale);
Vecd inlet_normal(0.1034, -0.0458, 0.9935);
Vecd inlet_buffer_translation = Vecd(1.5611, 5.8559, -30.8885) * length_scale + inlet_normal * 2.0 * resolution_ref;
RotationResult inlet_rotation_result = RotationCalculator(inlet_normal, standard_direction);
Rotation3d inlet_emitter_rotation(inlet_rotation_result.angle, inlet_rotation_result.axis);
Rotation3d inlet_disposer_rotation(inlet_rotation_result.angle + Pi, inlet_rotation_result.axis);

// outlet1 R=1.9416, (-2.6975, -0.4330, 21.7855), (-0.3160, -0.0009, 0.9488)
Real DW_up = 1.9416 * 2 * length_scale;
Vecd outlet_up_buffer_half = Vecd(2.0 * resolution_ref, 2.4 * length_scale, 2.4 * length_scale);
Vecd outlet_up_normal(-0.3160, -0.0009, 0.9488);
Vecd outlet_up_buffer_translation = Vecd(-2.6975, -0.4330, 21.7855) * length_scale - outlet_up_normal * 2.0 * resolution_ref;
RotationResult outlet_up_rotation_result = RotationCalculator(outlet_up_normal, standard_direction);
Rotation3d outlet_up_disposer_rotation(outlet_up_rotation_result.angle, outlet_up_rotation_result.axis);
Rotation3d outlet_up_emitter_rotation(outlet_up_rotation_result.angle + Pi, outlet_up_rotation_result.axis);

// outlet2 R=1.2760, (9.0465, 1.down552, 18.6363), (-0.0417, 0.0701, 0.9967)
Real DW_down = 1.2760 * 2 * length_scale;
Vecd outlet_down_buffer_half = Vecd(2.0 * resolution_ref, 1.5 * length_scale, 1.5 * length_scale);
Vecd outlet_down_normal(-0.0417, 0.0701, 0.9967);
Vecd outlet_down_buffer_translation = Vecd(9.0465, 1.02552, 18.6363) * length_scale - outlet_down_normal * 2.0 * resolution_ref;
RotationResult outlet_down_rotation_result = RotationCalculator(outlet_down_normal, standard_direction);
Rotation3d outlet_down_disposer_rotation(outlet_down_rotation_result.angle, outlet_down_rotation_result.axis);
Rotation3d outlet_down_emitter_rotation(outlet_down_rotation_result.angle + Pi, outlet_down_rotation_result.axis);
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1060;   /**< Reference density of fluid. */
Real U_f = 0.5;       /**< Characteristic velocity. */
Real U_max = 2 * U_f; /**< Characteristic velocity. */
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f * SMAX(Real(1), DW_in *DW_in / (DW_up * DW_up + DW_down * DW_down));
Real mu_f = 0.00355; /**< Dynamics viscosity. */
// Real Outlet_pressure = 0;  // for comparison with solely velocity inlet bc
// Real Outlet_pressure = 1.33e4;
Real Outlet_pressure = 2.666e3;
//----------------------------------------------------------------------
//	Global parameters on the solid properties (血管壁参数)
//----------------------------------------------------------------------
Real rho0_s = 1265;               /**< 血管壁的密度，单位：kg/m³ (千克每立方米) */
Real poisson = 0.45;              /**< 血管壁的泊松比，无量纲 (无单位) */
Real Youngs_modulus = 1e5;     /**< 血管壁的杨氏模量，单位：Pa (帕) */
Real physical_viscosity = 1000.0; /**< 血管壁的物理粘度，单位：Pa·s (帕·秒) */
//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        return p_;
    }
};

struct RightInflowPressure
{
    template <class BoundaryConditionType>
    RightInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        // Real run_time = GlobalStaticVariables::physical_time_;

        /*constant pressure*/
        Real pressure = Outlet_pressure;
        // return run_time < 0.5 ? 0.0: pressure;
        return pressure;
    }
};
//----------------------------------------------------------------------
// Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_, interval_;
    AlignedBoxShape &aligned_box_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(0.1), t_ref_(0.218), interval_(0.5),
          aligned_box_(boundary_condition.getAlignedBox()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;
        Real run_time = current_time; // 使用传入的 current_time 代替 GlobalStaticVariables::physical_time_
        int n = static_cast<int>(run_time / interval_);
        Real t_in_cycle = run_time - n * interval_;

        target_velocity[0] = t_in_cycle < t_ref_ ? 0.5 * sin(4 * Pi * (run_time + 0.0160236)) : u_ref_;
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

class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
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