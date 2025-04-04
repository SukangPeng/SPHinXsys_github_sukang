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
std::string vessel_wall_cut_bottom = "./input/carotid_cut_bottom.stl";
std::string vessel_wall_cut_up_left = "./input/carotid_cut_up_left.stl";
std::string vessel_fluid = "./input/carotid_fluid_geo.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Vec3d translation_water_block(0.0, 0.0, 0.0);                             /**< 水块的初始平移，单位：m (米) */
Vec3d translation_wall_boundary(0.0, 0.0, 0.0);                           /**< 血管壁的初始平移，单位：m (米) */
Real length_scale_vessel = 1.0;                                   /**< 长度比例因子，无量纲 (无单位) */
Real scaling = pow(10, -3);                                                  /**< 长度比例因子，无量纲 (无单位) */
Real dp_0 = 0.3 * scaling;                                /**< 初始参考粒子间距，单位：m (米) */
Real wall_resolution = dp_0;                                            /*thickness = 1.0 * shell_resolution*/
Real BW = dp_0 * 4.0;                                           /**< 发射器的参考大小，单位：m (米) */
Real diameter = 0.004;                                                    /**< 血管外径，单位：m (米) */
Vec3d domain_lower_bound(-8.0 * scaling, -6.0 * scaling, -35.0 * scaling);
Vec3d domain_upper_bound(18.0 * scaling, 15.0 * scaling, 30.0 * scaling); /**< 系统域的上边界，单位：m (米) */
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

// inlet R=2.9293, (1.5611, 5.8559, -30.8885), (0.1034, -0.0458, 0.9935)
Real DW_in = 2.9293 * 2 * scaling;
Vec3d inlet_half = Vec3d(2.0 * dp_0, 3.5 * scaling, 3.5 * scaling);
Vec3d inlet_half_cut = Vec3d(2.0 * dp_0, 6.0 * scaling, 6.0 * scaling);
Vec3d inlet_normal(0.1034, -0.0458, 0.9935);
Vec3d inlet_cut_translation = Vec3d(1.5611, 5.8559, -30.8885) * scaling - inlet_normal * 2.0 * dp_0;
Vec3d inlet_buffer_translation = Vec3d(1.5611, 5.8559, -30.8885) * scaling + inlet_normal * 2.0 * dp_0;
RotationResult inlet_rotation_result = RotationCalculator(inlet_normal, standard_direction);
Rotation3d inlet_emitter_rotation(inlet_rotation_result.angle, inlet_rotation_result.axis);
Rotation3d inlet_disposer_rotation(inlet_rotation_result.angle + Pi, inlet_rotation_result.axis);

// outlet1 R=1.9416, (-2.6975, -0.4330, 21.7855), (-0.3160, -0.0009, 0.9488)
Real DW_out_up = 1.9416 * 2 * scaling;
Vec3d outlet_up_half = Vec3d(2.0 * dp_0, 2.4 * scaling, 2.4 * scaling);
Vec3d outlet_up_half_cut = Vec3d(2.0 * dp_0, 5.0 * scaling, 5.0 * scaling);
Vec3d outlet_up_normal(-0.3160, -0.0009, 0.9488);
Vec3d outlet_up_cut_translation = Vec3d(-2.6975, -0.4330, 21.7855) * scaling + outlet_up_normal * 2.0 * dp_0;
Vec3d outlet_up_buffer_translation = Vec3d(-2.6975, -0.4330, 21.7855) * scaling - outlet_up_normal * 2.0 * dp_0;
RotationResult outlet_up_rotation_result = RotationCalculator(outlet_up_normal, standard_direction);
Rotation3d outlet_up_disposer_rotation(outlet_up_rotation_result.angle, outlet_up_rotation_result.axis);
Rotation3d outlet_up_emitter_rotation(outlet_up_rotation_result.angle + Pi, outlet_up_rotation_result.axis);

// outlet2 R=1.3261, (9.0220, 0.9750, 18.6389), (-0.0399, 0.0693, 0.9972)
Real DW_out_down = 1.3261 * 2 * scaling;
Vec3d outlet_down_half = Vec3d(2.0 * dp_0, 2.0 * scaling, 2.0 * scaling);
Vec3d outlet_down_half_cut = Vec3d(2.0 * dp_0, 4.0 * scaling, 4.0 * scaling);
Vec3d outlet_down_normal(-0.0399, 0.0693, 0.9972);
Vec3d outlet_down_cut_translation = Vec3d(9.0220, 0.9750, 18.6389) * scaling + outlet_down_normal * 2.0 * dp_0;
Vec3d outlet_down_buffer_translation = Vec3d(9.0220, 0.9750, 18.6389) * scaling - outlet_down_normal * 2.0 * dp_0;
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
Real c_f = 10.0 * U_f * SMAX(Real(1), DW_in *DW_in / (DW_out_up * DW_out_up + DW_out_down * DW_out_down));
Real mu_f = 0.00355;      /**< Dynamics viscosity. */
Real Outlet_pressure = 0.0;  // for comparison with solely velocity inlet bc
//Real Outlet_pressure = 2.666e3;
//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real curent_time)
    {
        return p;
    }
};

struct RightInflowPressure
{
    template <class BoundaryConditionType>
    RightInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real curent_time)
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
        add<TriangleMeshShapeSTL>(vessel_fluid, translation_water_block, scaling);
    }
};

class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<ExtrudeShape<TriangleMeshShapeSTL>>(4.0 * dp_0, vessel_fluid, translation_wall_boundary, scaling);
        subtract<TriangleMeshShapeSTL>(vessel_fluid, translation_water_block, scaling);

        subtract<AlignedBoxShape>(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_cut_translation)), inlet_half_cut);
        subtract<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_up_emitter_rotation), Vec3d(outlet_up_cut_translation)), outlet_up_half_cut);
        subtract<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_down_emitter_rotation), Vec3d(outlet_down_cut_translation)), outlet_down_half_cut);
    }
};
//----------------------------------------------------------------------
//	BoundaryGeometry.
//----------------------------------------------------------------------
class BoundaryGeometry : public BodyPartByParticle
{
  public:
    BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry() {};

  private:
    void tagManually(size_t index_i)
    {
        if (base_particles_.ParticlePositions()[index_i][2] < -30.8885 * scaling + 0.9935 * 4 * dp_0 
            || (base_particles_.ParticlePositions()[index_i][2] > 21.7855 * scaling - 0.9488 * 4 * dp_0 
                && base_particles_.ParticlePositions()[index_i][0] <= 0) 
            || (base_particles_.ParticlePositions()[index_i][2] > 18.6389 * scaling - 0.9972 * 4 * dp_0 
                && base_particles_.ParticlePositions()[index_i][0] > 0))
        {
            body_part_particles_.push_back(index_i);
        }
    };
};