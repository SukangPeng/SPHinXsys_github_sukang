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
//	RotationCalculator
//----------------------------------------------------------------------
struct RotationCalculator
{
    Vec3d rotation_axis;   // 旋转轴
    Real rotation_angle;   // 旋转角度（弧度）
    Mat3d rotation_matrix; // 旋转矩阵

    // 1. 从两个点计算旋转（不影响 inlet/outlet，但提供备用方法）
    RotationCalculator(const Vec3d &pointA, const Vec3d &pointB, const Vec3d &initial_direction)
    {
        Vec3d target_direction = (pointB - pointA).normalized(); // 目标方向向量
        initialize(target_direction, initial_direction);
    }

    // 2. 从方向向量计算旋转（用于 inlet/outlet）
    RotationCalculator(const Vec3d &target_direction, const Vec3d &initial_direction)
    {
        initialize(target_direction.normalized(), initial_direction.normalized());
    }

    // 统一初始化逻辑
    void initialize(const Vec3d &target, const Vec3d &initial)
    {
        // 计算旋转轴
        rotation_axis = initial.cross(target);

        // 计算旋转角度
        rotation_angle = std::acos(initial.dot(target));

        // 处理旋转轴数值稳定性问题
        if (rotation_axis.norm() < 1e-6)
        {
            if (initial.dot(target) < 0)
            {
                rotation_axis = Vec3d(1, 0, 0); // 取X轴作为默认旋转轴
                rotation_angle = M_PI;
            }
            else
            {
                rotation_axis = Vec3d(0, 0, 1); // 保持不变
                rotation_angle = 0;
            }
        }
        else
        {
            rotation_axis.normalize();
        }

        // 计算旋转矩阵
        rotation_matrix = Eigen::AngleAxis<Real>(rotation_angle, rotation_axis).toRotationMatrix();
    }

    // 打印旋转矩阵
    void printRotationMatrix() const
    {
        std::cout << "Rotation Matrix:\n"
                  << rotation_matrix << std::endl;
    }

    // 打印旋转轴
    void printRotationAxis() const
    {
        std::cout << "Rotation Axis: (" << rotation_axis.x() << ", " << rotation_axis.y() << ", " << rotation_axis.z() << ")" << std::endl;
    }

    //  打印旋转角度
    void printRotationAngle() const
    {
        std::cout << "Rotation Angle (in radians): " << rotation_angle << std::endl;
    }

    // 获取旋转矩阵
    Mat3d getRotationMatrix() const { return rotation_matrix; }

    // 获取旋转轴
    Vec3d getRotationAxis() const { return rotation_axis; }

    // 获取旋转角度
    Real getRotationAngle() const { return rotation_angle; }
};

// X轴标准方向
Vec3d x_Axis(1, 0, 0);

// inlet 参数 inlet R=2.9293, (1.5611, 5.8559, -30.8885), (0.1034, -0.0458, 0.9935)
Real DW_in = 2.9293 * 2 * scaling;
Vec3d inlet_half(2.0 * dp_0, 3.5 * scaling, 3.5 * scaling);
Vec3d inlet_half_cut(2.0 * dp_0, 6.0 * scaling, 6.0 * scaling);
Vec3d inlet_normal(0.1034, -0.0458, 0.9935);
Vec3d inlet_cut_translation = Vec3d(1.5611, 5.8559, -30.8885) * scaling - inlet_normal * 2.0 * dp_0;
Vec3d inlet_buffer_translation = Vec3d(1.5611, 5.8559, -30.8885) * scaling + inlet_normal * 2.0 * dp_0;

// 使用新的 RotationCalculator
RotationCalculator inlet_rotation_calculator(inlet_normal, x_Axis);
Rotation3d inlet_emitter_rotation(inlet_rotation_calculator.getRotationAngle(), inlet_rotation_calculator.getRotationAxis());
Rotation3d inlet_disposer_rotation(inlet_rotation_calculator.getRotationAngle() + M_PI, inlet_rotation_calculator.getRotationAxis());

// outlet1 参数 outlet1 R=1.9416, (-2.6975, -0.4330, 21.7855), (-0.3160, -0.0009, 0.9488)
Real DW_out_up = 1.9416 * 2 * scaling;
Vec3d outlet_large_half(2.0 * dp_0, 2.4 * scaling, 2.4 * scaling);
Vec3d outlet_large_cut_half(2.0 * dp_0, 5.0 * scaling, 5.0 * scaling);
Vec3d outlet_large_normal(-0.3160, -0.0009, 0.9488);
Vec3d outlet_large_cut_translation = Vec3d(-2.6975, -0.4330, 21.7855) * scaling + outlet_large_normal * 2.0 * dp_0;
Vec3d outlet_large_buffer_translation = Vec3d(-2.6975, -0.4330, 21.7855) * scaling - outlet_large_normal * 2.0 * dp_0;

// 使用新的 RotationCalculator
RotationCalculator outlet_large_rotation_calculator(outlet_large_normal, x_Axis);
Rotation3d outlet_large_disposer_rotation(outlet_large_rotation_calculator.getRotationAngle(), outlet_large_rotation_calculator.getRotationAxis());
Rotation3d outlet_large_emitter_rotation(outlet_large_rotation_calculator.getRotationAngle() + M_PI, outlet_large_rotation_calculator.getRotationAxis());

// outlet2 参数 outlet2 R=1.3261, (9.0220, 0.9750, 18.6389), (-0.0399, 0.0693, 0.9972)
Real DW_out_down = 1.3261 * 2 * scaling;
Vec3d outlet_small_half(2.0 * dp_0, 2.0 * scaling, 2.0 * scaling);
Vec3d outlet_small_cut_half(2.0 * dp_0, 4.0 * scaling, 4.0 * scaling);
Vec3d outlet_small_normal(-0.0399, 0.0693, 0.9972);
Vec3d outlet_small_cut_translation = Vec3d(9.0220, 0.9750, 18.6389) * scaling + outlet_small_normal * 2.0 * dp_0;
Vec3d outlet_small_buffer_translation = Vec3d(9.0220, 0.9750, 18.6389) * scaling - outlet_small_normal * 2.0 * dp_0;

// 使用新的 RotationCalculator
RotationCalculator outlet_small_rotation_calculator(outlet_small_normal, x_Axis);
Rotation3d outlet_small_disposer_rotation(outlet_small_rotation_calculator.getRotationAngle(), outlet_small_rotation_calculator.getRotationAxis());
Rotation3d outlet_small_emitter_rotation(outlet_small_rotation_calculator.getRotationAngle() + M_PI, outlet_small_rotation_calculator.getRotationAxis());
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
    AlignedBox &aligned_box_;

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
// AlignedBoxShape
//----------------------------------------------------------------------
/**
 * @class AlignedBoxShape
 * @brief 用于描述一个 bounding box，其中 upper bound 方向与形状上的平面法线对齐。
 */
class AlignedBoxShape : public TransformShape<GeometricShapeBox>
{
    const int alignment_axis_;

  public:
    /** 直接构造 AlignedBoxShape */
    template <typename... Args>
    explicit AlignedBoxShape(int upper_bound_axis, const Transform &transform, Args &&...args)
        : TransformShape<GeometricShapeBox>(transform, std::forward<Args>(args)...),
          alignment_axis_(upper_bound_axis) {}

    /** 从已有的形状构造 AlignedBoxShape */
    template <typename... Args>
    explicit AlignedBoxShape(int upper_bound_axis, const Shape &shape, Args &&...args)
        : TransformShape<GeometricShapeBox>(
              Transform(Vecd(0.5 * (shape.bounding_box_.second_ + shape.bounding_box_.first_))),
              0.5 * (shape.bounding_box_.second_ - shape.bounding_box_.first_), std::forward<Args>(args)...),
          alignment_axis_(upper_bound_axis) {}

    virtual ~AlignedBoxShape() {}

    Vecd HalfSize() { return halfsize_; }
    int AlignmentAxis() { return alignment_axis_; }

    bool checkInBounds(const Vecd &probe_point)
    {
        Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
        return position_in_frame[alignment_axis_] >= -halfsize_[alignment_axis_] &&
               position_in_frame[alignment_axis_] <= halfsize_[alignment_axis_];
    }

    bool checkUpperBound(const Vecd &probe_point)
    {
        Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
        return position_in_frame[alignment_axis_] > halfsize_[alignment_axis_];
    }

    bool checkLowerBound(const Vecd &probe_point)
    {
        Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
        return position_in_frame[alignment_axis_] < -halfsize_[alignment_axis_];
    }

    bool checkNearUpperBound(const Vecd &probe_point, Real threshold)
    {
        Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
        return ABS(position_in_frame[alignment_axis_] - halfsize_[alignment_axis_]) <= threshold;
    }

    bool checkNearLowerBound(const Vecd &probe_point, Real threshold)
    {
        Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
        return ABS(position_in_frame[alignment_axis_] + halfsize_[alignment_axis_]) <= threshold;
    }

    Vecd getUpperPeriodic(const Vecd &probe_point)
    {
        Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
        Vecd shift = Vecd::Zero();
        shift[alignment_axis_] -= 2.0 * halfsize_[alignment_axis_];
        return transform_.shiftFrameStationToBase(position_in_frame + shift);
    }

    Vecd getLowerPeriodic(const Vecd &probe_point)
    {
        Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
        Vecd shift = Vecd::Zero();
        shift[alignment_axis_] += 2.0 * halfsize_[alignment_axis_];
        return transform_.shiftFrameStationToBase(position_in_frame + shift);
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
        subtract<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_large_emitter_rotation), Vec3d(outlet_large_cut_translation)), outlet_large_cut_half);
        subtract<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_small_emitter_rotation), Vec3d(outlet_small_cut_translation)), outlet_small_cut_half);
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

