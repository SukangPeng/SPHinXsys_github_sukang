﻿/**
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
std::string stent_path = "./input/stent_test.stl";
//std::string stent_path = "./input/Stent.stl";
//std::string vessel_wall_path = "./input/vessel_wall.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Vec3d translation_vessel_wall(0.0, 0.0, 0.0); /**< 血管壁的初始平移，单位：m (米) */
Vec3d translation_stent(0.0, 0.0, 0.0);      /**< 支架的初始平移，单位：m (米) */
Real length_scale_stent = 1e-3;                 /**< 长度比例因子，无量纲 (无单位) */
Real length_scale = 1.0;                        /**< 长度比例因子，无量纲 (无单位) */
Real resolution_ref = 0.3 * length_scale_stent; /**< 初始参考粒子间距，单位：m (米) */
// resolution_ref = 0.0003;                                            /**< 初始参考粒子间距，单位：m (米) */
Real diameter = 0.004;                                                    /**< 血管外径，单位：m (米) */
Vec3d domain_lower_bound(-0.006, -0.002, -0.002);                         /**< 系统域的下边界，单位：m (米) */
Vec3d domain_upper_bound(0.006, 0.002, 0.002);                            /**< 系统域的上边界，单位：m (米) */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound); /**< 系统域的边界框，单位：m (米) */
Real full_length = 0.025;
//----------------------------------------------------------------------
//	Global parameters for vessel (血管壁参数)
//----------------------------------------------------------------------
/*Vessel Material*/
Real rho0_s_vessel = 1265;              /**< 血管壁的密度，单位：kg/m³ (千克每立方米) */
Real poisson_vessel = 0.45;             /**< 血管壁的泊松比，无量纲 (无单位) */
Real Youngs_modulus_vessel = 1e5;     /**< 血管壁的杨氏模量，单位：Pa (帕) */
Real physical_viscosity_vessel = 200.0; /**< 血管壁的物理粘度，单位：Pa·s (帕·秒) */
//----------------------------------------------------------------------
//	Global parameters for stent (血管壁参数)
//----------------------------------------------------------------------
/* Stent Material */
Real rho0_s_stent = 6450.0;            /**< 支架的密度，单位：kg/m³ (千克每立方米) */
Real poisson_stent = 0.33;             /**< 支架的泊松比，无量纲 (无单位) */
Real youngs_modulus_stent = 1e6;       /**< 支架的杨氏模量，单位：Pa (帕) */
Real physical_viscosity_stent = 100.0; /**< 支架的物理粘度，单位：Pa·s (帕·秒) */
//----------------------------------------------------------------------
//	Define SPH bodies.
//----------------------------------------------------------------------
class Stent : public ComplexShape
{
  public:
    explicit Stent(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(stent_path, translation_stent, length_scale_stent);
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
//	RadialForce.
//----------------------------------------------------------------------
class RadialForce
{
  protected:
    Vecd reference_acceleration_;
    int axis_; // xAxis, yAxis, or zAxis

  public:
    RadialForce(Real magnitude, int axis)
        : reference_acceleration_(Vecd::Zero()), axis_(axis)
    {
        // Set the magnitude along the specified axis.
        reference_acceleration_[axis_] = magnitude;
    }
    ~RadialForce() {}

    Vecd InducedAcceleration(const Vecd &position = Vecd::Zero(), Real physical_time = 0.0) const
    {
        // Calculate radial direction based on the axis.
        Vecd radial_direction = Vecd::Zero();
        if (axis_ == xAxis)
        {
            radial_direction[1] = position[1];
            radial_direction[2] = position[2];
        }
        else if (axis_ == yAxis)
        {
            radial_direction[0] = position[0];
            radial_direction[2] = position[2];
        }
        else if (axis_ == zAxis)
        {
            radial_direction[0] = position[0];
            radial_direction[1] = position[1];
        }

        // Normalize to get unit vector in radial direction.
        radial_direction /= (radial_direction.norm() + 1e-6); // Adding small value to avoid division by zero.

        // Return the radial acceleration based on the magnitude defined.
        return reference_acceleration_[axis_] * radial_direction;
    }
};

class StartupRadialForce : public RadialForce
{
    Real target_time_;

  public:
    StartupRadialForce(Real magnitude, int axis, Real target_time)
        : RadialForce(magnitude, axis), target_time_(target_time) {}
    ~StartupRadialForce() {}

    Vecd InducedAcceleration(const Vecd &position, Real physical_time) const
    {
        Real time_factor = physical_time / target_time_;
        Vecd acceleration = 0.5 * Pi * sin(Pi * time_factor) * RadialForce::InducedAcceleration(position);
        return time_factor < 1.0 ? acceleration : Vecd::Zero();
    }
};

class IncreaseToFullRadialForce : public RadialForce
{
    Real time_to_full_force_;

  public:
    IncreaseToFullRadialForce(Real magnitude, int axis, Real time_to_full_force)
        : RadialForce(magnitude, axis), time_to_full_force_(time_to_full_force) {}
    ~IncreaseToFullRadialForce() {}

    Vecd InducedAcceleration(const Vecd &position, Real physical_time) const
    {
        Real time_factor = physical_time / time_to_full_force_;
        Vecd full_acceleration = RadialForce::InducedAcceleration(position);
        return time_factor < 1.0 ? time_factor * full_acceleration : full_acceleration;
    }
};


//----------------------------------------------------------------------
//  RadialForceApplication.
//----------------------------------------------------------------------
template <class RadialForceType>
class RadialForceApplication : public ForcePrior
{
  protected:
    const RadialForceType radial_force_;
    Vecd *pos_;
    Real *mass_;
    Real *physical_time_;

  public:
    RadialForceApplication(SPHBody &sph_body, const RadialForceType &radial_force)
        : ForcePrior(sph_body, "RadialForce"), radial_force_(radial_force),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          mass_(particles_->registerStateVariable<Real>("Mass")),
          physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime")) {}
    virtual ~RadialForceApplication() {}

    void update(size_t index_i, Real dt = 0.0)
    {
        current_force_[index_i] =
            mass_[index_i] * radial_force_.InducedAcceleration(pos_[index_i], *physical_time_);
        ForcePrior::update(index_i, dt);
    }
};

//----------------------------------------------------------------------
// BoundingBox Calculation Functions
//----------------------------------------------------------------------
BoundingBox getRealTimeBoundingBox(BaseParticles &particles)
{
    // Initialize bounding box with extreme values
    Vecd lower_bound = Vecd::Constant(std::numeric_limits<Real>::max());
    Vecd upper_bound = Vecd::Constant(-std::numeric_limits<Real>::max());

    // Access particle position data directly
    const Vecd *positions = particles.ParticlePositions();

    // Iterate over all particle positions to update bounding box
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
    {
        const Vecd &particle_pos = positions[i];
        lower_bound = lower_bound.cwiseMin(particle_pos);
        upper_bound = upper_bound.cwiseMax(particle_pos);
    }

    return BoundingBox(lower_bound, upper_bound);
}

BoundingBox getAlignedBoundingBox(BaseParticles &particles, const Mat3d &rotation, const Vecd &translation)
{
    // Initialize bounding box in local coordinate system
    Vecd local_lower_bound = Vecd::Constant(std::numeric_limits<Real>::max());
    Vecd local_upper_bound = Vecd::Constant(-std::numeric_limits<Real>::max());

    // Access particle position data directly
    const Vecd *positions = particles.ParticlePositions();

    // Calculate inverse rotation matrix and translation vector
    Mat3d inv_rotation = rotation.transpose(); // Inverse rotation
    Vecd inv_translation = -translation;       // Inverse translation

    // Iterate over all particle positions
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
    {
        // Get global coordinate of the particle
        const Vecd &particle_pos = positions[i];

        // Apply inverse translation and inverse rotation to get local coordinates
        Vecd local_pos = inv_rotation * (particle_pos + inv_translation);

        // Update bounding box in the local coordinate system
        local_lower_bound = local_lower_bound.cwiseMin(local_pos);
        local_upper_bound = local_upper_bound.cwiseMax(local_pos);
    }

    // Return bounding box aligned with the stent
    return BoundingBox(local_lower_bound, local_upper_bound);
}

void printBoundingBoxAndDelta(const BoundingBox &bbox)
{
    Vecd lower_bound = bbox.first_;
    Vecd upper_bound = bbox.second_;
    Vecd delta = upper_bound - lower_bound;

    std::cout << "Current Bounding Box Lower Bound: ("
              << lower_bound[0] << ", " << lower_bound[1] << ", " << lower_bound[2] << ")\n";
    std::cout << "Current Bounding Box Upper Bound: ("
              << upper_bound[0] << ", " << upper_bound[1] << ", " << upper_bound[2] << ")\n";
    std::cout << lower_bound[0] << " to " << upper_bound[0] << " (delta: " << delta[0] << ")\n";
    std::cout << lower_bound[1] << " to " << upper_bound[1] << " (delta: " << delta[1] << ")\n";
    std::cout << lower_bound[2] << " to " << upper_bound[2] << " (delta: " << delta[2] << ")\n";
}


//----------------------------------------------------------------------
//	ReloadParticleRecordingToXml
//----------------------------------------------------------------------
/**
 * @class ReloadParticleRecordingToXml
 * @brief This class records the latest particle state in XML format, inheriting directly from BaseIO.
 * It writes the particle state to an XML file which can be used for reloading the particle state in future simulations.
 */
class ReloadParticleRecordingToXml : public BaseIO
{
  public:
    // 构造函数，传入需要记录状态的 SPHBody
    ReloadParticleRecordingToXml(SPHBody &sph_body)
        : BaseIO(sph_body.getSPHSystem()), sph_body_(sph_body), base_particles_(sph_body.getBaseParticles())
    {
        // 确保输出文件夹正确设置
        output_folder_ = io_environment_.output_folder_ + "/particle-reload";
        if (!fs::exists(output_folder_))
        {
            fs::create_directories(output_folder_);
        }
    }

    // 公共函数，将当前粒子状态写入 XML 文件
    void writeToFile(size_t iteration_step)
    {
        // 将迭代步骤转换为字符串，用于文件名
        std::string sequence = std::to_string(iteration_step);
        // 构造保存粒子重载数据的文件路径（XML格式）
        std::string filefullpath = output_folder_ + "/particle_reload_" + sph_body_.getName() + "_" + sequence + ".xml";

        // 如果文件已存在，删除旧文件
        if (fs::exists(filefullpath))
        {
            fs::remove(filefullpath);
        }

        // 打开输出文件流
        std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);

        // 开始写入 XML 结构
        out_file << "<?xml version=\"1.0\"?>\n";
        out_file << "<particles>\n";

        // 遍历所有真实粒子并写入数据
        size_t total_real_particles = base_particles_.TotalRealParticles();
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            Vecd position = base_particles_.ParticlePositions()[i];
            Real volume = base_particles_.VolumetricMeasures()[i];

            // 以单行格式写入每个粒子的属性
            out_file << "  <particle VolumetricMeasure=\"" << volume << "\" Position=\""
                     << position[0] << ", " << position[1];

            // 如果是三维项目，则添加第三个坐标
            if (position.size() == 3)
            {
                out_file << ", " << position[2];
            }

            out_file << "\"/>\n";
        }

        // 关闭 XML 结构
        out_file << "</particles>\n";
        out_file.close();

        // 调试：确认文件已写入
        std::cout << "Particle state for " << sph_body_.getName() << " written to " << filefullpath << std::endl;
    }

  private:
    SPHBody &sph_body_;             // 我们正在写入粒子状态的 SPHBody
    BaseParticles &base_particles_; // 该物体的粒子
    std::string output_folder_;     // 输出文件将被写入的文件夹
};