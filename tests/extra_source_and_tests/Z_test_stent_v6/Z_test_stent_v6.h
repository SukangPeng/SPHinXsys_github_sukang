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
std::string stent_path = "./input/stent_cylinder.stl";
std::string vessel_wall_path = "./input/vessel_wall_narrow.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Vec3d translation_vessel_wall(0.0, 0.0, 0.0);                                                               /**< 血管壁的初始平移，单位：m (米) */
Vec3d translation_stent(0.0, 0.0, 0.0);                                                                     /**< 支架的初始平移，单位：m (米) */
Real length_scale_stent = 1e-3;                                                                             /**< 长度比例因子，无量纲 (无单位) */
Real length_scale = 1.0;                                                                                    /**< 长度比例因子，无量纲 (无单位) */
Real resolution_ref = 0.3 * length_scale_stent;                                                             /**< 初始参考粒子间距，单位：m (米) */
Real diameter = 0.004;                                                                                      /**< 血管外径，单位：m (米) */
Vec3d domain_lower_bound(-16.0 * length_scale_stent, -2.5 * length_scale_stent, -2.5 * length_scale_stent); /**< 系统域的下边界，单位：m (米) */
Vec3d domain_upper_bound(16.0 * length_scale_stent, 2.5 * length_scale_stent, 2.5 * length_scale_stent);    /**< 系统域的上边界，单位：m (米) */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);                                   /**< 系统域的边界框，单位：m (米) */
Real full_length = 0.025;
//----------------------------------------------------------------------
//	Global parameters for vessel (血管壁参数)
//----------------------------------------------------------------------
/*Vessel Material*/
Real rho0_s_vessel = 1265;             /**< 血管壁的密度，单位：kg/m³ (千克每立方米) */
Real poisson_vessel = 0.45;            /**< 血管壁的泊松比，无量纲 (无单位) */
Real Youngs_modulus_vessel = 1e6;      /**< 血管壁的杨氏模量，单位：Pa (帕) */
Real physical_viscosity_vessel = 20.0; /**< 血管壁的物理粘度，单位：Pa·s (帕·秒) */
//----------------------------------------------------------------------
//	Global parameters for stent (血管壁参数)
//----------------------------------------------------------------------
/* Stent Material */
Real rho0_s_stent = 6450.0;           /**< 支架的密度，单位：kg/m³ (千克每立方米) */
Real poisson_stent = 0.33;            /**< 支架的泊松比，无量纲 (无单位) */
Real youngs_modulus_stent = 1e9;      /**< 支架的杨氏模量，单位：Pa (帕) */
Real physical_viscosity_stent = 40.0; /**< 支架的物理粘度，单位：Pa·s (帕·秒) */
//----------------------------------------------------------------------
//	Define SPH bodies.
//----------------------------------------------------------------------
class Stent : public ComplexShape
{
  public:
    explicit Stent(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(stent_path, translation_stent, length_scale);
    }
};

class VesselWall : public ComplexShape
{
  public:
    explicit VesselWall(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(vessel_wall_path, translation_vessel_wall, length_scale);
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

    Vecd GetCurrentForce() const
    {
        // Return the current force vector along the specified axis.
        Vecd force_vector = reference_acceleration_;
        return force_vector;
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
        Real sine_value = (time_factor < 1.0) ? sin(Pi / 2 * time_factor) : 1.0; // 只取上升部分，达到最大值后保持为 1.0

        // 使用相同的逻辑计算加速度
        Vecd acceleration = sine_value * RadialForce::InducedAcceleration(position);
        return acceleration;
    }

    Vecd GetCurrentForce(Real physical_time) const
    {
        Real time_factor = physical_time / target_time_;
        Real sine_value = (time_factor < 1.0) ? sin(Pi / 2 * time_factor) : 1.0; // 只取上升部分，达到最大值后保持为 1.0

        Real current_magnitude = sine_value * reference_acceleration_[axis_];

        Vecd force_vector = Vecd::Zero();
        force_vector[axis_] = current_magnitude;
        return force_vector;
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
//	RotationCalculator
//----------------------------------------------------------------------
struct RotationCalculator
{
    Vec3d rotation_axis;   // 旋转轴
    Real rotation_angle;   // 旋转角度（弧度）
    Mat3d rotation_matrix; // 旋转矩阵

    // 构造函数：根据两点和初始方向向量计算旋转矩阵、轴和角度
    RotationCalculator(const Vec3d &pointA, const Vec3d &pointB, const Vec3d &initial_direction)
    {
        Vec3d target_direction = (pointB - pointA).normalized();             // 目标方向向量
        Vec3d normalized_initial_direction = initial_direction.normalized(); // 规范化初始方向向量

        // 计算旋转轴
        rotation_axis = normalized_initial_direction.cross(target_direction).normalized();

        // 计算旋转角度
        rotation_angle = std::acos(normalized_initial_direction.dot(target_direction));

        // 构造旋转矩阵
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

    // 打印旋转角度
    void printRotationAngle() const
    {
        std::cout << "Rotation Angle (in radians): " << rotation_angle << std::endl;
    }

    // 获取旋转矩阵
    Mat3d getRotationMatrix() const
    {
        return rotation_matrix;
    }

    // 获取旋转轴
    Vec3d getRotationAxis() const
    {
        return rotation_axis;
    }

    // 获取旋转角度
    Real getRotationAngle() const
    {
        return rotation_angle;
    }
};

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

//----------------------------------------------------------------------
//	Define constrain class for stent translation and rotation.
//----------------------------------------------------------------------
/**
 * @class QuantityMassPosition
 * @brief Compute the mass-weighted position of a body
 */
template <typename DynamicsIdentifier = SPHBody>
class QuantityMassPosition : public QuantitySummation<Vecd, DynamicsIdentifier>
{
  protected:
    Real *mass_; // 指向粒子的质量数据

  public:
    explicit QuantityMassPosition(DynamicsIdentifier &identifier)
        : QuantitySummation<Vecd, DynamicsIdentifier>(identifier, "Position"),
          mass_(this->particles_->template getVariableDataByName<Real>("Mass"))
    {
        this->quantity_name_ = "MassWeightedPosition";
    };
    virtual ~QuantityMassPosition() {}

    // 重写 reduce 方法，用于计算质量加权的位置
    Vecd reduce(size_t index_i, Real dt = 0.0)
    {
        // 获取粒子的位置和质量
        Vecd position = this->variable_[index_i];
        Real mass = mass_[index_i];

        // 计算质量加权的位置
        return position * mass;
    }
};

Vecd computeAveragePosition(SolidBody &stent_body)
{
    // 获取 stent_body 的 BaseParticles 实例
    BaseParticles &particles = stent_body.getBaseParticles();

    // 获取粒子总数
    size_t total_particles = particles.TotalRealParticles();

    // 如果没有粒子，返回零向量
    if (total_particles == 0)
    {
        return Vecd::Zero();
    }

    // 获取粒子的位置向量数组
    Vecd *positions = particles.ParticlePositions();

    // 初始化总位置
    Vecd total_position = Vecd::Zero();

    // 遍历所有粒子并累加位置
    for (size_t i = 0; i < total_particles; ++i)
    {
        total_position += positions[i];
    }

    // 计算平均位置
    return total_position / static_cast<Real>(total_particles);
}

/**
 * @class QuantityMomentOfInertia
 * @brief Compute the moment of inertia of a body
 */
template <typename DynamicsIdentifier = SPHBody>
class QuantityMomentOfInertia : public QuantitySummation<Real, DynamicsIdentifier>
{
  protected:
    Vecd mass_center_; // 质心位置
    size_t p_1_, p_2_; // 惯性矩分量的索引
    Vecd *positions_;  // 粒子的位置数据

  public:
    explicit QuantityMomentOfInertia(DynamicsIdentifier &identifier, Vecd mass_center, size_t p_1, size_t p_2)
        : QuantitySummation<Real, DynamicsIdentifier>(identifier, "Mass"),
          mass_center_(mass_center), p_1_(p_1), p_2_(p_2),
          positions_(this->particles_->template getVariableDataByName<Vecd>("Position"))
    {
        this->quantity_name_ = "MomentOfInertia";

        // 检查 positions_ 是否初始化成功
        if (!positions_)
        {
            throw std::runtime_error("Error: Unable to initialize positions_. Check if 'Position' variable exists.");
        }
    }

    virtual ~QuantityMomentOfInertia() {}

    // 重写 reduce 方法，用于计算惯性矩分量
    Real reduce(size_t index_i, Real dt = 0.0)
    {
        Vecd relative_position = positions_[index_i] - mass_center_; // 访问粒子位置
        Real mass = this->variable_[index_i];                        // 获取粒子质量

        if (p_1_ == p_2_)
        {
            // 对角项惯性矩分量
            return mass * (relative_position.squaredNorm() - relative_position[p_1_] * relative_position[p_1_]);
        }
        else
        {
            // 非对角项惯性矩分量
            return -mass * relative_position[p_1_] * relative_position[p_2_];
        }
    }
};

/**
 * @class QuantityMomentOfMomentum
 * @brief Computes the moment of momentum (angular momentum) for a given SPHBody.
 */
class QuantityMomentOfMomentum : public QuantitySummation<Vecd, SPHBody>
{
  protected:
    Vecd *positions_;  // Pointer to particle positions
    Vecd *velocities_; // Pointer to particle velocities
    Real *masses_;     // Pointer to particle masses
    Vecd mass_center_; // Center of mass of the body

  public:
    /**
     * @brief Constructor for QuantityMomentOfMomentum
     * @param sph_body The SPHBody to compute the moment of momentum for.
     * @param mass_center The center of mass of the body.
     */
    explicit QuantityMomentOfMomentum(SPHBody &sph_body, Vecd mass_center)
        : QuantitySummation<Vecd, SPHBody>(sph_body, "Velocity"),
          positions_(sph_body.getBaseParticles().getVariableDataByName<Vecd>("Position")),
          velocities_(sph_body.getBaseParticles().getVariableDataByName<Vecd>("Velocity")),
          masses_(sph_body.getBaseParticles().getVariableDataByName<Real>("Mass")),
          mass_center_(mass_center)
    {
        this->quantity_name_ = "MomentOfMomentum";
    }

    /**
     * @brief Compute the contribution to the moment of momentum for a single particle.
     * @param index_i Index of the particle.
     * @param dt Time step size (not used in this calculation).
     * @return The moment of momentum contribution of the particle.
     */
    Vecd reduce(size_t index_i, Real dt = 0.0)
    {
        Vecd relative_position = positions_[index_i] - mass_center_;
        return masses_[index_i] * relative_position.cross(velocities_[index_i]);
    }
};

/**
 * @class Constrain3DSolidBodyTranslation
 * @brief Constrain the translation of a 3D solid body.
 */
class Constrain3DSolidBodyRotation : public LocalDynamics
{
  private:
    Vecd mass_center_;                                                          // 物体的质心
    Matd moment_of_inertia_;                                                    // 惯性张量
    Vecd angular_velocity_;                                                     // 当前角速度
    ReduceDynamics<QuantityMomentOfMomentum> compute_total_moment_of_momentum_; // 计算总角动量
    Vecd *positions_;                                                           // 粒子的位置数组
    Vecd *velocities_;                                                          // 粒子速度数组

  protected:
    // 在每个时间步开始时，计算角速度
    virtual void setupDynamics(Real dt = 0.0) override
    {
        angular_velocity_ = moment_of_inertia_.inverse() * compute_total_moment_of_momentum_.exec(dt);
    }

  public:
    // 构造函数，初始化变量
    explicit Constrain3DSolidBodyRotation(SPHBody &sph_body, const Vecd &mass_center, const Matd &inertia_tensor)
        : LocalDynamics(sph_body),
          mass_center_(mass_center),
          moment_of_inertia_(inertia_tensor),
          compute_total_moment_of_momentum_(sph_body, mass_center)
    {
        // 获取粒子数据
        BaseParticles &particles = sph_body.getBaseParticles();
        positions_ = particles.getVariableDataByName<Vecd>("Position");
        velocities_ = particles.getVariableDataByName<Vecd>("Velocity");

        // 检查是否正确初始化变量
        if (!positions_ || !velocities_)
        {
            throw std::runtime_error("Error: Unable to initialize 'Position' or 'Velocity' variable. Check if they are registered.");
        }
    }

    virtual ~Constrain3DSolidBodyRotation() {}

    // 对每个粒子，更新速度以阻止旋转
    void update(size_t index_i, Real dt = 0.0)
    {
        Vecd relative_position = positions_[index_i] - mass_center_;                       // 计算相对质心的位置
        Vecd linear_velocity_due_to_rotation = angular_velocity_.cross(relative_position); // 计算线速度
        velocities_[index_i] -= linear_velocity_due_to_rotation;                           // 调整速度以去除旋转分量
    }
};
