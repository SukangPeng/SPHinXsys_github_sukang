/**
 * @file     carotid_transient_state_simulation.h
 * @brief    Header file for transient-state blood flow simulation in the carotid artery
 * @details  Defines necessary classes and parameters for the SPH simulation.
 * @author   Sukang Peng
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
std::string vessel_fluid = "./input/carotid_fluid_geo.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Vec3d translation_water_block(0.0, 0.0, 0.0);                                             /**< Initial translation of the water block (fluid region), Unit: m */
Vec3d translation_wall_boundary(0.0, 0.0, 0.0);                                           /**< Initial translation of the vessel wall (solid boundary), Unit: m */
Real length_scale_vessel = 1.0;                                                           /**< Length scaling factor for the vessel, Dimensionless */
Real length_scale = pow(10, -3);                                                          /**< General length scaling factor, Dimensionless */
Real resolution_ref = 0.3 * length_scale;                                                 /**< Initial reference particle spacing, Unit: m */
Real BW = resolution_ref * 4.0;                                                           /**< Reference size of the emitter (boundary width), Unit: m */
Vec3d domain_lower_bound(-8.0 * length_scale, -6.0 * length_scale, -35.0 * length_scale); /**< Lower boundary of the simulation domain (bounding box), Unit: m */
Vec3d domain_upper_bound(18.0 * length_scale, 15.0 * length_scale, 30.0 * length_scale);  /**< Upper boundary of the simulation domain (bounding box), Unit: m */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);                 /**< Bounding box defining the simulation domain limits, Unit: m */
//----------------------------------------------------------------------
//	NormCalculator
//----------------------------------------------------------------------
// Compute normal using cross product method
Vecd computeNormalCrossProduct(const Vecd &pointA, const Vecd &pointB, const Vecd &pointC)
{
    Vecd v1 = pointB - pointA;
    Vecd v2 = pointC - pointA;
    Vecd normal = v1.cross(v2);
    return normal.normalized(); // Normalize the vector
}

// Compute normal using center point and direction point
Vecd computeNormalDirectional(const Vecd &centerPoint, const Vecd &directionPoint)
{
    Vecd normal = directionPoint - centerPoint;
    return normal.normalized(); // Normalize the vector
}
//----------------------------------------------------------------------
// RotationCalculator
//----------------------------------------------------------------------
struct RotationCalculator
{
    Vec3d rotation_axis;   // Rotation axis
    Real rotation_angle;   // Rotation angle (radians)
    Mat3d rotation_matrix; // Rotation matrix

    // 1. Compute rotation from two points (does not affect inlet/outlet but provides an alternative method)
    RotationCalculator(const Vec3d &pointA, const Vec3d &pointB, const Vec3d &initial_direction)
    {
        Vec3d target_direction = (pointB - pointA).normalized(); // Target direction vector
        initialize(target_direction, initial_direction);
    }

    // 2. Compute rotation from direction vectors (used for inlet/outlet)
    RotationCalculator(const Vec3d &target_direction, const Vec3d &initial_direction)
    {
        initialize(target_direction.normalized(), initial_direction.normalized());
    }

    // Unified initialization logic
    void initialize(const Vec3d &target, const Vec3d &initial)
    {
        // Compute rotation axis
        rotation_axis = initial.cross(target);

        // Compute rotation angle
        rotation_angle = std::acos(initial.dot(target));

        // Handle numerical stability of the rotation axis
        if (rotation_axis.norm() < 1e-6)
        {
            if (initial.dot(target) < 0)
            {
                rotation_axis = Vec3d(1, 0, 0); // Use X-axis as the default rotation axis
                rotation_angle = M_PI;
            }
            else
            {
                rotation_axis = Vec3d(0, 0, 1); // Keep unchanged
                rotation_angle = 0;
            }
        }
        else
        {
            rotation_axis.normalize();
        }

        // Compute rotation matrix
        rotation_matrix = Eigen::AngleAxis<Real>(rotation_angle, rotation_axis).toRotationMatrix();
    }

    // Print the rotation matrix
    void printRotationMatrix() const
    {
        std::cout << "Rotation Matrix:\n"
                  << rotation_matrix << std::endl;
    }

    // Print the rotation axis
    void printRotationAxis() const
    {
        std::cout << "Rotation Axis: (" << rotation_axis.x() << ", " << rotation_axis.y() << ", " << rotation_axis.z() << ")" << std::endl;
    }

    // Print the rotation angle
    void printRotationAngle() const
    {
        std::cout << "Rotation Angle (in radians): " << rotation_angle << std::endl;
    }

    // Get the rotation matrix
    Mat3d getRotationMatrix() const { return rotation_matrix; }

    // Get the rotation axis
    Vec3d getRotationAxis() const { return rotation_axis; }

    // Get the rotation angle
    Real getRotationAngle() const { return rotation_angle; }
};
//----------------------------------------------------------------------
// Buffer Emitter Disposer
//----------------------------------------------------------------------
// Standard direction of the X-axis
Vec3d x_Axis(1, 0, 0);

Real d_inlet = 6.0 * length_scale;
Real A_inlet = 30.7e-6;
Vec3d inlet_half(2.0 * resolution_ref, 3.5 * length_scale, 3.5 * length_scale);
Vec3d inlet_cut_half(2.0 * resolution_ref, 6.0 * length_scale, 6.0 * length_scale);
Vec3d inlet_fix_half(2.0 * resolution_ref, 6.0 * length_scale, 6.0 * length_scale);
Vec3d inlet_normal(0.1034, -0.0458, 0.9935);
Vec3d inlet_normal_flow(0.1034, -0.0458, 0.9935);
Vec3d inlet_cut_translation = Vec3d(1.56, 5.85, -30.9) * length_scale - inlet_normal * 2.0 * resolution_ref;
Vec3d inlet_buffer_translation = Vec3d(1.56, 5.85, -30.9) * length_scale + inlet_normal * 2.0 * resolution_ref;
Vec3d inlet_delete_translation = Vec3d(1.56, 5.85, -30.9) * length_scale - inlet_normal * 5.0 * resolution_ref;
Vec3d inlet_delete_half(5.0 * resolution_ref, 20.0 * length_scale, 20.0 * length_scale);
// RotationInlet
RotationCalculator inlet_rotation_calculator(inlet_normal, x_Axis);
Rotation3d inlet_emitter_rotation(inlet_rotation_calculator.getRotationAngle(), inlet_rotation_calculator.getRotationAxis());
Rotation3d inlet_disposer_rotation(inlet_rotation_calculator.getRotationAngle() + M_PI, inlet_rotation_calculator.getRotationAxis());

Real d_outlet_large = 4.0 * length_scale;
Real A_outlet_large = 15.05e-6;
Vec3d outlet_large_half(2.0 * resolution_ref, 2.4 * length_scale, 2.4 * length_scale);
Vec3d outlet_large_cut_half(2.0 * resolution_ref, 5.0 * length_scale, 5.0 * length_scale);
Vec3d outlet_large_fix_half(2.0 * resolution_ref, 5.0 * length_scale, 5.0 * length_scale);
Vec3d outlet_large_normal(-0.3160, -0.0009, 0.9488);
Vec3d outlet_large_normal_flow(-0.3160, -0.0009, 0.9488);
Vec3d outlet_large_cut_translation = Vec3d(-2.7, -0.4335, 21.78) * length_scale + outlet_large_normal * 2.0 * resolution_ref;
Vec3d outlet_large_buffer_translation = Vec3d(-2.7, -0.4335, 21.78) * length_scale - outlet_large_normal * 2.0 * resolution_ref;
Vec3d outlet_large_delete_translation = Vec3d(-2.7, -0.4335, 21.78) * length_scale + outlet_large_normal * 5.0 * resolution_ref;
Vec3d outlet_large_delete_half(5.0 * resolution_ref, 10.0 * length_scale, 10.0 * length_scale);
// RotationOutletLarge
RotationCalculator outlet_large_rotation_calculator(outlet_large_normal, x_Axis);
Rotation3d outlet_large_disposer_rotation(outlet_large_rotation_calculator.getRotationAngle(), outlet_large_rotation_calculator.getRotationAxis());
Rotation3d outlet_large_emitter_rotation(outlet_large_rotation_calculator.getRotationAngle() + M_PI, outlet_large_rotation_calculator.getRotationAxis());

Real d_outlet_small = 2.5 * length_scale;
Real A_outlet_small = 7.3e-6;
Vec3d outlet_small_half(2.0 * resolution_ref, 2.0 * length_scale, 2.0 * length_scale);
Vec3d outlet_small_cut_half(2.0 * resolution_ref, 4.0 * length_scale, 4.0 * length_scale);
Vec3d outlet_small_fix_half(2.0 * resolution_ref, 4.0 * length_scale, 4.0 * length_scale);
Vec3d outlet_small_normal(-0.0399, 0.0693, 0.9972);
Vec3d outlet_small_normal_flow(-0.0399, 0.0693, 1.1412);
Vec3d outlet_small_cut_translation = Vec3d(9.02, 0.976, 18.64) * length_scale + outlet_small_normal * 2.0 * resolution_ref;
Vec3d outlet_small_buffer_translation = Vec3d(9.02, 0.976, 18.64) * length_scale - outlet_small_normal * 2.0 * resolution_ref;
Vec3d outlet_small_delete_translation = Vec3d(9.02, 0.976, 18.64) * length_scale + outlet_small_normal * 5.0 * resolution_ref;
Vec3d outlet_small_delete_half(5.0 * resolution_ref, 4.5 * length_scale, 4.5 * length_scale);
// RotationOutletSmall
RotationCalculator outlet_small_rotation_calculator(outlet_small_normal, x_Axis);
Rotation3d outlet_small_disposer_rotation(outlet_small_rotation_calculator.getRotationAngle(), outlet_small_rotation_calculator.getRotationAxis());
Rotation3d outlet_small_emitter_rotation(outlet_small_rotation_calculator.getRotationAngle() + M_PI, outlet_small_rotation_calculator.getRotationAxis());
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1060;   /**< Reference density of fluid. */
Real U_f = 0.315;     /**< Characteristic velocity. */
Real U_max = 2 * U_f; /**< Characteristic velocity. */
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f * SMAX(Real(1), d_inlet *d_inlet / (d_outlet_large * d_outlet_large + d_outlet_small * d_outlet_small));
Real mu_f = 0.0035;         /**< Dynamics viscosity. */
Real Outlet_pressure = 0.0; // for comparison with solely velocity inlet bc
// Real Outlet_pressure = 2.666e3;
//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
struct InflowPressure
{
    template <class BoundaryConditionType>
    InflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real curent_time)
    {
        return p;
    }
};

struct OutletInflowPressure
{
    template <class BoundaryConditionType>
    OutletInflowPressure(BoundaryConditionType &boundary_condition) {}

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
// struct InflowVelocity
//{
//    Real u_ref_, t_ref_, interval_;
//    AlignedBox &aligned_box_;
//
//    template <class BoundaryConditionType>
//    InflowVelocity(BoundaryConditionType &boundary_condition)
//        : u_ref_(0.1), t_ref_(0.218), interval_(0.5),
//          aligned_box_(boundary_condition.getAlignedBox()) {}
//
//    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
//    {
//        Vecd target_velocity = velocity;
//        Real run_time = current_time; // 使用传入的 current_time 代替 GlobalStaticVariables::physical_time_
//        int n = static_cast<int>(run_time / interval_);
//        Real t_in_cycle = run_time - n * interval_;
//
//        target_velocity[0] = t_in_cycle < t_ref_ ? 0.5 * sin(4 * Pi * (run_time + 0.0160236)) : u_ref_;
//        return target_velocity;
//    }
//};
struct InflowVelocity
{
    std::vector<std::pair<Real, Real>> time_velocity_data; // 存储 (时间, 速度)
    AlignedBox & aligned_box_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType & boundary_condition)
        : aligned_box_(boundary_condition.getAlignedBox())
    {
        loadVelocityData("./input/V-inlet.csv"); // 读取 CSV 文件
    }

    // 读取 CSV 文件，存储时间和速度
    void loadVelocityData(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: Unable to open velocity file: " << filename << std::endl;
            return;
        }

        std::string line;
        Real time, velocity;
        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            ss >> time;
            if (ss.peek() == ',')
                ss.ignore();
            ss >> velocity;
            time_velocity_data.emplace_back(time, velocity);
        }
        file.close();

        if (time_velocity_data.empty())
        {
            std::cerr << "Error: No velocity data loaded!" << std::endl;
        }
    }

    // 线性插值计算速度
    Real interpolateVelocity(Real current_time)
    {
        if (time_velocity_data.empty())
            return 0.0; // 没有数据时返回 0

        // 如果当前时间小于最小时间，返回第一个速度
        if (current_time <= time_velocity_data.front().first)
            return time_velocity_data.front().second;

        // 如果当前时间超过最大时间，返回最后一个速度
        if (current_time >= time_velocity_data.back().first)
            return time_velocity_data.back().second;

        // 找到 current_time 所在的时间区间
        for (size_t i = 0; i < time_velocity_data.size() - 1; i++)
        {
            if (current_time >= time_velocity_data[i].first && current_time <= time_velocity_data[i + 1].first)
            {
                // 进行线性插值
                Real t1 = time_velocity_data[i].first;
                Real t2 = time_velocity_data[i + 1].first;
                Real v1 = time_velocity_data[i].second;
                Real v2 = time_velocity_data[i + 1].second;
                return v1 + (v2 - v1) * (current_time - t1) / (t2 - t1);
            }
        }

        return time_velocity_data.back().second; // 超过范围返回最后一个速度值
    }

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;
        target_velocity[0] = interpolateVelocity(current_time); // 使用插值速度
        return target_velocity;
    }
};

//----------------------------------------------------------------------
// AlignedBoxShape
//----------------------------------------------------------------------
/**
 * @class AlignedBoxShape
 * @brief Describes a bounding box where the upper bound direction aligns with the normal of the shape's plane.
 */
class AlignedBoxShape : public TransformShape<GeometricShapeBox>
{
    const int alignment_axis_;

  public:
    /** Directly construct an AlignedBoxShape */
    template <typename... Args>
    explicit AlignedBoxShape(int upper_bound_axis, const Transform &transform, Args &&...args)
        : TransformShape<GeometricShapeBox>(transform, std::forward<Args>(args)...),
          alignment_axis_(upper_bound_axis) {}

    /** Construct an AlignedBoxShape from an existing shape */
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
        add<TriangleMeshShapeSTL>(vessel_fluid, translation_water_block, length_scale);
    }
};

class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<ExtrudeShape<TriangleMeshShapeSTL>>(3.0 * resolution_ref, vessel_fluid, translation_wall_boundary, length_scale);

        subtract<TriangleMeshShapeSTL>(vessel_fluid, translation_water_block, length_scale);

        subtract<AlignedBoxShape>(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_cut_translation)), inlet_cut_half);
        subtract<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_large_emitter_rotation), Vec3d(outlet_large_cut_translation)), outlet_large_cut_half);
        subtract<AlignedBoxShape>(xAxis, Transform(Rotation3d(outlet_small_emitter_rotation), Vec3d(outlet_small_cut_translation)), outlet_small_cut_half);
    }
};
//----------------------------------------------------------------------
//	ConstrainGeometryAlignedBox.
//----------------------------------------------------------------------
class ConstrainGeometryAlignedBox : public BodyPartByParticle
{
  public:
    ConstrainGeometryAlignedBox(SPHBody &body, const std::string &body_part_name, const AlignedBox &aligned_box)
        : BodyPartByParticle(body, body_part_name), aligned_box_(aligned_box)
    {
        // Bind the custom tagging function
        TaggingParticleMethod tagging_particle_method = std::bind(&ConstrainGeometryAlignedBox::tagByAlignedBox, this, _1);
        tagParticles(tagging_particle_method);
    };

    virtual ~ConstrainGeometryAlignedBox() {};

  private:
    AlignedBox aligned_box_; // Bound `AlignedBox` used for particle selection

    void tagByAlignedBox(size_t index_i)
    {
        // Get the particle position
        Vecd particle_position = base_particles_.ParticlePositions()[index_i];

        // Use `AlignedBox`'s `checkInBounds()` to determine if the particle is inside the `AlignedBox`
        if (aligned_box_.checkInBounds(particle_position))
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
//----------------------------------------------------------------------
// DeleteParticlesInDetectionBox.
//----------------------------------------------------------------------
/**
 * @class DeleteParticlesInAlignedBoxByCell
 * @brief Deletes (switches to buffer particles) particles located within a specified detection region.
 *
 * This class inherits from BaseLocalDynamics<BodyPartByCell>,
 * uses the detection region (aligned_box) defined by BodyAlignedBoxByCell,
 * and calls particles_->switchToBufferParticle(index) for all particles within the region to handle deletion.
 */
class DeleteParticlesInAlignedBoxByCell : public BaseLocalDynamics<BodyPartByCell>
{
  public:
    explicit DeleteParticlesInAlignedBoxByCell(BodyAlignedBoxByCell &aligned_box_part)
        : BaseLocalDynamics<BodyPartByCell>(aligned_box_part),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          aligned_box_(aligned_box_part.getAlignedBox()) // In the new version, use getAlignedBox() to obtain the AlignedBox object
    {
    }

    virtual ~DeleteParticlesInAlignedBoxByCell() {}

    virtual void update(size_t index_i, Real dt = 0.0)
    {
        // Use lock_guard for simple locking to ensure memory safety in multithreading
        std::lock_guard<std::mutex> lock(mutex_);
        // If the particle is within the detection region and index_i is a real particle, switch it to a buffer particle
        if (aligned_box_.checkInBounds(pos_[index_i]) && index_i < particles_->TotalRealParticles())
        {
            particles_->switchToBufferParticle(index_i);
        }
    }

  protected:
    std::mutex mutex_;
    Vecd *pos_;               ///< Particle position data
    AlignedBox &aligned_box_; ///< AlignedBox object for the detection region
};
//----------------------------------------------------------------------
// FlowRateCalculator
//----------------------------------------------------------------------
/**
 * @class FlowRateCalculator
 * @brief Calculates the volumetric flow rate Q by summing, over all particles
 *        in BodyPartByCell, the product of the particle's velocity projection
 *        (onto a defined normal) and its volume, and then computes the
 *        volume-weighted average velocity v_avg using the given cross-sectional area.
 *        The results are written to output files.
 */
class FlowRateCalculator
    : public BaseLocalDynamicsReduce<ReduceSum<Real>, BodyPartByCell>
{
  public:
    using ReturnType = Real; // Reduction return type

    /**
     * Constructor.
     * @param aligned_box_part   The BodyAlignedBoxByCell object corresponding to the inlet/outlet region.
     * @param normal             The cross-sectional normal vector (which will be normalized).
     * @param cross_section_area The cross-sectional area.
     * @param prefix             The prefix used for output file names (e.g., "inlet", "outlet_large", "outlet_small").
     */
    explicit FlowRateCalculator(BodyAlignedBoxByCell &aligned_box_part,
                                const Vecd &normal,
                                Real cross_section_area,
                                const std::string &prefix)
        : BaseLocalDynamicsReduce<ReduceSum<Real>, BodyPartByCell>(aligned_box_part),
          normal_(normal),
          outlet_area_(cross_section_area),
          transient_flow_rate_(0.0),
          average_velocity_(0.0),
          name_prefix_(prefix)
    {
        // Retrieve velocity and volumetric measure data from BaseParticles.
        vel_ = this->particles_->getVariableDataByName<Vecd>("Velocity");
        Vol_ = this->particles_->getVariableDataByName<Real>("VolumetricMeasure");

        // Get the system physical time pointer.
        physical_time_ = this->sph_system_.getSystemVariableDataByName<Real>("PhysicalTime");
    }

    virtual ~FlowRateCalculator() {}

    /**
     * For each particle, return its contribution to the flow rate:
     * contribution = (velocity · normal_) * Vol.
     * The index here is determined by the reordered list of particles within the region.
     */
    virtual Real reduce(size_t index, Real dt = 0.0)
    {
        return vel_[index].dot(normal_) * Vol_[index];
    }

    /**
     * After the reduction, compute the volume-weighted average velocity:
     * v_avg = ∑(v · Vol) / ∑Vol,
     * then compute the flow rate Q = v_avg * outlet_area_, and write the results to files.
     *
     * Note: Instead of a simple iteration from 0 to SizeOfLoopRange(),
     *       we iterate over the list of cells in the region to ensure that
     *       only the volumes corresponding to particles within the region are included,
     *       avoiding counting particles outside the region.
     */
    virtual Real outputResult(Real reduced_value) override
    {
        Real total_volume = 0.0;
        // Accumulate the volume from the cell list corresponding to the region.
        auto &cell_lists = this->getDynamicsIdentifier().LoopRange();
        for (size_t c = 0; c < cell_lists.size(); ++c)
        {
            // Assume each cell_lists[c] is a container (e.g., vector<size_t>) storing global indices of particles in the region.
            for (size_t i = 0; i < cell_lists[c]->size(); ++i)
            {
                size_t global_index = (*cell_lists[c])[i];
                total_volume += Vol_[global_index];
            }
        }
        if (total_volume > 0.0)
        {
            average_velocity_ = reduced_value / total_volume;
        }
        else
        {
            average_velocity_ = 0.0;
        }
        // Compute the flow rate.
        transient_flow_rate_ = average_velocity_ * outlet_area_;

        // Output results to files.
        std::string output_folder = "./output";
        std::filesystem::create_directory(output_folder);

        // Flow rate file.
        std::string flow_rate_file = output_folder + "/flow_rate_" + name_prefix_ + ".dat";
        {
            std::ofstream flow_out(flow_rate_file, std::ios::app);
            if (flow_out.is_open())
            {
                flow_out << *physical_time_ << "\t" << transient_flow_rate_ << "\n";
            }
        }

        // Average velocity file.
        std::string velocity_file = output_folder + "/average_velocity_" + name_prefix_ + ".dat";
        {
            std::ofstream velocity_out(velocity_file, std::ios::app);
            if (velocity_out.is_open())
            {
                velocity_out << *physical_time_ << "\t" << average_velocity_ << "\n";
            }
        }

        return transient_flow_rate_;
    }

  private:
    // Normal vector of the cross-section (assumed normalized).
    Vecd normal_;
    // Cross-sectional area.
    Real outlet_area_;
    // Prefix for naming output files.
    std::string name_prefix_;

    // Intermediate variables during calculation.
    Real transient_flow_rate_;
    Real average_velocity_;

    // Pointers to particle data.
    Vecd *vel_;
    Real *Vol_;

    // Pointer to the system physical time.
    Real *physical_time_;
};
