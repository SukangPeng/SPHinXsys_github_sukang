/**
 * @file     coronary_artery_simulation.h
 * @brief    Header file for steady-state blood flow simulation in the coronary artery
 * @details  Defines necessary classes, parameters, and boundary conditions for the SPH simulation.
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
//std::string vessel_wall_path = "./input/vessel_wall_stenosis.stl";
std::string vessel_wall_path = "./input/vessel_wall.stl";
std::string vessel_fluid_path = "./input/vessel_fluid.stl";
//std::string vessel_fluid_path = "./input/vessel_fluid.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Vec3d translation_water_block(0.0, 0.0, 0.0);   /**< Initial translation of the water block, unit: meters (m) */
Vec3d translation_wall_boundary(0.0, 0.0, 0.0); /**< Initial translation of the vessel wall, unit: meters (m) */
// Real length_scale_vessel = pow(10, -3);                                   /**< Length scale factor, dimensionless (unitless) */
Real length_scale_vessel = 1.0;                                           /**< Length scale factor, dimensionless (unitless) */
Real length_scale = pow(10, -3);                                          /**< Length scale factor, dimensionless (unitless) */
Real resolution_ref = 0.3 * length_scale;                                 /**< Initial reference particle spacing, unit: meters (m) */
Vec3d domain_lower_bound = Vec3d(-1.0, -15.0, -35.0) * length_scale;      /**< Lower boundary of the system domain, unit: meters (m) */
Vec3d domain_upper_bound = Vec3d(45.0, 15.0, 35.0) * length_scale;        /**< Upper boundary of the system domain, unit: meters (m) */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound); /**< Defines the bounding box of the system domain */
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

// inlet Parameter: inlet d = 4mm  (0,0,0), (1,0,0)
Real d_inlet = 4.0 * length_scale;
Real A_inlet = 12.4e-6;
Vec3d inlet_half(2.0 * resolution_ref, 3.0 * length_scale, 3.0 * length_scale);
Vec3d inlet_fix_half(2.0 * resolution_ref, 7.0 * length_scale, 7.0 * length_scale);
Vec3d inlet_normal(1.0, 0, 0);
Vec3d inlet_buffer_translation = Vec3d(0, 0, 0) * length_scale + inlet_normal * 2.0 * resolution_ref;
RotationCalculator inlet_rotation_calculator(inlet_normal, x_Axis);
Rotation3d inlet_emitter_rotation(inlet_rotation_calculator.getRotationAngle(), inlet_rotation_calculator.getRotationAxis());
Rotation3d inlet_disposer_rotation(inlet_rotation_calculator.getRotationAngle() + M_PI, inlet_rotation_calculator.getRotationAxis());

// outlet_large Parameter: d=3.5, (24.26695907, -9.65881011, -3.31277237), (22.99372438, -8.15633227, -2.96429891)
Real d_outlet_large = 3.5 * length_scale;
Real A_outlet_large = 9.35e-6;
Vec3d outlet_large_half(2.0 * resolution_ref, 2.5 * length_scale, 2.5 * length_scale);
Vec3d outlet_large_fix_half(2.0 * resolution_ref, 5.5 * length_scale, 5.5 * length_scale);
Vecd outlet_large_centerPoint = Vecd(24.26695907, -9.65881011, -3.31277237) * length_scale;
Vecd outlet_large_directionPoint = Vecd(22.99372438, -8.15633227, -2.96429891) * length_scale;
Vecd outlet_large_normal = computeNormalDirectional(outlet_large_centerPoint, outlet_large_directionPoint);
Vec3d outlet_large_buffer_translation = Vec3d(24.26695907, -9.65881011, -3.31277237) * length_scale + outlet_large_normal * 2.0 * resolution_ref;
RotationCalculator outlet_large_rotation_calculator(outlet_large_normal, x_Axis);
Rotation3d outlet_large_emitter_rotation(outlet_large_rotation_calculator.getRotationAngle(), outlet_large_rotation_calculator.getRotationAxis());
Rotation3d outlet_large_disposer_rotation(outlet_large_rotation_calculator.getRotationAngle() + M_PI, outlet_large_rotation_calculator.getRotationAxis());

// outlet_middle Parameter: d = 3 , (38.75557092, 10.65711145, -25.28028013), (37.47065858, 9.48707165, -24.29032704) A=7.033589
Real d_outlet_middle = 3.0 * length_scale;
Real A_outlet_middle = 7.1e-6;
Vec3d outlet_middle_half(2.0 * resolution_ref, 2.5 * length_scale, 2.5 * length_scale);
Vec3d outlet_middle_fix_half(2.0 * resolution_ref, 4.0 * length_scale, 4.0 * length_scale);
Vecd outlet_middle_centerPoint = Vecd(38.75557092, 10.65711145, -25.28028013) * length_scale;
Vecd outlet_middle_directionPoint = Vecd(37.47065858, 9.48707165, -24.29032704) * length_scale;
Vecd outlet_middle_normal = computeNormalDirectional(outlet_middle_centerPoint, outlet_middle_directionPoint);
Vec3d outlet_middle_buffer_translation = Vec3d(38.75557092, 10.65711145, -25.28028013) * length_scale + outlet_large_normal * 2.0 * resolution_ref;
RotationCalculator outlet_middle_rotation_calculator(outlet_middle_normal, x_Axis);
Rotation3d outlet_middle_emitter_rotation(outlet_middle_rotation_calculator.getRotationAngle(), outlet_middle_rotation_calculator.getRotationAxis());
Rotation3d outlet_middle_disposer_rotation(outlet_middle_rotation_calculator.getRotationAngle() + M_PI, outlet_middle_rotation_calculator.getRotationAxis());

// outlet_small Parameter: d = 2.5  , (27.00875603, 12.17813813, -28.45272815), (26.23558344, 10.92904811, -27.09553315) A=5.021522
Real d_outlet_small = 2.5 * length_scale;
Real A_outlet_small = 5.2e-6;
Vec3d outlet_small_half(2.0 * resolution_ref, 2.0 * length_scale, 2.0 * length_scale);
Vec3d outlet_small_fix_half(2.0 * resolution_ref, 4.0 * length_scale, 4.0 * length_scale);
Vecd outlet_small_centerPoint = Vecd(27.00875603, 12.17813813, -28.45272815) * length_scale;
Vecd outlet_small_directionPoint = Vecd(26.23558344, 10.92904811, -27.09553315) * length_scale;
Vecd outlet_small_normal = computeNormalDirectional(outlet_small_centerPoint, outlet_small_directionPoint);
Vec3d outlet_small_buffer_translation = Vec3d(27.00875603, 12.17813813, -28.45272815) * length_scale + outlet_small_normal * 2.0 * resolution_ref;
RotationCalculator outlet_small_rotation_calculator(outlet_small_normal, x_Axis);
Rotation3d outlet_small_emitter_rotation(outlet_small_rotation_calculator.getRotationAngle(), outlet_small_rotation_calculator.getRotationAxis());
Rotation3d outlet_small_disposer_rotation(outlet_small_rotation_calculator.getRotationAngle() + M_PI, outlet_small_rotation_calculator.getRotationAxis());
//----------------------------------------------------------------------
// Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1060;   /**< Reference density of fluid, unit: kg/m³ (kilograms per cubic meter) */
Real U_f = 0.315;     /**< Characteristic velocity, unit: m/s (meters per second) */
Real U_max = 2 * U_f; /**< Maximum characteristic velocity, unit: m/s (meters per second) */
/**< Reference sound speed considering flow speed in narrow channels, unit: m/s (meters per second) */
Real c_f = 10.0 * U_f * SMAX(Real(1), d_inlet *d_inlet / (d_outlet_large * d_outlet_large + d_outlet_middle * d_outlet_middle + d_outlet_small * d_outlet_small));
Real mu_f = 0.0035;         /**< Dynamic viscosity, unit: Pa·s (Pascal-second) */
Real Outlet_pressure = 0.0; /**< Outlet pressure for comparison with solely velocity inlet boundary condition, unit: Pa (Pascal) */
// Real Outlet_pressure = 2.666e3;
//----------------------------------------------------------------------
// Global parameters on the solid properties (Vessel wall parameters)
//----------------------------------------------------------------------
Real rho0_s_vessel = 1080;              /**< Density of the vessel wall, unit: kg/m³ (kilograms per cubic meter) */
Real poisson_vessel = 0.45;             /**< Poisson's ratio of the vessel wall, dimensionless (unitless) */
Real Youngs_modulus_vessel = 1e6;       /**< Young's modulus of the vessel wall, unit: Pa (Pascal) */
Real physical_viscosity_vessel = 100.0; /**< Physical viscosity of the vessel wall, unit: Pa·s (Pascal-second) */
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
// Inflow velocity (Steady-state condition)
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_; // Constant reference velocity
    AlignedBox &aligned_box_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(0.315), aligned_box_(boundary_condition.getAlignedBox()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;
        target_velocity[0] = u_ref_; // Set constant velocity in x-direction
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
        add<TriangleMeshShapeSTL>(vessel_fluid_path, translation_water_block, length_scale_vessel);
    }
};
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(vessel_wall_path, translation_wall_boundary, length_scale_vessel);
    }
};
//----------------------------------------------------------------------
//	BoundaryGeometry.
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
//----------------------------------------------------------------------
// MaxStress
//----------------------------------------------------------------------
struct MaxStressResult
{
    Real stress;
    size_t index;
};

// Modified MaxStressLocalDynamics class, setting ReturnType as MaxStressResult
class MaxStressLocalDynamics
{
  public:
    using ReturnType = MaxStressResult;

    // To allow ReduceDynamics to call LoopRange() and other interfaces, set identifier_ as public
    SPHBody &identifier_;

    /**
     * Constructor.
     * @param body Target SPHBody object (e.g., wall_boundary), whose BaseParticles must have updated "VonMisesStress" data
     */
    explicit MaxStressLocalDynamics(SPHBody &body)
        : identifier_(body),
          output_folder_("./output"),
          name_prefix_("max_stress")
    {
        BaseParticles &particles = identifier_.getBaseParticles();
        stress_ = particles.getVariableDataByName<Real>("VonMisesStress");
        if (!stress_)
        {
            std::cerr << "Error: 'VonMisesStress' not found in BaseParticles! "
                      << "Ensure it has been created and updated." << std::endl;
            std::exit(1);
        }
        pos_ = particles.getVariableDataByName<Vecd>("Position");
        if (!pos_)
        {
            std::cerr << "Error: 'Position' not found in BaseParticles!" << std::endl;
            std::exit(1);
        }
        // Get system physical time pointer
        physical_time_ = identifier_.getSPHSystem().getSystemVariableDataByName<Real>("PhysicalTime");
    }

    // setupDynamics: Preprocessing before reduction, no operation needed here
    void setupDynamics(Real /*dt*/) {}

    // LoopRange: Returns the index range of all particles in wall_boundary
    const IndexRange &LoopRange() const
    {
        return identifier_.LoopRange();
    }

    /**
     * Returns a MaxStressResult for each particle,
     * where stress is the von Mises stress of the particle,
     * and index is the particle index.
     */
    ReturnType reduce(size_t index, Real /*dt*/ = 0.0)
    {
        return {stress_[index], index};
    }

    /**
     * Initial reduction value, set to the lowest possible value and an arbitrary index (set to 0 here).
     */
    ReturnType Reference()
    {
        return {std::numeric_limits<Real>::lowest(), 0};
    }

    /**
     * Returns the reduction operation as a lambda expression,
     * which compares two MaxStressResult structures and returns the one with the larger stress value.
     */
    auto getOperation()
    {
        return [](const MaxStressResult &a, const MaxStressResult &b) -> MaxStressResult
        {
            return (a.stress > b.stress) ? a : b;
        };
    }

    /**
     * Outputs the reduction result to a file and returns the reduction result.
     * Output file format:
     * First column: Time [s]; Second column: Maximum von Mises stress [Pa];
     * Third column: Particle coordinates (x, y, z), with components separated by spaces.
     */
    ReturnType outputResult(ReturnType result)
    {
        std::filesystem::create_directory(output_folder_);
        std::string file_name = output_folder_ + "/" + name_prefix_ + ".dat";
        std::ofstream out(file_name, std::ios::app);
        if (out.is_open())
        {
            // Output time, maximum stress, and the corresponding particle position
            out << std::setprecision(9) << *physical_time_ << "\t" << result.stress << "\t";
            // Assuming Vecd is a 3D vector, output its components
            out << pos_[result.index][0] << " " << pos_[result.index][1] << " " << pos_[result.index][2] << "\n";
            out.close();
        }
        else
        {
            std::cerr << "Error: Cannot open file " << file_name << " for writing." << std::endl;
        }
        return result;
    }

  private:
    Real *stress_;              ///< Pointer to "VonMisesStress" data
    Vecd *pos_;                 ///< Pointer to particle "Position" data
    Real *physical_time_;       ///< Pointer to system physical time data
    std::string output_folder_; ///< Output folder
    std::string name_prefix_;   ///< File name prefix (final file name: name_prefix_.dat)
};
using MaxStressCalculator = ReduceDynamics<MaxStressLocalDynamics, SPH::execution::ParallelPolicy>;

class PressureCalculator
    : public BaseLocalDynamicsReduce<ReduceSum<Real>, BodyPartByCell>
{
  public:
    using ReturnType = Real; // Reduction return type

    /**
     * Constructor.
     * @param aligned_box_part   The BodyAlignedBoxByCell object corresponding to the inlet/outlet region.
     * @param prefix             The prefix used for output file names (e.g., "inlet", "outlet_large", "outlet_small").
     */
    explicit PressureCalculator(BodyAlignedBoxByCell &aligned_box_part,
                                const std::string &prefix)
        : BaseLocalDynamicsReduce<ReduceSum<Real>, BodyPartByCell>(aligned_box_part),
          transient_pressure_(0.0),
          average_pressure_(0.0),
          name_prefix_(prefix)
    {
        // Retrieve pressure and volumetric measure data from BaseParticles.
        pressure_ = this->particles_->getVariableDataByName<Real>("Pressure");
        Vol_ = this->particles_->getVariableDataByName<Real>("VolumetricMeasure");

        // Get the system physical time pointer.
        physical_time_ = this->sph_system_.getSystemVariableDataByName<Real>("PhysicalTime");
    }

    virtual ~PressureCalculator() {}

    /**
     * For each particle, return its contribution to the total pressure sum:
     * contribution = Pressure * Volume
     */
    virtual Real reduce(size_t index, Real dt = 0.0)
    {
        return pressure_[index] * Vol_[index];
    }

    /**
     * After the reduction, compute the volume-weighted average pressure:
     * P_avg = ∑(P * V) / ∑V,
     * then write the results to files.
     */
    virtual Real outputResult(Real reduced_value) override
    {
        Real total_volume = 0.0;
        auto &cell_lists = this->getDynamicsIdentifier().LoopRange();
        for (size_t c = 0; c < cell_lists.size(); ++c)
        {
            for (size_t i = 0; i < cell_lists[c]->size(); ++i)
            {
                size_t global_index = (*cell_lists[c])[i];
                total_volume += Vol_[global_index];
            }
        }
        if (total_volume > 0.0)
        {
            average_pressure_ = reduced_value / total_volume;
        }
        else
        {
            average_pressure_ = 0.0;
        }

        // Output results to files.
        std::string output_folder = "./output";
        std::filesystem::create_directory(output_folder);

        // Pressure file.
        std::string pressure_file = output_folder + "/pressure_" + name_prefix_ + ".dat";
        {
            std::ofstream pressure_out(pressure_file, std::ios::app);
            if (pressure_out.is_open())
            {
                pressure_out << *physical_time_ << "\t" << average_pressure_ << "\n";
            }
        }

        return average_pressure_;
    }

  private:
    // Prefix for naming output files.
    std::string name_prefix_;

    // Intermediate variables during calculation.
    Real transient_pressure_;
    Real average_pressure_;

    // Pointers to particle data.
    Real *pressure_;
    Real *Vol_;

    // Pointer to the system physical time.
    Real *physical_time_;
};
