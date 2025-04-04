
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
std::string stent_path = "./input/stent.stl";
std::string vessel_wall_path = "./input/vessel_wall_stenosis.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Vec3d translation_stent(0.0, 0.0, 0.0);                                   /**< Initial translation of the stent, unit: meters (m) */
Vec3d translation_wall_boundary(0.0, 0.0, 0.0);                           /**< Initial translation of the vessel wall, unit: meters (m) */
Real length_scale_stent = pow(10, -3);                                    /**< Length scale factor, dimensionless (unitless) */
Real length_scale_vessel = 1.0;                                   /**< Length scale factor, dimensionless (unitless) */
Real length_scale = pow(10, -3);                                          /**< Length scale factor, dimensionless (unitless) */
Real resolution_ref = 0.2 * length_scale;                                 /**< Initial reference particle spacing, unit: meters (m) */
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
Real A_inlet = 12.67e-6;
Vec3d inlet_half(2.0 * resolution_ref, 3.0 * length_scale, 3.0 * length_scale);
Vec3d inlet_fix_half(2.0 * resolution_ref, 7.0 * length_scale, 7.0 * length_scale);
Vec3d inlet_normal(1.0, 0, 0);
Vec3d inlet_buffer_translation = Vec3d(0, 0, 0) * length_scale + inlet_normal * 2.0 * resolution_ref;
RotationCalculator inlet_rotation_calculator(inlet_normal, x_Axis);
Rotation3d inlet_emitter_rotation(inlet_rotation_calculator.getRotationAngle(), inlet_rotation_calculator.getRotationAxis());
Rotation3d inlet_disposer_rotation(inlet_rotation_calculator.getRotationAngle() + M_PI, inlet_rotation_calculator.getRotationAxis());

// outlet_large Parameter: d=3.5, (24.26695907, -9.65881011, -3.31277237), (22.99372438, -8.15633227, -2.96429891)
Real d_outlet_large = 3.5 * length_scale;
Real A_outlet_large = 9.295576e-6;
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
Real A_outlet_middle = 7.033589e-6;
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
Real A_outlet_small = 5.021522e-6;
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
// Global parameters on the solid properties (Vessel wall parameters)
//----------------------------------------------------------------------
Real rho0_s_vessel = 1080;              /**< Density of the vessel wall, unit: kg/m³ (kilograms per cubic meter) */
Real poisson_vessel = 0.45;             /**< Poisson's ratio of the vessel wall, dimensionless (unitless) */
Real Youngs_modulus_vessel = 1e6;       /**< Young's modulus of the vessel wall, unit: Pa (Pascal) */
Real physical_viscosity_vessel = 100.0; /**< Physical viscosity of the vessel wall, unit: Pa·s (Pascal-second) */
//----------------------------------------------------------------------
//	Global parameters for stent (Stent parameters)
//----------------------------------------------------------------------
Real rho0_s_stent = 6450.0;
Real poisson_stent = 0.3;
Real youngs_modulus_stent = 2e10;
Real physical_viscosity_stent = 100.0;

Real yield_stress_stent = 8e7;
Real hardening_modulus_stent = 5e8;
Real saturation_flow_stress_stent = 2.5e8;
Real saturation_exponent_stent = 4.0;
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
class Stent : public ComplexShape
{
  public:
    explicit Stent(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(stent_path, translation_stent, length_scale_stent);
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
//----------------------------------------------------------------------
// BoundingBox Calculation Functions
//----------------------------------------------------------------------
class BoundingBoxCalculator
{
  public:
    /**
     * @brief Obtain the Bounding Box in the global coordinate system and print information
     */
    static BoundingBox getGlobalBoundingBox(BaseParticles &particles, bool print_info = true)
    {
        Vecd lower_bound = particles.ParticlePositions()[0]; // Initialize to the position of the first particle
        Vecd upper_bound = lower_bound;

        const Vecd *positions = particles.ParticlePositions();
        for (size_t i = 1; i < particles.TotalRealParticles(); ++i)
        {
            lower_bound = lower_bound.cwiseMin(positions[i]);
            upper_bound = upper_bound.cwiseMax(positions[i]);
        }

        Vecd bbox_size = upper_bound - lower_bound;

        if (print_info)
        {
            printBoundingBoxInfo(lower_bound, upper_bound, bbox_size, "Global");
        }

        return BoundingBox(lower_bound, upper_bound);
    }

    /**
     * @brief Obtain the Bounding Box in the stent's own coordinate system
     * @details Requires the rotation matrix and the stent's center point
     */
    static BoundingBox getLocalBoundingBox(BaseParticles &particles, const Mat3d &rotation_matrix, const Vecd &center, bool print_info = true)
    {
        Mat3d inv_rotation = rotation_matrix.transpose(); // Compute the inverse rotation matrix
        Vecd inv_translation = -center;                   // Translate to the local coordinate system

        Vecd local_lower_bound = Vecd::Constant(std::numeric_limits<Real>::max());
        Vecd local_upper_bound = Vecd::Constant(-std::numeric_limits<Real>::max());

        const Vecd *positions = particles.ParticlePositions();
        for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        {
            Vecd local_pos = inv_rotation * (positions[i] + inv_translation);
            local_lower_bound = local_lower_bound.cwiseMin(local_pos);
            local_upper_bound = local_upper_bound.cwiseMax(local_pos);
        }

        Vecd bbox_size = local_upper_bound - local_lower_bound;

        if (print_info)
        {
            printBoundingBoxInfo(local_lower_bound, local_upper_bound, bbox_size, "Local (Stent)");
        }

        return BoundingBox(local_lower_bound, local_upper_bound);
    }

    /**
     * @brief Print Bounding Box information
     */
    static void printBoundingBoxInfo(const Vecd &lower_bound, const Vecd &upper_bound, const Vecd &size, const std::string &type)
    {
        std::cout << "========================================\n";
        std::cout << "         " << type << " Bounding Box\n";
        std::cout << "----------------------------------------\n";
        std::cout << "Lower Bound:  (" << lower_bound[0] << ", " << lower_bound[1] << ", " << lower_bound[2] << ")\n";
        std::cout << "Upper Bound:  (" << upper_bound[0] << ", " << upper_bound[1] << ", " << upper_bound[2] << ")\n";
        std::cout << "----------------------------------------\n";
        std::cout << "Bounding Box Size:\n";
        std::cout << "  X-axis: [" << lower_bound[0] << ", " << upper_bound[0] << "]  (Δx: " << size[0] << ")\n";
        std::cout << "  Y-axis: [" << lower_bound[1] << ", " << upper_bound[1] << "]  (Δy: " << size[1] << ")\n";
        std::cout << "  Z-axis: [" << lower_bound[2] << ", " << upper_bound[2] << "]  (Δz: " << size[2] << ")\n";
        std::cout << "========================================\n";
    }
};
//----------------------------------------------------------------------
// Modified RadialForce: Supports Arbitrary Direction
//----------------------------------------------------------------------
class RadialForce
{
  protected:
    Real magnitude_;        // Magnitude of the radial force
    Vecd center_;           // Center of the support
    Mat3d rotation_matrix_; // Rotation matrix (global -> local)
    Vecd axis_direction_;   // Axis direction of the support (global)

  public:
    // **Constructor: input radial force magnitude and endpoints of the support**
    RadialForce(Real magnitude, const Vecd &pointA, const Vecd &pointB)
        : magnitude_(magnitude), center_((pointA + pointB) * 0.5)
    {
        // **Calculate the main axis direction of the support**
        axis_direction_ = (pointB - pointA).normalized();

        // **Construct the rotation matrix such that the local X axis aligns with the support**
        Vecd local_x = axis_direction_;
        Vecd local_z = Vecd(0, 0, 1);              // Assume the Z axis as the reference vector
        if (std::abs(local_x.dot(local_z)) > 0.99) // Prevent parallel with Z axis
        {
            local_z = Vecd(0, 1, 0);
        }
        Vecd local_y = local_z.cross(local_x).normalized();
        local_z = local_x.cross(local_y).normalized();
        rotation_matrix_ << local_x, local_y, local_z;
    }

    // **Calculate the radial direction in the local coordinate system and apply the force**
    Vecd InducedAcceleration(const Vecd &global_position, Real physical_time = 0.0) const
    {
        // **Convert to the local coordinate system**
        Vecd local_position = rotation_matrix_.transpose() * (global_position - center_);

        // **Calculate the radial component in the local YZ plane**
        Vecd radial_direction = local_position;
        radial_direction[0] = 0.0; // **No force applied in the X direction**
        radial_direction.normalize();

        // **Apply force: return the radial force in global coordinates**
        Vecd global_force = rotation_matrix_ * (magnitude_ * radial_direction);

        return global_force;
    }

    // **Return the total current radial force being applied**
    Vecd GetCurrentForce() const
    {
        Vecd local_force = Vecd(0.0, magnitude_, magnitude_); // Force in local YZ directions
        Vecd global_force = rotation_matrix_ * local_force;
        return global_force;
    }

    Mat3d getRotationMatrix() const
    {
        return rotation_matrix_;
    }
};

class FastStartRadialForce : public RadialForce
{
    Real target_time_; // Time to reach maximum force
    Real growth_rate_; // Control growth rate

  public:
    FastStartRadialForce(Real magnitude, const Vecd &pointA, const Vecd &pointB, Real target_time, Real growth_rate = 5.0)
        : RadialForce(magnitude, pointA, pointB), target_time_(target_time), growth_rate_(growth_rate) {}
    ~FastStartRadialForce() {}

    /** Calculate the current force amplification factor */
    Real GetScalingFactor(Real physical_time) const
    {
        if (physical_time < target_time_)
        {
            return 1.0 - exp(-growth_rate_ * (physical_time / target_time_));
        }
        return 1.0; // After reaching target time, maintain maximum force
    }

    /** Calculate the acceleration acting on the particle */
    Vecd InducedAcceleration(const Vecd &global_position, Real physical_time) const
    {
        Real scaling_factor = GetScalingFactor(physical_time);
        return scaling_factor * RadialForce::InducedAcceleration(global_position, physical_time);
    }

    /** Get the current total force being applied */
    Vecd GetCurrentForce(Real physical_time) const
    {
        Real scaling_factor = GetScalingFactor(physical_time);
        return scaling_factor * RadialForce::GetCurrentForce();
    }
};

class StartupRadialForce : public RadialForce
{
    Real target_time_; // Time period for force increase

  public:
    StartupRadialForce(Real magnitude, const Vecd &pointA, const Vecd &pointB, Real target_time)
        : RadialForce(magnitude, pointA, pointB), target_time_(target_time) {}
    ~StartupRadialForce() {}

    /** Calculate the current force amplification factor */
    Real GetScalingFactor(Real physical_time) const
    {
        if (physical_time < target_time_)
        {
            return sin(PI / 2 * (physical_time / target_time_)); // Smooth growth
        }
        return 1.0; // After reaching target time, maintain maximum force
    }

    /** Calculate the acceleration acting on the particle */
    Vecd InducedAcceleration(const Vecd &global_position, Real physical_time) const
    {
        Real scaling_factor = GetScalingFactor(physical_time);
        return scaling_factor * RadialForce::InducedAcceleration(global_position, physical_time);
    }

    /** Get the current total force being applied */
    Vecd GetCurrentForce(Real physical_time) const
    {
        Real scaling_factor = GetScalingFactor(physical_time);
        return scaling_factor * RadialForce::GetCurrentForce();
    }
};

class IncreasingRadialForce : public RadialForce
{
    Real target_time_; // Target time; within this time, linearly increase to maximum force

  public:
    IncreasingRadialForce(Real max_magnitude, const Vecd &pointA, const Vecd &pointB, Real target_time)
        : RadialForce(max_magnitude, pointA, pointB), target_time_(target_time) {}

    // **Modify InducedAcceleration method, so the force increases over 0.02s**
    Vecd InducedAcceleration(const Vecd &global_position, Real physical_time = 0.0) const
    {
        Real scaling_factor = (physical_time < target_time_) ? (physical_time / target_time_) : 1.0;
        return scaling_factor * RadialForce::InducedAcceleration(global_position, physical_time);
    }

    // **Return the current magnitude of the radial force**
    Vecd GetCurrentForce(Real physical_time) const
    {
        Real scaling_factor = (physical_time < target_time_) ? (physical_time / target_time_) : 1.0;
        return scaling_factor * RadialForce::GetCurrentForce();
    }
};
//----------------------------------------------------------------------
//  RadialForceApplication.
//----------------------------------------------------------------------
template <class RadialForceType>
class RadialForceApplication : public ForcePrior
{
  protected:
    RadialForceType &radial_force_;
    Vecd *pos_;
    Real *mass_;
    Real *physical_time_;

    // ✅ **New variable: for storing debug information of selected particles**
    std::vector<size_t> sampled_indices_; // Stores indices of selected particles

  public:
    RadialForceApplication(SPHBody &sph_body, RadialForceType &radial_force)
        : ForcePrior(sph_body, "RadialForce"), radial_force_(radial_force),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          mass_(particles_->registerStateVariable<Real>("Mass")),
          physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime"))
    {
        // ✅ **Randomly select some particles for debugging**
        size_t total_particles = particles_->TotalRealParticles();
        for (size_t i = 0; i < total_particles; i += std::max(1, (int)(total_particles / 10))) // Select one particle every 10%
        {
            sampled_indices_.push_back(i);
        }
    }

    virtual ~RadialForceApplication() {}

    void update(size_t index_i, Real dt = 0.0)
    {
        Vecd global_position = pos_[index_i];
        Vecd applied_force = mass_[index_i] * radial_force_.InducedAcceleration(global_position, *physical_time_);

        current_force_[index_i] = applied_force;
        ForcePrior::update(index_i, dt);
    }

    /**
     * @brief ✅ **Calculate and print the total radial force applied to the support**
     */
    void printAppliedForce() const
    {
        Vecd total_force = Vecd::Zero();

        for (size_t i = 0; i < particles_->TotalRealParticles(); ++i)
        {
            total_force += mass_[i] * radial_force_.InducedAcceleration(pos_[i], *physical_time_);
        }

        // ✅ **Use radial_force_ to obtain the rotation matrix**
        Mat3d rotation_matrix = radial_force_.getRotationMatrix();
        Vecd local_force = rotation_matrix.transpose() * total_force;

        std::cout << "========================================\n";
        std::cout << "         Applied Radial Force\n";
        std::cout << "----------------------------------------\n";
        std::cout << "Total Force (Stent Frame):\n";
        std::cout << "  X: " << local_force[0] << " N\n";
        std::cout << "  Y: " << local_force[1] << " N\n";
        std::cout << "  Z: " << local_force[2] << " N\n";
        std::cout << "  |F|: " << local_force.norm() << " N (Magnitude)\n";
        std::cout << "========================================\n\n";
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
    // Constructor, pass in the SPHBody whose state needs to be recorded
    ReloadParticleRecordingToXml(SPHBody &sph_body)
        : BaseIO(sph_body.getSPHSystem()), sph_body_(sph_body), base_particles_(sph_body.getBaseParticles())
    {
        // Ensure the output folder is correctly set
        output_folder_ = io_environment_.output_folder_ + "/particle-reload";
        if (!fs::exists(output_folder_))
        {
            fs::create_directories(output_folder_);
        }
    }

    // Public function: write the current particle state to an XML file
    void writeToFile(size_t iteration_step)
    {
        // Convert iteration step to string for file naming
        std::string sequence = std::to_string(iteration_step);
        // Construct the file path for saving the particle reload data (in XML format)
        std::string filefullpath = output_folder_ + "/particle_reload_" + sph_body_.getName() + "_" + sequence + ".xml";

        // If the file already exists, remove the old file
        if (fs::exists(filefullpath))
        {
            fs::remove(filefullpath);
        }

        // Open the output file stream
        std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);

        // Begin writing the XML structure
        out_file << "<?xml version=\"1.0\"?>\n";
        out_file << "<particles>\n";

        // Iterate over all real particles and write the data
        size_t total_real_particles = base_particles_.TotalRealParticles();
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            Vecd position = base_particles_.ParticlePositions()[i];
            Real volume = base_particles_.VolumetricMeasures()[i];

            // Write each particle's attributes in a single line
            out_file << "  <particle VolumetricMeasure=\"" << volume << "\" Position=\""
                     << position[0] << ", " << position[1];

            // If it is a 3D project, add the third coordinate
            if (position.size() == 3)
            {
                out_file << ", " << position[2];
            }

            out_file << "\"/>\n";
        }

        // Close the XML structure
        out_file << "</particles>\n";
        out_file.close();

        // Debug: Confirm that the file has been written
        std::cout << "Particle state for " << sph_body_.getName() << " written to " << filefullpath << std::endl;
    }

  private:
    SPHBody &sph_body_;             // The SPHBody whose particle state we are writing
    BaseParticles &base_particles_; // Particles of the body
    std::string output_folder_;     // Folder where the output file will be written
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
    Real *mass_; // Pointer to particle mass data

  public:
    explicit QuantityMassPosition(DynamicsIdentifier &identifier)
        : QuantitySummation<Vecd, DynamicsIdentifier>(identifier, "Position"),
          mass_(this->particles_->template getVariableDataByName<Real>("Mass"))
    {
        this->quantity_name_ = "MassWeightedPosition";
    };
    virtual ~QuantityMassPosition() {}

    // Override the reduce method to compute the mass-weighted position
    Vecd reduce(size_t index_i, Real dt = 0.0)
    {
        // Get the particle's position and mass
        Vecd position = this->variable_[index_i];
        Real mass = mass_[index_i];

        // Compute the mass-weighted position
        return position * mass;
    }
};

Vecd computeAveragePosition(SolidBody &stent_body)
{
    // Get the BaseParticles instance from stent_body
    BaseParticles &particles = stent_body.getBaseParticles();

    // Get the total number of particles
    size_t total_particles = particles.TotalRealParticles();

    // If there are no particles, return the zero vector
    if (total_particles == 0)
    {
        return Vecd::Zero();
    }

    // Get the array of particle positions
    Vecd *positions = particles.ParticlePositions();

    // Initialize the total position
    Vecd total_position = Vecd::Zero();

    // Iterate over all particles and accumulate the positions
    for (size_t i = 0; i < total_particles; ++i)
    {
        total_position += positions[i];
    }

    // Calculate the average position
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
    Vecd mass_center_; // Center of mass position
    size_t p_1_, p_2_; // Indices for the inertia matrix components
    Vecd *positions_;  // Particle position data

  public:
    explicit QuantityMomentOfInertia(DynamicsIdentifier &identifier, Vecd mass_center, size_t p_1, size_t p_2)
        : QuantitySummation<Real, DynamicsIdentifier>(identifier, "Mass"),
          mass_center_(mass_center), p_1_(p_1), p_2_(p_2),
          positions_(this->particles_->template getVariableDataByName<Vecd>("Position"))
    {
        this->quantity_name_ = "MomentOfInertia";

        // Check whether positions_ is successfully initialized
        if (!positions_)
        {
            throw std::runtime_error("Error: Unable to initialize positions_. Check if 'Position' variable exists.");
        }
    }

    virtual ~QuantityMomentOfInertia() {}

    // Override the reduce method to compute the inertia matrix component
    Real reduce(size_t index_i, Real dt = 0.0)
    {
        Vecd relative_position = positions_[index_i] - mass_center_; // Access the particle position
        Real mass = this->variable_[index_i];                        // Get the particle mass

        if (p_1_ == p_2_)
        {
            // Diagonal inertia component
            return mass * (relative_position.squaredNorm() - relative_position[p_1_] * relative_position[p_1_]);
        }
        else
        {
            // Off-diagonal inertia component
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
 * @class Constrain3DSolidBodyRotation
 * @brief Constrain the rotation of a 3D solid body.
 */
class Constrain3DSolidBodyRotation : public LocalDynamics
{
  private:
    Vecd mass_center_;                                                          // Center of mass of the body
    Matd moment_of_inertia_;                                                    // Inertia tensor
    Vecd angular_velocity_;                                                     // Current angular velocity
    ReduceDynamics<QuantityMomentOfMomentum> compute_total_moment_of_momentum_; // Calculate total angular momentum
    Vecd *positions_;                                                           // Array of particle positions
    Vecd *velocities_;                                                          // Array of particle velocities

  protected:
    // At the start of each time step, calculate the angular velocity
    virtual void setupDynamics(Real dt = 0.0) override
    {
        angular_velocity_ = moment_of_inertia_.inverse() * compute_total_moment_of_momentum_.exec(dt);
    }

  public:
    // Constructor, initialize variables
    explicit Constrain3DSolidBodyRotation(SPHBody &sph_body, const Vecd &mass_center, const Matd &inertia_tensor)
        : LocalDynamics(sph_body),
          mass_center_(mass_center),
          moment_of_inertia_(inertia_tensor),
          compute_total_moment_of_momentum_(sph_body, mass_center)
    {
        // Obtain particle data
        BaseParticles &particles = sph_body.getBaseParticles();
        positions_ = particles.getVariableDataByName<Vecd>("Position");
        velocities_ = particles.getVariableDataByName<Vecd>("Velocity");

        // Check if the variables are correctly initialized
        if (!positions_ || !velocities_)
        {
            throw std::runtime_error("Error: Unable to initialize 'Position' or 'Velocity' variable. Check if they are registered.");
        }
    }

    virtual ~Constrain3DSolidBodyRotation() {}

    // For each particle, update the velocity to constrain rotation
    void update(size_t index_i, Real dt = 0.0)
    {
        Vecd relative_position = positions_[index_i] - mass_center_;                       // Compute the position relative to the center of mass
        Vecd linear_velocity_due_to_rotation = angular_velocity_.cross(relative_position); // Calculate linear velocity due to rotation
        velocities_[index_i] -= linear_velocity_due_to_rotation;                           // Adjust velocity to remove the rotational component
    }
};

//----------------------------------------------------------------------
//	CustomContactFactorSummation.
//----------------------------------------------------------------------
class CustomContactFactorSummation : public solid_dynamics::RepulsionFactorSummation<Contact<>>
{
  public:
    explicit CustomContactFactorSummation(SurfaceContactRelation &contact_relation)
        : RepulsionFactorSummation<Contact<>>(contact_relation) {}

    // Add a method to modify offset_W_ij_
    void scaleOffsetW(Real factor)
    {
        for (size_t k = 0; k < offset_W_ij_.size(); ++k)
        {
            offset_W_ij_[k] *= factor;
        }
    }
};

//----------------------------------------------------------------------
//	CustomContactForce.
//----------------------------------------------------------------------
class CustomContactForce : public solid_dynamics::RepulsionForce<SPH::Contact<>>
{
  public:
    explicit CustomContactForce(BaseContactRelation &solid_body_contact_relation)
        : RepulsionForce<SPH::Contact<>>(solid_body_contact_relation) {}

    void modifyContactStiffness(Real factor)
    {
        for (size_t k = 0; k < contact_stiffness_ave_.size(); ++k)
        {
            contact_stiffness_ave_[k] *= factor;
        }
    }

    void modifyForceK(Real softening)
    {
        size_t total_particles = particles_->TotalRealParticles(); // Get the number of particles

        for (size_t index_i = 0; index_i < total_particles; ++index_i)
        {
            repulsion_force_[index_i] = Vecd::Zero(); // Ensure initialization
            for (size_t k = 0; k < contact_configuration_.size(); ++k)
            {
                Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
                for (size_t n = 0; n < contact_neighborhood.current_size_; ++n)
                {
                    size_t index_j = contact_neighborhood.j_[n];
                    Vecd e_ij = contact_neighborhood.e_ij_[n];

                    Real sigma_star = 0.5 * (repulsion_factor_[index_i] + contact_repulsion_factor_[k][index_j]);
                    repulsion_force_[index_i] -= 2.0 * sigma_star * e_ij *
                                                 (contact_neighborhood.dW_ij_[n] + softening) * contact_Vol_[k][index_j];
                }
            }
        }
    }
};