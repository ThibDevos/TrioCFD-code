#ifndef PARTICLES_PROPERTIES
#define PARTICLES_PROPERTIES

#include <TRUSTTabs_forward.h>
/* Small structure to regroup particle properties that are usually passed as arguments*/
struct particle_properties
{
  particle_properties(DoubleTab& particles_position_, DoubleTab& particles_velocity_,
                      DoubleTab& particles_rot_velocity_, DoubleTab& inertia_tensor_, DoubleTab& volume_)
    : position(particles_position_), velocity(particles_velocity_), rot_velocity(particles_rot_velocity_),
      inertia_tensor(inertia_tensor_), volume(volume_) { ; }

  particle_properties(DoubleTab const& particles_position_, DoubleTab const& particles_velocity_,
                      DoubleTab const& particles_rot_velocity_, DoubleTab const& inertia_tensor_, DoubleTab const& volume_)
    : position(particles_position_), velocity(particles_velocity_), rot_velocity(particles_rot_velocity_),
      inertia_tensor(inertia_tensor_), volume(volume_) { ; }

  DoubleTab const& position;
  DoubleTab const& velocity;
  DoubleTab const& rot_velocity;
  DoubleTab const& inertia_tensor;
  DoubleTab const& volume;
};

#endif