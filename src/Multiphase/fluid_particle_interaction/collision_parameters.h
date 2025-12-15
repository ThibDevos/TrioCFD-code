#ifndef COLLISION_PARAMETERS
#define COLLISION_PARAMETERS
#include <TRUSTTabFT.h>
#include <limits>
struct collision_parameters
{
  int particle_i =-1;
  int particle_j =-1;
  int i_closest =-1;
  int j_closest =-1;
  bool i_facet = false;
  bool j_facet = false;
  bool i_j_are_close = false;
  bool part_part_collision = false;
  DoubleTab dX;
  double distance = std::numeric_limits<double>::max();
  collision_parameters(int dimension) {dX.resize(dimension);}
  collision_parameters() {dX.resize(3);}
};
#endif