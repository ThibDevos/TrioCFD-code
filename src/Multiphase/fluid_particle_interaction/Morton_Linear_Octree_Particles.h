#ifndef MORTON_LINEAR_OCTREE_PARTICLES
#define MORTON_LINEAR_OCTREE_PARTICLES
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <collision_parameters.h>
#include <TRUSTTabFT.h>
#include <fstream>

class Morton_Linear_Octree_Particles
{
public:
  Morton_Linear_Octree_Particles(const DoubleTab& sommets, double h,   const ArrOfInt& compo_connexes_sommets);
  void float_coordinates_to_int(DoubleTab const& sommets, std::vector<std::array<uint32_t,3>>& int_coords, double h);
  void int_coordinates_to_Morton(std::vector<std::array<uint32_t,3>>& int_coords);
  void build_leaves(const ArrOfInt& compo_connexes_sommets);
  void compute_neighbours();

  void find_closest(const DoubleTab& sommets, const ArrOfInt& compo_connexes_sommets, std::vector<std::vector<collision_parameters>>& param, DoubleVect const& origin, DoubleVect const& domain_dimensions);


  //Methods relative to Morton code
  uint64_t expand_bits(uint32_t x);
  uint32_t compact_bits(uint64_t x);
  uint64_t encodeMorton(uint32_t x, uint32_t y, uint32_t z);
  void DecodeMorton(uint64_t code, uint32_t& x, uint32_t& y, uint32_t& z);
  // void build();

  void print_indices(const ArrOfInt& compo_connexes_sommets);
  void which_leaf(int s);

private:
  static constexpr int dimension = 3; //This has been implemented for 3D XXX
  static constexpr uint32_t L = 21; //each coordinate is represented by 21 bits (so the Morton code is 3*21=63 bits)
  int l=0;  //the cube will be devided in 2^l*2^l*2^l cubes (called leaves) in which the vertices are

  uint32_t grid_limits[3];

  struct leaf
  {
    uint64_t Morton_code; //coordinates of the leaf
    std::vector<int> vertex_indices; //indices of the vertices in it (indices correspond to the indices of the mesh vertices)
    std::vector<int> closest_vertex_indices; //closest_vertex_indices[i] is the index of the closest vertex of vertex_indices[i]
    std::vector<leaf*> neighbours; //pointers to the neighbour leaves (26 in 3d)
    int compo = -1; //compo is the index of the particle of the vertices inside leaf. Is equal to -1 if it is a mix
    int boundary = 0; // boundary \in [0,63], such that, each bit correspond to a boundary and is 1 if the leaf is close to it.
    //  The bits correspond, from right to left to x_-, y_-, z_-, x_+, y_+, z_+
  };

  std::vector<leaf> leaves;
  std::vector<uint64_t> Morton_code; // Morton code is a way to express in 1D, 3D coordinates (integer)
  std::vector<int> indices; //Morton_code will be sorted. indices helps keeping a relation of indices between Morton_code and mesh.sommets().
  // indices[i] is the original index of Morton_code[i]. Is identity before Morton_code being sorted
};

#endif