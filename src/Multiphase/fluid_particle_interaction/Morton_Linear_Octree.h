#ifndef MORTON_LINEAR_OCTREE
#define MORTON_LINEAR_OCTREE
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <TRUSTTabFT.h>

class Morton_Linear_Octree
{
public:
  Morton_Linear_Octree(const DoubleTab& sommets, double h);
  void float_coordinates_to_int(DoubleTab const& sommets, std::vector<std::array<uint32_t,3>>& int_coords, double h);
  void int_coordinates_to_Morton(std::vector<std::array<uint32_t,3>>& int_coords);
  void build_leaves();
  void compute_neighbours();




  //Methods relative to Morton code
  uint64_t expand_bits(uint32_t x);
  uint32_t compact_bits(uint64_t x);
  uint64_t encodeMorton(uint32_t x, uint32_t y, uint32_t z);
  void DecodeMorton(uint64_t code, uint32_t& x, uint32_t& y, uint32_t& z);
  // void build();

private:
  static constexpr int dimension = 3; //This has been implemented for 3D XXX
  static constexpr uint32_t L = 21; //each coordinate is represented by 21 bits (so the Morton code is 3*21=63 bits)
  int l=5;  //the cube will be devided in 2^l*2^l*2^l cubes (called leaves) in which the vertices are

  struct leaf
  {
    uint64_t Morton_code; //coordinates of the leaf
    std::vector<int> vertex_indices; //indices of the vertices in it (indices correspond to the indices of the mesh vertices)
    std::vector<int> closest_vertex_indices; //closest_vertex_indices[i] is the index of the closest vertex of vertex_indices[i]
    std::vector<leaf*> neighbours; //pointers to the neighbour leaves (26 in 3d)
  };

  std::vector<leaf> leaves;
  std::vector<uint64_t> Morton_code; // Morton code is a way to express in 1D, 3D coordinates (integer)
  std::vector<int> indices; //Morton_code will be sorted. indices helps keeping a relation of indices between Morton_code and mesh.sommets().
  // indices[i] is the original index of Morton_code[i]. Is identity before Morton_code being sorted
};

#endif