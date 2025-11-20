#include <Morton_Linear_Octree_Particles.h>


Morton_Linear_Octree_Particles::Morton_Linear_Octree_Particles(const DoubleTab& sommets, double h,  const ArrOfInt& compo_connexes_sommets)
{
  int nb_sommets=sommets.dimension(0);
  indices.resize(nb_sommets);
  Morton_code.resize(nb_sommets);
  std::iota(indices.begin(), indices.end(), 0);
  std::vector<std::array<uint32_t,3>> int_coords(nb_sommets);

  //transform the floating points coordiantes of the vertices to Morton code
  float_coordinates_to_int(sommets, int_coords, h);
  int_coordinates_to_Morton(int_coords);

  //Sort Morton_code (and indices to keep trace of the original indices)
  std::sort(indices.begin(), indices.end(),
  [&](int a, int b) { return Morton_code[a] < Morton_code[b]; });
  std::sort(Morton_code.begin(), Morton_code.end());

  build_leaves(compo_connexes_sommets);
  compute_neighbours();

}

//Transform floating points coordinates to integer coordinates
void Morton_Linear_Octree_Particles::float_coordinates_to_int(DoubleTab const& sommets, std::vector<std::array<uint32_t,3>>& int_coords, double h)
{
  std::cout<<"float_coordinates_to_int"<<std::endl;
  constexpr uint32_t max_coord = ((1u<<L) - 1);
  double min[3];
  min[0] = std::numeric_limits<double>::max();
  min[1]=min[0];
  min[2] = min[0];
  double max[3];
  max[0] = std::numeric_limits<double>::lowest();
  max[1]=max[0];
  max[2] = max[0];
  for(int i=0; i<sommets.dimension(0); ++i)
    {
      for(int d=0; d<dimension; ++d)
        {
          if(sommets(i,d)<min[d]) min[d] = sommets(i,d);
          if(sommets(i,d)>max[d]) max[d] = sommets(i,d);
        }
    }

  for(int i=0; i<sommets.dimension(0); ++i)
    {
      for(int d=0; d<dimension; ++d)
        {
          int_coords[i][d] = static_cast<uint32_t>(std::clamp((sommets(i,d) - min[d]) / (max[d] - min[d]),0.,1.) * max_coord);
        }
    }

  // use max and min to determine l
  double Max = std::max( (max[0]-min[0]), std::max( (max[1]-min[1]), (max[2]-min[2]) ));
  l = (int)std::floor(std::log(Max/h) / std::log(2));
  if(l<0)l=0;
  std::cout<<"l = "<<l<<std::endl;
}

void Morton_Linear_Octree_Particles::int_coordinates_to_Morton(std::vector<std::array<uint32_t,3>>& int_coords)
{
  std::cout<<"int_coordinates_to_Morton"<<std::endl;
  int i=0;
  for(auto &c : int_coords)
    {
      Morton_code[i] = encodeMorton(c[0], c[1], c[2]);
      ++i;
    }
}

/*
This methods builds the leaves.
A leaf as a Morton code written using l<=L bits.
A vertex, is in a leaf if its l first bits (from the left) are the same (example 1011 is in leaf 10)
The first leaf is created from the first vertex (by taking its l first bits). Then,
if the next vertex as the same first l bits of the leaf, we store it. Otherwise, we browsed all vertices which
are in the leaf. We determine then the next leaf, using the first l bits of the first vertex that is not in the previous leaf.

The l first bits are taken using   >> (3*(L-l))
*/
void Morton_Linear_Octree_Particles::build_leaves(const ArrOfInt& compo_connexes_sommets)
{
  std::cout<<"build_leaves"<<std::endl;
  uint64_t leaf_code = Morton_code[0] >> (3*(L-l)); //Morton code of the first leaf
  leaf current_leaf;
  current_leaf.Morton_code = leaf_code;
  current_leaf.vertex_indices.push_back(indices[0]); //indices[0] is the index of the vertex in "sommets" corresponding to Morton_code[0]
  int current_compo = compo_connexes_sommets(indices[0]);
  current_leaf.compo = current_compo;
  for(size_t i=1; i<Morton_code.size(); ++i)
    {
      uint64_t next_vertex_code = Morton_code[i] >>(3*(L-l));
      if(next_vertex_code==leaf_code) //Morton_code[i] is still in current_leaf
        {
          current_leaf.vertex_indices.push_back(indices[i]);
          if( (compo_connexes_sommets(indices[i]) != current_leaf.compo ) ) { current_leaf.compo = -1;}
        }
      else //Morton_code[i] is not in current_leaf. We create a new leaf from it
        {
          current_leaf.closest_vertex_indices.resize(current_leaf.vertex_indices.size());
          leaves.push_back(current_leaf);

          current_leaf = leaf(); //reset current leaf //XXX might not be the best way to do it

          current_compo = compo_connexes_sommets(indices[i]);
          current_leaf.Morton_code = next_vertex_code;
          current_leaf.vertex_indices.push_back(indices[i]);
          current_leaf.compo = current_compo;
          leaf_code = next_vertex_code;
        }
    }
  // we add the last leaf
  current_leaf.closest_vertex_indices.resize(current_leaf.vertex_indices.size());
  leaves.push_back(current_leaf);
}


void Morton_Linear_Octree_Particles::compute_neighbours()
{
  std::cout<<"compute_neighbours"<<std::endl;

  // int32_t max_coord = (1 << (L-l)) - 1; //for bound checking
  std::map<uint64_t, size_t> morton_to_idx;
  for(size_t i=0; i<leaves.size(); ++i) {morton_to_idx[leaves[i].Morton_code] = i;}
  uint32_t x, y, z;
  for(auto &current_leaf : leaves)
    {
      DecodeMorton(current_leaf.Morton_code, x, y, z); //need to go back to interger coordinates
      //each neighbours has the same coordinates with +/- 1 on the different coordinates
      for (int dx = -1; dx <= 1; ++dx)
        {
          for (int dy = -1; dy <= 1; ++dy)
            {
              for (int dz = -1; dz <= 1; ++dz)
                {
                  if (dx == 0 && dy == 0 && dz == 0) continue; // this corresponds to itself
                  // if ( x + dx < 0 ||  y + dy < 0 ||  z + dz < 0 ||  x + dx > max_coord ||  y + dy > max_coord || z + dz > max_coord) continue; //bounds
                  uint32_t nx = x + dx;
                  uint32_t ny = y + dy;
                  uint32_t nz = z + dz;
                  uint64_t n_code = encodeMorton(nx, ny, nz); //recode Morton
                  auto it = morton_to_idx.find(n_code); //search for the Morton code of the neighbour in the list
                  if (it != morton_to_idx.end()) // n_code might not exist if there is no vertex in this zone
                    {
                      current_leaf.neighbours.push_back(&leaves[it->second]);
                    }
                }
            }
        }
    }
  std::cout<<"finished neighbour"<<std::endl;
}

void Morton_Linear_Octree_Particles::print_indices(const ArrOfInt& compo_connexes_sommets)
{
  std::ofstream f;
  for(auto & current_leaf : leaves )
    {
      for(int j=0; j<(int)current_leaf.vertex_indices.size(); ++j)
        {
          std::string path;
          path = "MLO_P_";
          int index_j = current_leaf.vertex_indices[j];
          int compo_j = compo_connexes_sommets(index_j);
          path += std::to_string(compo_j);
          path += ".txt";
          f.open(path, std::ios::app);
          f<<index_j<<"\n";
          f.close();
        }
    }

}

void Morton_Linear_Octree_Particles::find_closest(const DoubleTab& sommets, const ArrOfInt& compo_connexes_sommets, std::vector<std::vector<collision_parameters>>& param)
{
  DoubleTab dX_loc(dimension);
  DoubleTab ci(dimension);
  DoubleTab cj(dimension);
  for (auto &current_leaf : leaves)
    {
      for (int i = 0; i < (int)current_leaf.vertex_indices.size(); ++i)
        {
          current_leaf.closest_vertex_indices[i] = -1;
          int index_i = current_leaf.vertex_indices[i];
          if(index_i==1429)std::cout<<"got 1429\n";
          int compo_i = compo_connexes_sommets(index_i);
          ci[0] = sommets(index_i, 0);
          ci[1] = sommets(index_i, 1);
          ci[2] = sommets(index_i, 2);
          if (current_leaf.compo == -1) // we have to look inside the leaf
            {
              for (int j = 0; j < (int)current_leaf.vertex_indices.size(); ++j) // loop on the other points of the leaf
                {
                  if (i == j)
                    continue;

                  int index_j = current_leaf.vertex_indices[j];
                  if(index_j==2315)std::cout<<"got 2315\n";
                  cj[0] = sommets(index_j, 0);
                  cj[1] = sommets(index_j, 1);
                  cj[2] = sommets(index_j, 2);

                  for (int d = 0; d < dimension; ++d)
                    {
                      dX_loc(d) = ci(d) - cj(d);
                    }
                  double d = sqrt(local_carre_norme_vect(dX_loc));
                  int compo_j = compo_connexes_sommets(index_j);
                  if(compo_j==compo_i) continue;
                  bool swaped = false;
                  if (compo_j < compo_i)
                    {
                      std::swap(compo_j, compo_i);
                      std::swap(index_i, index_j);
                      swaped = true;
                    }
                  assert(compo_i<compo_j);
                  if (d < (param[compo_i][compo_j].distance ))
                    {
                      std::cout<<"new couple for MLO  "<<index_i<<" "<<index_j<<" in compos "<<compo_i<<" and "<<compo_j<<std::endl;
                      current_leaf.closest_vertex_indices[i] = index_j;
                      param[compo_i][compo_j].i_closest = index_i;
                      param[compo_i][compo_j].j_closest = index_j;
                      param[compo_i][compo_j].i_j_are_close = true;
                      param[compo_i][compo_j].part_part_collision = true;
                      param[compo_i][compo_j].distance = d;
                      for (int k = 0; k < dimension; ++k)
                        {
                          param[compo_i][compo_j].dX(k) = dX_loc(k);
                        }
                    }
                  if(swaped)
                    {
                      std::swap(compo_j, compo_i);
                      std::swap(index_i, index_j);
                    }
                }
            }
          for (auto neighbour_leaf : current_leaf.neighbours)
            {
              // if (neighbour_leaf->compo == current_leaf.compo)
              //   continue;
              for (int j = 0; j < (int)neighbour_leaf->vertex_indices.size(); ++j)
                {
                  int index_j = neighbour_leaf->vertex_indices[j];
                  if(index_j==2315)std::cout<<"got 2315\n";
                  cj[0] = sommets(index_j, 0);
                  cj[1] = sommets(index_j, 1);
                  cj[2] = sommets(index_j, 2);
                  for (int d = 0; d < dimension; ++d)
                    {
                      dX_loc(d) = ci(d) - cj(d);
                    }
                  double d = sqrt(local_carre_norme_vect(dX_loc));
                  int compo_j = compo_connexes_sommets(index_j);
                  if(compo_j==compo_i) continue;
                  bool swaped = false;
                  if (compo_j < compo_i)
                    {
                      std::swap(compo_j, compo_i);
                      std::swap(index_i, index_j);
                      swaped = true;
                    }
                  if(index_i==1429 && index_j==2315) {std::cout<<"got them "<<d<<"    "<<param[compo_i][compo_j].distance<<"\n";}
                  assert(compo_i<compo_j);
                  if (d < (param[compo_i][compo_j].distance ))
                    {
                      std::cout<<"new couple for MLO in neighbour  "<<index_i<<" "<<index_j<<" in compos "<<compo_i<<" and "<<compo_j<<std::endl;
                      current_leaf.closest_vertex_indices[i] = index_j;
                      param[compo_i][compo_j].i_closest = index_i;
                      param[compo_i][compo_j].j_closest = index_j;
                      param[compo_i][compo_j].i_j_are_close = true;
                      param[compo_i][compo_j].part_part_collision = true;
                      param[compo_i][compo_j].distance = d;
                      for (int k = 0; k < dimension; ++k)
                        {
                          param[compo_i][compo_j].dX(k) = dX_loc(k);
                        }
                    }
                  if(swaped)
                    {
                      std::swap(compo_j, compo_i);
                      std::swap(index_i, index_j);
                    }
                }
            }
        }
    }
}

//x is written using 21 bits (all 9 bits forming the 32bits int are 0).
//This function returns the same number but with two 0 between each bit
//ex : x =1101 -> y = 1001000001
// see https://stackoverflow.com/questions/18529057/produce-interleaving-bit-patterns-morton-keys-for-3d-coordinates-for-32-bit
// see https://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/
uint64_t Morton_Linear_Octree_Particles::expand_bits(uint32_t x)
{
  uint64_t y = x & 0x1fffff;
  y = (y | y << 32) & 0x1f00000000ffff; // shift left 32 bits, OR with self, and 00011111000000000000000000000000000000001111111111111111
  y = (y | y << 16) & 0x1f0000ff0000ff; // shift left 32 bits, OR with self, and 00011111000000000000000011111111000000000000000011111111
  y = (y | y << 8) & 0x100f00f00f00f00f; // shift left 32 bits, OR with self, and 0001000000001111000000001111000000001111000000001111000000000000
  y = (y | y << 4) & 0x10c30c30c30c30c3; // shift left 32 bits, OR with self, and 0001000011000011000011000011000011000011000011000011000100000000
  y = (y | y << 2) & 0x1249249249249249;
  return y;
}

//inverse of expand_bits.
uint32_t Morton_Linear_Octree_Particles::compact_bits(uint64_t x)
{
  x &= 0x1249249249249249;
  x = (x ^ (x >> 2)) & 0x10c30c30c30c30c3;
  x = (x ^ (x >> 4)) & 0x100f00f00f00f00f;
  x = (x ^ (x >> 8)) & 0x1f0000ff0000ff;
  x = (x ^ (x >> 16)) & 0x1f00000000ffff;
  x = (x ^ (x >> 32)) & 0x1fffff;
  return static_cast<uint32_t>(x);
}

//from integer coordinates to Morton code
uint64_t Morton_Linear_Octree_Particles::encodeMorton(uint32_t x, uint32_t y, uint32_t z)
{
  uint64_t a=0;
  a|= expand_bits(x) | expand_bits(y)<<1 | expand_bits(z)<<2;
  return a;
}

//from a Morton code to integer coordinates
void Morton_Linear_Octree_Particles::DecodeMorton(uint64_t code, uint32_t& x, uint32_t& y, uint32_t& z)
{
  x = compact_bits(code >> 0);
  y = compact_bits(code >> 1);
  z = compact_bits(code >> 2);
}