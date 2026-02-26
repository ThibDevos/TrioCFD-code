#include <Z_Curve_Link_Cell.h>
#include <bitset>

// Sort 'codes' and synchronize 'indices'
void radix_sort_indices(std::vector<uint64_t>& codes, std::vector<int>& indices)
{
  size_t n = codes.size();
  if (n == 0) return;

  std::vector<uint64_t> temp_codes(n);
  std::vector<int> temp_indices(n);

  for (int shift = 0; shift < 64; shift += 8)
    {
      size_t count[257] = {0};

      for (size_t i = 0; i < n; i++)
        {
          size_t bucket = (codes[i] >> shift) & 0xFF;
          count[bucket + 1]++;
        }

      for (size_t i = 0; i < 256; i++)
        {
          count[i + 1] += count[i];
        }

      for (size_t i = 0; i < n; i++)
        {
          size_t bucket = (codes[i] >> shift) & 0xFF;
          size_t pos = count[bucket]++;
          temp_codes[pos] = codes[i];
          temp_indices[pos] = indices[i];
        }

      codes.swap(temp_codes);
      indices.swap(temp_indices);
    }
}

Z_Curve_Link_Cell::Z_Curve_Link_Cell() {}
Z_Curve_Link_Cell::Z_Curve_Link_Cell(const DoubleTab& sommets, double h,  const ArrOfInt& compo_connexes_sommets)
{
  int nb_sommets=sommets.dimension(0);
  indices.resize(nb_sommets);
  Morton_code.resize(nb_sommets);
  if(nb_sommets==0) return;
  std::iota(indices.begin(), indices.end(), 0);
  std::vector<std::array<uint32_t,3>> int_coords(nb_sommets);

  //transform the floating points coordiantes of the vertices to Morton code
  float_coordinates_to_int(sommets, int_coords, h);
  int_coordinates_to_Morton(int_coords);


  radix_sort_indices(Morton_code, indices);

  build_leaves(compo_connexes_sommets);
  compute_neighbours();
}

void Z_Curve_Link_Cell::build(const DoubleTab& sommets, double h,  const ArrOfInt& compo_connexes_sommets)
{
  int nb_sommets=sommets.dimension(0);
  indices.resize(nb_sommets);
  Morton_code.resize(nb_sommets);
  if(nb_sommets==0) return;
  std::iota(indices.begin(), indices.end(), 0);
  std::vector<std::array<uint32_t,3>> int_coords(nb_sommets);

  //transform the floating points coordiantes of the vertices to Morton code
  float_coordinates_to_int(sommets, int_coords, h);
  int_coordinates_to_Morton(int_coords);

  radix_sort_indices(Morton_code, indices);

  build_leaves(compo_connexes_sommets);
  compute_neighbours();
}

//Transform floating points coordinates to integer coordinates
void Z_Curve_Link_Cell::float_coordinates_to_int(DoubleTab const& sommets, std::vector<std::array<uint32_t,3>>& int_coords, double h)
{
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

  //Max length
  double extent[3];
  extent[0] = max[0] - min[0];
  extent[1] = max[1] - min[1];
  extent[2] = max[2] - min[2];

  double Max = std::max(extent[0], std::max(extent[1], extent[2]));
  if(Max == 0.0) Max = 1.0;

  for(int i=0; i<sommets.dimension(0); ++i)
    {
      for(int d=0; d<dimension; ++d)
        {
          int_coords[i][d] = static_cast<uint32_t>(std::clamp( (sommets(i,d) - min[d]) / (Max), 0.,1.) * max_coord);
        }
    }
  uint32_t max_global_index = (1u << l) - 1;
  // On projette l'extension physique sur la grille d'indices
  // Note: extent/Max est <= 1.0.
  grid_limits[0] = static_cast<uint32_t>( (extent[0]/Max) * max_global_index );
  grid_limits[1] = static_cast<uint32_t>( (extent[1]/Max) * max_global_index );
  grid_limits[2] = static_cast<uint32_t>( (extent[2]/Max) * max_global_index );
  // use max and min to determine l
  l = (int)std::floor(std::log(Max/(h)) / std::log(2));
  std::cout<<"l= "<<l<<std::endl;
  if(l<0)l=0;
}

void Z_Curve_Link_Cell::int_coordinates_to_Morton(std::vector<std::array<uint32_t,3>>& int_coords)
{
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
void Z_Curve_Link_Cell::build_leaves(const ArrOfInt& compo_connexes_sommets)
{
  grid_limits[0] = 0;
  grid_limits[1] = 0;
  grid_limits[2] = 0;
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
          uint32_t x, y, z;
          DecodeMorton(current_leaf.Morton_code, x, y, z);
          if(x>grid_limits[0]) grid_limits[0] = x;
          if(y>grid_limits[1]) grid_limits[1] = y;
          if(z>grid_limits[2]) grid_limits[2] = z;
          current_leaf.closest_vertex_indices.resize(current_leaf.vertex_indices.size());
          leaves.push_back(current_leaf);

          current_leaf = leaf(); //reset current leaf

          current_compo = compo_connexes_sommets(indices[i]);
          current_leaf.Morton_code = next_vertex_code;
          current_leaf.vertex_indices.push_back(indices[i]);
          current_leaf.compo = current_compo;
          leaf_code = next_vertex_code;
        }
    }

  // we add the last leaf
  uint32_t x, y, z;
  DecodeMorton(current_leaf.Morton_code, x, y, z);
  if(x>grid_limits[0]) grid_limits[0] = x;
  if(y>grid_limits[1]) grid_limits[1] = y;
  if(z>grid_limits[2]) grid_limits[2] = z;
  current_leaf.closest_vertex_indices.resize(current_leaf.vertex_indices.size());
  leaves.push_back(current_leaf);
}


void Z_Curve_Link_Cell::compute_neighbours()
{
  std::ofstream f;
  f.open("unordered_map.txt", std::ios::app);
  std::chrono::steady_clock::time_point begin;
  std::chrono::steady_clock::time_point end;
  auto time =  std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
  std::map<uint64_t, size_t> morton_to_idx;
  for(size_t i=0; i<leaves.size(); ++i) {morton_to_idx[leaves[i].Morton_code] = i;}
  uint32_t x, y, z;
  int nb_n = 0;
  for (auto it = leaves.begin(); it != leaves.end(); ++it)
    {
      auto& current_leaf = *it;
      DecodeMorton(current_leaf.Morton_code, x, y, z); //need to go back to interger coordinates
      //set if the leaf is near a boundary or not
      int boundary = 0;
      boundary |=  (x==0);
      boundary <<= 1;
      boundary |=  (y==0);
      boundary <<= 1;
      boundary |=  (z==0);
      boundary <<= 1;
      boundary |=  (x==grid_limits[0]);
      boundary <<= 1;
      boundary |=  (y==grid_limits[1]);
      boundary <<= 1;
      boundary |=  (z==grid_limits[2]);

      current_leaf.boundary = boundary;

      //each neighbours has the same coordinates with +/- 1 on the different coordinates
      for (int dx = -1; dx <= 1; ++dx)
        {
          for (int dy = -1; dy <= 1; ++dy)
            {
              for (int dz = -1; dz <= 1; ++dz)
                {
                  if (dx == 0 && dy == 0 && dz == 0) continue; // this corresponds to itself
                  if ( (int)(x + dx) < 0 ||  (int)(y + dy) < 0 ||  (int)(z + dz) < 0 ||  x + dx > grid_limits[0] ||  y + dy > grid_limits[1] || z + dz > grid_limits[2]) continue; //bounds
                  uint32_t nx = x + dx;
                  uint32_t ny = y + dy;
                  uint32_t nz = z + dz;
                  uint64_t n_code = encodeMorton(nx, ny, nz); //recode Morton
                  if(n_code<=current_leaf.Morton_code) continue;
                  leaf search_leaf;
                  search_leaf.Morton_code = n_code;
                  begin = std::chrono::steady_clock::now();
                  auto found_it = morton_to_idx.find(n_code);

		  //auto found_it = std::lower_bound(leaves.begin(), leaves.end(), search_leaf, [](const leaf& a, const leaf&b) {return a.Morton_code<b.Morton_code;}); //we use the fact that the leaves are ordered
                  end = std::chrono::steady_clock::now();
                  time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
                  //if(found_it!=leaves.end() && found_it->Morton_code == n_code)
                  //  {
                  //    current_leaf.neighbours.push_back(&leaves[found_it - leaves.begin()]);
                  //  }
		  if (found_it != morton_to_idx.end()) // n_code might not exist if there is no vertex in this zone
                  {
                    current_leaf.neighbours.push_back(&leaves[found_it->second]);
                  }
                }
            }
        }
      nb_n += (int)current_leaf.neighbours.size();
    }
  f<<time<<" "<<time/leaves.size()<<"\n";
  f.close();
}

void Z_Curve_Link_Cell::print_indices(const ArrOfInt& compo_connexes_sommets)
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

void Z_Curve_Link_Cell::which_leaf(int s)
{
  bool got_it = false;
  uint64_t leaf_code;
  for (auto &current_leaf : leaves)
    {
      for (int i = 0; i < (int)current_leaf.vertex_indices.size(); ++i)
        {
          if(current_leaf.vertex_indices[i] == s)
            {
              got_it =true;
              leaf_code = current_leaf.Morton_code;
              std::cout<<"neighbours : ";
              for(auto & n : current_leaf.neighbours)
                {
                  std::cout<<n->Morton_code<<" ";
                }
              std::cout<<std::endl;
            }
        }
    }
  if(got_it) {std::cout<<"Leaf for "<<s<<" : "<<leaf_code<<std::endl;}
  else {std::cout<<"not found"<<std::endl;}
}

void Z_Curve_Link_Cell::find_closest(const DoubleTab& sommets, const ArrOfInt& compo_connexes_sommets, std::vector<std::vector<collision_parameters>>& param,
                                     DoubleVect const& origin, DoubleVect const& domain_dimensions)
{
  DoubleTab dX_loc(dimension);
  for(auto & current_leaf : leaves)
    {
      std::vector<int> vertices;
      std::vector<double> coordinates;
      //we first load all vertices indices and brut force on the list
      //If current_leaf.compo != -1, then we don't need to load vertices in current_leaf, and only need to load vertices of the neighbour leaves that have neighbour->compo != current_leaf.compo

      for(int i=0; i<(int)current_leaf.vertex_indices.size(); ++i)
        {
          vertices.push_back(current_leaf.vertex_indices[i]);
          for(int d=0; d<dimension; ++d) {coordinates.push_back(sommets(current_leaf.vertex_indices[i],d));}
        }
      for(auto neighbour_leaf : current_leaf.neighbours)
        {
          //if(neighbour_leaf->compo == current_leaf.compo && current_leaf.compo !=-1) continue;
          for(int i=0; i<(int)neighbour_leaf->vertex_indices.size(); ++i)
            {
              vertices.push_back(neighbour_leaf->vertex_indices[i]);
              for(int d=0; d<dimension; ++d) {coordinates.push_back(sommets(neighbour_leaf->vertex_indices[i],d));}
            }
        }
      for(int i=0; i<(int)current_leaf.vertex_indices.size(); ++i)
        {
          int index_i = vertices[i];
          int compo_i = compo_connexes_sommets(index_i);
          int boundary_mask = current_leaf.boundary;
          int nb_particles_tot = (int)param[compo_i].size();
          //check boundaries
          // /!\ As for the rest of the Collision Models at this time, it is for rectangle domain with faces parallel to the canonical basis
	  if(boundary_mask!=0)
            {
		    int ind_wall = 0;
              if(boundary_mask & 32)
                {
                  double d = std::fabs(coordinates[3*i + 0] - origin[0]);
                  ind_wall = nb_particles_tot - 6;

                  if(d<param[compo_i][ind_wall].distance)
                    {
                      param[compo_i][ind_wall].i_closest = index_i;
                      param[compo_i][ind_wall].i_j_are_close = true;
                      param[compo_i][ind_wall].distance = d;
                      param[compo_i][ind_wall].dX(0)= d;
                    }
                }
              if(boundary_mask & 16)
                {
                  double d = std::fabs(coordinates[3*i + 1]- origin[1]);
                  ind_wall = nb_particles_tot -5 ;
                  if(d<param[compo_i][ind_wall].distance)
                    {
                      param[compo_i][ind_wall].i_closest = index_i;
                      param[compo_i][ind_wall].i_j_are_close = true;
                      param[compo_i][ind_wall].distance = d;
                      param[compo_i][ind_wall].dX(1)= d;
                    }
                }
              if(boundary_mask & 8)
                {
                  double d = std::fabs(coordinates[3*i + 2] - origin[2]);
                  ind_wall = nb_particles_tot - 4;
                  if(d<param[compo_i][ind_wall].distance)
                    {
                      param[compo_i][ind_wall].i_closest = index_i;
                      param[compo_i][ind_wall].i_j_are_close = true;
                      param[compo_i][ind_wall].distance = d;
                      param[compo_i][ind_wall].dX(2)= d;
                    }
                }
              if(boundary_mask & 4)
                {
                  // std::cout<<"dom dim : "<<-origin[0] - domain_dimensions[0]<<std::endl;
                  double d = std::fabs(coordinates[3*i + 0] - origin[0] - domain_dimensions[0]);
                  // std::cout<<d<<std::endl;
                  ind_wall = nb_particles_tot - 3;
                  if(d<param[compo_i][ind_wall].distance)
                    {
                      param[compo_i][ind_wall].i_closest = index_i;
                      param[compo_i][ind_wall].i_j_are_close = true;
                      param[compo_i][ind_wall].distance = d;
                      param[compo_i][ind_wall].dX(0)= d;
                    }
                }
              if(boundary_mask & 2)
                {
                  double d = std::fabs(coordinates[3*i + 1] - origin[1] - domain_dimensions[1]);
                  ind_wall = nb_particles_tot - 2;
                  if(d<param[compo_i][ind_wall].distance)
                    {
                      param[compo_i][ind_wall].i_closest = index_i;
                      param[compo_i][ind_wall].i_j_are_close = true;
                      param[compo_i][ind_wall].distance = d;
                      param[compo_i][ind_wall].dX(1)= d;
                    }
                }
              if(boundary_mask & 1)
                {
                  double d = std::fabs(coordinates[3*i + 2] - origin[2] - domain_dimensions[2]);
                  ind_wall = nb_particles_tot - 1;
                  if(d<param[compo_i][ind_wall].distance)
                    {
                      param[compo_i][ind_wall].i_closest = index_i;
                      param[compo_i][ind_wall].i_j_are_close = true;
                      param[compo_i][ind_wall].distance = d;
                      param[compo_i][ind_wall].dX(2)= d;
                    }
                }
            }

          int j = current_leaf.compo == -1 ? i+1 : (int)current_leaf.vertex_indices.size() +1;
          for(j=i+1; j<(int)vertices.size(); ++j)
            {

              int index_j = vertices[j];
              int compo_j = compo_connexes_sommets(index_j);
              if(compo_i==compo_j) continue;


              for (int d = 0; d < dimension; ++d)
                {
                  dX_loc(d) = coordinates[3*i + d] - coordinates[3*j + d];
                }
              double d = (local_carre_norme_vect(dX_loc));
              bool swaped = false;
              if (compo_j < compo_i)
                {
                  std::swap(compo_j, compo_i);
                  std::swap(index_i, index_j);
                  swaped = true;
                }
              assert(compo_i<compo_j);
              if (d < (param[compo_i][compo_j-compo_i-1].distance ))
                {
                  current_leaf.closest_vertex_indices[i] = index_j;
                  param[compo_i][compo_j-compo_i-1].i_closest = index_i;
                  param[compo_i][compo_j-compo_i-1].j_closest = index_j;
                  param[compo_i][compo_j-compo_i-1].i_j_are_close = true;
                  param[compo_i][compo_j-compo_i-1].part_part_collision = true;
                  param[compo_i][compo_j-compo_i-1].distance = d;
                  for (int k = 0; k < dimension; ++k)
                    {
                      param[compo_i][compo_j-compo_i-1].dX(k) = dX_loc(k) * (swaped ? -1 : 1);
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



//x is written using 21 bits (all 9 bits forming the 32bits int are 0).
//This function returns the same number but with two 0 between each bit
//ex : x =1101 -> y = 1001000001
// see https://stackoverflow.com/questions/18529057/produce-interleaving-bit-patterns-morton-keys-for-3d-coordinates-for-32-bit
// see https://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/
uint64_t Z_Curve_Link_Cell::expand_bits(uint32_t x)
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
uint32_t Z_Curve_Link_Cell::compact_bits(uint64_t x)
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
uint64_t Z_Curve_Link_Cell::encodeMorton(uint32_t x, uint32_t y, uint32_t z)
{
  uint64_t a=0;
  a|= expand_bits(x) | expand_bits(y)<<1 | expand_bits(z)<<2;
  return a;
}

//from a Morton code to integer coordinates
void Z_Curve_Link_Cell::DecodeMorton(uint64_t code, uint32_t& x, uint32_t& y, uint32_t& z)
{
  x = compact_bits(code >> 0);
  y = compact_bits(code >> 1);
  z = compact_bits(code >> 2);
}
