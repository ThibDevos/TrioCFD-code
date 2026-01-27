/****************************************************************************
* Copyright (c) 2025, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include <Collision_Model_FT_ellipsoid.h>
#include <Probleme_FT_Disc_gen.h>
#include <Solid_Particle_ellipsoid.h>


Implemente_instanciable_sans_constructeur(Collision_Model_FT_ellipsoid,"Collision_Model_FT_ellipsoid",Collision_Model_FT_base);

Collision_Model_FT_ellipsoid::Collision_Model_FT_ellipsoid()
{
}
int Collision_Model_FT_ellipsoid::lire_motcle_non_standard(const Motcle& word, Entree& is)
{
  if(word=="collision_detection")
    {
      Motcles words;
      words.add("deepest_points");
      words.add("closest_point");
      Motcle secondword;
      is >> secondword;
      const int r = words.search(secondword);
      switch(r)
        {
        case 0 :
          collision_detection_ = Collision_detection::DEEPEST;
          break;
        case 1 :
          collision_detection_ = Collision_detection::CLOSEST;
          break;
        default:
          Cerr << "Error " << words << "was expected whereas " << secondword <<
               " has been found."<< finl;
          Process::exit();
        }
    }
  else if(word=="collision_normal")
    {
      Motcles words;
      words.add("normal_i");
      words.add("normal_average_ij");
      Motcle secondword;
      is >> secondword;
      const int r = words.search(secondword);
      switch(r)
        {
        case 0 :
          collision_normal_ = Collision_normal::NORMAL_I;
          break;
        case 1 :
          collision_normal_ = Collision_normal::NORMAL_IJ;
          break;
        default:
          Cerr << "Error " << words << "was expected whereas " << secondword <<
               " has been found."<< finl;
          Process::exit();
        }
    }
  else if(word=="detection_option")
    {
      Motcles words;
      words.add("naive");
      words.add("octree");
      words.add("zlc");
      Motcle secondword;
      is >> secondword;
      const int r = words.search(secondword);
      switch(r)
        {
        case 0 :
          detection_option = Detection_Option::NAIVE;
          break;
        case 1 :
          detection_option = Detection_Option::OCTREE;
          break;
        case 2 :
          detection_option = Detection_Option::ZLC;
          break;
        default:
          Cerr << "Error " << words << "was expected whereas " << secondword <<
               " has been found."<< finl;
          Process::exit();
        }
    }
  else if(word == "fichier_debug")// XXX debug file
    {
      is >> fichier_debug;
      std::transform(fichier_debug.getString().begin(), fichier_debug.getString().end(), fichier_debug.getString().begin(), ::tolower);
    }
  else
    {
      return Collision_Model_FT_base::lire_motcle_non_standard(word, is);
    }
  return 1;

}
Entree& Collision_Model_FT_ellipsoid::readOn (Entree& is)
{
  Collision_Model_FT_base::readOn(is);
  return is;
}

Sortie& Collision_Model_FT_ellipsoid::printOn(Sortie& os) const
{
  Cerr << "Error::printOn is not implemented." << finl;
  Process::exit();
  return os;
}

void Collision_Model_FT_ellipsoid::deepest_points(IntLists const& compo_sommets, IntLists const& sommets_facets, Maillage_FT_Disc const& mesh, collision_parameters& param,
                                                  DoubleTab const& positions, bool check_cg, DoubleTab& dX)
{
  const DoubleTab& sommets = mesh.sommets();
  DoubleTab dX_min(dimension);
  double dist_cg = 0.;
  int i_som_closest = -1;
  int j_som_closest = -1;

  for (int d = 0; d < dimension; ++d)
    dX_min(d) = 1000.;
  double dX_min_norm_i = sqrt(local_carre_norme_vect(dX_min));

  for (int i_som = 0; i_som < compo_sommets[param.particle_i].size(); ++i_som)
    {
      int i_global = compo_sommets[param.particle_i][i_som];
      for (int d = 0; d < dimension; ++d)
        {
          dX(d) = sommets(i_global, d) - positions(param.particle_j, d);
        }
      dist_cg = sqrt(local_carre_norme_vect(dX));
      if (dist_cg < dX_min_norm_i)
        {
          for (int d = 0; d < dimension; ++d)
            {
              dX_min(d) = dX(d);
            }
          dX_min_norm_i = dist_cg;
          i_som_closest = i_global;
        }
    }
  assert(i_som_closest >= 0);

  for (int d = 0; d < dimension; ++d)
    dX_min(d) = 1000.;
  double dX_min_norm_j = sqrt(local_carre_norme_vect(dX_min));

  for (int j_som = 0; j_som < compo_sommets[param.particle_j].size(); ++j_som)
    {
      int j_global = compo_sommets[param.particle_j][j_som];
      for (int d = 0; d < dimension; ++d)
        {
          dX(d) = sommets(j_global, d) - positions(param.particle_i, d);
        }
      dist_cg = sqrt(local_carre_norme_vect(dX));
      if (dist_cg < dX_min_norm_j)
        {
          for (int d = 0; d < dimension; ++d)
            {
              dX_min(d) = dX(d);
            }
          dX_min_norm_j = dist_cg;
          j_som_closest = j_global;
        }
    }
  assert(j_som_closest >= 0);

  int fi_closest = -1;
  int fj_closest = -1;
  DoubleTab const& cg_fa7 = mesh.get_gravity_center_fa7();
  if (check_cg)
    {
      // check if a center of gravity of the incident facets of i_som_closest is deeper in the particle than the i_som_closest
      for (int fi = 0; fi < sommets_facets[i_som_closest].size(); ++fi)
        {
          int fa7 = sommets_facets[i_som_closest][fi];
          for (int d = 0; d < dimension; ++d)
            {
              dX(d) = cg_fa7(fa7, d) - positions(param.particle_j, d);
            }
          dist_cg = sqrt(local_carre_norme_vect(dX));
          if (dist_cg < dX_min_norm_i)
            {
              for (int d = 0; d < dimension; ++d)
                {
                  dX_min(d) = dX(d);
                }
              dX_min_norm_i = dist_cg;
              fi_closest = fa7;
            }
        }
      for (int fj = 0; fj < sommets_facets[j_som_closest].size(); ++fj)
        {
          int fa7 = sommets_facets[j_som_closest][fj];
          for (int d = 0; d < dimension; ++d)
            {
              dX(d) = cg_fa7(fa7, d) - positions(param.particle_i, d);
            }
          dist_cg = sqrt(local_carre_norme_vect(dX));
          if (dist_cg < dX_min_norm_j)
            {
              for (int d = 0; d < dimension; ++d)
                {
                  dX_min(d) = dX(d);
                }
              dX_min_norm_j = dist_cg;
              fj_closest = fa7;
            }
        }
    }

  DoubleTab deepest_i(dimension);
  DoubleTab deepest_j(dimension);

  if (fi_closest < 0)
    {
      param.i_facet = false;
      param.i_closest = i_som_closest;
      for (int d = 0; d < dimension; ++d)
        {
          deepest_i(d) = sommets(i_som_closest, d);
        }
    }
  else
    {
      param.i_facet = true;
      param.i_closest = fi_closest;
      for (int d = 0; d < dimension; ++d)
        {
          deepest_i(d) = cg_fa7(fi_closest, d);
        }
    }
  if (fj_closest < 0)
    {
      param.j_facet = false;
      param.j_closest = j_som_closest;
      for (int d = 0; d < dimension; ++d)
        {
          deepest_j(d) = sommets(j_som_closest, d);
        }
    }
  else
    {
      param.j_facet = true;
      param.j_closest = fj_closest;
      for (int d = 0; d < dimension; ++d)
        {
          deepest_j(d) = cg_fa7(fj_closest, d);
        }
    }

  for (int d = 0; d < dimension; ++d)
    {
      dX(d) = deepest_j(d) - deepest_i(d);
    }
}

void Collision_Model_FT_ellipsoid::closest_nodes(IntLists const& compo_sommets, Maillage_FT_Disc const& mesh, collision_parameters& param,
                                                 DoubleTab& dX)
{
  auto sommets = mesh.sommets();
  DoubleTab dX_min(dimension);
  double dist_cg = 0.;
  param.i_closest = -1;
  param.j_closest = -1;
  for (int d = 0; d < dimension; ++d)
    dX_min(d) = 1000.;


  double dX_min_norm = sqrt(local_carre_norme_vect(dX_min));
  for (int i = 0; i < compo_sommets[param.particle_i].size(); ++i)
    {
      int i_global = compo_sommets[param.particle_i][i];
      for (int j = 0; j < compo_sommets[param.particle_j].size(); ++j)
        {
          int j_global = compo_sommets[param.particle_j][j];
          for (int d = 0; d < dimension; ++d)
            {
              dX(d) = sommets(i_global, d) - sommets(j_global, d);
            }
          dist_cg = sqrt(local_carre_norme_vect(dX));
          if (dist_cg < dX_min_norm)
            {
              for (int d = 0; d < dimension; ++d)
                {
                  dX_min(d) = dX(d);
                }
              dX_min_norm = dist_cg;
              param.i_closest = i_global;
              param.j_closest = j_global;
              param.i_j_are_close  = true;
            }
        }
    }

  for (int d = 0; d < dimension; ++d)
    {
      dX(d) = dX_min(d);
    }
}
void Collision_Model_FT_ellipsoid::closest_nodes(IntLists const& compo_sommets, const ArrOfInt& compo_connexes_fa7, const IntLists& compo_connexe_facets, Maillage_FT_Disc const& mesh, collision_parameters& param,
                                                 DoubleTab& dX,  const Octree_Double& octree)
{
  DoubleTab dX_min(dimension);
  auto sommets = mesh.sommets();
  auto facets = mesh.facettes();
  ArrOfInt liste_facettes;
  double distmax = mesh.get_global_mesh_size();
  double dist_cg=0;
  double dX_min_norm=std::numeric_limits<double>::max();
  DoubleTab coord(dimension);
  for (int i = 0; i < compo_sommets[param.particle_i].size(); ++i) //loop on all the vertices of the compo i
    {
      int i_global = compo_sommets[param.particle_i][i];
      for(int d=0; d<dimension; ++d)
        coord(d) = sommets(i_global, d);
      octree.search_elements_box(coord[0] - distmax, coord[1] - distmax, coord[2] - distmax, coord[0] + distmax,
                                 coord[1] + distmax, coord[2] + distmax, liste_facettes); // get the list of facets that are in a box centered at coord and of size distmax*2


      const int nliste = liste_facettes.size_array();
      for(int j=0; j<nliste; ++j) //loop on the vertices
        {
          int num_facette = 0;
          if(detection_option != Detection_Option::NAIVE)
            {
              num_facette = liste_facettes[j];
              if(compo_connexes_fa7(num_facette) != param.particle_j) continue; //if facet does not belong to compo j, we don't need to check
            }
          else
            {
              num_facette = compo_connexe_facets[param.particle_j][liste_facettes[j]];
            }
          for(int v=0; v<3; ++v) //loop on the vertices of the facet
            {
              for(int d=0; d<dimension; ++d)
                {
                  dX(d) = coord[d] - sommets(facets(num_facette,v),d);
                }
              dist_cg = sqrt(local_carre_norme_vect(dX));
              if (dist_cg < dX_min_norm)
                {
                  std::cout<<"new couple "<<i_global<<" "<<facets(num_facette,v)<<std::endl;
                  for (int d = 0; d < dimension; ++d)
                    {
                      dX_min(d) = dX(d);
                    }
                  dX_min_norm = dist_cg;
                  param.i_closest = i_global;
                  param.j_closest = facets(num_facette,v);
                  param.i_j_are_close  = true;
                }
            }
        }
    }
  for (int d = 0; d < dimension; ++d)
    {
      dX(d) = dX_min(d);
    }
}
void Collision_Model_FT_ellipsoid::normal_i(Maillage_FT_Disc const& mesh, IntLists const& sommets_facets, collision_parameters& param, DoubleTab& n)
{
  auto nfacets = mesh.get_normale_facettes();
  if(param.i_facet)
    {
      for(int d = 0; d<dimension; ++d) {n(d) = nfacets(param.i_closest,d);}
    }
  else //the deepest point is the node -> the normal at the node is the average of the normal of the incident facets
    {
      for(int fi = 0; fi<sommets_facets[param.i_closest].size(); ++fi)
        {
          int fa7 = sommets_facets[param.i_closest][fi];
          for(int d = 0; d<dimension; ++d) {n(d) += nfacets(fa7,d);}
        }
    }
  double nn = sqrt(local_carre_norme_vect(n));
  for(int d = 0; d<dimension; ++d) {n(d)/=nn;}
}

void Collision_Model_FT_ellipsoid::normal_average_ij(Maillage_FT_Disc const& mesh,IntLists const& sommets_facets, collision_parameters& param, DoubleTab& n)
{
  auto nfacets = mesh.get_normale_facettes();
  DoubleTab ni(dimension);
  DoubleTab nj(dimension);

  if(param.i_facet)
    {
      for(int d = 0; d<dimension; ++d) {ni(d) = nfacets(param.i_closest,d);}
    }
  else //the deepest point is the node -> the normal at the node is the average of the normal of the incident facets
    {
      for(int fi = 0; fi<sommets_facets[param.i_closest].size(); ++fi)
        {
          int fa7 = sommets_facets[param.i_closest][fi];
          for(int d = 0; d<dimension; ++d) {ni(d) += nfacets(fa7,d);}
        }
    }
  if(param.j_facet)
    {
      for(int d = 0; d<dimension; ++d) {nj(d) = nfacets(param.j_closest,d);}
    }
  else //the deepest point is the node -> the normal at the node is the average of the normal of the incident facets
    {
      for(int fi = 0; fi<sommets_facets[param.j_closest].size(); ++fi)
        {
          int fa7 = sommets_facets[param.j_closest][fi];
          for(int d = 0; d<dimension; ++d) {nj(d) += nfacets(fa7,d);}
        }
    }

  double nni = sqrt(local_carre_norme_vect(ni));
  double nnj = sqrt(local_carre_norme_vect(nj));
  for(int d = 0; d<dimension; ++d) {ni(d)/=nni; nj(d)/=nnj;}
  for (int d = 0; d < dimension; ++d)
    {
      n(d) = (ni(d) - nj(d));
    }
  double nn = sqrt(local_carre_norme_vect(n));
  for (int d = 0; d < dimension; ++d)
    {
      n(d) /= nn;
    }
}

void Collision_Model_FT_ellipsoid::collision_point(Maillage_FT_Disc const& mesh, int i_closest, int j_closest, bool i_fa7, bool j_fa7, DoubleTab& xc)
{
  IntTab facettes= mesh.facettes();
  DoubleTab sommets = mesh.sommets();
  DoubleTab const& cg_fa7 = mesh.get_gravity_center_fa7();
  if(i_fa7)
    {
      for(int d = 0; d<dimension; ++d) {xc(d) = cg_fa7(facettes(i_closest),d)/2.;}
    }
  else
    {
      for(int d = 0; d<dimension; ++d) {xc(d) = sommets(i_closest,d)/2.;}
    }
  if(j_fa7)
    {
      for(int d = 0; d<dimension; ++d) {xc(d) += cg_fa7(facettes(j_closest),d)/2.;}
    }
  else
    {
      for(int d = 0; d<dimension; ++d) {xc(d) += sommets(j_closest,d)/2.;}
    }
}

//used
void Collision_Model_FT_ellipsoid::compute_dX(DoubleTab& dX, collision_parameters& param, const DoubleTab& particles_position,
                                              const bool is_particle_particle_collision, const IntLists& compo_sommets,
                                              const IntLists& sommets_facets, const Maillage_FT_Disc& mesh)
{
  const DoubleTab& sommets = mesh.sommets();
  if(is_particle_particle_collision)
    {
      bool check_cg = true;
      if(collision_detection_==Collision_detection::CLOSEST)
        {
          closest_nodes(compo_sommets, mesh, param ,dX); //check_cg not implemented
        }

      if(collision_detection_==Collision_detection::DEEPEST)
        {
          deepest_points(compo_sommets,sommets_facets,mesh, param, particles_position, check_cg, dX);
        }
    }
  else
    {
      if(sommets.dimension(0)==0) {return;} //for Procs that don't have any vertex
      int ind_wall = param.particle_j - nb_particles_tot_;
      int ori = ind_wall < dimension ? ind_wall : ind_wall - dimension;
      double dX_min =  std::numeric_limits<double>::max();
      for(int i_som = 0; i_som < compo_sommets[param.particle_i].size(); ++i_som) //recherche le sommet le plus proche du mur
        {
          dX(ori) = std::fabs(sommets(compo_sommets[param.particle_i][i_som],ori) -  (origin_(ori) + (ind_wall>2)*domain_dimensions_(ori)));
          if(dX(ori) < dX_min )
            {
              param.i_closest = compo_sommets[param.particle_i][i_som];
              param.i_j_are_close  = true;
            }
          dX_min = dX(ori) <= dX_min ? dX(ori) : dX_min;
        }
      dX(ori) = dX_min;
    }
}
void Collision_Model_FT_ellipsoid::compute_dX_octree(DoubleTab& dX, collision_parameters& param, const DoubleTab& particles_position,
                                                     const bool is_particle_particle_collision, const IntLists& compo_sommets,
                                                     const IntLists& sommets_facets, const ArrOfInt& compo_connexes_fa7, const IntLists& compo_connexe_facets,  const Maillage_FT_Disc& mesh, const Octree_Double& octree)
{
  const DoubleTab& sommets = mesh.sommets();
  if(is_particle_particle_collision)
    {
      bool check_cg = true;
      if(collision_detection_==Collision_detection::CLOSEST)
        {
          closest_nodes(compo_sommets,compo_connexes_fa7, compo_connexe_facets,mesh, param ,dX,  octree); //check_cg not implemented
        }

      if(collision_detection_==Collision_detection::DEEPEST)
        {
          deepest_points(compo_sommets,sommets_facets,mesh, param, particles_position, check_cg, dX);
        }
    }
  else
    {
      if(sommets.dimension(0)==0) {return;} //for Procs that don't have any vertex
      int ind_wall = param.particle_j - nb_particles_tot_;
      int ori = ind_wall < dimension ? ind_wall : ind_wall - dimension;
      double dX_min =  std::numeric_limits<double>::max();
      for(int i_som = 0; i_som < compo_sommets[param.particle_i].size(); ++i_som) //recherche le sommet le plus proche du mur
        {
          dX(ori) = std::fabs(sommets(compo_sommets[param.particle_i][i_som],ori) -  (origin_(ori) + (ind_wall>2)*domain_dimensions_(ori)));
          if(dX(ori) < dX_min ) param.i_closest = i_som;
          dX_min = dX(ori) <= dX_min ? dX(ori) : dX_min;
        }
      dX(ori) = dX_min;
    }
}

void Collision_Model_FT_ellipsoid::compute_norm(DoubleTab& norm, collision_parameters& param, const IntLists& sommets_facets, const Maillage_FT_Disc& mesh)
{

  if(collision_normal_==Collision_normal::NORMAL_I)
    {
      normal_i(mesh, sommets_facets, param, norm);
    }
  if(collision_normal_==Collision_normal::NORMAL_IJ)
    {
      normal_average_ij(mesh, sommets_facets, param, norm);
    }
}

void Collision_Model_FT_ellipsoid::compute_dU_cp(DoubleTab& dU, DoubleTab& cp, collision_parameters& param, bool is_particle_particle_collision,
                                                 const IntLists& compo_sommets,const particle_properties& part_prop, const Maillage_FT_Disc& mesh)
{
  const DoubleTab& sommets = mesh.sommets();
  if(is_particle_particle_collision)
    {
      collision_point(mesh, param.i_closest, param.j_closest, param.i_facet, param.j_facet, cp);
      DoubleTab velocity_node_i(dimension), r_i(dimension);
      DoubleTab velocity_node_j(dimension), r_j(dimension);

      for(int d=0; d<dimension; ++d)
        {
          r_i(d) = cp(d) - part_prop.position(param.particle_i, d);
          r_j(d) = cp(d) - part_prop.position(param.particle_j, d);
        }

      for(int d=0; d<dimension; ++d)
        {
          velocity_node_i(d) = part_prop.velocity(param.particle_i,d) + part_prop.rot_velocity(param.particle_i,(d+1)%3) * r_i((d+2)%3) - part_prop.rot_velocity(param.particle_i,(d+2)%3) * r_i((d+1)%3);
          velocity_node_j(d) = part_prop.velocity(param.particle_j,d) + part_prop.rot_velocity(param.particle_j,(d+1)%3) * r_j((d+2)%3) - part_prop.rot_velocity(param.particle_j,(d+2)%3) * r_j((d+1)%3);
        }
      for (int d = 0; d < dimension; d++)
        {
          dU(d) = velocity_node_i(d) - velocity_node_j(d);
        }
    }
  else
    {
      DoubleTab r_i(dimension);
      for(int d=0; d<dimension; ++d)
        {
          cp(d) = sommets(param.i_closest,d);
          r_i(d) =cp(d) - part_prop.position(param.particle_i, d);
        }
      for (int d = 0; d < dimension; d++)
        dU(d) = part_prop.velocity(param.particle_i,d) + (part_prop.rot_velocity(param.particle_i,(d+1)%3) * r_i((d+2)%3) - part_prop.rot_velocity(param.particle_i,(d+2)%3) * r_i((d+1)%3));
    }

}

void Collision_Model_FT_ellipsoid::compute_dX_boundary(collision_parameters& param, IntLists const& compo_sommets, const Maillage_FT_Disc& mesh)
{
  auto sommets = mesh.sommets();
  if(sommets.dimension(0)==0) {return;} //for Procs that don't have any vertex
  int ind_wall = param.particle_j - nb_particles_tot_;
  int ori = ind_wall < dimension ? ind_wall : ind_wall - dimension;
  double dX_min =  std::numeric_limits<double>::max();
  DoubleTab dX_loc(dimension);
  for(int i_som = 0; i_som < compo_sommets[param.particle_i].size(); ++i_som) //recherche le sommet le plus proche du mur
    {
      dX_loc(ori) = std::fabs(sommets(compo_sommets[param.particle_i][i_som],ori) -  (origin_(ori) + (ind_wall>2)*domain_dimensions_(ori)));
      if(dX_loc(ori) < dX_min )
        {
          param.i_closest = compo_sommets[param.particle_i][i_som];
          dX_min = dX_loc(ori);
          param.i_j_are_close = true;
        }
    }
  param.dX(0) = dX_min;
}

//used
void Collision_Model_FT_ellipsoid::detect_collision(int part_i, std::vector<collision_parameters>& col_param, int nb_part_j, int start_j, const Octree_Double& octree, IntLists const& compo_sommets,
                                                    const ArrOfInt& compo_connexes_fa7, const ArrOfInt& compo_connexes_sommets, const Maillage_FT_Disc& mesh)
{
  //initialization
  DoubleTab dX_min_norm(nb_part_j);
  std::map<int, int> mapping_partj_idloc;
  for (int loc_compo_j =start_j; loc_compo_j < nb_part_j; loc_compo_j++)
    {
      col_param[loc_compo_j].particle_i = part_i;
      col_param[loc_compo_j].particle_j = get_particle_j(part_i,loc_compo_j);
      col_param[loc_compo_j].part_part_collision = col_param[loc_compo_j].particle_j < nb_particles_tot_;
      if(!col_param[loc_compo_j].part_part_collision) {compute_dX_boundary(col_param[loc_compo_j], compo_sommets, mesh); }
      dX_min_norm(loc_compo_j) = std::numeric_limits<double>::max();
      mapping_partj_idloc[col_param[loc_compo_j].particle_j] = loc_compo_j;
    }
  auto sommets = mesh.sommets();
  auto facets = mesh.facettes();
  DoubleTab dX_loc(dimension);
  DoubleTab dX_min(nb_part_j, dimension);
  DoubleTab coord(dimension);
  ArrOfInt liste_facettes;
  ArrOfInt liste_sommets;
  int num_sommet = -1;
  int num_compo_j = -1;
  double dist = 0.;
  double distmax =1.5e-4;
  for (int i = 0; i < compo_sommets[part_i].size(); ++i) //loop on all the vertices of the compo i
    {


      int i_global = compo_sommets[part_i][i];
      for (int d = 0; d < dimension; ++d)
        coord(d) = sommets(i_global, d);
      octree.search_elements_box(coord[0] - distmax, coord[1] - distmax, coord[2] - distmax, coord[0] + distmax,
                                 coord[1] + distmax, coord[2] + distmax, liste_sommets); // get the list of vertices that are in a box centered at coord and of size distmax*2

      for (int s = 0; s < liste_sommets.size_array(); ++s)
        {
          num_sommet = liste_sommets[s];
          num_compo_j = compo_connexes_sommets(num_sommet);
          int loc_compo_j = mapping_partj_idloc[num_compo_j];
          if (loc_compo_j < start_j)
            continue;
          for (int d = 0; d < dimension; ++d)
            {
              dX_loc(d) = coord[d] - sommets(num_sommet, d);
            }
          dist = sqrt(local_carre_norme_vect(dX_loc));
          if (dist < dX_min_norm[loc_compo_j])
            {
              for (int d = 0; d < dimension; ++d)
                {
                  dX_min(loc_compo_j, d) = dX_loc(d);
                }
              dX_min_norm[loc_compo_j] = dist;

              col_param[loc_compo_j].i_closest = i_global;
              col_param[loc_compo_j].j_closest = num_sommet;
              col_param[loc_compo_j].i_j_are_close = true;
            }
        }
    }

  for(int loc_compo_j = start_j; loc_compo_j<nb_part_j; ++loc_compo_j)
    {
      if(!col_param[loc_compo_j].part_part_collision) continue;
      for(int d=0; d<dimension; ++d)
        {
          col_param[loc_compo_j].dX(d) = dX_min(loc_compo_j, d);
        }
    }
}

void Collision_Model_FT_ellipsoid::compute_lagrangian_contact_forces(const Fluide_Diphasique& two_phase_fluid,
                                                                     const DoubleTab& particles_position,
                                                                     const DoubleTab& particles_velocity,
                                                                     const DoubleTab& particles_rot_velocity,
                                                                     const double& deltat_simu,
                                                                     const Maillage_FT_Disc& mesh) {Cerr<<"Do nothing\n";}

void Collision_Model_FT_ellipsoid::compute_lagrangian_contact_forces(const Fluide_Diphasique& two_phase_fluid,
                                                                     const particle_properties& part_prop,
                                                                     const double& deltat_simu,
                                                                     const Maillage_FT_Disc& mesh)
{
  double t=part_prop.t;
  const int& id_fluid_phase= two_phase_fluid.get_id_fluid_phase();
  const int& id_solid_phase=1-id_fluid_phase;
  const auto& solid_particle=ref_cast(Solid_Particle_ellipsoid,two_phase_fluid.fluide_phase(id_solid_phase));
  const auto& incompressible_fluid=ref_cast(Fluide_Incompressible,
                                            two_phase_fluid.fluide_phase(id_fluid_phase));
  const double& solid_density = solid_particle.masse_volumique().valeurs()(0, 0);
  const double& fluid_density = incompressible_fluid.masse_volumique().valeurs()(0, 0);
  const double& fluid_viscosity  = fluid_density
                                   * incompressible_fluid.viscosite_cinematique().valeurs()(0, 0);
  // const double& long_radius=solid_particle.get_max_radius();
  const double& volume=solid_particle.get_volume();
  const double& density=solid_particle.get_mass()/volume;
  auto particles_position = part_prop.position;
  auto particles_velocity = part_prop.velocity;
  auto particles_rot_velocity = part_prop.rot_velocity;
  auto particles_inertia_tensor = part_prop.inertia_tensor;
  auto particles_volume = part_prop.volume;

  auto facets = mesh.facettes();
  mesh.compute_mesh_size();

  const double& e_dry=solid_particle.get_e_dry();
  const double min_threshold=1e-10;
  double friction_coef=0.15; //XXX to pass in data file
  DoubleTab dX(dimension), dU(dimension), norm(dimension), tang(dimension), collision_point(dimension);
  DoubleTab dUt(dimension); //tangential velocity
  lagrangian_contact_forces_=0;
  lagrangian_contact_moments_=0;
  collision_number_=0;
  particles_collision_number_=0;

  /*--------------all these functions should be refactored as search_connex_components_local_FT and compute_global_connex_components_FT are called in each one*/
  test_connex_compo(mesh);
  IntLists compo_sommets; //compo_sommet[i] contient les indices des sommets (dans le proc) composant la compo i
  connec_compo_sommets(mesh, compo_sommets);
  IntLists sommets_facets; //compo_sommet[i] contient les indices des sommets composant la compo i
  connec_sommets_fa7(mesh, sommets_facets);
  const int& nb_fa7 = mesh.nb_facettes();
  ArrOfInt compo_connexes_fa7(nb_fa7); //compo_connexes_fa7(fa7) donne l'indice de la compo (particule) contenant la facette fa7
  int n = search_connex_components_local_FT(mesh, compo_connexes_fa7);
  compute_global_connex_components_FT(mesh, compo_connexes_fa7, n);
  IntLists compo_connexe_facets;
  connec_compo_facettes(mesh, compo_connexe_facets); //compo_connexe_facets[i] contient les indices des facettes composant la compo i
  ArrOfInt compo_connecs_sommets(mesh.sommets().dimension(0));
  compo_connexe_sommets(mesh, compo_connexes_fa7, compo_connecs_sommets);


  std::vector<collision_parameters> collisions_param;
  Octree_Double octree;
  Z_Curve_Link_Cell ZLC;
  std::vector<std::vector<collision_parameters>> collisions_param_zlc;
  collisions_param_zlc.resize(nb_real_particles_);
  for (int ind_particle_i = 0; ind_particle_i < nb_real_particles_; ind_particle_i++) collisions_param_zlc[ind_particle_i].resize(nb_real_particles_ + 2*dimension);


  //XXX check time
  std::chrono::steady_clock::time_point begin_octree;
  std::chrono::steady_clock::time_point end_octree;
  std::chrono::steady_clock::time_point begin_MLO;
  std::chrono::steady_clock::time_point end_MLO;
  std::chrono::steady_clock::time_point begin_naive;
  std::chrono::steady_clock::time_point end_naive;
  auto time_octree =  std::chrono::duration_cast<std::chrono::nanoseconds>(end_octree - begin_octree).count();
  auto construction_octree = std::chrono::duration_cast<std::chrono::nanoseconds>(end_MLO - begin_MLO).count();
  auto closest_octree = std::chrono::duration_cast<std::chrono::nanoseconds>(end_MLO - begin_MLO).count();
  auto time_mlo = std::chrono::duration_cast<std::chrono::nanoseconds>(end_MLO - begin_MLO).count();
  auto construction_mlo = std::chrono::duration_cast<std::chrono::nanoseconds>(end_MLO - begin_MLO).count();
  auto closest_mlo = std::chrono::duration_cast<std::chrono::nanoseconds>(end_MLO - begin_MLO).count();
  auto time_naive = std::chrono::duration_cast<std::chrono::nanoseconds>(end_naive - begin_naive).count();
  ///

  if(detection_option==Detection_Option::ZLC)
    {
      begin_MLO = std::chrono::steady_clock::now();
      ZLC.build(mesh.sommets(), link_cell_size * solid_particle.get_max_radius(), compo_connecs_sommets);
      end_MLO = std::chrono::steady_clock::now();
      construction_mlo += std::chrono::duration_cast<std::chrono::nanoseconds>(end_MLO - begin_MLO).count();
      //Contact detection already done here
      begin_MLO = std::chrono::steady_clock::now();
      ZLC.find_closest(mesh.sommets(),compo_connecs_sommets, collisions_param_zlc, origin_, domain_dimensions_ );
      end_MLO = std::chrono::steady_clock::now();
      closest_mlo += std::chrono::duration_cast<std::chrono::nanoseconds>(end_MLO - begin_MLO).count();
    }

  //XXX We could also detect all closest nodes here for octree and naive
  // and check only collisions in the next loop
  else if(detection_option==Detection_Option::OCTREE)
    {
      begin_octree = std::chrono::steady_clock::now();
      octree.build_nodes(mesh.sommets(), 0.,0.);
      end_octree = std::chrono::steady_clock::now();
      construction_octree+= std::chrono::duration_cast<std::chrono::nanoseconds>(end_octree - begin_octree).count();
    }



  /*******************************END DEBUG******************************* */

  for (int ind_particle_i = 0; ind_particle_i < nb_real_particles_; ind_particle_i++)
    {
      int particle_i=get_particle_i(ind_particle_i);
      int nb_particles_j=get_nb_particles_j(ind_particle_i);
      int ind_start_part_j=get_ind_start_particles_j(ind_particle_i);

      collisions_param.resize(nb_particles_j);

      if(detection_option==Detection_Option::OCTREE)
        {
          begin_octree = std::chrono::steady_clock::now();
          detect_collision(particle_i, collisions_param, nb_particles_j, ind_start_part_j, octree, compo_sommets, compo_connexes_fa7, compo_connecs_sommets, mesh);
          end_octree = std::chrono::steady_clock::now();
          closest_octree+= std::chrono::duration_cast<std::chrono::nanoseconds>(end_octree - begin_octree).count();
        }

      Matrice_Dense Ji(dimension, dimension);
      DoubleTab ri(dimension);
      DoubleTab Omega_i(dimension);
      DoubleTab tangential_force_contact(dimension);

      for(int i=0; i<dimension; ++i)
        {
          for(int j=0; j<dimension; ++j)
            {
              if(std::fabs(particles_inertia_tensor(particle_i, i, j))>1e-15)
                Ji(i,j) = particles_inertia_tensor(particle_i, i, j)*density; //particles_inertia_tensor is computed without density in Transport_Interfaces_FT_Disc
            }
        }

      for (int ind_particle_j =ind_start_part_j; ind_particle_j < nb_particles_j; ind_particle_j++)
        {
          dX = 0.;
          dU = 0;
          norm = 0.;
          double dist_between_particles = std::numeric_limits<double>::max();
          int particle_j=get_particle_j(ind_particle_i,ind_particle_j);
          int is_particle_particle_collision = particle_j < nb_particles_tot_;

          if(detection_option==Detection_Option::ZLC)
            {
              collisions_param[ind_particle_j] = collisions_param_zlc[ind_particle_i][ind_particle_j];
            }
          collisions_param[ind_particle_j].part_part_collision = particle_j < nb_particles_tot_;
          collisions_param[ind_particle_j].particle_i = particle_i;
          collisions_param[ind_particle_j].particle_j = particle_j;




          if(!( (is_particle_particle_collision && compo_sommets[particle_j].size()==0) || compo_sommets[particle_i].size()==0))  //particle_j has no vertex in the local proc, so we continue
            {

              if(detection_option==Detection_Option::NAIVE)
                {
                  std::cout<<"use naive"<<std::endl;
                  begin_naive = std::chrono::steady_clock::now();
                  compute_dX(collisions_param[ind_particle_j].dX, collisions_param[ind_particle_j], particles_position, is_particle_particle_collision, compo_sommets, sommets_facets, mesh);
                  end_naive = std::chrono::steady_clock::now();
                  time_naive+= std::chrono::duration_cast<std::chrono::nanoseconds>(end_naive - begin_naive).count();
                }

              if(is_particle_particle_collision)
                {
                  std::cout<<"is part part col"<<std::endl;
                  if(collisions_param[ind_particle_j].i_j_are_close)
                    {
                      compute_norm(norm, collisions_param[ind_particle_j], sommets_facets, mesh);
                    }
                }
              else
                {
                  int ind_wall = particle_j - nb_particles_tot_;
                  int ori = ind_wall < dimension ? ind_wall : ind_wall - dimension;
                  norm(ori) = (ind_wall < 3 ? 1 : -1); //La normale correspond à la normale à la paroi (seuls les parallélépipèdes sont considérés)

                }


              if(is_particle_particle_collision)
                {
                  if(collisions_param[ind_particle_j].i_j_are_close)
                    {
                      if(closest_indices(ind_particle_i,ind_particle_j)==-1 || closest_indices(ind_particle_j,ind_particle_i)==-1)
                        {
                          int i = collisions_param[ind_particle_j].i_closest;
                          int j = collisions_param[ind_particle_j].j_closest;
                          for(int d=0; d<dimension; ++d)
                            {
                              dX(d) = mesh.sommets()(i,d) - mesh.sommets()(j,d);
                            }
                          dist_between_particles = -local_prodscal(dX,norm);
                          if(dist_between_particles<=0) //we have a collision so we store the indices
                            {
                              closest_indices(ind_particle_j,ind_particle_i)=collisions_param[ind_particle_j].i_closest;
                              closest_indices(ind_particle_i,ind_particle_j)=collisions_param[ind_particle_j].j_closest;
                            }
                        }
                      else //we already have a pair from previous time step
                        {
                          int i = closest_indices(ind_particle_j,ind_particle_i);
                          int j = closest_indices(ind_particle_i,ind_particle_j);
                          for(int d=0; d<dimension; ++d)
                            {
                              dX(d) = mesh.sommets()(i,d) - mesh.sommets()(j,d);
                            }
                          dist_between_particles = -local_prodscal(dX,norm);
                          if(dist_between_particles>0) //we are not on a collision so we uptdate
                            {
                              closest_indices(ind_particle_i,ind_particle_j)=collisions_param[ind_particle_j].j_closest;
                              closest_indices(ind_particle_j,ind_particle_i)=collisions_param[ind_particle_j].i_closest;
                              i = closest_indices(ind_particle_j,ind_particle_i);
                              j = closest_indices(ind_particle_i,ind_particle_j);
                              for(int d=0; d<dimension; ++d)
                                {
                                  dX(d) = mesh.sommets()(i,d) - mesh.sommets()(j,d);
                                }
                            }
                          else //we are still in the collision
                            {
                              collisions_param[ind_particle_j].j_closest = closest_indices(ind_particle_i,ind_particle_j);
                              collisions_param[ind_particle_j].i_closest = closest_indices(ind_particle_j,ind_particle_i);
                            }
                        }


                      dist_between_particles = -local_prodscal(dX,norm); //normal penetration distance
                      if(dist_between_particles>0) //no collision
                        {
                          closest_indices(ind_particle_i,ind_particle_j)=-1;
                          closest_indices(ind_particle_j,ind_particle_i)=-1;
                        }
                    }
                  else
                    {
                      closest_indices(ind_particle_i,ind_particle_j)=-1;
                      closest_indices(ind_particle_j,ind_particle_i)=-1;
                    }
                }
              else
                {
                  if(collisions_param[ind_particle_j].i_j_are_close)
                    {
                      double dist_gravity_center = sqrt(local_carre_norme_vect(collisions_param[ind_particle_j].dX));
                      dist_between_particles = dist_gravity_center - activation_distance_ ;
                    }
                }
            }
          std::ofstream ffn, fft, fm, fu, fd, ft;
          std::string path;
          path = fichier_debug + "normal_force_" + std::to_string(particle_i)+"_P"+std::to_string(Process::me())+".txt";
          ffn.open(path, std::ios::app);
          path = fichier_debug + "tangent_force_" + std::to_string(particle_i)+"_P"+std::to_string(Process::me())+".txt";
          fft.open(path, std::ios::app);
          path = fichier_debug + "moment_" + std::to_string(particle_i)+"_P"+std::to_string(Process::me())+".txt";
          fm.open(path, std::ios::app);
          path = fichier_debug + "vitesse_impact_" + std::to_string(particle_i)+"_P"+std::to_string(Process::me())+".txt";
          fu.open(path, std::ios::app);
          path = fichier_debug + "distance_" + std::to_string(particle_i)+"_P"+std::to_string(Process::me())+".txt";
          fd.open(path, std::ios::app);

          // Check if the current proc is the one that has to compute the force
          DoubleTab dist(1); //can use mp_min_for_each_item only with TRUSTArray
          dist(0) = dist_between_particles;
          mp_min_for_each_item(dist);
          if(dist(0)!=dist_between_particles) continue;



          F_now_(particle_i, particle_j) = 0;
          double max_dist = 0.;
          if (dist_between_particles <= 0) // contact
            {
              compute_dU_cp(dU, collision_point, collisions_param[ind_particle_j], is_particle_particle_collision, compo_sommets, part_prop, mesh);
              std::cout<<"Computed du cp"<<std::endl;
              std::ofstream fcol;
              path = fichier_debug + "is_collision_" + std::to_string(particle_i)+"_P"+std::to_string(Process::me())+".txt";
              fcol.open(path, std::ios::app);
              fcol<<collisions_param[ind_particle_j].i_closest<<" "<<collisions_param[ind_particle_j].j_closest<<std::endl;
              max_dist = std::max(max_dist, -dist_between_particles);

              if(is_particle_particle_collision) {dist_between_particles*=-1;}
              add_collision(particle_i,particle_j,is_particle_particle_collision);

              double dU_scal_norm = local_prodscal(dU,norm);
              DoubleTab dUn(dimension);
              for (int d = 0; d < dimension; d++)
                dUn(d) = dU_scal_norm * norm(d);

              const double impact_velocity = sqrt(local_carre_norme_vect(dUn));

              F_now_(particle_i, particle_j) = 1;
              int is_start_of_collision = F_now_(particle_i, particle_j) >
                                          F_old_(particle_i, particle_j); // We need to know
              // if this is the first time step of the collision to compute the impact velocity
              double a = solid_particle.get_x_radius();
              double b = solid_particle.get_y_radius();
              double c = solid_particle.get_z_radius();
              double effective_radius = a*b*c * (pow(collision_point(0)/(a*a),2) + pow(collision_point(1)/(b*b),2) + pow(collision_point(2)/(c*c),2));

              // const double effective_radius = is_particle_particle_collision ? solid_particle.get_equivalent_radius()/2 :
              //                                 solid_particle.get_equivalent_radius();
              const double impact_Stokes = solid_density * 2 * effective_radius * impact_velocity /
                                           (9 * fluid_viscosity);
              if (is_start_of_collision)
                {std::cout<<"Compute e_eff"<<std::endl; e_eff_(particle_i,particle_j)=e_dry*compute_ewet_legendre(impact_Stokes);}
              std::cout<<"e_eff = " <<e_eff_(particle_i,particle_j)<<" e_dry = "<<e_dry<<endl;
              std::cout<<solid_density<<" "<<effective_radius<<" "<<impact_velocity<<" "<<fluid_viscosity<<std::endl;
              std::cout<<impact_Stokes<<" "<<compute_ewet_legendre(impact_Stokes)<<std::endl;;
              DoubleTab force_contact=compute_contact_force(
                                        dist_between_particles,
                                        norm,
                                        dUn,
                                        particle_i,
                                        particle_j,
                                        dU_scal_norm<=0,
                                        is_particle_particle_collision);

              ffn<<"Proc "<<Process::me()<<" "<<t<<" "<<force_contact(0)<<" "<<force_contact(1)<<" "<<force_contact(2)<<"\n";
              fu<<"Proc "<<Process::me()<<" "<<t<<" "<<dUn(0)<<" "<<dUn(1)<<" "<<dUn(2)<<" "<<impact_velocity<<"\n";
              fd<<"Proc "<<Process::me()<<" "<<t<<" "<<dist_between_particles<<" "<<norm(0)<<" "<<norm(1)<<" "<<norm(2)<<"\n";
              if(is_particle_particle_collision)
                {

                  std::ofstream fdd;
                  path = fichier_debug + "distance.txt";
                  fdd.open(path, std::ios::app);
                  fdd<<dist_between_particles<<std::endl;
                }
              for(int d=0; d<dimension; ++d)
                {
                  dUt(d) = dU(d) - dUn(d);
                }
              double dUt_norm = 0.;
              dUt_norm = sqrt(local_carre_norme_vect(dUt));
              std::cout<<"dut_norm: "<<dUt_norm<<std::endl;
              for(int d=0; d<dimension; ++d)
                {
                  tang(d) = dUt(d)/dUt_norm;
                }
              double tangential_displacement = sqrt(local_carre_norme_vect(dUt)) * deltat_simu;
              if (is_start_of_collision)
                e_eff_t(particle_i,particle_j)=0.39*compute_ewet_legendre(impact_Stokes); //XXX 0.39 comes from Jain 2019. Has to be a .data parameter

              if(is_start_of_collision)
                {
                  tangential_force_contact = 0.;
                }
              compute_tangential_contact_force(tangential_displacement, tang, e_eff_t(particle_i,particle_j), friction_coef, force_contact, dU_scal_norm<=0, is_particle_particle_collision,tangential_force_contact);
              for(int d=0; d<dimension; ++d)
                {
                  force_contact(d) += tangential_force_contact(d);
                }
              fft<<"Proc "<<Process::me()<<" "<<t<<" "<<tangential_force_contact(0)<<" "<<tangential_force_contact(1)<<" "<<tangential_force_contact(2)<<"\n";



              Matrice_Dense Jj(dimension, dimension);
              Matrice_Dense Jj_temp(dimension, dimension);
              DoubleTab rj(dimension);
              DoubleTab Omega_j(dimension);

              for(int d=0; d<dimension; ++d)
                {
                  ri(d) = collision_point(d) - particles_position(particle_i,d) ;
                  Omega_i(d) = particles_rot_velocity(particle_i,d);
                  if(!is_particle_particle_collision)
                    continue;
                  rj(d) = collision_point(d) - particles_position(particle_j,d);
                  Omega_j(d) = particles_rot_velocity(particle_j,d);
                }
              DoubleTab moment_contact_i=compute_contact_moment(Ji,force_contact, ri, Omega_i);
              DoubleTab moment_contact_j(dimension);
              if(is_particle_particle_collision)
                {
                  for(int i=0; i<dimension; ++i)
                    {
                      for(int j=0; j<dimension; ++j)
                        {
                          Jj(i,j) = particles_inertia_tensor(particle_j, i, j)*density; //particles_inertia_tensor is computed without density in Transport_Interfaces_FT_Disc;
                        }
                    }
                  moment_contact_j=compute_contact_moment(Jj,force_contact, rj, Omega_j);
                }
              fm<<"Proc "<<Process::me()<<" "<<t<<" "<<moment_contact_i(0)<<" "<<moment_contact_i(1)<<" "<<moment_contact_i(2)<<"\n";

              ffn.close();
              fft.close();
              fm.close();
              for (int d = 0; d < dimension; d++)
                {
                  lagrangian_contact_forces_(particle_i, d) += fabs(force_contact(d)) <=
                                                               min_threshold ? 0 : force_contact(d) / volume;
                  lagrangian_contact_moments_(particle_i, d) += fabs(moment_contact_i(d)) <=
                                                                min_threshold ? 0 : moment_contact_i(d) * density;
                  if (!is_particle_particle_collision)
                    continue; // wall collision, no force to apply on the wall
                  lagrangian_contact_forces_(particle_j, d) -= fabs(force_contact(d)) <=
                                                               min_threshold ? 0 :  force_contact(d) / volume;
                  lagrangian_contact_moments_(particle_j, d) -= fabs(moment_contact_j(d)) <=
                                                                min_threshold ? 0 : moment_contact_j(d) * density;//XXX assumption: all particles have the same density
                }

              F_old_(particle_i, particle_j) = F_now_(particle_i, particle_j);
            }


        }
      //XXX debug file, energy
      DoubleTab v(dimension), om(dimension);
      DoubleTab InertiaOmega(dimension);
      for(int d=0; d<dimension; ++d) {v(d) = particles_velocity(particle_i,d); om(d) = particles_rot_velocity(particle_i,d);}
      Ji.ajouter_multvect_(om,InertiaOmega);
      double energy_potentielle = 9.81 * solid_particle.get_mass() * particles_position(particle_i,1);
      double energy_cinetique = 0.5 * solid_particle.get_mass() * local_carre_norme_vect(v);
      double energy_rotation = 0.5 * local_prodscal(InertiaOmega,om);
      double energy = energy_potentielle + energy_cinetique + energy_rotation;
      std::string path;
      std::fstream f;
      path = fichier_debug + "energy_" + std::to_string(Process::me())+".txt";
      f.open(path, std::ios::app);
      f<<t<<" "<<energy<<" "<<energy_potentielle<<" "<<energy_cinetique<<" "<<energy_rotation<<"\n";


    }
  auto tot_const_mlo =Process::mp_sum(construction_mlo);
  auto tot_closest_mlo =Process::mp_sum(closest_mlo);
  auto tot_const_octrcee =Process::mp_sum(construction_octree);
  auto tot_closest_octree =Process::mp_sum(closest_octree);
  std::string path;
  path = fichier_debug + "time_P"+std::to_string(Process::me())+".txt";
  std::ofstream ft;
  ft.open(path, std::ios::app);
  ft<<t<<" "<<time_octree<<" "<<time_mlo<<" "<<time_naive<<std::endl;
  ft.close();

  if(Process::me()==0)
    {
      path = fichier_debug + "mlo.txt";
      ft.open(path, std::ios::app);

      ft<<t<<" "<<tot_const_mlo<<" "<<tot_closest_mlo<<" "<<tot_const_mlo+tot_closest_mlo<<std::endl;
      ft.close();
      path = fichier_debug + "octree.txt";
      ft.open(path, std::ios::app);
      ft<<t<<" "<<tot_const_octrcee<<" "<<tot_closest_octree<<" "<<tot_const_octrcee+tot_closest_octree<<std::endl;
      ft.close();
    }

  mp_sum_for_each_item(lagrangian_contact_forces_);
  mp_sum_for_each_item(lagrangian_contact_moments_);
  mp_max_for_each_item(F_old_);
  mp_max_for_each_item(F_now_);
  mp_max_for_each_item(e_eff_);
  mp_sum_for_each_item(particles_collision_number_);
  collision_number_=Process::check_int_overflow(Process::mp_sum(collision_number_));


}

//XXX Inertia not a ref to debug only
DoubleTab Collision_Model_FT_ellipsoid::compute_contact_moment(Matrice_Dense Inertia, DoubleTab const& force, DoubleTab const& r, DoubleTab const& Omega)
// DoubleTab Collision_Model_FT_ellipsoid::compute_contact_moment(Matrice_Dense& Inertia, DoubleTab const& force, DoubleTab const& r, DoubleTab const& Omega)
{
  DoubleTab contact_moment(dimension);
  DoubleTab temp(dimension);
  DoubleTab InertiaOmega(dimension);
  Inertia.ajouter_multvect_(Omega,InertiaOmega);
  Inertia.inverse(); //Inertia is now its inverse
  for(int d = 0; d<dimension; ++d)
    {
      temp[d] = (r[(d+1)%3] * force[(d+2)%3] - r[(d+2)%3] * force[(d+1)%3]);
    }
  Inertia.ajouter_multvect_(temp,contact_moment);
  return contact_moment;
}

void Collision_Model_FT_ellipsoid::discretize_contact_forces_eulerian_field(
  const DoubleTab& volumic_phase_indicator_function,
  const Domaine_VF& domain_vf,
  const IntTab& particles_eulerian_id_number,
  const DoubleTab& particles_position,
  DoubleTab& contact_force_source_term)
{}
void Collision_Model_FT_ellipsoid::discretize_contact_forces_eulerian_field(
  const DoubleTab& volumic_phase_indicator_function,
  const Domaine_VF& domain_vf,
  const IntTab& particles_eulerian_id_number,
  const DoubleTab& particles_position,
  DoubleTab& contact_force_source_term,
  const particle_properties& part_prop)
{
  static int t=0;
  const DoubleTab Omega = part_prop.rot_velocity;
  DoubleTab Omega_i(dimension);
  const DoubleTab v = part_prop.velocity;
  const DoubleVect& interlaced_volumes=domain_vf.volumes_entrelaces();
  const DoubleTab& cg_faces=domain_vf.xv();
  const int nb_faces=interlaced_volumes.size_array();
  const IntVect& orientation = domain_vf.orientation();
  const IntTab& face_voisins=domain_vf.face_voisins();
  double max_force = 0.;
  double max_moment = 0.;
  DoubleTab Omega_r(dimension);
  DoubleTab OmegaJOmgea(dimension);
  DoubleTab inertia(dimension);
  auto particles_inertia_tensor = part_prop.inertia_tensor;
  Matrice_Dense Ji(dimension, dimension);
  // bool collision = false; // if no collision, we don't apply centrifugal force
  for (int face=0; face<nb_faces; face++)
    {
      const int left_elem=face_voisins(face,0);
      const int right_elem=face_voisins(face,1);
      int id_left = left_elem != -1 ? particles_eulerian_id_number(left_elem) : -1;
      int id_right = right_elem != -1 ? particles_eulerian_id_number(right_elem) : -1;
      const int id_number=std::max(id_left,id_right);
      if (id_number!=-1)
        {
          // if(lagrangian_contact_forces_(id_number,0)!=0 || lagrangian_contact_forces_(id_number,1)!=0 || lagrangian_contact_forces_(id_number,2)!=0) {collision = true;}
          const int ori=orientation(face);

          for(int i=0; i<dimension; ++i)
            {
              Omega_i(i) = Omega(id_number, i);
              inertia(i) = 0.;
              for(int j=0; j<dimension; ++j)
                {
                  if(std::fabs(particles_inertia_tensor(id_number, i, j))>1e-15)
                    Ji(i,j) = particles_inertia_tensor(id_number, i, j)*part_prop.density; //particles_inertia_tensor is computed without density in Transport_Interfaces_FT_Disc
                }
            }
          DoubleTab InertiaOmega(dimension);
          Ji.ajouter_multvect_(Omega_i,InertiaOmega);
          for(int d=0; d<dimension; ++d)
            {
              Omega_r(d) = Omega_i((d+1)%3) * (cg_faces(face,(ori+2)%3) - particles_position(id_number,(ori+2)%3)) - Omega_i((d+2)%3) * (cg_faces(face,(ori+1)%3) - particles_position(id_number,(ori+1)%3) );
              OmegaJOmgea(d) = Omega_i((d+1)%3) * InertiaOmega((d+2)%3) - Omega_i((d+2)%3) * InertiaOmega((d+1)%3);
            }
          Ji.inverse();//Ji is now its inverse
          Ji.ajouter_multvect_(OmegaJOmgea,inertia);
          //Omega_r = \omega \times r
          //inertia = J^{-1} \omega \times J \omega


          contact_force_source_term(face)=(1-volumic_phase_indicator_function(face))
                                          *interlaced_volumes(face)*(
                                            lagrangian_contact_forces_(id_number,ori)
                                            + (lagrangian_contact_moments_(id_number,(ori+1)%3) * (cg_faces(face,(ori+2)%3) - particles_position(id_number,(ori+2)%3)) -
                                               lagrangian_contact_moments_(id_number,(ori+2)%3) * (cg_faces(face,(ori+1)%3) - particles_position(id_number,(ori+1)%3) ))
                                            -  (inertia((ori+1)%3) * (cg_faces(face,(ori+2)%3) - particles_position(id_number,(ori+2)%3)) -
                                                inertia((ori+2)%3) * (cg_faces(face,(ori+1)%3) - particles_position(id_number,(ori+1)%3) )) +
                                            /*collision **/ part_prop.density * ( Omega_i((ori+1)%3) * Omega_r((ori+2)%3) - Omega_i((ori+2)%3) * Omega_r((ori+1)%3))
                                          );
        }
    }
  std::ofstream f,g;
  f.open("check_energy_files/max_force.txt", std::ios::app);
  g.open("check_energy_files/max_moment.txt", std::ios::app);
  f<<t<<" "<<max_force<<"\n";
  g<<t++<<" "<<max_moment<<"\n";
}

void Collision_Model_FT_ellipsoid::test_connex_compo(Maillage_FT_Disc const& maillage)
{
  const int& nb_fa7 = maillage.nb_facettes();

  ArrOfInt compo_connexes_fa7(nb_fa7); //compo_connexes_fa7(fa7) donne l'indice de la compo (particule) contenant la facette fa7
  int n = search_connex_components_local_FT(maillage, compo_connexes_fa7);
  int nb_compo_tot=compute_global_connex_components_FT(maillage, compo_connexes_fa7, n);
  std::ofstream f;
  std::string path;
  path = "test_connex_compo/Proc_"+std::to_string(Process::me())+".txt";
  f.open(path);
  f<<"nb_fa7 = "<<nb_fa7<<std::endl;
  f<<"n = "<<n<<std::endl;
  f<<"nb_compo_tot = "<<nb_compo_tot<<std::endl;

}