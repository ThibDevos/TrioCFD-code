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
                                                 DoubleTab& dX, bool chek_cg)
{
  auto sommets = mesh.sommets();
  DoubleTab dX_min(dimension);
  double dist_cg = 0.;
  param.i_closest = -1;
  param.j_closest = -1;
  for (int d = 0; d < dimension; ++d)
    dX_min(d) = 1000.;

  if (chek_cg)
    {
    }
  else
    {
      std::cout<<"["<<Process::me()<<"] : "<<"Particule "<<param.particle_i <<" compare with particle "<<param.particle_j<<std::endl;
      std::cout<<"["<<Process::me()<<"] : "<<compo_sommets[param.particle_i].size()<<" vs "<<compo_sommets[param.particle_j].size()<<std::endl;
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
          closest_nodes(compo_sommets, mesh, param ,dX, false); //check_cg not implemented
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
          cp(d) = sommets(compo_sommets[param.particle_i][param.i_closest],d);
          r_i(d) =cp(d) - part_prop.position(param.particle_i, d);
        }
      for (int d = 0; d < dimension; d++)
        dU(d) = part_prop.velocity(param.particle_i,d) + (part_prop.rot_velocity(param.particle_i,(d+1)%3) * r_i((d+2)%3) - part_prop.rot_velocity(param.particle_i,(d+2)%3) * r_i((d+1)%3));
    }

}

// void Collision_Model_FT_ellipsoid::compute_dX_dU_normal(DoubleTab& dX, DoubleTab& dU, DoubleTab& norm, DoubleTab& cp, int particle,
//                                                         int neighbor, const DoubleTab& particles_position, const
//                                                         DoubleTab& particles_velocity, DoubleTab const& particles_rot_velocity, const bool is_particle_particle_collision,
//                                                         const IntLists& compo_sommets,const IntLists& sommets_facets, const Maillage_FT_Disc& mesh )
// {
//   const DoubleTab& sommets = mesh.sommets();
//   if (is_particle_particle_collision)
//     {
//       bool check_cg = true;
//       int i_closest, j_closest;
//       bool i_facet = false, j_facet=false;
//       if(collision_detection_==Collision_detection::CLOSEST)
//         {
//           closest_nodes(compo_sommets, mesh, particle, neighbor, particles_position ,dX, i_closest, j_closest, false, i_facet, j_facet); //check_cg not implemented
//         }

//       if(collision_detection_==Collision_detection::DEEPEST)
//         {
//           deepest_points(compo_sommets,sommets_facets,mesh,particle, neighbor, particles_position, check_cg, dX, i_closest, j_closest, i_facet, j_facet);
//         }
//       if(collision_normal_==Collision_normal::NORMAL_I)
//         {
//           normal_i(mesh, sommets_facets, i_facet, i_closest, norm);
//         }
//       if(collision_normal_==Collision_normal::NORMAL_IJ)
//         {
//           normal_average_ij(mesh, sommets_facets, i_facet, i_closest, j_facet, j_closest, norm);
//         }
//       collision_point(mesh, i_closest, j_closest, i_facet, j_facet, cp);
//       DoubleTab const& cg_fa7 = mesh.get_gravity_center_fa7();
//       DoubleTab velocity_node_i(dimension), r_i(dimension);
//       DoubleTab velocity_node_j(dimension), r_j(dimension);
//       if(i_facet)
//         {
//           for(int d=0; d<dimension; ++d)
//             {
//               r_i(d) = cg_fa7(i_closest,d) - particles_position(particle, d);
//             }
//         }
//       else
//         {
//           for(int d=0; d<dimension; ++d)
//             {
//               r_i(d) = sommets(i_closest,d) - particles_position(particle, d);
//             }
//         }
//       if(j_facet)
//         {
//           for(int d=0; d<dimension; ++d)
//             {
//               r_j(d) = cg_fa7(j_closest,d) - particles_position(neighbor, d);
//             }
//         }
//       else
//         {
//           for(int d=0; d<dimension; ++d)
//             {
//               r_j(d) = sommets(j_closest,d) - particles_position(neighbor, d);
//             }
//         }
//       for(int d=0; d<dimension; ++d)
//         {
//           velocity_node_i(d) = particles_velocity(particle,d) + particles_rot_velocity(particle,(d+1)%3) * r_i((d+2)%3) - particles_rot_velocity(particle,(d+2)%3) * r_i((d+1)%3);
//           velocity_node_j(d) = particles_velocity(neighbor,d) + particles_rot_velocity(neighbor,(d+1)%3) * r_j((d+2)%3) - particles_rot_velocity(neighbor,(d+2)%3) * r_j((d+1)%3);
//         }
//       for (int d = 0; d < dimension; d++)
//         {
//           dU(d) = velocity_node_i(d) - velocity_node_j(d);
//         }
//     }
//   else
//     {
//       if(sommets.dimension(0)==0) {return;} //for Procs that don't have any vertex
//       int ind_wall = neighbor - nb_particles_tot_;
//       int ori = ind_wall < dimension ? ind_wall : ind_wall - dimension;
//       double dX_min =  std::numeric_limits<double>::max();
//       int i_som_closest=-1;
//       for(int i_som = 0; i_som < compo_sommets[particle].size(); ++i_som) //recherche le sommet le plus proche du mur
//         {
//           dX(ori) = std::fabs(sommets(compo_sommets[particle][i_som],ori) -  (origin_(ori) + (ind_wall>2)*domain_dimensions_(ori)));
//           if(dX(ori) < dX_min ) i_som_closest = i_som;
//           dX_min = dX(ori) <= dX_min ? dX(ori) : dX_min;
//         }
//       dX(ori) = dX_min;
//       DoubleTab r_i(dimension);
//       for(int d=0; d<dimension; ++d)
//         {
//           cp(d) = sommets(compo_sommets[particle][i_som_closest],d);
//           r_i(d) =cp(d) - particles_position(particle, d);
//         }

//       for (int d = 0; d < dimension; d++)
//         dU(d) = particles_velocity(particle,d) + (particles_rot_velocity(particle,(d+1)%3) * r_i((d+2)%3) - particles_rot_velocity(particle,(d+2)%3) * r_i((d+1)%3));

//       norm(ori) = (ind_wall < 3 ? 1 : -1); //La normale correspond à la normale à la paroi (seuls les parallélépipèdes sont considérés)


//     }
// }
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


  const double& e_dry=solid_particle.get_e_dry();
  const double min_threshold=1e-10;
  double friction_coef=0.15; //XXX to pass in data file
  DoubleTab dX(dimension), dU(dimension), norm(dimension), tang(dimension), collision_point(dimension);
  DoubleTab dUt(dimension); //tangential velocity
  lagrangian_contact_forces_=0;
  lagrangian_contact_moments_=0;
  collision_number_=0;
  particles_collision_number_=0;

  test_connex_compo(mesh);
  IntLists compo_sommets; //compo_sommet[i] contient les indices des sommets (dans le proc) composant la compo i
  connec_compo_sommets(mesh, compo_sommets);
  IntLists sommets_facets; //compo_sommet[i] contient les indices des sommets composant la compo i
  connec_sommets_fa7(mesh, sommets_facets);
  IntLists compo_connexe_facets; //compo_connexes_fa7(fa7) donne l'indice de la compo (particule) contenant la facette fa7
  connec_compo_facettes(mesh, compo_connexe_facets);
  for (int ind_particle_i = 0; ind_particle_i < nb_real_particles_; ind_particle_i++)
    {
      int particle_i=get_particle_i(ind_particle_i);
      int nb_particles_j=get_nb_particles_j(ind_particle_i);
      int ind_start_part_j=get_ind_start_particles_j(ind_particle_i);
      std::cout<<"Je suis le proc "<<Process::me()<<" et j'ai "<<nb_particles_j<<" nb_particles_j et je demarre a "<<ind_start_part_j<<std::endl;
      // if(compo_sommets[ind_particle_i].size()==0)continue;

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
          std::cout<<"                    Je suis le proc "<<Process::me()<<" et j "<<particle_j<<" part_j"<<std::endl;
          int is_particle_particle_collision = particle_j < nb_particles_tot_;
          collision_parameters param;
          if(!( (is_particle_particle_collision && compo_sommets[particle_j].size()==0) || compo_sommets[particle_i].size()==0))  //particle_j has no vertex in the local proc, so we continue
            {
              param.particle_i = particle_i;
              param.particle_j = particle_j;


              compute_dX(dX, param, particles_position, is_particle_particle_collision, compo_sommets, sommets_facets, mesh);
              if(is_particle_particle_collision)
                {
                  compute_norm(norm, param, sommets_facets, mesh);
                }
              else
                {
                  int ind_wall = particle_j - nb_particles_tot_;
                  int ori = ind_wall < dimension ? ind_wall : ind_wall - dimension;
                  norm(ori) = (ind_wall < 3 ? 1 : -1); //La normale correspond à la normale à la paroi (seuls les parallélépipèdes sont considérés)
                }


              if(is_particle_particle_collision)
                {
                  dist_between_particles = -local_prodscal(dX,norm); //normal penetration distance
                }
              else
                {
                  double dist_gravity_center = sqrt(local_carre_norme_vect(dX));//project on normal ?? XXX
                  dist_between_particles = dist_gravity_center - activation_distance_ ;
                }
            }
          // Check if the current proc is the one that has to compute the force
          DoubleTab dist(1); //can use mp_min_for_each_item only with TRUSTArray
          dist(0) = dist_between_particles;
          mp_min_for_each_item(dist);
          // std::cout<<"["<<Process::me()<<"] : "<<"min  : "<<dist(0)<<" my distance : "<<dist_between_particles<<std::endl;
          if(dist(0)!=dist_between_particles) continue;


          std::cout<<"["<<Process::me()<<"] : "<<"finished compute dx du normal "<<is_particle_particle_collision<<std::endl;
          std::ofstream ffn, fft, fm, fu, fd;
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
          std::ofstream f_cp;
          // XXX debug files
          // std::ofstream f;
          // path = fichier_debug + "_closest_" + std::to_string(particle_i)+"_P"+std::to_string(Process::me())+".txt";
          // f_cp.open(path, std::ios::app);
          // if(particle_j - nb_particles_tot_==1) //ground
          //   f_cp<<t<<" "<<collision_point(0)<<" "<<collision_point(1)<<" "<<collision_point(2)<<"\n";
          // f_cp.close();



          F_now_(particle_i, particle_j) = 0;
          double max_dist = 0.;
          if (dist_between_particles <= 0) // contact
            {

              compute_dU_cp(dU, collision_point, param, is_particle_particle_collision, compo_sommets, part_prop, mesh);
              std::ofstream fcol;
              path = fichier_debug + "is_collision_" + std::to_string(particle_i)+"_P"+std::to_string(Process::me())+".txt";
              fcol.open(path, std::ios::app);
              fcol<<t<<" ["<<Process::me()<<"] : "<<"Collision detected !!!!"<<std::endl;
              max_dist = std::max(max_dist, -dist_between_particles);

              if(is_particle_particle_collision) {dist_between_particles*=-1;}
              add_collision(particle_i,particle_j,is_particle_particle_collision);

              // double dX_scal_dU = local_prodscal(dX,dU) / (dist_gravity_center>0 ? dist_gravity_center : 1);
              double dU_scal_norm = local_prodscal(dU,norm);
              DoubleTab dUn(dimension);
              for (int d = 0; d < dimension; d++)
                dUn(d) = dU_scal_norm * norm(d);

              const double impact_velocity = sqrt(local_carre_norme_vect(dUn));

              F_now_(particle_i, particle_j) = 1;
              int is_start_of_collision = F_now_(particle_i, particle_j) >
                                          F_old_(particle_i, particle_j); // We need to know
              // if this is the first time step of the collision to compute the impact velocity

              // DoubleTab next_dX(dimension);
              // for (int d = 0; d < dimension; d++)
              //   next_dX(d) = dX(d) + deltat_simu * dU(d);
              const double effective_radius = is_particle_particle_collision ? solid_particle.get_equivalent_radius()/2 :
                                              solid_particle.get_equivalent_radius();
              const double impact_Stokes = solid_density * 2 * effective_radius * impact_velocity /
                                           (9 * fluid_viscosity);
              if (is_start_of_collision)
                e_eff_(particle_i,particle_j)=e_dry*compute_ewet_legendre(impact_Stokes);
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
              for(int d=0; d<dimension; ++d)
                {
                  dUt(d) = dU(d) - dUn(d);
                }
              double dUt_norm = 0.;
              dUt_norm = sqrt(local_carre_norme_vect(dUt));
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
          if(particle_j - nb_particles_tot_==1) //ground
            fd<<"Proc "<<Process::me()<<" "<<t<<" "<<max_dist<<"\n";

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
  // if (detection_method_==Detection_method::LC_VERLET)
  //   {
  mp_sum_for_each_item(lagrangian_contact_forces_);
  mp_sum_for_each_item(lagrangian_contact_moments_);
  mp_max_for_each_item(F_old_);
  mp_max_for_each_item(F_now_);
  mp_max_for_each_item(e_eff_);
  mp_sum_for_each_item(particles_collision_number_);
  collision_number_=Process::check_int_overflow(Process::mp_sum(collision_number_));

  // output_orientation(particles_position,compo_sommets, mesh,particles_position.dimension(0), t);
  // }
  // t+=deltat_simu;
}

void Collision_Model_FT_ellipsoid::output_orientation(const DoubleTab& particles_position, IntLists const& compo_sommets, Maillage_FT_Disc const& mesh, int nb_compo, double t)
{
  auto const& sommets = mesh.sommets();
  for(int compo = 0; compo<nb_compo; ++compo)
    {
      std::cout<<"compo : "<<compo<<std::endl;
      static int further1=-1, further2=-1;
      double d_further1=std::numeric_limits<double>::min(), d_further2=std::numeric_limits<double>::min();
      static int closest1=-1, closest2=-1;
      double d_closest1=std::numeric_limits<double>::max(), d_closest2=std::numeric_limits<double>::max();
      DoubleTab X(dimension);
      DoubleTab xs(dimension);
      DoubleTab xf1(dimension), xf2(dimension), xc1(dimension), xc2(dimension);
      for(int d=0; d<dimension; ++d) {X(d) = 0.*particles_position(compo, d);}
      DoubleTab coord_f1(dimension +1);
      DoubleTab coord_f2(dimension +1);
      DoubleTab coord_c1(dimension +1);
      DoubleTab coord_c2(dimension +1);

      further1 = 0;
      closest1 = 0;
      for(int d=0; d<dimension; ++d) {xf1(d) = sommets(0, d) - X(d); xc1(d) = sommets(0, d) - X(d);}
      d_further1 = sqrt(local_prodscal(xs, xs));
      d_closest1 = sqrt(local_prodscal(xs, xs));
      for (int i = 1; i < sommets.dimension(0); ++i)
        // for (int i = 0; i < compo_sommets[compo].size(); ++i)
        {
          // int i_global = compo_sommets[compo][i];
          int i_global = i;
          for(int d=0; d<dimension; ++d) {xs(d) = sommets(i_global, d) - X(d);}

          double dist = sqrt(local_prodscal(xs, xs));
          if(dist>d_further2 && local_prodscal(xs,xf1) >= 0)
            {
              further2 = i_global;
              d_further2=dist;
              for(int d=0; d<dimension; ++d) {xf2(d) = xs(d);}
            }
          if (dist>d_further1 && i_global != further2)
            {
              further1 = i_global;
              d_further1=dist;
              for(int d=0; d<dimension; ++d) {xf1(d) = xs(d);}
            }
          else if (dist<d_closest2 && local_prodscal(xs,xc1) >= 0)
            {
              closest2 = i_global;
              d_closest2=dist;
              for(int d=0; d<dimension; ++d) {xc2(d) = xs(d);}
            }
          else if (dist<d_closest1 /*&& local_prodscal(xs,xc2) <= 0*/)
            {
              closest1 = i_global;
              d_closest1=dist;
              for(int d=0; d<dimension; ++d) {xc1(d) = xs(d);}
            }
        }
      std::cout<<"d closest 1 : "<<d_closest1<<std::endl;
      std::cout<<"d closest 2 : "<<d_closest2<<std::endl;
      coord_f1(0) = d_further1;
      coord_f2(0) = d_further2;
      coord_c1(0) = d_closest1;
      coord_c2(0) = d_closest2;

      for(int d=0; d<dimension; ++d)
        {
          coord_f1(d+1) = sommets(further1, d);
          coord_f2(d+1) = sommets(further2, d);
          coord_c1(d+1) = sommets(closest1, d);
          coord_c2(d+1) = sommets(closest2, d);
        }
      if(Process::is_parallel())
        {
          if(Process::me()!=0)
            {
              envoyer(coord_f1, Process::me(), 0, 100+1);
              envoyer(coord_f2, Process::me(), 0, 100+2);
              envoyer(coord_c1, Process::me(), 0, 100+3);
              envoyer(coord_c2, Process::me(), 0, 100+4);
            }
          else
            {
              int nb_proc = Process::nproc();
              DoubleTab coord_f1_list(nb_proc, dimension);
              DoubleTab coord_f2_list(nb_proc, dimension);
              DoubleTab coord_c1_list(nb_proc, dimension);
              DoubleTab coord_c2_list(nb_proc, dimension);
              double global_d_f1 = std::numeric_limits<double>::min();
              double global_d_f2 = std::numeric_limits<double>::min();
              double global_d_c1 = std::numeric_limits<double>::max();
              double global_d_c2 = std::numeric_limits<double>::max();
              int max_f1=-1, max_f2=-1, min_c1=-1, min_c2=-1;
              for(int i=0; i<nb_proc; ++i)
                {
                  DoubleTab temp_f1(dimension+1), temp_f2(dimension+1), temp_c1(dimension+1), temp_c2(dimension+1);
                  recevoir(temp_f1, i, 0, 100+1);
                  recevoir(temp_f2, i, 0, 100+2);
                  recevoir(temp_c1, i, 0, 100+3);
                  recevoir(temp_c2, i, 0, 100+4);

                  for(int d=0; d<dimension; ++d)
                    {
                      coord_f1_list(i,d) = temp_f1(d+1);
                      coord_f2_list(i,d) = temp_f2(d+1);
                      coord_c1_list(i,d) = temp_c1(d+1);
                      coord_c2_list(i,d) = temp_c2(d+1);
                    }
                  if(temp_f1(0)>global_d_f1) {max_f1 = i;}
                  if(temp_f2(0)>global_d_f2) {max_f2 = i;}
                  if(temp_c1(0)<global_d_c1) {min_c1 = i;}
                  if(temp_c2(0)<global_d_c2) {min_c2 = i;}

                }
              DoubleTab global_f1(dimension), global_f2(dimension), global_c1(dimension), global_c2(dimension);
              DoubleTab long_axe(dimension), small_axe(dimension);
              for(int d=0; d<dimension; ++d)
                {
                  global_f1(d) = coord_f1_list(max_f1, d);
                  global_f2(d) = coord_f2_list(max_f2, d);
                  long_axe(d) = global_f1(d) - global_f2(d);

                  global_c1(d) = coord_c1_list(min_c1, d);
                  global_c2(d) = coord_c2_list(min_c2, d);
                  small_axe(d) = global_c1(d) - global_c2(d);
                }

              double angle_long_axe = acos(long_axe(1)/local_prodscal(long_axe, long_axe));//angle wrt the y axis
              double angle_small_axe = acos(long_axe(0)/local_prodscal(small_axe, small_axe));//angle wrt the x axis
              std::string path;
              std::fstream f;
              path = fichier_debug + "orientation_points.txt";
              f.open(path, std::ios::app);
              f<<t<<" ";
              f.close();
              path = fichier_debug + "orientation.txt";
              f.open(path, std::ios::app);
              f<<t<<" "<<angle_long_axe<<" "<<angle_small_axe<<"\n";
            }
        }
      else
        {
          std::cout<<"got here"<<std::endl;
          DoubleTab global_f1(dimension), global_f2(dimension), global_c1(dimension), global_c2(dimension);
          DoubleTab long_axe(dimension), small_axe(dimension);
          for(int d=0; d<dimension; ++d)
            {
              global_f1(d) = coord_f1(d+1);
              global_f2(d) = coord_f2(d+1);
              long_axe(d) = global_f1(d) - global_f2(d);
              // long_axe(d) = xf1(d);

              global_c1(d) = coord_c1(d+1);
              global_c2(d) = coord_c2(d+1);
              small_axe(d) = global_c1(d) - global_c2(d);
              // small_axe(d) = xc1(d);
            }
          std::cout<<global_f1(0)<<" "<<global_f1(1)<<" "<<global_f1(2)<<std::endl;
          std::cout<<global_f2(0)<<" "<<global_f2(1)<<" "<<global_f2(2)<<std::endl;

          double angle_long_axe = acos(long_axe(1)/sqrt(local_prodscal(long_axe, long_axe)));//angle wrt the y axis
          std::cout<<"angle_long_axe "<<angle_long_axe<<std::endl;
          double angle_small_axe = acos(small_axe(3)/sqrt(local_prodscal(small_axe, small_axe)));//angle wrt the z axis
          std::cout<<"angle_small_axe "<<angle_small_axe<<std::endl;
          std::string path;
          std::fstream f;
          path = fichier_debug + "orientation_points.txt";
          f.open(path, std::ios::app);
          f<<t<<" ";
          f.close();
          path = fichier_debug + "orientation.txt";
          f.open(path, std::ios::app);
          f<<t<<" "<<angle_long_axe<<" "<<angle_small_axe<<"\n";
        }
    }

}
//XXX Inertia not a ref because .inverse() modify the object
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
      temp[d] = (r[(d+1)%3] * force[(d+2)%3] - r[(d+2)%3] * force[(d+1)%3]) -
                (Omega[(d+1)%3] * InertiaOmega[(d+2)%3] - Omega[(d+2)%3] * InertiaOmega[(d+1)%3]);
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


          contact_force_source_term(face)=(1-volumic_phase_indicator_function(face))
                                          *interlaced_volumes(face)*(lagrangian_contact_forces_(id_number,ori)
                                                                     + (lagrangian_contact_moments_(id_number,(ori+1)%3) * (cg_faces(face,(ori+2)%3) - particles_position(id_number,(ori+2)%3)) -
                                                                        lagrangian_contact_moments_(id_number,(ori+2)%3) * (cg_faces(face,(ori+1)%3) - particles_position(id_number,(ori+1)%3) )) +
                                                                     /*collision **/ part_prop.density * ( Omega_i((ori+1)%3) * Omega_r((ori+2)%3) - Omega_i((ori+2)%3) * Omega_r((ori+1)%3) +
                                                                                                           inertia(ori)));
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