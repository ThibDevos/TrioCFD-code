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

void Collision_Model_FT_ellipsoid::deepest_points(IntLists const& compo_sommets, IntLists const& sommets_facets, Maillage_FT_Disc const& mesh, int compo_i, int compo_j,
                                                  DoubleTab const& positions, bool check_cg, DoubleTab& dX, int& i_closest, int& j_closest, bool& i_fa7, bool& j_fa7)
{
  const DoubleTab& sommets = mesh.sommets();
  DoubleTab dX_min(dimension);
  double dist_cg = 0.;
  int i_som_closest = -1;
  int j_som_closest = -1;

  for (int d = 0; d < dimension; ++d)
    dX_min(d) = 1000.;
  double dX_min_norm_i = sqrt(local_carre_norme_vect(dX_min));

  for (int i_som = 0; i_som < compo_sommets[compo_i].size(); ++i_som)
    {
      int i_global = compo_sommets[compo_i][i_som];
      for (int d = 0; d < dimension; ++d)
        {
          dX(d) = sommets(i_global, d) - positions(compo_j, d);
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

  for (int j_som = 0; j_som < compo_sommets[compo_j].size(); ++j_som)
    {
      int j_global = compo_sommets[compo_j][j_som];
      for (int d = 0; d < dimension; ++d)
        {
          dX(d) = sommets(j_global, d) - positions(compo_i, d);
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
              dX(d) = cg_fa7(fa7, d) - positions(compo_j, d);
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
              dX(d) = cg_fa7(fa7, d) - positions(compo_i, d);
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
      i_fa7 = false;
      i_closest = i_som_closest;
      for (int d = 0; d < dimension; ++d)
        {
          deepest_i(d) = sommets(i_som_closest, d);
        }
    }
  else
    {
      i_fa7 = true;
      i_closest = fi_closest;
      for (int d = 0; d < dimension; ++d)
        {
          deepest_i(d) = cg_fa7(fi_closest, d);
        }
    }
  if (fj_closest < 0)
    {
      j_fa7 = false;
      j_closest = j_som_closest;
      for (int d = 0; d < dimension; ++d)
        {
          deepest_j(d) = sommets(j_som_closest, d);
        }
    }
  else
    {
      j_fa7 = true;
      j_closest = fj_closest;
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

void Collision_Model_FT_ellipsoid::closest_nodes(IntLists const& compo_sommets, Maillage_FT_Disc const& mesh, int particle, int neighbor,
                                                 DoubleTab const& positions, DoubleTab& dX, int& i_closest, int& j_closest, bool chek_cg, bool i_facet, bool j_facet)
{
  auto sommets = mesh.sommets();
  DoubleTab dX_min(dimension);
  double dist_cg = 0.;
  i_closest = -1;
  j_closest = -1;
  for (int d = 0; d < dimension; ++d)
    dX_min(d) = 1000.;

  if (chek_cg)
    {
    }
  else
    {
      double dX_min_norm = sqrt(local_carre_norme_vect(dX_min));
      for (int i = 0; i < compo_sommets[particle].size(); ++i)
        {
          int i_global = compo_sommets[particle][i];
          for (int j = 0; j < compo_sommets[neighbor].size(); ++j)
            {
              int j_global = compo_sommets[neighbor][j];
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
                  i_closest = i_global;
                  j_closest = j_global;
                }
            }
        }
    }
  for (int d = 0; d < dimension; ++d)
    {
      dX(d) = -dX_min(d);
    }
  assert(i_closest >= 0);
  assert(j_closest >= 0);
}

void Collision_Model_FT_ellipsoid::normal_i(Maillage_FT_Disc const& mesh, IntLists const& sommets_facets, bool i_facet, int i_closest, DoubleTab& n)
{
  auto nfacets = mesh.get_normale_facettes();
  if(i_facet)
    {
      for(int d = 0; d<dimension; ++d) {n(d) = nfacets(i_closest,d);}
    }
  else //the deepest point is the node -> the normal at the node is the average of the normal of the incident facets
    {
      for(int fi = 0; fi<sommets_facets[i_closest].size(); ++fi)
        {
          int fa7 = sommets_facets[i_closest][fi];
          for(int d = 0; d<dimension; ++d) {n(d) += nfacets(fa7,d);}
        }
    }
  double nn = sqrt(local_carre_norme_vect(n));
  for(int d = 0; d<dimension; ++d) {n(d)/=nn;}
}

void Collision_Model_FT_ellipsoid::normal_average_ij(Maillage_FT_Disc const& mesh,IntLists const& sommets_facets, bool i_facet, int i_closest, bool j_facet, int j_closest, DoubleTab& n)
{
  auto nfacets = mesh.get_normale_facettes();
  DoubleTab ni(dimension);
  DoubleTab nj(dimension);

  if(i_facet)
    {
      std::cout<<"a\n";
      for(int d = 0; d<dimension; ++d) {ni(d) = nfacets(i_closest,d);}
      std::cout<<"aa\n";
    }
  else //the deepest point is the node -> the normal at the node is the average of the normal of the incident facets
    {
      std::cout<<"b\n";
      for(int fi = 0; fi<sommets_facets[i_closest].size(); ++fi)
        {
          std::cout<<fi<<"\n";
          int fa7 = sommets_facets[i_closest][fi];
          std::cout<<fi<<"\n";
          for(int d = 0; d<dimension; ++d) {ni(d) += nfacets(fa7,d);}
        }
      std::cout<<"bb\n";
    }
  if(j_facet)
    {
      std::cout<<"c\n";
      for(int d = 0; d<dimension; ++d) {nj(d) = nfacets(j_closest,d);}
      std::cout<<"cc\n";
    }
  else //the deepest point is the node -> the normal at the node is the average of the normal of the incident facets
    {
      std::cout<<"d\n";
      for(int fi = 0; fi<sommets_facets[j_closest].size(); ++fi)
        {
          int fa7 = sommets_facets[j_closest][fi];
          for(int d = 0; d<dimension; ++d) {nj(d) += nfacets(fa7,d);}
        }
      std::cout<<"dd\n";
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
  std::cout<<"leave normal_ij\n";
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

void Collision_Model_FT_ellipsoid::compute_dX_dU_normal(DoubleTab& dX, DoubleTab& dU, DoubleTab& norm, DoubleTab& cp, int particle,
                                                        int neighbor, const DoubleTab& particles_position, const
                                                        DoubleTab& particles_velocity, const bool is_particle_particle_collision,
                                                        const IntLists& compo_sommets,const IntLists& sommets_facets, const Maillage_FT_Disc& mesh )
{
  std::cout<<"==========DX_DU=============="<<is_particle_particle_collision<<"\n";
  const DoubleTab& sommets = mesh.sommets();
  if (is_particle_particle_collision)
    {
      bool check_cg = true;
      int i_closest, j_closest;
      bool i_facet, j_facet;
      if(collision_detection_==Collision_detection::CLOSEST)
        {
          std::cout<<"================closest\n";
          closest_nodes(compo_sommets, mesh, particle, neighbor, particles_position ,dX, i_closest, j_closest, false, i_facet, j_facet); //check_cg not implemented
        }
      if(collision_detection_==Collision_detection::DEEPEST)
        {
          std::cout<<"================deppest\n";
          deepest_points(compo_sommets,sommets_facets,mesh,particle, neighbor, particles_position, check_cg, dX, i_closest, j_closest, i_facet, j_facet);
        }
      if(collision_normal_==Collision_normal::NORMAL_I)
        {
          std::cout<<"================normal_i\n";
          normal_i(mesh, sommets_facets, i_facet, i_closest, norm);
        }
      if(collision_normal_==Collision_normal::NORMAL_IJ)
        {
          std::cout<<"================normal_ij\n";
          normal_average_ij(mesh, sommets_facets, i_facet, i_closest, j_facet, j_closest, norm);
        }
      for (int d = 0; d < dimension; d++)
        {
          dU(d) = particles_velocity(particle, d) - particles_velocity(neighbor, d);
        }
      collision_point(mesh, i_closest, j_closest, i_facet, j_facet, cp);
    }
  else
    {

      int ind_wall = neighbor - nb_particles_tot_;
      int ori = ind_wall < dimension ? ind_wall : ind_wall - dimension;
      double dX_min =  std::numeric_limits<double>::max();
      int i_som_closest=-1;
      for(int i_som = 0; i_som < compo_sommets[particle].size(); ++i_som) //recherche le sommet le plus proche du mur
        {
          dX(ori) = std::fabs(sommets(compo_sommets[particle][i_som],ori) -  (origin_(ind_wall) + (ind_wall>2)*domain_dimensions_(ori)));
          if(dX(ori) < dX_min ) i_som_closest = i_som;
          dX_min = dX(ori) <= dX_min ? dX(ori) : dX_min;
        }
      dX(ori) = dX_min;
      ofstream f;
      f.open("dist.txt", std::ios::app);
      if(ind_wall==1 || ind_wall == 4)
        f<<ind_wall<<" "<<dX(0)<<" "<<dX(1)<<" "<<dX(2)<<" "<<fictive_wall_coordinates_(ind_wall)<<"\n";
      for (int d = 0; d < dimension; d++)
        dU(d) = particles_velocity(particle, d);

      norm(ori) = (ind_wall < 3 ? 1 : -1); //La normale correspond à la normale à la paroi (seuls les parallélépipèdes sont considérés)

      for (int d = 0; d < dimension; d++)
        cp(d) = sommets(compo_sommets[particle][i_som_closest],d);
    }
}

//Only for 3x3 matrices
double compute_det_M33(const Matrice_Dense& M)
{
  return M(0,0) * M(1,1) * M(2,2) + M(1,0) * M(2,1) * M(0,2) + M(0,1) * M(1,2) * M(2,0)
         - ( M(2,0) * M(1,1) * M(0,2) + M(1,0) * M(0,1) * M(2,2) + M(0,0) * M(2,1) * M(1,2) );

}

void Collision_Model_FT_ellipsoid::compute_inertia_tensor(const Maillage_FT_Disc& mesh, int ind_particle_i,
                                                          const DoubleTab& particles_position, Matrice_Dense& J)
{


  const DoubleTab& sommets = mesh.sommets();
  const IntTab& facets = mesh.facettes();
  IntLists compo_connexe_facets;
  connec_compo_facettes(mesh, compo_connexe_facets);
  int nb_facets = compo_connexe_facets[ind_particle_i].size();


  J.clean();

  DoubleTab b_K(dimension);

  //Prepare Gauss quadrature
  // double a = (5.-std::sqrt(5.))/20.;
  DoubleTab points(5,3);
  DoubleTab poids(5);
  // points(0,0) = a;
  // points(0,1) = a;
  // points(0,2) = 1-3*a;

  // points(1,0) = a;
  // points(1,1) = a;
  // points(1,2) = a;

  // points(2,0) = 1-3*a;
  // points(2,1) = a;
  // points(2,2) = a;

  // points(3,0) = a;
  // points(3,1) = 1-3*a;
  // points(3,2) = a;

  points(0,0) = 0.25;
  points(0,1) = 0.25;
  points(0,2) = 0.25;
  poids(0) = -4./30.;

  points(1,0) = 1./6.;
  points(1,1) = 1./6.;
  points(1,2) = 1./6.;
  poids(1) = 9./120.;

  points(2,0) = 1./2.;
  points(2,1) = 1./6.;
  points(2,2) = 1./6.;
  poids(2) = 9./120.;

  points(3,0) = 1./6.;
  points(3,1) = 1./2.;
  points(3,2) = 1./6.;
  poids(3) = 9./120.;

  points(4,0) = 1./6.;
  points(4,1) = 1./6.;
  points(4,2) = 1./2.;
  poids(4) = 9./120.;


  for(int d=0; d<dimension; ++d) {b_K(d) = particles_position(ind_particle_i,d);}
  std::cout<<"position "<<b_K(0)<<" "<<b_K(1)<<" "<<b_K(2)<<"\n";
  for(int f = 0; f<nb_facets; ++f)
    {
      int f_global = compo_connexe_facets[ind_particle_i][f];
      Matrice_Dense T_K(dimension, dimension);
      Matrice_Dense J_loc(dimension, dimension);
      J_loc.clean();
      for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
            {
              T_K(i, j) = sommets(facets(f_global, j), i) - b_K(i);
            }
        }
      double det_T_K = std::fabs(compute_det_M33(T_K));
      // det_T_K /= det_T_K; //XXX
      DoubleTab x_q(dimension);
      DoubleTab point(dimension);
      for(int q=0; q<5; ++q)
        {
          for(int d=0; d<dimension; ++d) {point(d) = points(q,d); x_q(d)=b_K(d);}
          T_K.ajouter_multvect_(point,x_q);
          // x_q += b_K;

          for(int i=0; i<dimension; ++i)
            {
              J_loc(i,i) += det_T_K * ((x_q((i+1)%3) - b_K((i+1)%3) ) * (x_q((i+1)%3) - b_K((i+1)%3)) + (x_q((i+2)%3) - b_K((i+2)%3)) * (x_q((i+2)%3) - b_K((i+2)%3))) * poids(q);
              for(int j=i+1; j<dimension; ++j)
                {
                  J_loc(i,j) += -det_T_K * (x_q(i) - b_K(i)) * (x_q(j)- b_K(j))* poids(q);
                  J_loc(j,i) = J_loc(i,j);
                }
            }
        }
      //J = J + J_loc;
      for(int i=0; i<dimension; ++i)
        {
          for(int j=0; j<dimension; ++j)
            {
              J(i,j) += J_loc(i,j);
            }
        }
    }
}

void Collision_Model_FT_ellipsoid::compute_lagrangian_contact_forces(const Fluide_Diphasique& two_phase_fluid,
                                                                     const DoubleTab& particles_position,
                                                                     const DoubleTab& particles_velocity,
                                                                     const DoubleTab& particles_rot_velocity,
                                                                     const double& deltat_simu,
                                                                     const Maillage_FT_Disc& mesh)
{
  static double t=0.;
  const int& id_fluid_phase= two_phase_fluid.get_id_fluid_phase();
  const int& id_solid_phase=1-id_fluid_phase;
  const auto& solid_particle=ref_cast(Solid_Particle_ellipsoid,two_phase_fluid.fluide_phase(id_solid_phase));
  const auto& incompressible_fluid=ref_cast(Fluide_Incompressible,
                                            two_phase_fluid.fluide_phase(id_fluid_phase));
  const double& solid_density = solid_particle.masse_volumique().valeurs()(0, 0);
  const double& fluid_density = incompressible_fluid.masse_volumique().valeurs()(0, 0);
  const double& fluid_viscosity  = fluid_density
                                   * incompressible_fluid.viscosite_cinematique().valeurs()(0, 0);
  const double& long_radius=solid_particle.get_max_radius();
  const double& volume=solid_particle.get_volume();
  const double& e_dry=solid_particle.get_e_dry();
  const double min_threshold=1e-10;
  DoubleTab dX(dimension), dU(dimension), norm(dimension), collision_point(dimension);
  lagrangian_contact_forces_=0;
  lagrangian_contact_moments_=0;
  collision_number_=0;
  particles_collision_number_=0;


  IntLists compo_sommets; //compo_sommet[i] contient les indices des sommets composant la compo i
  connec_compo_sommets(mesh, compo_sommets);
  IntLists sommets_facets; //compo_sommet[i] contient les indices des sommets composant la compo i
  connec_sommets_fa7(mesh, sommets_facets);

  for (int ind_particle_i = 0; ind_particle_i < nb_real_particles_; ind_particle_i++)
    {
      int particle_i=get_particle_i(ind_particle_i);
      int nb_particles_j=get_nb_particles_j(ind_particle_i);
      int ind_start_part_j=get_ind_start_particles_j(ind_particle_i);
      for (int ind_particle_j =ind_start_part_j; ind_particle_j < nb_particles_j; ind_particle_j++)
        {
          dX = 0;
          dU = 0;
          norm = 0;
          int particle_j=get_particle_j(ind_particle_i,ind_particle_j);
          int is_particle_particle_collision = particle_j < nb_particles_tot_;
          compute_dX_dU_normal(dX, dU, norm, collision_point, particle_i, particle_j, particles_position,
                               particles_velocity, is_particle_particle_collision, compo_sommets, sommets_facets, mesh);
          double dist_gravity_center = sqrt(local_carre_norme_vect(dX));
          double dist_between_particles = 0.;
          if(is_particle_particle_collision)
            {
              dist_between_particles = local_prodscal(dX,norm);
            }
          else
            {
              dist_between_particles = dist_gravity_center - activation_distance_ ;
            }
          F_now_(particle_i, particle_j) = 0;
          Matrice_Dense J(dimension, dimension);

          if (dist_between_particles <= 0) // contact
            {
              if(is_particle_particle_collision) {dist_between_particles*=-1;}
              add_collision(particle_i,particle_j,is_particle_particle_collision);


              double dX_scal_dU = local_prodscal(dX,dU) / (dist_gravity_center>0 ? dist_gravity_center : 1);
              DoubleTab dUn(dimension);
              for (int d = 0; d < dimension; d++)
                dUn(d) = dX_scal_dU * norm(d);

              const double impact_velocity = sqrt(local_carre_norme_vect(dUn));

              F_now_(particle_i, particle_j) = 1;
              int is_start_of_collision = F_now_(particle_i, particle_j) >
                                          F_old_(particle_i, particle_j); // We need to know
              // if this is the first time step of the collision to compute the impact velocity

              DoubleTab next_dX(dimension);
              for (int d = 0; d < dimension; d++)
                next_dX(d) = dX(d) + deltat_simu * dU(d);
              const double effective_radius = is_particle_particle_collision ? long_radius/2 :
                                              long_radius;
              const double impact_Stokes = solid_density * 2 * effective_radius * impact_velocity /
                                           (9 * fluid_viscosity);
              if (is_start_of_collision)
                e_eff_(particle_i,particle_j)=e_dry*compute_ewet_legendre(impact_Stokes)/compute_ewet_legendre(impact_Stokes);
              DoubleTab force_contact=compute_contact_force(
                                        dist_between_particles,
                                        norm,
                                        dUn,
                                        particle_i,
                                        particle_j,
                                        dX_scal_dU<=0,
                                        is_particle_particle_collision);

              // XXX test with real tensor for sphere
              for(int i=0; i<dimension; ++i)
                {
                  J(i,i) = 0.4*(1.5e-3)*1.5e-3;
                  for(int j=i+1; j<dimension; ++j)
                    {
                      J(i,j)=0.;
                      J(j,i)=0.;
                    }
                }
              DoubleTab r(dimension);
              r(0) = particles_position(particle_i,0) - collision_point(0);
              r(1) = particles_position(particle_i,1) - collision_point(1);
              r(2) = particles_position(particle_i,2) - collision_point(2);
              DoubleTab moment_contact=compute_contact_moment(J,force_contact, r, particles_rot_velocity);
              std::cout<<force_contact(0)<<" "<<force_contact(1)<<" "<<force_contact(2)<<"\n";
              std::cout<<moment_contact(0)<<" "<<moment_contact(1)<<" "<<moment_contact(2)<<"\n";
              // exit();
              for (int d = 0; d < dimension; d++)
                {
                  lagrangian_contact_forces_(particle_i, d) += fabs(force_contact(d)) <=
                                                               min_threshold ? 0 : force_contact(d) / volume;
                  lagrangian_contact_moments_(particle_i, d) += fabs(moment_contact(d)) <=
                                                                min_threshold ? 0 : moment_contact(d) / volume;
                  if (!is_particle_particle_collision)
                    continue; // wall collision, no force to apply on the wall
                  lagrangian_contact_forces_(particle_j, d) -= fabs(force_contact(d)) <=
                                                               min_threshold ? 0 :  force_contact(d) / volume;
                  lagrangian_contact_moments_(particle_j, d) -= fabs(moment_contact(d)) <=
                                                                min_threshold ? 0 : moment_contact(d) / volume;
                }

              F_old_(particle_i, particle_j) = F_now_(particle_i, particle_j);
            }
        }
    }

  if (detection_method_==Detection_method::LC_VERLET)
    {
      mp_sum_for_each_item(lagrangian_contact_forces_);
      mp_sum_for_each_item(lagrangian_contact_moments_);
      mp_max_for_each_item(F_old_);
      mp_max_for_each_item(F_now_);
      mp_max_for_each_item(e_eff_);
      mp_sum_for_each_item(particles_collision_number_);
      collision_number_=Process::check_int_overflow(Process::mp_sum(collision_number_));
    }
  t+=deltat_simu;
}

DoubleTab Collision_Model_FT_ellipsoid::compute_contact_moment(Matrice_Dense& Inertia, DoubleTab const& force, DoubleTab const& r, DoubleTab const& Omega)
{
  DoubleTab contact_moment(dimension);
  DoubleTab temp(dimension);
  DoubleTab InertiaOmega(dimension);
  Inertia.ajouter_multvect_(Omega,InertiaOmega);
  std::cout<<Inertia(0,0)<<" "<<Inertia(0,1)<<" "<<Inertia(0,2)<<"\n";
  std::cout<<Inertia(1,0)<<" "<<Inertia(1,1)<<" "<<Inertia(1,2)<<"\n";
  std::cout<<Inertia(2,0)<<" "<<Inertia(2,1)<<" "<<Inertia(2,2)<<"\n";
  std::cout<<"=========================\n";
  Inertia.inverse(); //Inertia is now its inverse
  for(int d = 0; d<dimension; ++d)
    {
      temp[d] = (r[(d+1)%3] * force[(d+2)%3] - r[(d+2)%3] * force[(d+1)%3]) -
                (Omega[(d+1)%3] * InertiaOmega[(d+2)%3] - Omega[(d+2)%3] * InertiaOmega[(d+1)%3]);
    }
  std::cout<<Inertia(0,0)<<" "<<Inertia(0,1)<<" "<<Inertia(0,2)<<"\n";
  std::cout<<Inertia(1,0)<<" "<<Inertia(1,1)<<" "<<Inertia(1,2)<<"\n";
  std::cout<<Inertia(2,0)<<" "<<Inertia(2,1)<<" "<<Inertia(2,2)<<"\n";
  std::cout<<"=========================\n";
  std::cout<<r(0)<<" "<<r(1)<<" "<<r(2)<<"\n";
  std::cout<<InertiaOmega(0)<<" "<<InertiaOmega(1)<<" "<<InertiaOmega(2)<<"\n";
  std::cout<<Omega(0)<<" "<<Omega(1)<<" "<<Omega(2)<<"\n";
  std::cout<<temp(0)<<" "<<temp(1)<<" "<<temp(2)<<"\n";
  Inertia.ajouter_multvect_(temp,contact_moment);
  std::cout<<contact_moment(0)<<" "<<contact_moment(1)<<" "<<contact_moment(2)<<"\n";
  std::cout<<"---------------------------\n";
  return contact_moment;
}

void Collision_Model_FT_ellipsoid::discretize_contact_forces_eulerian_field(
  const DoubleTab& volumic_phase_indicator_function,
  const Domaine_VF& domain_vf,
  const IntTab& particles_eulerian_id_number,
  const DoubleTab& particles_position,
  DoubleTab& contact_force_source_term)
{
  const DoubleVect& interlaced_volumes=domain_vf.volumes_entrelaces();
  const DoubleTab& cg_faces=domain_vf.xp();
  const int nb_faces=interlaced_volumes.size_array();
  const IntVect& orientation = domain_vf.orientation();
  const IntTab& face_voisins=domain_vf.face_voisins();
  for (int face=0; face<nb_faces; face++)
    {
      const int left_elem=face_voisins(face,0);
      const int right_elem=face_voisins(face,1);
      int id_left = left_elem != -1 ? particles_eulerian_id_number(left_elem) : -1;
      int id_right = right_elem != -1 ? particles_eulerian_id_number(right_elem) : -1;
      const int id_number=std::max(id_left,id_right);
      if (id_number!=-1)
        {
          const int ori=orientation(face);
          contact_force_source_term(face)=(1-volumic_phase_indicator_function(face))
                                          *interlaced_volumes(face)*(lagrangian_contact_forces_(id_number,ori)
                                                                     + (lagrangian_contact_moments_(id_number,(ori+1)%3) * (particles_position(id_number,(ori+2)%3) - cg_faces(face,(ori+2)%3)) -
                                                                        lagrangian_contact_moments_(id_number,(ori+2)%3) * (particles_position(id_number,(ori+1)%3) - cg_faces(face,(ori+1)%3)))) ;
        }
    }
}

