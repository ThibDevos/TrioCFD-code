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

#ifndef Modele_Collision_FT_ellipsoid_included
#define Modele_Collision_FT_ellipsoid_included

#include <Collision_Model_FT_base.h>
#include <Connex_components_FT.h>
#include <Matrice_Dense.h>
#include <fstream>

/*! @brief : class Collision_Model_FT
 *
 *  Description: This class enables to compute solid-solid
 *  interactions for fpi module under the framework of
 *  soft-sphere collision model. Under this framework,
 *  multiple collisions can occurs at the same time (ie a
 *  particle can collide with 2 or more particles). The
 *  collision is spread out on multiple time steps. A slight
 *  overlap (less than the mesh grid size) occurs during the
 *  process.
 */

class Collision_Model_FT_ellipsoid : public Collision_Model_FT_base
{
  Declare_instanciable_sans_constructeur(Collision_Model_FT_ellipsoid);

public:
  void test_connex_compo(Maillage_FT_Disc const& maillage);
  Collision_Model_FT_ellipsoid();
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int preparer_calcul(const Domaine_VDF& domain_vdf,
                      const int nb_particles_tot,
                      const Navier_Stokes_FT_Disc& ns,
                      const Transport_Interfaces_FT_Disc& eq_transport,
                      const Schema_Comm& schema_comm_FT);
  void compute_fictive_wall_coordinates(const double& radius);
  void  collision_point(Maillage_FT_Disc const& maillage, int i_closest, int j_closest, bool i_fa7, bool j_fa7, DoubleTab& xc);
  void deepest_points(IntLists const& compo_sommets, IntLists const& sommets_facets, Maillage_FT_Disc const& mesh, int particle, int neighbor,
                      DoubleTab const& positions, bool check_cg, DoubleTab& dX, int& i_closest, int& j_closest, bool& i_facet, bool& j_facet);
  void closest_nodes(IntLists const& compo_sommets, Maillage_FT_Disc const& mesh, int particle, int neighbor,
                     DoubleTab const& positions, DoubleTab& dX, int& i_closest, int& j_closest, bool chek_cg, bool i_facet, bool j_facet);
  void normal_i(Maillage_FT_Disc const& mesh, IntLists const& sommets_facets, bool i_facet, int i_closest, DoubleTab& n);
  void normal_average_ij(Maillage_FT_Disc const& mesh,IntLists const& sommets_facets, bool i_facet, int i_closest, bool j_facet, int j_closest, DoubleTab& n);
  void compute_dX_dU_normal(DoubleTab& dX, DoubleTab& dU, DoubleTab& norm, DoubleTab& cp, const int particle,
                            const int neighbor, const DoubleTab& particles_position, const DoubleTab& particles_velocity, const DoubleTab& particles_rot_velocity, const bool is_particle_particle_collision,
                            const IntLists& compo_sommets, const IntLists& sommets_facets, const Maillage_FT_Disc& mesh);



  DoubleTab compute_contact_moment(Matrice_Dense Inertia, DoubleTab const& force, DoubleTab const& contact_point, DoubleTab const& Omega);
  virtual void compute_lagrangian_contact_forces(const Fluide_Diphasique& two_phase_fluid,
                                                 const particle_properties& part_prop,
                                                 const double& deltat_simu,
                                                 const Maillage_FT_Disc& mesh) override;
  void compute_lagrangian_contact_forces(const Fluide_Diphasique& two_phase_fluid,
                                         const DoubleTab& particles_position,
                                         const DoubleTab& particles_velocity,
                                         const DoubleTab& particles_rot_velocity,
                                         const double& deltat_simu,
                                         const Maillage_FT_Disc& mesh) override;

  void discretize_contact_forces_eulerian_field(const DoubleTab& volumic_phase_indicator_function,
                                                const Domaine_VF& domain_vf,
                                                const IntTab& particles_eulerian_id_number,
                                                const DoubleTab& particles_position,
                                                DoubleTab& contact_force_source_term) override;


private:
  enum class Collision_detection {DEEPEST, CLOSEST};
  Collision_detection collision_detection_ = Collision_detection::DEEPEST;
  enum class Collision_normal {NORMAL_I, NORMAL_IJ};
  Collision_normal collision_normal_ = Collision_normal::NORMAL_IJ;

  // XXX debug
  Motcle fichier_debug;
};

#endif

