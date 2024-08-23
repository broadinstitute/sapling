#ifndef SAPLING_EVO_MODEL_H_
#define SAPLING_EVO_MODEL_H_

#include "sequence.h"

namespace sapling {

struct Site_evo_model {
  double mu{};
  Seq_vector<double> pi_a{};
  Seq_matrix<double> q_ab{};

  auto q_a(Real_seq_letter a) const -> double { return -q_ab[a][a]; }
};

struct Global_evo_model {
  Site_vector<double> nu_l;
  Site_evo_model site_evo_model;

  auto pi_a(Real_seq_letter a) const -> double {
    return site_evo_model.pi_a[a];
  }

  auto Q_l_a(Site_index l, Real_seq_letter a) const -> double {
    return site_evo_model.mu * nu_l.at(l) * site_evo_model.q_a(a);
  }

  auto Q_l_ab(Site_index l, Real_seq_letter a, Real_seq_letter b) const -> double {
    return site_evo_model.mu * nu_l.at(l) * site_evo_model.q_ab[a][b];
  }
};

auto make_single_partition_global_evo_model(Site_index num_sites) -> Global_evo_model;

}  // namespace sapling

#endif // SAPLING_EVO_MODEL_H_
