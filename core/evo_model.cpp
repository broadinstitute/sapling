#include "evo_model.h"

#include "absl/log/check.h"

namespace sapling {

auto make_single_partition_global_evo_model(Site_index num_sites) -> Global_evo_model {
  return Global_evo_model{
    Site_vector<double>(num_sites, 1.0), // nu_l = 1.0 for all l
    Site_evo_model{}                     // Site_evo_model for the one and only partition
  };
}

}  // namespace sapling
