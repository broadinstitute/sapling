#ifndef SAPLING_COAL_SIM_H_
#define SAPLING_COAL_SIM_H_

#include <span>

#include "absl/random/bit_gen_ref.h"

#include "phylo_tree.h"
#include "pop_model.h"

namespace sapling {

auto coal_sim(const Pop_model& pop_model, std::span<const double> tip_times, absl::BitGenRef rng) -> Phylo_tree;

}  // namespace sapling

#endif // SAPLING_COAL_SIM_H_
