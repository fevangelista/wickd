#include <algorithm>

#define FMT_HEADER_ONLY
#include <fmt/core.h> // for fmt::format

#include "helpers/combinatorics.h"
#include "helpers/helpers.h"
#include "helpers/stl_utils.hpp"

#include "contraction.h"
#include "operator.h"
#include "operator_product.h"

#include "wick_theorem.h"

#define PRINT(detail, code)                                                    \
  if (print_ >= detail) {                                                      \
    code                                                                       \
  }

using namespace std;

namespace {

void generate_bounded_allocations_backtrack(
    const std::vector<int> &capacities,
    const std::vector<int> &suffix_capacities, int pos, int remaining,
    std::vector<int> &allocation, std::vector<std::vector<int>> &allocations) {
  if (remaining < 0 || remaining > suffix_capacities[pos]) {
    return;
  }
  if (pos == static_cast<int>(capacities.size())) {
    if (remaining == 0) {
      allocations.push_back(allocation);
    }
    return;
  }

  const int max_count = std::min(capacities[pos], remaining);
  for (int count = 0; count <= max_count; ++count) {
    allocation[pos] = count;
    generate_bounded_allocations_backtrack(capacities, suffix_capacities,
                                           pos + 1, remaining - count,
                                           allocation, allocations);
  }
  allocation[pos] = 0;
}

std::vector<std::vector<int>>
generate_bounded_allocations(int total, const std::vector<int> &capacities) {
  std::vector<int> suffix_capacities(capacities.size() + 1, 0);
  for (int pos = static_cast<int>(capacities.size()) - 1; pos >= 0; --pos) {
    suffix_capacities[pos] = suffix_capacities[pos + 1] + capacities[pos];
  }

  std::vector<std::vector<int>> allocations;
  std::vector<int> allocation(capacities.size(), 0);
  generate_bounded_allocations_backtrack(capacities, suffix_capacities, 0,
                                         total, allocation, allocations);
  return allocations;
}

} // namespace

std::vector<ElementaryContraction>
WickTheorem::generate_elementary_contractions(const OperatorProduct &ops) {
  PRINT(PrintLevel::Summary,
        std::cout << "\n- Step 1. Generating elementary contractions"
                  << std::endl;)

  int nops = ops.size();

  // a vector that will hold all the contractions
  std::vector<ElementaryContraction> contr_vec;

  // Fermionic general spaces that may participate in mixed-space cumulants.
  // The existing per-space generator remains responsible for contractions
  // confined to a single general space.
  std::vector<int> general_spaces;

  PRINT(
      PrintLevel::Summary, cout << "\n  Operator   Space   Cre.   Ann.";
      cout << "\n  ------------------------------";
      for (int op = 0; op < nops; ++op) {
        for (int s = 0; s < orbital_subspaces->num_spaces(); s++) {
          cout << "\n      " << op << "        " << orbital_subspaces->label(s)
               << "      " << ops[op].cre(s) << "      " << ops[op].ann(s);
        }
      };
      cout << "\n";)

  // loop over orbital spaces
  for (int s = 0; s < orbital_subspaces->num_spaces(); s++) {
    PRINT(PrintLevel::Summary, std::cout
                                   << "\n  Elementary contractions for space "
                                   << orbital_subspaces->label(s) << ": ";)

    // differentiate between various types of spaces
    SpaceType space_type = orbital_subspaces->space_type(s);

    // 1. Pairwise contractions 1 cre + 1 ann operator:
    // ┌───┐
    // a^+ a
    if (space_type == SpaceType::Occupied) {
      elementary_contractions_occupied(ops, s, contr_vec);
    }

    // 2. Pairwise contractions 1 ann + 1 cre operator:
    // ┌───┐
    // a   a^+
    if (space_type == SpaceType::Unoccupied) {
      elementary_contractions_unoccupied(ops, s, contr_vec);
    }

    // 3. 2k-legged contractions (k >= 1) of k cre + k ann operators:
    // ┌───┬───┬───┐
    // a^+ a   a   a^+
    if (space_type == SpaceType::General) {
      elementary_contractions_general(ops, s, contr_vec);
      if (orbital_subspaces->field_type(s) == FieldType::Fermion) {
        general_spaces.push_back(s);
      }
    }
  }

  if (general_spaces.size() > 1) {
    // Stage mixed-space contractions separately until graph canonicalization
    // handles two-leg contractions whose legs belong to different spaces.
    std::vector<ElementaryContraction> mixed_contr_vec;
    elementary_contractions_general_mixed(ops, general_spaces, mixed_contr_vec);
    PRINT(PrintLevel::Summary,
          cout << "\n  Mixed-space elementary contractions staged: "
               << mixed_contr_vec.size() << "\n";)
  }
  return contr_vec;
}

void WickTheorem::elementary_contractions_occupied(
    const OperatorProduct &ops, int s,
    std::vector<ElementaryContraction> &contr_vec) {
  int nops = ops.size();
  for (int l = 0; l < nops; l++) {             // loop over creation (left)
    for (int r = l + 1; r < nops; r++) {       // loop over annihilation (right)
      if (ops[l].cre(s) * ops[r].ann(s) > 0) { // is contraction viable?
        std::vector<GraphMatrix> new_contr(nops);
        new_contr[l].set_cre(s, 1);
        new_contr[r].set_ann(s, 1);
        contr_vec.push_back(new_contr);
        PRINT(PrintLevel::Summary,
              cout << fmt::format("\n    {:5d}:", contr_vec.size());
              PRINT_ELEMENTS(new_contr, " "););
      }
    }
  }
}

void WickTheorem::elementary_contractions_unoccupied(
    const OperatorProduct &ops, int s,
    std::vector<ElementaryContraction> &contr_vec) {
  int nops = ops.size();
  for (int l = 0; l < nops; l++) {             // loop over annihilation (left)
    for (int r = l + 1; r < nops; r++) {       // loop over creation (right)
      if (ops[l].ann(s) * ops[r].cre(s) > 0) { // is contraction viable?
        std::vector<GraphMatrix> new_contr(nops);
        new_contr[l].set_ann(s, 1);
        new_contr[r].set_cre(s, 1);
        contr_vec.push_back(new_contr);
        PRINT(PrintLevel::Summary,
              cout << fmt::format("\n    {:5d}:", contr_vec.size());
              PRINT_ELEMENTS(new_contr, " "););
      }
    }
  }
}

void WickTheorem::elementary_contractions_general(
    const OperatorProduct &ops, int s,
    std::vector<ElementaryContraction> &contr_vec) {
  const int nops = ops.size();
  // compute the largest possible cumulant for this space
  int sumcre = 0;
  int sumann = 0;
  for (int A = 0; A < nops; A++) {
    sumcre += ops[A].cre(s);
    sumann += ops[A].ann(s);
  }
  // the number of legs is limited by the smallest of number of cre/ann
  // operators and the maximum cumulant level allowed
  const int max_half_legs = std::min(std::min(sumcre, sumann), maxcumulant_);

  // in this algorithm we loop over all possible lengths of half-leg
  // contractions, partition this number into integers, permute these integers
  // to generate elementary contractions, and test if these are valid
  for (int half_legs = 1; half_legs <= max_half_legs; half_legs++) {
    PRINT(PrintLevel::Summary,
          cout << "\n    " << 2 * half_legs << "-legs contractions";)
    // create partitions of the number of half legs into at most nops numbers.
    // For half_legs = 2 and nops = 2, half_legs_part = [[2],[1,1]]
    auto half_legs_part = integer_partitions(half_legs, nops);
    // create lists of leg partitionings among all operators that are
    // compatible with the number of creation and annihilation operators
    //
    // these vectors store the number of cre/ann operators contracted per
    // operator
    std::vector<std::vector<int>> cre_legs_vec, ann_legs_vec;
    for (const auto part : half_legs_part) {
      // here we copy the partition and permute it (with added zeros, which
      // signify no contraction)
      std::vector<int> perm(nops, 0);
      std::copy(part.begin(), part.end(), perm.begin());
      std::sort(perm.begin(), perm.end());
      do {
        // check if compatible with creation/annihilation operators
        bool cre_compatible = true;
        bool ann_compatible = true;
        for (int A = 0; A < nops; A++) {
          if (ops[A].cre(s) < perm[A]) {
            cre_compatible = false;
          }
          if (ops[A].ann(s) < perm[A]) {
            ann_compatible = false;
          }
        }
        if (cre_compatible) {
          cre_legs_vec.push_back(perm);
        }
        if (ann_compatible) {
          ann_legs_vec.push_back(perm);
        }
      } while (std::next_permutation(perm.begin(), perm.end()));
    }

    // combine the creation and annihilation operators
    for (const auto cre_legs : cre_legs_vec) {
      for (const auto ann_legs : ann_legs_vec) {
        // count number of operators contracted
        int nops_contracted = 0;
        for (int A = 0; A < nops; A++) {
          nops_contracted += (cre_legs[A] + ann_legs[A] > 0);
        }
        // exclude operators that have legs only on one operator
        if (nops_contracted < 2)
          continue;
        // for a vector of GraphMatrix objects that represent this
        // contraction
        std::vector<GraphMatrix> new_contr(nops);
        for (int A = 0; A < nops; A++) {
          new_contr[A].set_cre(s, cre_legs[A]);
          new_contr[A].set_ann(s, ann_legs[A]);
        }
        contr_vec.push_back(new_contr);

        PRINT(PrintLevel::Summary,
              cout << fmt::format("\n    {:5d}:", contr_vec.size());
              PRINT_ELEMENTS(new_contr, " "););
      }
    }
  }
}

void WickTheorem::elementary_contractions_general_mixed(
    const OperatorProduct &ops, const std::vector<int> &spaces,
    std::vector<ElementaryContraction> &contr_vec) {
  const int nops = static_cast<int>(ops.size());
  const int nspaces = static_cast<int>(spaces.size());

  // Flatten the (operator, space) pairs into a list of slots. These vectors
  // store the maximum number of creation/annihilation legs that can be taken
  // from each slot. The slot corresponding to operator A and spaces[p] is
  // stored at A * nspaces + p.
  std::vector<int> cre_capacities;
  std::vector<int> ann_capacities;
  int total_cre = 0;
  int total_ann = 0;
  for (int A = 0; A < nops; ++A) {
    for (int s : spaces) {
      const int ncre = ops[A].cre(s);
      const int nann = ops[A].ann(s);
      cre_capacities.push_back(ncre);
      ann_capacities.push_back(nann);
      total_cre += ncre;
      total_ann += nann;
    }
  }

  // The number of half legs is limited by the total number of available
  // creation and annihilation operators across all general spaces and by the
  // maximum cumulant level allowed.
  const int max_half_legs =
      std::min(std::min(total_cre, total_ann), maxcumulant_);

  PRINT(PrintLevel::Summary,
        cout << "\n  Elementary contractions spanning general spaces:";)

  // Loop over all possible cumulant ranks. A contraction with half_legs = k
  // contains k creation and k annihilation operators and produces a 2k-legged
  // density/cumulant tensor.
  for (int half_legs = 1; half_legs <= max_half_legs; ++half_legs) {
    // Generate every distribution of half_legs creation/annihilation operators
    // over the flattened slots that is compatible with the available legs.
    // Bounded allocations generate each distribution once, without padding
    // and permuting integer partitions.
    const auto cre_allocations =
        generate_bounded_allocations(half_legs, cre_capacities);
    const auto ann_allocations =
        generate_bounded_allocations(half_legs, ann_capacities);

    // Combine the creation and annihilation distributions.
    for (const auto &cre : cre_allocations) {
      for (const auto &ann : ann_allocations) {
        // Count the operator blocks and orbital spaces touched by this
        // contraction.
        std::vector<bool> touched_ops(nops, false);
        std::vector<bool> touched_spaces(nspaces, false);
        for (int A = 0; A < nops; ++A) {
          for (int p = 0; p < nspaces; ++p) {
            const int pos = A * nspaces + p;
            if (cre[pos] + ann[pos] > 0) {
              touched_ops[A] = true;
              touched_spaces[p] = true;
            }
          }
        }

        const int num_touched_ops =
            std::count(touched_ops.begin(), touched_ops.end(), true);
        const int num_touched_spaces =
            std::count(touched_spaces.begin(), touched_spaces.end(), true);

        // As for the one-space generator, exclude contractions internal to a
        // single normal-ordered operator block. Also require at least two
        // spaces so this generator does not duplicate contractions produced
        // by elementary_contractions_general.
        if (num_touched_ops < 2 || num_touched_spaces < 2) {
          continue;
        }

        // Convert the flattened allocations into one GraphMatrix per operator
        // block. Unlike the one-space case, a GraphMatrix may have nonzero
        // entries in several general spaces.
        std::vector<GraphMatrix> new_contr(nops);
        for (int A = 0; A < nops; ++A) {
          for (int p = 0; p < nspaces; ++p) {
            const int pos = A * nspaces + p;
            new_contr[A].set_cre(spaces[p], cre[pos]);
            new_contr[A].set_ann(spaces[p], ann[pos]);
          }
        }
        contr_vec.push_back(new_contr);

        PRINT(PrintLevel::Summary,
              cout << fmt::format("\n    {:5d}:", contr_vec.size());
              PRINT_ELEMENTS(new_contr, " "););
      }
    }
  }
}
