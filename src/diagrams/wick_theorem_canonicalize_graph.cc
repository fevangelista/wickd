#include <algorithm>
#include <iomanip>
#include <iostream>

#include "helpers/combinatorics.h"
#include "helpers/helpers.h"
#include "helpers/orbital_space.h"
#include "helpers/stl_utils.hpp"

#include "contraction.h"
#include "operator.h"
#include "operator_expression.h"

#include "../algebra/expression.h"
#include "../algebra/sqoperator.h"
#include "../algebra/tensor.h"
#include "../algebra/term.h"

#include "wick_theorem.h"

#define PRINT(detail, code)                                                    \
  if (print_ >= detail) {                                                      \
    code                                                                       \
  }

void print_key(std::tuple<int, int, bool, int> key, int n);
void print_contraction(const OperatorProduct &ops,
                       const std::vector<Tensor> &tensors,
                       const std::vector<std::vector<bool>> &bit_map_vec,
                       const std::vector<SQOperator> &sqops,
                       const std::vector<int> sign_order);

void print_contraction_graph(const OperatorProduct &ops,
                             const CompositeContraction &contractions,
                             const std::vector<int> ops_perm,
                             const std::vector<int> contr_perm);

std::string contraction_signature(const OperatorProduct &ops,
                                  const CompositeContraction &contractions,
                                  const std::vector<int> &ops_perm,
                                  const std::vector<int> &contr_perm);

using namespace std;

std::pair<bool, scalar_t>
is_ops_permutation_valid(const OperatorProduct &ops, std::vector<int> ops_perm,
                         const std::vector<std::vector<bool>> &permutable) {
  // ok, so we are given a permutation of some operators and a matrix that tells
  // us which operators/contractions commute. Let's find out if this permutation
  // is consistent with the allowed permutations.

  // Here we reorder the operators one by one and in the process check if
  // this permutation is allowed. We also keep track of the sign, which
  // depends on the number of second quantized operators permuted
  int sign = 1;
  int n = ops_perm.size();
  for (int i = 0; i < n - 1; i++) {
    // Last i elements are already in place
    for (int j = 0; j < n - i - 1; j++) {
      if (ops_perm[j + 1] < ops_perm[j]) {
        // permutable
        if (permutable[ops_perm[j + 1]][ops_perm[j]]) {
          // the parity is given by the product of number of operators
          sign *= 1 - 2 * ((ops[ops_perm[j + 1]].num_ops() *
                            ops[ops_perm[j]].num_ops()) %
                           2);
          swap(ops_perm[j], ops_perm[j + 1]);
        } else {
          // this permutation is not valid
          return std::pair(false, scalar_t(1));
        }
      }
    }
  }
  return std::pair(true, scalar_t(sign));
}

bool do_contractions_commute(int i, int j, const OperatorProduct &ops,
                             const CompositeContraction &contractions) {
  // contractions commute if rearranging two operators does not change the final
  // result
  // loop over all elementary contractions
  for (const auto &el_contr : contractions) {
    // Check if this is a single contraction. Its creation and annihilation
    // legs may belong to different general spaces, so determine the role of
    // each operator across all spaces before comparing them.
    if (el_contr.num_ops() == 2) {
      int i_cre = 0;
      int i_ann = 0;
      int j_cre = 0;
      int j_ann = 0;
      for (int s = 0; s < orbital_subspaces->num_spaces(); ++s) {
        i_cre += el_contr[i].cre(s);
        i_ann += el_contr[i].ann(s);
        j_cre += el_contr[j].cre(s);
        j_ann += el_contr[j].ann(s);
      }
      if (i_cre * j_ann > 0 || i_ann * j_cre > 0) {
        return false;
      }
    }
  }
  return true;
}

bool graph_less(const std::pair<int, int> &l, const std::pair<int, int> &r,
                const std::vector<std::vector<int>> &ops_perms,
                std::vector<std::vector<int>> &con_perms,
                const OperatorProduct &ops,
                const CompositeContraction &contractions) {
  // here we are given the indices of the permutations of the operators and we
  // have to determine if permutation l is less than r

  const auto &l_ops_perm = ops_perms[l.first];
  const auto &r_ops_perm = ops_perms[r.first];
  const auto &l_con_perm = con_perms[l.second];
  const auto &r_con_perm = con_perms[r.second];

  // first compare the operators
  // consider three cases:
  // 1.  lops <  rops (return true)
  // 2.  lops == rops (must check the contractions)
  // 3.  lops >  rops (return false)
  const int no = l_ops_perm.size();
  for (int i = 0; i < no; i++) {
    if (ops[l_ops_perm[i]] < ops[r_ops_perm[i]]) {
      return true;
    }
    if (ops[r_ops_perm[i]] < ops[l_ops_perm[i]]) {
      return false;
    }
  }

  // if we made it here lops == rops
  const int nc = l_con_perm.size();
  // loop over all the contractions permuted
  for (int j = 0; j < nc; j++) {
    const auto &l_con = contractions[l_con_perm[j]];
    const auto &r_con = contractions[r_con_perm[j]];
    for (int i = 0; i < no; i++) {
      if (l_con[l_ops_perm[i]] < r_con[r_ops_perm[i]]) {
        return false;
      }
      if (r_con[r_ops_perm[i]] < l_con[l_ops_perm[i]]) {
        return true;
      }
    }
  }
  return false;
}

std::tuple<OperatorProduct, CompositeContraction, scalar_t>
WickTheorem::canonicalize_contraction_graph(
    const OperatorProduct &ops, const CompositeContraction &contractions) {

  PRINT(PrintLevel::Detailed,
        cout << "  Graph of the contraction to canonicalize:" << endl;
        const auto input_ops_perm = iota_vector<int>(ops.size());
        const auto input_contr_perm = iota_vector<int>(contractions.size());
        print_contraction_graph(ops, contractions, input_ops_perm,
                                input_contr_perm););

  const int nops = ops.size();

  // create a matrix that tells us if we can permute the position of two
  // operators
  std::vector<std::vector<bool>> commutable(nops,
                                            std::vector<bool>(nops, false));
  for (int i = 0; i < nops; i++) {
    for (int j = 0; j < nops; j++) {
      // check if commuting operators i and j changes the contraction
      if (do_contractions_commute(i, j, ops, contractions)) {
        commutable[i][j] = true;
      }
    }
  }

  PRINT(
      PrintLevel::Detailed, cout << "\n  Commutable operator matrix:" << endl;
      for (const auto row : commutable) {
        cout << "    ";
        PRINT_ELEMENTS(row);
        cout << endl;
      };
      cout << endl;);

  // 1. Generate all possible permutation of operators
  // these vectors store the allowed permutations
  std::vector<std::vector<int>> ops_perms;
  std::vector<scalar_t> ops_perms_sign;
  {
    // Loop over all permutations of operators
    std::vector<int> ops_perm(ops.size());
    std::iota(ops_perm.begin(), ops_perm.end(), 0);
    do {
      if (const auto [is_valid, sign] =
              is_ops_permutation_valid(ops, ops_perm, commutable);
          is_valid) {
        ops_perms.push_back(ops_perm);
        ops_perms_sign.push_back(sign);
        PRINT(PrintLevel::Detailed, cout << "  Operator permutation: ";
              PRINT_ELEMENTS(ops_perm);
              cout << ", sign = " << sign.repr() << " (allowed)" << endl;);
      } else {
        PRINT(PrintLevel::Detailed, cout << "  Operator permutation: ";
              PRINT_ELEMENTS(ops_perm); cout << " (not allowed)" << endl;);
      }
    } while (std::next_permutation(ops_perm.begin(), ops_perm.end()));
  }
  PRINT(PrintLevel::Detailed, cout << "\n  Found " << ops_perms.size()
                                   << " valid operator permutations\n"
                                   << endl;);

  std::vector<std::vector<int>> con_perms;
  {
    // Loop over all permutations of contractions
    std::vector<int> con_perm(contractions.size());
    std::iota(con_perm.begin(), con_perm.end(), 0);
    do {
      PRINT(PrintLevel::Detailed, cout << "  Contraction permutation: ";
            PRINT_ELEMENTS(con_perm); cout << endl;);
      con_perms.push_back(con_perm);
    } while (std::next_permutation(con_perm.begin(), con_perm.end()));
  }
  PRINT(PrintLevel::Detailed, cout << "\n  Found " << con_perms.size()
                                   << " valid contraction permutations\n"
                                   << endl;);

  // Decoupled two-pass search
  // graph_less is lexicographic: operators first, then contractions.
  // We exploit this to avoid forming the full nops! × c! Cartesian product.

  // Pass 1: find the minimum operator sequence among all valid op permutations.
  // Cost: O(|ops_perms| × nops).
  auto ops_seq_less = [&](int a, int b) -> bool {
    for (int i = 0; i < nops; i++) {
      if (ops[ops_perms[a][i]] < ops[ops_perms[b][i]])
        return true;
      if (ops[ops_perms[b][i]] < ops[ops_perms[a][i]])
        return false;
    }
    return false; // equal
  };

  int min_ops_idx = 0;
  for (int a = 1; a < static_cast<int>(ops_perms.size()); a++) {
    if (ops_seq_less(a, min_ops_idx))
      min_ops_idx = a;
  }

  // Pass 2: collect the tie set - all op permutations with the same minimum
  // operator sequence.  k == 1 whenever the operators are distinguishable,
  // which is the common case.
  std::vector<int> tied_ops;
  for (int a = 0; a < static_cast<int>(ops_perms.size()); a++) {
    if ((not ops_seq_less(a, min_ops_idx)) and
        (not ops_seq_less(min_ops_idx, a))) {
      tied_ops.push_back(a);
    }
  }

  PRINT(PrintLevel::Detailed, cout << "  Tie set size: " << tied_ops.size()
                                   << " op permutation(s)\n"
                                   << endl;);

  // Pass 3: among the |tied_ops| × c! pairs scan for the overall minimum.
  std::pair<int, int> best_graph{tied_ops[0], 0};
  for (int oi : tied_ops) {
    for (int ci = 0; ci < static_cast<int>(con_perms.size()); ci++) {
      std::pair<int, int> candidate{oi, ci};
      if (graph_less(candidate, best_graph, ops_perms, con_perms, ops,
                     contractions))
        best_graph = candidate;
    }
  }

  const auto [canonical_ops_perm_idx, canonical_contr_perm_idx] = best_graph;

  // Get the sign associated with rearranging the operators
  const scalar_t canonical_sign = ops_perms_sign[canonical_ops_perm_idx];

  // Get the canonical order of the operators
  OperatorProduct canonical_ops;
  for (int o : ops_perms[canonical_ops_perm_idx]) {
    canonical_ops.push_back(ops[o]);
  }

  // Get the canonical order of the contractions
  // permute the order and operator upon a contraction acts
  CompositeContraction canonical_contr;
  for (int c : con_perms[canonical_contr_perm_idx]) {
    std::vector<GraphMatrix> permuted_contr;
    for (int o : ops_perms[canonical_ops_perm_idx]) {
      permuted_contr.push_back(contractions[c][o]);
    }
    canonical_contr.push_back(permuted_contr);
  }

  PRINT(PrintLevel::Detailed,
        cout << "\n  Canonical form of the contraction:" << endl;
        cout << "    Sign = " << canonical_sign.repr() << endl;
        cout << "    Operator permutation: ";
        PRINT_ELEMENTS(ops_perms[canonical_ops_perm_idx]); cout << endl;
        cout << "    Contraction permutation: ";
        PRINT_ELEMENTS(con_perms[canonical_contr_perm_idx]); cout << endl;
        cout << "    Graph of the canonical contraction:" << endl;
        print_contraction_graph(ops, contractions,
                                ops_perms[canonical_ops_perm_idx],
                                con_perms[canonical_contr_perm_idx]);
        cout << endl;);

  return std::make_tuple(canonical_ops, canonical_contr, canonical_sign);
}
