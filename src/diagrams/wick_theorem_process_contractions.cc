#include <algorithm>
#include <future>
#include <iomanip>
#include <iostream>

#include "helpers/unordered_dense.h"

#define FMT_HEADER_ONLY
#include <fmt/core.h> // for fmt::format

#include "helpers/combinatorics.h"
#include "helpers/helpers.h"
#include "helpers/orbital_space.h"
#include "helpers/stl_utils.hpp"
#include "helpers/timer.hpp"

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
void print_key(OpKey key, int n);
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

namespace {

CompositeContraction
canonical_port_order(const CompositeContraction &contractions) {
  CompositeContraction ordered_contractions(contractions);

  // Explicit operator legs are assigned using offsets accumulated while
  // elementary contractions are evaluated. Sort the contractions first so
  // that this port assignment does not depend on their incoming order. This
  // is the same descending lexicographic order selected by graph_less for a
  // fixed operator order.
  std::sort(ordered_contractions.begin(), ordered_contractions.end(),
            [](const ElementaryContraction &lhs,
               const ElementaryContraction &rhs) { return rhs < lhs; });
  return ordered_contractions;
}

} // namespace

Expression WickTheorem::process_contractions(scalar_t factor,
                                             const OperatorProduct &ops,
                                             const int minrank,
                                             const int maxrank) {
  PRINT(PrintLevel::Summary,
        std::cout << "\n- Step 3. Processing contractions" << std::endl;)

  Expression result;

  // loop over all contractions and keep those with the correct rank
  // contraction_vec stores a list of elementary contractions appearing
  // in a term
  int nprocessed = 0;
  int ops_rank = ops.num_ops();
  for (const auto &contraction_vec : contractions_) {
    int contr_rank = 0;
    for (int c : contraction_vec) {
      contr_rank += elementary_contractions_[c].num_ops();
    }
    int term_rank = ops_rank - contr_rank;

    if ((term_rank >= minrank) and (term_rank <= maxrank)) {
      nprocessed++;
      PRINT(PrintLevel::Basic,
            cout << "\n\n  Contraction: " << nprocessed
                 << "  Operator rank: " << ops_rank - contr_rank << endl;)

      CompositeContraction contraction;
      for (int c : contraction_vec) {
        contraction.push_back(elementary_contractions_[c]);
      }

      timer tc;
      const auto [best_ops, best_contractions, sign] =
          do_canonicalize_graph_
              ? canonicalize_contraction_graph(ops, contraction)
              : std::make_tuple(ops, contraction, scalar_t(1));
      timers_["canonicalize_contraction_graph"] += tc.get();

      timer te;
      std::pair<SymbolicTerm, scalar_t> term_factor =
          evaluate_contraction(best_ops, best_contractions, factor * sign);
      timers_["evaluate_contraction"] += te.get();

      SymbolicTerm &term = term_factor.first;
      scalar_t canonicalize_factor = term.canonicalize();
      canonicalize_factor *= term_factor.second;
      result += {term, canonicalize_factor};
      // result.add(
      // std::make_pair(term, term_factor.second * canonicalize_factor));

      PRINT(PrintLevel::Summary,
            Term t(term_factor.second * canonicalize_factor, term);
            cout << "\n    term: " << t << endl;)
    }
  }
  if (nprocessed == 0) {
    PRINT(PrintLevel::Summary, std::cout << "\n  No contractions generated\n"
                                         << std::endl;)
  }
  return result;
}

Expression WickTheorem::process_contractions_threads(scalar_t factor,
                                                     const OperatorProduct &ops,
                                                     const int minrank,
                                                     const int maxrank) {
  PRINT(PrintLevel::Summary,
        std::cout << "\n- Step 3. Processing contractions" << std::endl;)

  // Vector to store futures of each async call
  std::vector<std::future<Expression>> futures;
  int ops_rank = ops.num_ops();

  // Launch async tasks for each contraction
  for (const auto &contraction_vec : contractions_) {
    futures.push_back(
        std::async(std::launch::async, [this, &ops, &contraction_vec, factor,
                                        minrank, maxrank, ops_rank] {
          return process_single_contraction(ops, contraction_vec, minrank,
                                            maxrank, factor, ops_rank);
        }));
  }

  // Combine the results of all futures into a single Expression
  Expression result;

  // Use std::accumulate to sum all the results
  for (auto &fut : futures) {
    result += fut.get(); // Get the result from the future and add to the result
  }

  // Check if any contractions were processed
  if (futures.empty()) {
    PRINT(PrintLevel::Summary, std::cout << "\n  No contractions generated\n"
                                         << std::endl;)
  }

  return result;
}

Expression WickTheorem::process_single_contraction(
    const OperatorProduct &ops, const std::vector<int> &contraction_vec,
    const int minrank, const int maxrank, const scalar_t factor,
    const int ops_rank) {
  Expression result;
  int contr_rank = 0;
  for (int c : contraction_vec) {
    contr_rank += elementary_contractions_[c].num_ops();
  }
  int term_rank = ops_rank - contr_rank;

  if ((term_rank >= minrank) and (term_rank <= maxrank)) {

    CompositeContraction contraction;
    for (int c : contraction_vec) {
      contraction.push_back(elementary_contractions_[c]);
    }

    timer tc;
    const auto [best_ops, best_contractions, sign] =
        do_canonicalize_graph_
            ? canonicalize_contraction_graph(ops, contraction)
            : std::make_tuple(ops, contraction, scalar_t(1));
    timers_["canonicalize_contraction_graph"] += tc.get();

    timer te;
    std::pair<SymbolicTerm, scalar_t> term_factor =
        evaluate_contraction(best_ops, best_contractions, factor * sign);
    timers_["evaluate_contraction"] += te.get();

    SymbolicTerm &term = term_factor.first;
    scalar_t canonicalize_factor = term.canonicalize();
    canonicalize_factor *= term_factor.second;
    result += {term, canonicalize_factor};
    // result += (std::make_pair(term, canonicalize_factor));
    // // result.add(
    // std::make_pair(term, term_factor.second * canonicalize_factor));

    PRINT(PrintLevel::Summary,
          Term t(term_factor.second * canonicalize_factor, term);
          cout << "\n    term: " << t << endl;)
  }
  return result;
}

std::string contraction_signature(const OperatorProduct &ops,
                                  const CompositeContraction &contractions,
                                  const std::vector<int> &ops_perm,
                                  const std::vector<int> &contr_perm) {

  std::string s;
  int nops = ops.size();
  for (int i = 0; i < nops; i++) {
    s += ops[ops_perm[i]].label();
    s += signature(ops[ops_perm[i]].graph_matrix());
  }

  // 2. Compare contractions
  int ncontr = contractions.size();
  for (int i = 0; i < ncontr; i++) {
    for (int j = 0; j < nops; j++) {
      s += signature(contractions[contr_perm[i]][ops_perm[j]]);
    }
  }
  return s;
}

std::pair<SymbolicTerm, scalar_t>
WickTheorem::evaluate_contraction(const OperatorProduct &ops,
                                  const CompositeContraction &contractions,
                                  scalar_t factor) {
  // Canonicalize the order in which contractions consume the otherwise
  // indistinguishable legs of each operator tensor. Graph canonicalization
  // already produces this order, while the noncanonical path may not.
  const auto ordered_contractions = canonical_port_order(contractions);

  // 1. Get the Tensor objects and SQOperator vector corresponding to the
  // uncontracted term
  auto tensors_sqops_op_map = contraction_tensors_sqops(ops);
  std::vector<Tensor> &tensors = std::get<0>(tensors_sqops_op_map);
  std::vector<SQOperator> &sqops = std::get<1>(tensors_sqops_op_map);

  // this map takes the operator index (op), orbital space (s), the sqop type
  // (true = cre, false = ann), and an index and maps it to the operators as
  // they are stored in a vector:
  //
  //  ankerl::unordered_dense::map<OpKey, int> op_map;
  //  key = {v, s, creation, i}
  //
  ankerl::unordered_dense::map<OpKey, int> &op_map =
      std::get<2>(tensors_sqops_op_map);

  // 2. Apply the contractions to the second quantized operators and add new
  // tensors (density matrices, cumulants)

  // counts of how many second quantized operators are not contracted
  std::vector<GraphMatrix> ops_offset(ops.size());

  // a counter to keep track of the positions assigned to operators
  int sorted_position = 0;

  // a counter that keeps track of the number of sq operators contracted
  int nsqops_contracted = 0;

  // a sign factor that keeps into account a negative sign introduced by
  // sorting the operators in a unoccupied-unoccupied contraction
  int unoccupied_sign = 1;

  index_map_t pair_contraction_reindex_map;

  // vector to store the order of operators
  std::vector<int> sign_order(sqops.size(), -1);
  std::vector<std::vector<bool>> bit_map_vec;

  // Loop over elementary contractions
  for (const ElementaryContraction &contraction : ordered_contractions) {
    // a bit array to keep track of which operators are contracted
    std::vector<bool> bit_map(sqops.size(), false);

    // Find the rank and spaces of this contraction
    int rank = contraction.num_ops();
    const auto spaces = contraction.spaces_in_elementary_contraction();
    if (spaces.empty()) {
      throw std::runtime_error(
          "WickTheorem::evaluate_contraction: empty elementary contraction");
    }
    nsqops_contracted += rank;

    // find the position of the creation operators
    std::vector<int> pos_cre_sqops =
        elements_vec_to_pos(contraction, ops_offset, op_map, true);
    // find the position of the annihilation operators
    std::vector<int> pos_ann_sqops =
        elements_vec_to_pos(contraction, ops_offset, op_map, false);

    if (static_cast<int>(pos_cre_sqops.size() + pos_ann_sqops.size()) != rank) {
      throw std::runtime_error(
          "WickTheorem::evaluate_contraction: mapped operator count does not "
          "match the elementary contraction rank");
    }

    // mark the creation operators contracted and their order
    for (int c : pos_cre_sqops) {
      bit_map[c] = true;
      sign_order[c] = sorted_position;
      sorted_position += 1;
    }
    // mark the annihilation operators contracted and their order
    for (int a : pos_ann_sqops) {
      bit_map[a] = true;
      sign_order[a] = sorted_position;
      sorted_position += 1;
    }

    // Occupied and unoccupied contractions act in exactly one space. A
    // contraction may span several spaces only when all of them are general.
    SpaceType dmstruc = orbital_subspaces->space_type(spaces[0]);
    if (spaces.size() > 1) {
      const bool all_general =
          std::all_of(spaces.begin(), spaces.end(), [](int s) {
            return orbital_subspaces->space_type(s) == SpaceType::General;
          });
      if (not all_general) {
        throw std::runtime_error(
            "WickTheorem::evaluate_contraction: a multi-space contraction "
            "may contain only general spaces");
      }
      dmstruc = SpaceType::General;
    }

    // Pairwise contractions creation-annihilation:
    // ________
    // |      |
    // a^+(i) a(j) = delta(i,j)

    if (dmstruc == SpaceType::Occupied) {
      // Reindex the annihilator (j) to the creator (i)
      Index cre_index = sqops[pos_cre_sqops[0]].index();
      Index ann_index = sqops[pos_ann_sqops[0]].index();
      pair_contraction_reindex_map[ann_index] = cre_index;
    }

    // Pairwise contractions creation-annihilation:
    // ______
    // |    |
    // a(i) a^+(j) = delta(i,j)

    if (dmstruc == SpaceType::Unoccupied) {
      // Reindex the creator (j) to the annihilator (i)
      Index cre_index = sqops[pos_cre_sqops[0]].index();
      Index ann_index = sqops[pos_ann_sqops[0]].index();
      pair_contraction_reindex_map[cre_index] = ann_index;
      unoccupied_sign *= -1; // this factor is to compensate for the fact that
                             // we order operator in a canonical form in which
                             // annihilators are to the left of creation
                             // operators
    }

    // 2k-legged contractions (k >= 2) of k creation and k annihilation
    // operators:
    // _____________
    // |   |   |   |
    // a^+ a   a   a^+

    if (dmstruc == SpaceType::General) {
      std::vector<Index> lower;
      std::vector<Index> upper;
      // collect indices of creation operators for the upper indices
      for (int c : pos_cre_sqops) {
        upper.push_back(sqops[c].index());
      }
      // collect indices of annihilation operators for the lower indices
      for (int a : pos_ann_sqops) {
        lower.push_back(sqops[a].index());
      }
      // reverse the lower indices
      std::reverse(lower.begin(), lower.end());
      // prepare the label
      std::string label;
      if (rank == 2) {
        if (pos_cre_sqops[0] < pos_ann_sqops[0]) {
          label = "gamma1";
        } else {
          label = "eta1";
          unoccupied_sign *= -1; // this factor is to compensate for the fact
                                 // that we order operator in a canonical form
                                 // in which annihilators are to the left of
                                 // creation operators
        }
      } else {
        label = "lambda" + std::to_string(rank / 2);
      }
      // add the cumulant to the list of tensors
      tensors.push_back(
          Tensor(label, lower, upper, SymmetryType::Antisymmetric));
    }
    bit_map_vec.push_back(bit_map);
  }

  // assign an order to the uncontracted operators
  // creation operators come before annihilation operators
  for (SQOperatorType type :
       {SQOperatorType::Creation, SQOperatorType::Annihilation}) {
    for (int s = 0; s < orbital_subspaces->num_spaces(); s++) {
      for (int i = 0; i < sqops.size(); i++) {
        if ((sign_order[i] == -1) and (sqops[i].index().space() == s) and
            (sqops[i].type() == type)) {
          sign_order[i] = sorted_position;
          sorted_position += 1;
        }
      }
    }
  }

  PRINT(PrintLevel::Basic,
        print_contraction(ops, tensors, bit_map_vec, sqops, sign_order);)

  int sign = unoccupied_sign * permutation_sign(sign_order);

  PRINT(PrintLevel::All, PRINT_ELEMENTS(sign_order, "\n  positions: "););

  std::vector<std::pair<int, SQOperator>> sorted_sqops;
  sorted_position = 0;
  for (const auto &sqop : sqops) {
    sorted_sqops.push_back(std::make_pair(sign_order[sorted_position], sqop));
    sorted_position += 1;
  }
  std::sort(sorted_sqops.begin(), sorted_sqops.end());

  sqops.clear();
  for (int i = nsqops_contracted; i < sorted_sqops.size(); ++i) {
    sqops.push_back(sorted_sqops[i].second);
  }

  // find the combinatorial factor associated with this contraction
  scalar_t comb_factor = combinatorial_factor(ops, ordered_contractions);

  SymbolicTerm term;

  for (const auto &tensor : tensors) {
    term.add(tensor);
  }
  for (const auto &sqop : sqops) {
    term.add(sqop);
  }
  for (const auto &op : ops) {
    factor *= op.factor();
  }

  term.reindex(pair_contraction_reindex_map);

  PRINT(PrintLevel::Summary,
        cout << fmt::format("  sign =                 {:d}", sign) << endl;
        cout << fmt::format("  factor =               {:s}", factor.repr())
             << endl;
        cout << fmt::format("  combinatorial factor = {:s}", comb_factor.repr())
             << endl;);

  return std::make_pair(term, sign * factor * comb_factor);
}

static inline OpKey make_op_key(int v, int s, bool creation, int i) {
  return {static_cast<uint16_t>(v), static_cast<uint16_t>(s),
          static_cast<uint16_t>(creation ? 1 : 0), static_cast<uint16_t>(i)};
}

std::tuple<std::vector<Tensor>, std::vector<SQOperator>,
           ankerl::unordered_dense::map<OpKey, int>>
WickTheorem::contraction_tensors_sqops(const OperatorProduct &ops) {

  std::vector<SQOperator> sqops;
  std::vector<Tensor> tensors;
  ankerl::unordered_dense::map<OpKey, int> op_map;

  index_counter ic(orbital_subspaces->num_spaces());

  // Loop over all operators
  int n = 0;
  for (int o = 0, maxo = ops.size(); o < maxo; o++) {
    const auto &op = ops[o];
    // Loop over creation operators (lower indices)
    std::vector<Index> lower;
    for (int s = 0; s < orbital_subspaces->num_spaces(); s++) {
      for (int c = 0; c < op.cre(s); c++) {
        Index idx(s, ic.next_index(s)); // get next available index
        sqops.push_back(SQOperator(SQOperatorType::Creation, idx));
        lower.push_back(idx);
        auto key = make_op_key(o, s, true, c);
        op_map[key] = n;
        PRINT(PrintLevel::All, print_key(key, n););
        n += 1;
      }
    }

    // Loop over annihilation operators (upper indices)
    // the annihilation operators are layed out in a reversed order (hence the
    // need to reverse the upper indices of the tensor, see below)
    std::vector<Index> upper;
    for (int s = orbital_subspaces->num_spaces() - 1; s >= 0; s--) {
      for (int a = op.ann(s) - 1; a >= 0; a--) {
        Index idx(s, ic.next_index(s)); // get next available index
        sqops.push_back(SQOperator(SQOperatorType::Annihilation, idx));
        upper.push_back(idx);
        auto key = make_op_key(o, s, false, a);
        op_map[key] = n;
        PRINT(PrintLevel::All, print_key(key, n););
        n += 1;
      }
    }

    // reverse the order of the upper indices
    std::reverse(upper.begin(), upper.end());
    tensors.push_back(
        Tensor(op.label(), lower, upper, SymmetryType::Antisymmetric));
  }
  return make_tuple(tensors, sqops, op_map);
}

std::vector<int> WickTheorem::elements_vec_to_pos(
    const ElementaryContraction &elements_vec,
    std::vector<GraphMatrix> &ops_offset,
    ankerl::unordered_dense::map<OpKey, int> &op_map, bool creation) {

  std::vector<int> result;

  const auto spaces = elements_vec.spaces_in_elementary_contraction();

  PRINT(PrintLevel::All, cout << "\n  GraphMatrix to position:" << endl;);

  // Loop over all graph matrices. Creation operators are laid out with spaces
  // in ascending order, while annihilation operators use descending space
  // order. This is the same ordering used by contraction_tensors_sqops.
  for (int v = 0; v < elements_vec.size(); v++) {
    const auto &graph_matrix = elements_vec[v];
    for (int p = 0; p < static_cast<int>(spaces.size()); ++p) {
      const int space_pos =
          creation ? p : static_cast<int>(spaces.size()) - p - 1;
      const int s = spaces[space_pos];
      const int nops = creation ? graph_matrix.cre(s) : graph_matrix.ann(s);

      // Assign the explicit operator positions for this block and space. The
      // offset records legs already consumed by earlier contractions.
      const int ops_off =
          creation ? ops_offset[v].cre(s) : ops_offset[v].ann(s);
      for (int i = 0; i < nops; i++) {
        auto key = make_op_key(v, s, creation, ops_off + i);
        if (op_map.count(key) == 0) {
          PRINT(PrintLevel::All, print_key(key, -1););
          throw std::runtime_error(
              "WickTheorem::elements_vec_to_pos: operator position not found");
        }
        const int sqop_pos = op_map[key];
        result.push_back(sqop_pos);
        PRINT(PrintLevel::All, print_key(key, sqop_pos););
      }

      // Update the offset independently for every operator and space.
      if (creation) {
        ops_offset[v].set_cre(s, ops_off + nops);
      } else {
        ops_offset[v].set_ann(s, ops_off + nops);
      }
    }
  }
  return result;
}

scalar_t
WickTheorem::combinatorial_factor(const OperatorProduct &ops,
                                  const CompositeContraction &contractions) {

  scalar_t factor = 1;

  // stores the offset for each uncontracted operator
  std::vector<GraphMatrix> free_graph_matrix;
  for (const auto &op : ops) {
    free_graph_matrix.push_back(op.graph_matrix());
  }

  // for each contraction find the combinatorial factor associated to
  // permutations of contracted indices
  for (const auto &contraction : contractions) {
    for (int v = 0; v < contraction.size(); v++) {
      const auto &graph_matrix = contraction[v];
      for (int s = 0; s < orbital_subspaces->num_spaces(); s++) {
        const auto &[kcre, kann] = graph_matrix.elements(s);
        const auto &[ncre, nann] = free_graph_matrix[v].elements(s);
        factor *= binomial(ncre, kcre);
        factor *= binomial(nann, kann);
      }
      free_graph_matrix[v] -= graph_matrix;
    }
  }

  // This last factor accounts for permutations of contractions with
  // multiplicity higher than one
  std::map<ElementaryContraction, int> contraction_count;
  for (const auto &contraction : contractions) {
    contraction_count[contraction] += 1;
  }
  for (const auto &kv : contraction_count) {
    factor /= factorial(kv.second);
  }
  return factor;
}
