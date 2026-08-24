#include <algorithm>
#include <cassert>
#include <numeric>

#include "contraction.h"

#include "helpers/orbital_space.h"

int ElementaryContraction::num_ops() const {
  return std::accumulate(elements_.begin(), elements_.end(), 0,
                         [&](int a, const auto &b) { return a + b.num_ops(); });
}

std::vector<int>
ElementaryContraction::spaces_in_elementary_contraction() const {
  std::vector<int> vec;
  for (int s = 0; s < orbital_subspaces->num_spaces(); ++s) {
    for (const auto &graph_matrix : elements_) {
      if (graph_matrix.ann(s) + graph_matrix.cre(s) > 0) {
        vec.push_back(s);
        break;
      }
    }
  }
  return vec;
}

bool canonical_contraction_less(const ElementaryContraction &lhs,
                                const std::vector<int> &lhs_operator_order,
                                const ElementaryContraction &rhs,
                                const std::vector<int> &rhs_operator_order) {
  assert(lhs_operator_order.size() == rhs_operator_order.size());
  assert(lhs.size() == lhs_operator_order.size());
  assert(rhs.size() == rhs_operator_order.size());

  for (std::size_t i = 0; i < lhs_operator_order.size(); ++i) {
    const auto &lhs_graph = lhs[lhs_operator_order[i]];
    const auto &rhs_graph = rhs[rhs_operator_order[i]];
    if (lhs_graph < rhs_graph) {
      return false;
    }
    if (rhs_graph < lhs_graph) {
      return true;
    }
  }
  return false;
}

CompositeContraction
canonical_contraction_order(const CompositeContraction &contractions) {
  CompositeContraction ordered_contractions(contractions);
  if (ordered_contractions.size() < 2) {
    return ordered_contractions;
  }

  std::vector<int> operator_order(ordered_contractions[0].size());
  std::iota(operator_order.begin(), operator_order.end(), 0);
  std::stable_sort(
      ordered_contractions.begin(), ordered_contractions.end(),
      [&](const ElementaryContraction &lhs, const ElementaryContraction &rhs) {
        return canonical_contraction_less(lhs, operator_order, rhs,
                                          operator_order);
      });
  return ordered_contractions;
}
