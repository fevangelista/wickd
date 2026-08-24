#pragma once

#include <string>
#include <vector>

#include "../wickd-def.h"
#include "graph_matrix.h"
#include "helpers/product.hpp"

/// A class to represent an elementary contraction
class ElementaryContraction : public Product<GraphMatrix> {
public:
  /// Constructor
  ElementaryContraction() : Product<GraphMatrix>() {}

  /// Constructor. Set number of creation and annihilation operators
  ElementaryContraction(const std::vector<GraphMatrix> &graph_matrix)
      : Product<GraphMatrix>(graph_matrix) {}

  /// The number of second quantization operator contracted
  int num_ops() const;

  /// Return the sorted unique spaces touched by this contraction
  std::vector<int> spaces_in_elementary_contraction() const;
};

/// A class to represent an elementary contraction
class CompositeContraction : public Product<ElementaryContraction> {
public:
  /// Constructor
  CompositeContraction() : Product<ElementaryContraction>() {}
};

/// Compare two elementary contractions in canonical port-assignment order.
/// The operator-order vectors map each comparison position to an operator in
/// the corresponding contraction.
bool canonical_contraction_less(const ElementaryContraction &lhs,
                                const std::vector<int> &lhs_operator_order,
                                const ElementaryContraction &rhs,
                                const std::vector<int> &rhs_operator_order);

/// Return contractions sorted in canonical port-assignment order for their
/// current operator order.
CompositeContraction
canonical_contraction_order(const CompositeContraction &contractions);
