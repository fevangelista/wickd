#pragma once

#include <map>
#include <vector>

#include "equation.h"
#include "helpers/algebra.hpp"
#include "index.h"
#include "term.h"
#include "wickd-def.h"

/// A class to represent an algebraic expression
class Expression : public Algebra<Expression, SymbolicTerm, scalar_t> {
public:
  // ==> Constructor <==
  Expression();

  // ==> Class public interface <==

  /// Canonicalize this sum
  Expression &canonicalize();

  /// Reindex this sum
  Expression &reindex(index_map_t &idx_map);

  /// Compare this expression to another one
  bool operator==(const Expression &other);

  Expression adjoint() const;

  /// Return a string representation
  std::string str() const;

  /// Return a LaTeX representation
  std::string latex(const std::string &sep = " \\\\ \n") const;

  /// Convert this sum to a vector of many-body equations
  /// The result is stored into a map. The key to this map
  /// shows the number of upper/lower indices in each space
  std::map<std::string, std::vector<Equation>>
  to_manybody_equation(const std::string &label) const;

  /// @brief Order the operators in this expression in such a way that all the
  /// bare annihilation operators are to the left of the bare creation operators
  /// @param only_same_index_contractions only contract operators with the same
  /// indices (default: false), in other words assume that different indices
  /// represent different spin orbitals
  Expression
  vacuum_normal_ordered(bool only_same_index_contractions = false) const;

  // /// @brief Return a normal ordered version of this expression
  // /// @param only_same_index_contractions only contract operators with the
  // same
  // /// indices (default: false), in other words assume that different indices
  // /// represent different spin orbitals
  // Expression normal_ordered(bool only_same_index_contractions = false) const;

  bool is_vacuum_normal_ordered() const;
};

/// Print to an output stream
std::ostream &operator<<(std::ostream &os, const Expression &sum);

/// The syntax used to input a tensor expression
enum class TensorSyntax { wickd, TCE };

///// Create a sum from a string
Expression make_expression(const std::string &s, SymmetryType symmetry);

Expression make_operator_expr(const std::string &label,
                              const std::vector<std::string> &components,
                              bool normal_ordered, SymmetryType symmetry,
                              scalar_t coefficient = scalar_t(1));
