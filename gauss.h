#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <optional>

namespace gauss {

/**
 * Solves linear system Ax = b using Gaussian elimination
 * @param A The coefficient matrix
 * @param b The right-hand side vector
 * @return The solution vector or nullopt if system has no unique solution
 */
std::optional<Eigen::VectorXd> solve(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);

/**
 * Reads linear system from CSV file
 * @param filename Path to CSV file
 * @return Pair of matrix A and vector b
 */
std::pair<Eigen::MatrixXd, Eigen::VectorXd> readSystemFromCSV(const std::string& filename);

/**
 * Writes solution to CSV file
 * @param filename Path to output CSV file
 * @param solution Solution vector
 */
void writeSolutionToCSV(const std::string& filename, const Eigen::VectorXd& solution);

/**
 * Generates random linear system with known solution
 * @param size Size of the system (number of equations and variables)
 * @param seed Random seed for reproducibility
 * @return Pair of matrix A and vector b
 */
std::pair<Eigen::MatrixXd, Eigen::VectorXd> generateRandomSystem(int size, unsigned int seed);

} // namespace gauss 