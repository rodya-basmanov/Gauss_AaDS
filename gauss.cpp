#include "gauss.h"
#include "lazycsv.hpp"
#include <fstream>
#include <random>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <string>

namespace gauss {

std::optional<Eigen::VectorXd> solve(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
    // Make a copy of the input matrices to avoid modifying them
    Eigen::MatrixXd augmented(A.rows(), A.cols() + 1);
    augmented.leftCols(A.cols()) = A;
    augmented.rightCols(1) = b;
    
    const double eps = 1e-10; // Threshold for zero check
    const int n = A.rows();
    
    // Forward elimination with partial pivoting
    for (int i = 0; i < n; i++) {
        // Find pivot (the row with the largest absolute value in the current column)
        int pivot_row = i;
        double pivot_value = std::abs(augmented(i, i));
        
        for (int j = i + 1; j < n; j++) {
            double val = std::abs(augmented(j, i));
            if (val > pivot_value) {
                pivot_value = val;
                pivot_row = j;
            }
        }
        
        // Check if matrix is singular
        if (pivot_value < eps) {
            return std::nullopt; // No unique solution
        }
        
        // Swap rows if needed
        if (pivot_row != i) {
            augmented.row(i).swap(augmented.row(pivot_row));
        }
        
        // Eliminate below the pivot
        for (int j = i + 1; j < n; j++) {
            const double factor = augmented(j, i) / augmented(i, i);
            // Using vectorized row operations (hardware acceleration)
            augmented.row(j) -= factor * augmented.row(i);
        }
    }
    
    // Back substitution
    Eigen::VectorXd x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = augmented(i, n);
        for (int j = i + 1; j < n; j++) {
            sum -= augmented(i, j) * x(j);
        }
        if (std::abs(augmented(i, i)) < eps) {
            return std::nullopt; // No unique solution
        }
        x(i) = sum / augmented(i, i);
    }
    
    return x;
}

std::pair<Eigen::MatrixXd, Eigen::VectorXd> readSystemFromCSV(const std::string& filename) {
    lazycsv::parser parser(filename);
    
    // Get dimensions
    int rows = 0;
    int cols = 0;
    
    // Count rows and determine columns from the first row
    for (const auto& row : parser) {
        if (rows == 0) {
            // Count cells in the first row to determine matrix width
            int cell_count = 0;
            for (const auto& _ : row) {
                cell_count++;
            }
            cols = cell_count;
        }
        rows++;
    }
    
    // Last column is the b vector
    const int n = cols - 1;
    
    // Initialize matrices
    Eigen::MatrixXd A(rows, n);
    Eigen::VectorXd b(rows);
    
    // Load data
    int row_idx = 0;
    for (const auto& row : parser) {
        int col_idx = 0;
        for (const auto& cell : row) {
            // Convert string_view to string before using std::stod
            std::string cell_str(cell.raw());
            double value = std::stod(cell_str);
            
            if (col_idx < n) {
                A(row_idx, col_idx) = value;
            } else {
                b(row_idx) = value;
            }
            col_idx++;
        }
        row_idx++;
    }
    
    return {A, b};
}

void writeSolutionToCSV(const std::string& filename, const Eigen::VectorXd& solution) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for writing: " + filename);
    }
    
    for (int i = 0; i < solution.size(); i++) {
        file << "x" << i << "," << solution(i) << std::endl;
    }
    
    file.close();
}

std::pair<Eigen::MatrixXd, Eigen::VectorXd> generateRandomSystem(int size, unsigned int seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    
    // Generate a random solution vector
    Eigen::VectorXd x_true(size);
    for (int i = 0; i < size; i++) {
        x_true(i) = dist(gen);
    }
    
    // Generate a random matrix with good conditioning
    Eigen::MatrixXd A(size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            A(i, j) = dist(gen);
        }
    }
    
    // Ensure diagonal dominance for better conditioning
    for (int i = 0; i < size; i++) {
        double row_sum = 0;
        for (int j = 0; j < size; j++) {
            if (i != j) {
                row_sum += std::abs(A(i, j));
            }
        }
        // Make diagonal element larger than the sum of other elements
        A(i, i) = dist(gen) + row_sum + 1.0;
    }
    
    // Calculate right-hand side b = A*x_true
    Eigen::VectorXd b = A * x_true;
    
    return {A, b};
}

} // namespace gauss 