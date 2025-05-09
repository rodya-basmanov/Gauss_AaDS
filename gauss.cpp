#include "gauss.h"
#include "lazycsv.hpp"
#include <fstream>
#include <random>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace gauss {

std::optional<Eigen::VectorXd> solve(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
    // We'll use Eigen's built-in solver for better accuracy
    Eigen::VectorXd x;
    
    // A direct solver from Eigen for dense matrices
    x = A.colPivHouseholderQr().solve(b);
    
    // Check if the solution is valid
    double relative_error = (A * x - b).norm() / b.norm();
    if (relative_error > 1e-8) {
        // If the relative error is too large, check the rank
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
        if (qr.rank() < A.cols()) {
            // The matrix is rank-deficient, no unique solution
            return std::nullopt;
        }
    }
    
    return x;
}

std::pair<Eigen::MatrixXd, Eigen::VectorXd> readSystemFromCSV(const std::string& filename) {
    // First read all data into memory
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for reading: " + filename);
    }
    
    std::string line;
    int cols = 0;
    
    while (std::getline(file, line)) {
        std::vector<double> row_data;
        std::string cell;
        std::stringstream lineStream(line);
        
        while (std::getline(lineStream, cell, ',')) {
            row_data.push_back(std::stod(cell));
        }
        
        if (data.empty()) {
            cols = row_data.size();
        } else if (row_data.size() != cols) {
            throw std::runtime_error("Inconsistent number of columns in CSV");
        }
        
        data.push_back(row_data);
    }
    
    // Now create the matrix and vector
    int rows = data.size();
    int n = cols - 1; // Last column is b
    
    Eigen::MatrixXd A(rows, n);
    Eigen::VectorXd b(rows);
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < n; j++) {
            A(i, j) = data[i][j];
        }
        b(i) = data[i][n];
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