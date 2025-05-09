#include "gauss.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <fstream>
#include <random>

// Test solving a simple 2x2 system
TEST(GaussSolverTest, SimpleSystem) {
    Eigen::MatrixXd A(2, 2);
    A << 2, 1,
         1, 3;
    
    Eigen::VectorXd b(2);
    b << 5, 7;
    
    auto solution = gauss::solve(A, b);
    
    ASSERT_TRUE(solution.has_value());
    EXPECT_NEAR((*solution)(0), 2.0, 1e-10);
    EXPECT_NEAR((*solution)(1), 5.0/3.0, 1e-10);
    
    // Verify Ax = b
    Eigen::VectorXd expected_b = A * (*solution);
    for (int i = 0; i < b.size(); i++) {
        EXPECT_NEAR(expected_b(i), b(i), 1e-10);
    }
}

// Test solving a larger system
TEST(GaussSolverTest, LargerSystem) {
    int size = 10;
    auto [A, expected_x] = gauss::generateRandomSystem(size, 42);
    Eigen::VectorXd b = A * expected_x;
    
    auto solution = gauss::solve(A, b);
    
    ASSERT_TRUE(solution.has_value());
    
    // Verify solution is close to expected
    for (int i = 0; i < size; i++) {
        EXPECT_NEAR((*solution)(i), expected_x(i), 1e-8);
    }
}

// Test with a singular matrix
TEST(GaussSolverTest, SingularMatrix) {
    Eigen::MatrixXd A(3, 3);
    A << 1, 2, 3,
         2, 4, 6,
         3, 6, 9;
    
    Eigen::VectorXd b(3);
    b << 6, 12, 18;
    
    auto solution = gauss::solve(A, b);
    
    EXPECT_FALSE(solution.has_value());
}

// Test CSV file reading and writing
TEST(GaussCSVTest, ReadWriteCSV) {
    // Create a test CSV file
    std::string test_input = "test_input.csv";
    std::string test_output = "test_output.csv";
    
    // Create a simple 3x3 system
    std::ofstream input_file(test_input);
    input_file << "1,2,3,14\n"
               << "4,5,6,32\n"
               << "7,8,10,58\n";
    input_file.close();
    
    // Read the system
    auto [A, b] = gauss::readSystemFromCSV(test_input);
    
    // Check dimensions
    EXPECT_EQ(A.rows(), 3);
    EXPECT_EQ(A.cols(), 3);
    EXPECT_EQ(b.size(), 3);
    
    // Check values
    EXPECT_EQ(A(0, 0), 1);
    EXPECT_EQ(A(0, 1), 2);
    EXPECT_EQ(A(0, 2), 3);
    EXPECT_EQ(A(1, 0), 4);
    EXPECT_EQ(A(1, 1), 5);
    EXPECT_EQ(A(1, 2), 6);
    EXPECT_EQ(A(2, 0), 7);
    EXPECT_EQ(A(2, 1), 8);
    EXPECT_EQ(A(2, 2), 10);
    
    EXPECT_EQ(b(0), 14);
    EXPECT_EQ(b(1), 32);
    EXPECT_EQ(b(2), 58);
    
    // Solve the system
    auto solution = gauss::solve(A, b);
    ASSERT_TRUE(solution.has_value());
    
    // Write the solution
    gauss::writeSolutionToCSV(test_output, *solution);
    
    // Clean up
    std::remove(test_input.c_str());
    std::remove(test_output.c_str());
}

// Test random system generation
TEST(GaussRandomTest, GenerateRandomSystem) {
    int size = 5;
    unsigned int seed = 123;
    
    auto [A, b] = gauss::generateRandomSystem(size, seed);
    
    // Check dimensions
    EXPECT_EQ(A.rows(), size);
    EXPECT_EQ(A.cols(), size);
    EXPECT_EQ(b.size(), size);
    
    // Solve the system
    auto solution = gauss::solve(A, b);
    
    // Should have a unique solution
    ASSERT_TRUE(solution.has_value());
    
    // Verify Ax = b
    Eigen::VectorXd calculated_b = A * (*solution);
    for (int i = 0; i < b.size(); i++) {
        EXPECT_NEAR(calculated_b(i), b(i), 1e-8);
    }
}

// Test performance with a large system
TEST(GaussPerformanceTest, LargeSystem) {
    int size = 100;
    auto [A, expected_x] = gauss::generateRandomSystem(size, 12345);
    Eigen::VectorXd b = A * expected_x;
    
    auto solution = gauss::solve(A, b);
    
    ASSERT_TRUE(solution.has_value());
    
    // Verify solution
    Eigen::VectorXd residual = A * (*solution) - b;
    double error = residual.norm() / b.norm();
    
    EXPECT_LT(error, 1e-8);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
} 