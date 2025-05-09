#include "gauss.h"
#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include <stdexcept>

void printUsage() {
    std::cout << "Usage:" << std::endl;
    std::cout << "  ./main [input_file] [output_file]" << std::endl;
    std::cout << "  ./main --generate [size] [seed] [output_file]" << std::endl;
}

int main(int argc, char* argv[]) {
    try {
        if (argc == 3) {
            // Solve system from file
            std::string input_file = argv[1];
            std::string output_file = argv[2];
            
            auto start = std::chrono::high_resolution_clock::now();
            
            // Read system
            auto [A, b] = gauss::readSystemFromCSV(input_file);
            
            std::cout << "System size: " << A.rows() << "x" << A.cols() << std::endl;
            
            // Solve system
            auto solution = gauss::solve(A, b);
            
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            
            if (solution) {
                std::cout << "System solved successfully in " << elapsed.count() << " seconds" << std::endl;
                
                // Write solution to file
                gauss::writeSolutionToCSV(output_file, *solution);
                std::cout << "Solution written to " << output_file << std::endl;
            } else {
                std::cout << "System has no unique solution" << std::endl;
                return 1;
            }
        } 
        else if (argc == 5 && std::string(argv[1]) == "--generate") {
            // Generate random system
            int size = std::stoi(argv[2]);
            unsigned int seed = std::stoul(argv[3]);
            std::string output_file = argv[4];
            
            std::cout << "Generating random system of size " << size << " with seed " << seed << std::endl;
            
            auto [A, b] = gauss::generateRandomSystem(size, seed);
            
            // Solve system to verify it has a unique solution
            auto solution = gauss::solve(A, b);
            
            if (solution) {
                // Create augmented matrix [A|b] for the output
                Eigen::MatrixXd augmented(size, size + 1);
                augmented.leftCols(size) = A;
                augmented.rightCols(1) = b;
                
                // Write system to CSV file
                std::ofstream file(output_file);
                if (!file.is_open()) {
                    throw std::runtime_error("Failed to open file for writing: " + output_file);
                }
                
                for (int i = 0; i < size; i++) {
                    for (int j = 0; j < size + 1; j++) {
                        file << augmented(i, j);
                        if (j < size) {
                            file << ",";
                        }
                    }
                    file << std::endl;
                }
                
                file.close();
                std::cout << "Random system written to " << output_file << std::endl;
            } else {
                std::cout << "Generated system has no unique solution. Try a different seed." << std::endl;
                return 1;
            }
        } 
        else {
            printUsage();
            return 1;
        }
        
        return 0;
    } 
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
} 