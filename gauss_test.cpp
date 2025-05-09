#include "gauss.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <fstream>
#include <random>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

// Создаем класс для тестирования с общей функциональностью
class GaussTest : public ::testing::Test {
protected:
    // Функция для проверки наличия файла
    bool fileExists(const std::string& filename) {
        std::ifstream file(filename);
        return file.good();
    }
    
    // Функция для удаления файла, если он существует
    void removeFileIfExists(const std::string& filename) {
        if (fileExists(filename)) {
            std::remove(filename.c_str());
        }
    }
    
    // Функция для чтения CSV-файла с решением
    std::vector<double> readSolutionCSV(const std::string& filename) {
        std::vector<double> solution;
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            return solution;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            size_t pos = line.find(',');
            if (pos != std::string::npos) {
                std::string value = line.substr(pos + 1);
                solution.push_back(std::stod(value));
            }
        }
        
        return solution;
    }
    
    // Функция для создания системы с известным решением
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> createSystemWithKnownSolution(int size, const Eigen::VectorXd& expected_solution) {
        // Создаем диагонально доминирующую матрицу для обеспечения устойчивости
        Eigen::MatrixXd A = Eigen::MatrixXd::Random(size, size);
        
        // Обеспечиваем диагональное доминирование
        for (int i = 0; i < size; i++) {
            double row_sum = 0;
            for (int j = 0; j < size; j++) {
                if (i != j) {
                    row_sum += std::abs(A(i, j));
                }
            }
            A(i, i) = row_sum + 1.0; // Делаем диагональный элемент больше суммы остальных
        }
        
        // Создаем правую часть b = A * expected_solution
        Eigen::VectorXd b = A * expected_solution;
        
        return {A, b};
    }
};

// Тест на простую систему 2x2
TEST_F(GaussTest, SimpleSystem) {
    // Создаем систему с известным решением
    Eigen::VectorXd expected_solution(2);
    expected_solution << 1.0, 2.0;
    
    auto [A, b] = createSystemWithKnownSolution(2, expected_solution);
    
    // Решаем систему
    auto solution = gauss::solve(A, b);
    
    // Проверяем наличие решения
    ASSERT_TRUE(solution.has_value());
    
    // Проверяем правильность решения
    for (int i = 0; i < 2; i++) {
        EXPECT_NEAR((*solution)(i), expected_solution(i), 1e-8);
    }
    
    // Проверяем, что A * solution = b
    Eigen::VectorXd computed_b = A * (*solution);
    for (int i = 0; i < b.size(); i++) {
        EXPECT_NEAR(computed_b(i), b(i), 1e-8);
    }
}

// Тест на вырожденную матрицу
TEST_F(GaussTest, SingularMatrix) {
    // Создаем вырожденную матрицу (ранг < n)
    Eigen::MatrixXd A(3, 3);
    A << 1, 2, 3,
         2, 4, 6,
         3, 6, 9;  // Третья строка = (строка 1) + (строка 2)
    
    Eigen::VectorXd b(3);
    b << 1, 2, 3;  // Не коллинеарный вектор
    
    // Решаем систему
    auto solution = gauss::solve(A, b);
    
    // Проверяем отсутствие решения (матрица вырожденная)
    EXPECT_FALSE(solution.has_value());
}

// Тест на большую систему со случайными числами
TEST_F(GaussTest, RandomLargeSystem) {
    // Создаем генератор с фиксированным начальным значением
    unsigned int seed = 42;
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distribution(-10.0, 10.0);
    
    // Размер системы
    const int size = 20;
    
    // Генерируем ожидаемое решение
    Eigen::VectorXd expected_solution(size);
    for (int i = 0; i < size; i++) {
        expected_solution(i) = distribution(generator);
    }
    
    // Создаем систему с известным решением
    auto [A, b] = createSystemWithKnownSolution(size, expected_solution);
    
    // Решаем систему
    auto solution = gauss::solve(A, b);
    
    // Проверяем наличие решения
    ASSERT_TRUE(solution.has_value());
    
    // Проверяем правильность решения (с меньшей точностью для больших систем)
    for (int i = 0; i < size; i++) {
        EXPECT_NEAR((*solution)(i), expected_solution(i), 1e-6);
    }
}

// Тест на чтение и запись CSV
TEST_F(GaussTest, CSVReadWrite) {
    // Имена файлов для тестирования
    std::string input_filename = "test_input.csv";
    std::string output_filename = "test_output.csv";
    
    // Убеждаемся, что файлы не существуют
    removeFileIfExists(input_filename);
    removeFileIfExists(output_filename);
    
    // Создаем тестовую систему
    const int size = 3;
    Eigen::MatrixXd A(size, size);
    A << 4, 1, 2,
         3, 5, 1,
         1, 2, 3;
    
    Eigen::VectorXd expected_solution(size);
    expected_solution << 1.0, 2.0, 3.0;
    
    Eigen::VectorXd b = A * expected_solution;
    
    // Записываем систему в CSV-файл
    std::ofstream file(input_filename);
    ASSERT_TRUE(file.is_open());
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            file << A(i, j);
            if (j < size - 1) {
                file << ",";
            }
        }
        file << "," << b(i) << std::endl;
    }
    file.close();
    
    // Читаем систему из CSV
    auto [A_read, b_read] = gauss::readSystemFromCSV(input_filename);
    
    // Проверяем корректность чтения
    ASSERT_EQ(A_read.rows(), size);
    ASSERT_EQ(A_read.cols(), size);
    ASSERT_EQ(b_read.size(), size);
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            EXPECT_DOUBLE_EQ(A_read(i, j), A(i, j));
        }
        EXPECT_DOUBLE_EQ(b_read(i), b(i));
    }
    
    // Решаем систему
    auto solution = gauss::solve(A_read, b_read);
    ASSERT_TRUE(solution.has_value());
    
    // Записываем решение в CSV
    gauss::writeSolutionToCSV(output_filename, *solution);
    
    // Проверяем существование выходного файла
    ASSERT_TRUE(fileExists(output_filename));
    
    // Читаем решение из CSV
    auto solution_read = readSolutionCSV(output_filename);
    ASSERT_EQ(solution_read.size(), size);
    
    // Проверяем корректность записанного решения
    for (int i = 0; i < size; i++) {
        EXPECT_NEAR(solution_read[i], expected_solution(i), 1e-8);
    }
    
    // Очищаем после теста
    removeFileIfExists(input_filename);
    removeFileIfExists(output_filename);
}

// Тест на генерацию случайной системы
TEST_F(GaussTest, RandomSystemGeneration) {
    // Фиксированное начальное значение для воспроизводимости
    unsigned int seed = 123;
    const int size = 5;
    
    // Генерируем случайную систему
    auto [A, b] = gauss::generateRandomSystem(size, seed);
    
    // Проверяем размерности
    EXPECT_EQ(A.rows(), size);
    EXPECT_EQ(A.cols(), size);
    EXPECT_EQ(b.size(), size);
    
    // Решаем сгенерированную систему
    auto solution = gauss::solve(A, b);
    
    // Проверяем наличие решения
    ASSERT_TRUE(solution.has_value());
    
    // Проверяем, что A * solution = b
    Eigen::VectorXd computed_b = A * (*solution);
    double relative_error = (computed_b - b).norm() / b.norm();
    EXPECT_LT(relative_error, 1e-10);
}

// Тест производительности на больших системах
TEST_F(GaussTest, PerformanceLargeSystem) {
    // Фиксированное начальное значение для воспроизводимости
    unsigned int seed = 12345;
    const int size = 100;
    
    // Генерируем случайную систему
    auto [A, b] = gauss::generateRandomSystem(size, seed);
    
    // Засекаем время
    auto start = std::chrono::high_resolution_clock::now();
    
    // Решаем сгенерированную систему
    auto solution = gauss::solve(A, b);
    
    // Замеряем время
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    // Выводим информацию о производительности
    std::cout << "Solving system of size " << size << " took " << elapsed.count() << " seconds" << std::endl;
    
    // Проверяем наличие решения
    ASSERT_TRUE(solution.has_value());
    
    // Проверяем точность решения
    Eigen::VectorXd computed_b = A * (*solution);
    double relative_error = (computed_b - b).norm() / b.norm();
    EXPECT_LT(relative_error, 1e-8);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
} 