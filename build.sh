#!/bin/bash

set -e

# Сборка программы
make main

# запуск программы
./main

# если все прошло успешно то выводит All right, а если нет то ничего не выводит, и мы понимаем, что сработало set -e
echo "Build and execution completed successfully!" 