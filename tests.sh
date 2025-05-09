#!/bin/bash

set -e

# Сборка программы
make tests

# запуск программы
./tests

# если все прошло успешно то выводит All right, а если нет то ничего не выводит, и мы понимаем, что сработало set -e
echo "All tests passed successfully!" 