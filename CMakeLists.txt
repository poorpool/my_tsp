cmake_minimum_required(VERSION 2.8)

set(CMAKE_C_STANDARD 99)

# 生成可执行文件
add_executable(greedy greedy.c)
add_executable(data_gen generator.c)

# 添加编译选项和宏定义
target_compile_options(greedy PRIVATE -std=c99)
target_compile_definitions(greedy PRIVATE _GNU_SOURCE)
target_compile_options(data_gen PRIVATE -std=c99)
target_compile_definitions(data_gen PRIVATE _GNU_SOURCE)

