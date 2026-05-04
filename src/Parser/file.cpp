/**
 * @file file.cpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief Source code file for file::file class
 * @version 0.1
 * @date 2021-08-19
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "file.hpp"

#include <filesystem>


namespace uepm {
namespace file {

/**
 * @brief Construct a new file::file object.
 *
 * The function will check if the file exist or not.
 *
 * @param filename
 */
file::file(const std::string &filename) : m_file_path(filename) {
    std::filesystem::path file_path(m_file_path);
    if (!std::filesystem::exists(file_path)) {
        throw std::invalid_argument("The file does not exists : " + filename);
    }
    if (!std::filesystem::is_regular_file(file_path)) {
        throw std::invalid_argument("The path is not a regular file : " + filename);
    }
}

file::~file() {}

}  // namespace file

}  // namespace uepm