/**
 * @file amc_run_manifest.hpp
 * @brief Utilities for exporting reproducible AMC run metadata.
 */

#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace uepm::amc {

class simulation_manifest {
 public:
    void add(std::string_view section, std::string_view key, std::string_view value);
    void add(std::string_view section, std::string_view key, const char* value);
    void add(std::string_view section, std::string_view key, bool value);
    void add(std::string_view section, std::string_view key, std::size_t value);
    void add(std::string_view section, std::string_view key, int value);
    void add(std::string_view section, std::string_view key, double value);

    void write(const std::filesystem::path& filename) const;

 private:
    struct entry {
        std::string section;
        std::string key;
        std::string value;
    };

    std::vector<entry> m_entries;
};

std::string current_utc_timestamp();
std::string command_line_from_arguments(int argc, const char* const* argv);
std::string amc_project_version();
std::string amc_build_type();
std::string amc_compiler();

}  // namespace uepm::amc
