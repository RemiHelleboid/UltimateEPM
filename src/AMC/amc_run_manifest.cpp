/**
 * @file amc_run_manifest.cpp
 * @brief Utilities for exporting reproducible AMC run metadata.
 */

#include "amc_run_manifest.hpp"

#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace uepm::amc {

namespace {

std::string escape_value(std::string_view value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (const char character : value) {
        switch (character) {
            case '\\':
                escaped += "\\\\";
                break;
            case '\n':
                escaped += "\\n";
                break;
            case '\r':
                escaped += "\\r";
                break;
            default:
                escaped.push_back(character);
                break;
        }
    }
    return escaped;
}

std::string format_double(double value) {
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return stream.str();
}

std::string quote_command_argument(std::string_view argument) {
    if (argument.find_first_of(" \t\"'\\") == std::string_view::npos) {
        return std::string(argument);
    }

    std::string quoted = "\"";
    for (const char character : argument) {
        if (character == '"' || character == '\\') {
            quoted.push_back('\\');
        }
        quoted.push_back(character);
    }
    quoted.push_back('"');
    return quoted;
}

}  // namespace

void simulation_manifest::add(std::string_view section, std::string_view key, std::string_view value) {
    m_entries.push_back(entry{std::string(section), std::string(key), escape_value(value)});
}

void simulation_manifest::add(std::string_view section, std::string_view key, const char* value) {
    add(section, key, std::string_view(value));
}

void simulation_manifest::add(std::string_view section, std::string_view key, bool value) {
    add(section, key, value ? "true" : "false");
}

void simulation_manifest::add(std::string_view section, std::string_view key, std::size_t value) {
    add(section, key, std::to_string(value));
}

void simulation_manifest::add(std::string_view section, std::string_view key, int value) {
    add(section, key, std::to_string(value));
}

void simulation_manifest::add(std::string_view section, std::string_view key, double value) {
    add(section, key, format_double(value));
}

void simulation_manifest::write(const std::filesystem::path& filename) const {
    std::ofstream stream(filename);
    if (!stream.is_open()) {
        throw std::runtime_error("Could not open simulation manifest for writing: " + filename.string());
    }

    stream << "# UltimateEPM AMC simulation manifest\n";
    stream << "# format_version = 1\n";

    std::vector<std::string_view> sections;
    for (const auto& manifest_entry : m_entries) {
        bool section_seen = false;
        for (const auto section : sections) {
            if (section == manifest_entry.section) {
                section_seen = true;
                break;
            }
        }
        if (!section_seen) {
            sections.push_back(manifest_entry.section);
        }
    }

    for (const auto section : sections) {
        stream << "\n[" << section << "]\n";
        for (const auto& manifest_entry : m_entries) {
            if (manifest_entry.section == section) {
                stream << manifest_entry.key << " = " << manifest_entry.value << '\n';
            }
        }
    }
}

std::string current_utc_timestamp() {
    const auto now    = std::chrono::system_clock::now();
    const auto time_t = std::chrono::system_clock::to_time_t(now);
    std::tm    utc_time{};
#if defined(_WIN32)
    gmtime_s(&utc_time, &time_t);
#else
    gmtime_r(&time_t, &utc_time);
#endif

    std::ostringstream stream;
    stream << std::put_time(&utc_time, "%Y-%m-%dT%H:%M:%SZ");
    return stream.str();
}

std::string command_line_from_arguments(int argc, const char* const* argv) {
    std::string command_line;
    for (int argument_index = 0; argument_index < argc; ++argument_index) {
        if (!command_line.empty()) {
            command_line.push_back(' ');
        }
        command_line += quote_command_argument(argv[argument_index]);
    }
    return command_line;
}

std::string amc_project_version() { return UEPM_PROJECT_VERSION; }

std::string amc_build_type() { return UEPM_BUILD_TYPE; }

std::string amc_compiler() { return UEPM_COMPILER; }

}  // namespace uepm::amc
