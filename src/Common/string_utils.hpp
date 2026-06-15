/**
 * @file string_utils.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-10
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace uepm::utils {

bool starts_with_string(const std::string& text, const std::string& prefix) {
    return text.size() >= prefix.size() && text.compare(0, prefix.size(), prefix) == 0;
}

std::vector<std::string> make_mesh_stem_prefixes(std::string stem) {
    std::vector<std::string> prefixes;

    while (!stem.empty()) {
        prefixes.push_back(stem);

        const std::size_t pos = stem.rfind('_');
        if (pos == std::string::npos) {
            break;
        }

        stem = stem.substr(0, pos);
    }

    return prefixes;
}

bool looks_like_kstar_file(const std::filesystem::path& path) {
    const std::string filename = path.filename().string();

    return filename.find("_kstar_ibz_to_bz") != std::string::npos && path.extension() == ".txt";
}

std::filesystem::path find_kstar_file_for_mesh(const std::filesystem::path& mesh_path,
                                               const std::string&           maybe_kstar_filename) {
    namespace fs = std::filesystem;

    const fs::path mesh_dir = mesh_path.has_parent_path() ? mesh_path.parent_path() : fs::path(".");

    const std::string mesh_stem    = mesh_path.stem().string();
    const std::string kstar_suffix = "_kstar_ibz_to_bz.txt";

    /*
     * Only use the explicit argument if it really looks like a k-star file.
     * This prevents accidentally reading the .msh file as text.
     */
    if (!maybe_kstar_filename.empty()) {
        const fs::path explicit_path(maybe_kstar_filename);

        if (looks_like_kstar_file(explicit_path)) {
            if (fs::exists(explicit_path) && fs::is_regular_file(explicit_path)) {
                return explicit_path;
            }

            throw std::runtime_error("find_kstar_file_for_mesh: explicit k-star file does not exist: " +
                                     explicit_path.string());
        }
    }

    if (!fs::exists(mesh_dir) || !fs::is_directory(mesh_dir)) {
        throw std::runtime_error("find_kstar_file_for_mesh: mesh directory does not exist: " + mesh_dir.string());
    }

    const auto prefixes = make_mesh_stem_prefixes(mesh_stem);

    for (const auto& prefix : prefixes) {
        std::vector<fs::path> candidates;

        for (const auto& entry : fs::directory_iterator(mesh_dir)) {
            if (!entry.is_regular_file()) {
                continue;
            }

            const fs::path    candidate_path = entry.path();
            const std::string filename       = candidate_path.filename().string();

            if (!starts_with_string(filename, prefix)) {
                continue;
            }

            if (filename.find(kstar_suffix) == std::string::npos) {
                continue;
            }

            candidates.push_back(candidate_path);
        }

        std::sort(candidates.begin(), candidates.end());

        if (candidates.size() == 1) {
            return candidates.front();
        }

        if (candidates.size() > 1) {
            std::string message =
                "find_kstar_file_for_mesh: multiple k-star candidates found for prefix '" + prefix + "':";

            for (const auto& candidate : candidates) {
                message += "\n  - " + candidate.string();
            }

            message += "\nPlease pass the desired k-star file explicitly.";

            throw std::runtime_error(message);
        }
    }

    throw std::runtime_error("find_kstar_file_for_mesh: no k-star file found next to mesh " + mesh_path.string() +
                             " using prefixes derived from mesh stem '" + mesh_stem + "'");
}
}  // namespace uepm::utils