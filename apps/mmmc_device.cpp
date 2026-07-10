/**
 * @file mmmc_device.cpp
 * @brief Self-consistent Mixed Method Monte Carlo device simulation.
 */

#include <fmt/core.h>
#include <tclap/CmdLine.h>
#include <tclap/MultiArg.h>

#include <filesystem>
#include <string>
#include <vector>

#include "mmmc_device_config.hpp"
#include "mmmc_device_runner.hpp"

namespace {

std::vector<std::string> normalize_set_arguments(int argc, const char** argv) {
    std::vector<std::string> arguments;
    arguments.reserve(static_cast<std::size_t>(argc));
    for (int index = 0; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--set" && index + 2 < argc) {
            const std::string path = argv[index + 1];
            if (path.find('=') == std::string::npos) {
                arguments.push_back(argument);
                arguments.push_back(path + "=" + argv[index + 2]);
                index += 2;
                continue;
            }
        }
        arguments.push_back(argument);
    }
    return arguments;
}

}  // namespace

int main(int argc, const char** argv) {
    try {
        const std::string        original_command_line = uepm::MMMC::mmmc_command_line_from_arguments(argc, argv);
        std::vector<std::string> normalized_arguments  = normalize_set_arguments(argc, argv);
        std::vector<const char*> normalized_argv;
        normalized_argv.reserve(normalized_arguments.size());
        for (const auto& argument : normalized_arguments) {
            normalized_argv.push_back(argument.c_str());
        }

        TCLAP::CmdLine               cmd("Self-consistent Mixed Method Monte Carlo device simulation.", ' ', "1.0");
        TCLAP::ValueArg<std::string> arg_config("c",
                                                "config",
                                                "YAML simulation configuration file.",
                                                false,
                                                "",
                                                "path");
        TCLAP::MultiArg<std::string> arg_overrides(
            "",
            "set",
            "Override one YAML value: path.to.setting=value or path.to.setting value.",
            false,
            "assignment");
        TCLAP::ValueArg<std::string> arg_write_basic_config("",
                                                            "write-config",
                                                            "Write a complete YAML configuration file and exit.",
                                                            false,
                                                            "",
                                                            "path");

        cmd.add(arg_config);
        cmd.add(arg_overrides);
        cmd.add(arg_write_basic_config);
        cmd.parse(static_cast<int>(normalized_argv.size()), normalized_argv.data());

        if (!arg_write_basic_config.getValue().empty()) {
            const std::filesystem::path output_file = arg_write_basic_config.getValue();
            uepm::MMMC::write_basic_device_mmmc_config(output_file);
            fmt::print("Wrote complete MMMC configuration to {}\n", output_file.string());
            return 0;
        }
        if (arg_config.getValue().empty()) {
            throw std::invalid_argument("--config is required unless --write-config is used.");
        }

        auto run_config         = uepm::MMMC::load_device_mmmc_config(arg_config.getValue(), arg_overrides.getValue());
        run_config.command_line = original_command_line;
        uepm::MMMC::run_self_consistent_device_mmmc_simulation(run_config);
        return 0;
    } catch (const TCLAP::ArgException& error) {
        fmt::print(stderr, "Argument error: {} for argument {}\n", error.error(), error.argId());
        return 1;
    } catch (const std::exception& error) {
        fmt::print(stderr, "Error: {}\n", error.what());
        return 1;
    }
}
