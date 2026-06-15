/**
 * @file device.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-18
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "device.hpp"

#include <algorithm>
#include <filesystem>
#include <regex>

#include "contact.hpp"

namespace uepm::device {

void device::add_contact(const std::string &contact_name, const mesh::bbox contact_box, const double ohmic_resistance) {
    double         null_ohmic_resistance{ohmic_resistance};
    device_contact new_contact{contact_name, contact_box, null_ohmic_resistance};
    m_list_contacts.push_back(new_contact);
}

void device::add_contact(const std::string  &contact_name,
                         const mesh::vector3 corner1,
                         const mesh::vector3 corner2,
                         const double        ohmic_resistance) {
    double         null_ohmic_resistance{ohmic_resistance};
    device_contact new_contact{contact_name, corner1, corner2, null_ohmic_resistance};
    m_list_contacts.push_back(new_contact);
}

bool device::check_crossing_contact(const mesh::vector3 &point_A, const mesh::vector3 &point_B) const {
    return std::any_of(m_list_contacts.begin(), m_list_contacts.end(), [&](const device_contact &device_contact) {
        return device_contact.line_intersect_contact(point_A, point_B);
    });
}

bool device::check_enters_contact(const mesh::vector3 &point) {
    // return std::any_of(m_list_contacts.begin(), m_list_contacts.end(),
    //                    [&](const device_contact &device_contact) { return device_contact.point_in_contact_box(point);
    //                    });
    for (auto &device_contact : m_list_contacts) {
        if (device_contact.point_in_contact_box(point)) {
            device_contact.add_contact_current(1.0);
            return true;
        }
    }
    return false;
}

bool device::check_enters_contact(mesh::sp_element crossing_element) const {
    return std::any_of(m_list_contacts.begin(), m_list_contacts.end(), [&](const device_contact &device_contact) {
        return device_contact.is_element_in_contact(crossing_element);
    });
}

std::string get_node_number_from_filename(const std::string &filename) {
    std::regex  pattern_node_number{"n(\\d+)_"};
    std::smatch regex_result;
    std::regex_search(filename, regex_result, pattern_node_number);
    return regex_result.str();
}

std::optional<std::string> get_voltage_from_filename(const std::string &filename) {
    std::regex  pattern_node_number{"n\\d+_(\\d+.\\d+)"};
    std::smatch regex_result;
    std::regex_search(filename, regex_result, pattern_node_number);
    if (regex_result.size() > 1) {
        return regex_result[1].str();
    } else {
        return {};
    };
}

device_metadata::device_metadata(const std::string &device_file_name) {
    std::filesystem::path device_path(device_file_name);
    m_name                 = device_path.stem().string();
    m_node_number          = get_node_number_from_filename(m_name);
    auto voltage_extracted = get_voltage_from_filename(device_file_name);
    if (voltage_extracted.has_value()) {
        m_voltage_bias = std::stod(voltage_extracted.value());
    } else {
        m_voltage_bias = std::numeric_limits<double>::quiet_NaN();
    }
}

}  // namespace uepm::device