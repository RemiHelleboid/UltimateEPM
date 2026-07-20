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
#include <cmath>
#include <filesystem>
#include <limits>
#include <regex>

#include "contact.hpp"

namespace uepm::device {
namespace {

std::optional<double> segment_box_entry_fraction(const mesh::bbox    &box,
                                                 const mesh::vector3 &start,
                                                 const mesh::vector3 &end) {
    if (box.is_inside(start) || box.is_inside(end)) {
        if (box.is_inside(start)) {
            return 0.0;
        }
    }

    double t_min = 0.0;
    double t_max = 1.0;

    const auto update_axis = [&](double start_value, double end_value, double min_value, double max_value) {
        const double delta = end_value - start_value;
        if (std::abs(delta) <= std::numeric_limits<double>::epsilon()) {
            return start_value >= min_value && start_value <= max_value;
        }

        double t1 = (min_value - start_value) / delta;
        double t2 = (max_value - start_value) / delta;
        if (t1 > t2) {
            std::swap(t1, t2);
        }

        t_min = std::max(t_min, t1);
        t_max = std::min(t_max, t2);
        return t_min <= t_max;
    };

    if (update_axis(start.x(), end.x(), box.get_x_min(), box.get_x_max()) &&
        update_axis(start.y(), end.y(), box.get_y_min(), box.get_y_max()) &&
        update_axis(start.z(), end.z(), box.get_z_min(), box.get_z_max())) {
        return t_min;
    }
    return std::nullopt;
}

}  // namespace

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
        return segment_box_entry_fraction(device_contact.get_contact_box(), point_A, point_B).has_value();
    });
}

std::optional<contact_crossing> device::find_first_contact_crossing(const mesh::vector3 &point_A,
                                                                    const mesh::vector3 &point_B) const {
    std::optional<contact_crossing> first_crossing;
    for (std::size_t contact_index = 0; contact_index < m_list_contacts.size(); ++contact_index) {
        const auto &contact        = m_list_contacts[contact_index];
        const auto  entry_fraction = segment_box_entry_fraction(contact.get_contact_box(), point_A, point_B);
        if (entry_fraction.has_value() &&
            (!first_crossing.has_value() || *entry_fraction < first_crossing->segment_fraction)) {
            first_crossing = contact_crossing{contact_index, contact.get_contact_name(), *entry_fraction};
        }
    }
    return first_crossing;
}

bool device::check_enters_contact(const mesh::vector3 &point) {
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
