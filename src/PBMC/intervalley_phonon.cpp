/**
 * @file intervalley_phonon.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-04-13
 *
 *
 */

#include "intervalley_phonon.hpp"

#include "pbmc_material_model.hpp"

namespace uepm::PBMC {

std::vector<intervalley_phonon_branch> make_silicon_intervalley_phonon_branches() {
    return load_pbmc_material_model().m_electron_intervalley_transitions;
}

}  // namespace uepm::PBMC
