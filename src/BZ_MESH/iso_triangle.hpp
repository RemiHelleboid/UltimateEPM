/**
 * @file iso_triangle.hpp
 * @author Rémi Helleboid (remi.helleboid@st.com)
 * @brief
 * @version 0.1
 * @date 2022-07-16
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <array>
#include <random>
#include <vector>

#include "vector.hpp"

namespace uepm::mesh_bz {

/**
 * @brief This class represent the iso energy triangle within a Tetra.
 *
 */
class IsoTriangle {
 private:
    /**
     * @brief The positions of the 3 vertices of the triangle.
     *
     */
    std::array<vector3, 3> m_list_vertices{};

    /**
     * @brief Iso energy the triangle was constructed with.
     *
     */
    double m_iso_energy;

    /**
     * @brief List of the edges vectors of the tetrahedra.
     * Stored as follows: [v01, v02, v12]
     * where for example v12 = m_list_vertices[1] - m_list_vertices[2].
     *
     */
    std::array<vector3, 3> m_list_edges{};

    /**
     * @brief Signed volume of the tetrahedra.
     * The sign depends on the "orientation" of the tetrahedra.
     *
     */
    double m_surface;

 public:
    IsoTriangle() = delete;

    /**
     * @brief Construct a new Iso Triangle object from the list of vertices and the iso energy.
     *
     * @param VtxA
     * @param VtxB
     * @param VtxC
     * @param iso_energy
     */
    IsoTriangle(const vector3& VtxA, const vector3& VtxB, const vector3& VtxC, double iso_energy)
        : m_list_vertices{VtxA, VtxB, VtxC},
          m_iso_energy(iso_energy),
          m_list_edges{VtxA - VtxB, VtxA - VtxC, VtxB - VtxC},
          m_surface{0.5 * (cross_product(m_list_edges[1], m_list_edges[2]).norm())} {}

    /**
     * @brief Get the iso energy of the triangle.
     *
     * @return double
     */
    double get_iso_energy() const { return m_iso_energy; }

    /**
     * @brief Get the surface of the triangle.
     *
     * @return double
     */
    double get_signed_surface() const { return m_surface; }

    /**
     * @brief Draw a random point in the triangle. Uniform distribution.
     * See https://math.stackexchange.com/questions/538458/how-to-sample-points-on-a-triangle-surface-in-3d
     *
     * @param gen
     * @return vector3
     */
    vector3 draw_random_uniform_point_in_triangle(std::mt19937& gen) const {
        std::uniform_real_distribution<> dis_double(0.0, 1.0);
        double                           random_a = dis_double(gen);
        double                           random_b = dis_double(gen);
        double                           lambda_1 = (1.0 - sqrt(random_a));
        double                           lambda_2 = sqrt(random_a) * (1.0 - random_b);
        double                           lambda_3 = sqrt(random_a) * random_b;
        //   return get_center();
        return m_list_vertices[0] * lambda_1 + m_list_vertices[1] * lambda_2 + m_list_vertices[2] * lambda_3;
    }
};

}  // namespace uepm::mesh_bz