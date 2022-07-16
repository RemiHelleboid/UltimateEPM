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

#include "vector.hpp"


namespace bz_mesh {


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
    double m_signed_surface;

 public:
    IsoTriangle() = delete;
    IsoTriangle(const vector3& VtxA, const vector3& VtxB, const vector3& VtxC, double iso_energy):
        m_list_vertices{VtxA, VtxB, VtxC},
        m_iso_energy(iso_energy),
        m_list_edges{VtxA-VtxB, VtxA-VtxC, VtxB-VtxC},
        m_signed_surface{0.5 * (cross_product(m_list_edges[0], m_list_edges[1]).norm())} {}

    double get_iso_energy() const {return m_iso_energy;}
    double get_signed_surface() const {return m_signed_surface;}

};

}  // namespace bz_mesh