/**
 * @file mesh_tetra.cpp
 * @author Rémi Helleboid (remi.helleboid@st.com)
 * @brief
 * @version 0.1
 * @date 2022-07-14
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "mesh_tetra.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <vector>

#include "bbox_mesh.hpp"
#include "iso_triangle.hpp"
#include "physical_constants.hpp"

namespace uepm::mesh_bz {
namespace {

double geometry_scale(const std::array<vector3, 6>& edges) {
    return std::max({edges[0].norm() * edges[1].norm() * edges[2].norm(), 1.0});
}

bool nearly_same_point(const vector3& lhs, const vector3& rhs, double tolerance) {
    return (lhs - rhs).norm() <= tolerance;
}

void add_unique_point(std::vector<vector3>& points, const vector3& point, double tolerance) {
    if (std::none_of(points.begin(), points.end(), [&](const vector3& existing) {
            return nearly_same_point(existing, point, tolerance);
        })) {
        points.push_back(point);
    }
}

void add_unique_point(IsoEnergyPolygon& polygon, const vector3& point, double tolerance) {
    for (std::size_t index = 0; index < polygon.size; ++index) {
        if (nearly_same_point(polygon.points[index], point, tolerance)) {
            return;
        }
    }
    if (polygon.size < polygon.points.size()) {
        polygon.points[polygon.size++] = point;
    }
}

IsoEnergyPolygon compute_band_iso_energy_polygon_impl(const std::array<Vertex*, 4>& vertices,
                                                      const bbox_mesh&              bounding_box,
                                                      const std::array<double, 4>& energies,
                                                      double                        iso_energy) {
    const auto   minmax           = std::minmax_element(energies.begin(), energies.end());
    const double energy_scale     = std::max({std::abs(*minmax.first), std::abs(*minmax.second), 1.0});
    const double energy_tolerance = 1e-12 * energy_scale;
    if (iso_energy < *minmax.first - energy_tolerance || iso_energy > *minmax.second + energy_tolerance ||
        *minmax.second - *minmax.first <= energy_tolerance) {
        return {};
    }

    constexpr std::array<std::array<std::size_t, 2>, 6> edge_vertices = {
        {{{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}, {{1, 3}}, {{2, 3}}}};
    const double point_tolerance = 1e-12 * std::max(bounding_box.get_diagonal_size(), 1.0);
    IsoEnergyPolygon intersections;

    for (const auto& edge : edge_vertices) {
        const std::size_t i    = edge[0];
        const std::size_t j    = edge[1];
        const double      di   = energies[i] - iso_energy;
        const double      dj   = energies[j] - iso_energy;
        const bool        i_on = std::abs(di) <= energy_tolerance;
        const bool        j_on = std::abs(dj) <= energy_tolerance;

        if (i_on) {
            add_unique_point(intersections, vertices[i]->get_position(), point_tolerance);
        }
        if (j_on) {
            add_unique_point(intersections, vertices[j]->get_position(), point_tolerance);
        }
        if (!i_on && !j_on && ((di < 0.0) != (dj < 0.0))) {
            const double fraction = (iso_energy - energies[i]) / (energies[j] - energies[i]);
            const vector3 point =
                (1.0 - fraction) * vertices[i]->get_position() + fraction * vertices[j]->get_position();
            add_unique_point(intersections, point, point_tolerance);
        }
    }
    return intersections;
}

double triangle_area_fast(const vector3& a, const vector3& b, const vector3& c) noexcept {
    return 0.5 * cross_product(b - a, c - a).norm();
}

double order_iso_polygon_and_compute_area(IsoEnergyPolygon& polygon) noexcept {
    if (polygon.size == 3) {
        return triangle_area_fast(polygon.points[0], polygon.points[1], polygon.points[2]);
    }
    if (polygon.size != 4) {
        return 0.0;
    }

    // Fix point 0 and test the three distinct cyclic orders of the other
    // points. A convex planar quadrilateral has the largest shoelace area in
    // either of its two cyclic orientations; crossed orders have less area.
    constexpr std::array<std::array<std::size_t, 4>, 3> orders = {
        std::array<std::size_t, 4>{0, 1, 2, 3},
        std::array<std::size_t, 4>{0, 1, 3, 2},
        std::array<std::size_t, 4>{0, 2, 1, 3},
    };

    double                     maximum_area = 0.0;
    std::array<std::size_t, 4> best_order   = orders[0];
    for (const auto& order : orders) {
        vector3 area_vector{};
        for (std::size_t index = 0; index < order.size(); ++index) {
            area_vector += cross_product(polygon.points[order[index]],
                                         polygon.points[order[(index + 1) % order.size()]]);
        }
        const double area = 0.5 * area_vector.norm();
        if (area > maximum_area) {
            maximum_area = area;
            best_order   = order;
        }
    }
    const auto unordered_points = polygon.points;
    for (std::size_t index = 0; index < best_order.size(); ++index) {
        polygon.points[index] = unordered_points[best_order[index]];
    }
    return maximum_area;
}

}  // namespace

bbox_mesh Tetra::compute_bounding_box() const {
    std::array<double, 4> coordinates_x;
    std::array<double, 4> coordinates_y;
    std::array<double, 4> coordinates_z;
    for (int i = 0; i < 4; ++i) {
        coordinates_x[i] = m_list_vertices[i]->get_position().x();
        coordinates_y[i] = m_list_vertices[i]->get_position().y();
        coordinates_z[i] = m_list_vertices[i]->get_position().z();
    }
    auto min_max_x = std::minmax_element(coordinates_x.begin(), coordinates_x.end());
    auto min_max_y = std::minmax_element(coordinates_y.begin(), coordinates_y.end());
    auto min_max_z = std::minmax_element(coordinates_z.begin(), coordinates_z.end());
    return bbox_mesh(*min_max_x.first,
                     *min_max_x.second,
                     *min_max_y.first,
                     *min_max_y.second,
                     *min_max_z.first,
                     *min_max_z.second);
}

const bbox_mesh& Tetra::get_bounding_box() const { return m_bbox; }

/**
 * @brief Construct a new Tetra by passing directly the array of the four pointers to the vertices.
 *
 * @param list_vertices
 */
Tetra::Tetra(std::size_t index, const std::array<Vertex*, 4>& list_vertices)
    : m_index(index),
      m_list_vertices(list_vertices),
      m_nb_bands{0} {
    if (std::any_of(m_list_vertices.begin(), m_list_vertices.end(), [](const Vertex* vertex) {
            return vertex == nullptr;
        })) {
        throw std::invalid_argument("A tetrahedron cannot contain null vertex pointers.");
    }
    m_nb_bands = m_list_vertices[0]->get_number_bands();
    if (std::any_of(m_list_vertices.begin(), m_list_vertices.end(), [&](const Vertex* vertex) {
            return vertex->get_number_bands() != m_nb_bands;
        })) {
        throw std::invalid_argument("All tetrahedron vertices must contain the same number of bands.");
    }
    m_list_edges[0] = compute_edge(1, 0);
    m_list_edges[1] = compute_edge(2, 0);
    m_list_edges[2] = compute_edge(3, 0);
    m_list_edges[3] = compute_edge(2, 1);
    m_list_edges[4] = compute_edge(3, 1);
    m_list_edges[5] = compute_edge(3, 2);
    m_signed_volume = compute_signed_volume();
    m_bbox          = compute_bounding_box();
    m_barycenter    = compute_barycenter();

    // const double vol_threshold = 1e-12;
    // if (m_signed_volume <= vol_threshold) {
    //     std::cerr << "Nul volume !! " << std::endl;
    // }
}

vector3 Tetra::compute_barycenter() const {
    return (m_list_vertices[0]->get_position() + m_list_vertices[1]->get_position() +
            m_list_vertices[2]->get_position() + m_list_vertices[3]->get_position()) /
           4.0;
}

/**
 * @brief Compute the gradient of the energy within the tetrahedron.
 *
 * @param values_at_vertices
 * @return vector3
 */
vector3 Tetra::compute_gradient_at_tetra(const array4d& values_at_vertices) const {
    const vector3 a = m_list_edges[0];
    const vector3 b = m_list_edges[1];
    const vector3 c = m_list_edges[2];

    const double du1 = values_at_vertices[1] - values_at_vertices[0];
    const double du2 = values_at_vertices[2] - values_at_vertices[0];
    const double du3 = values_at_vertices[3] - values_at_vertices[0];

    const double     det = dot(a, cross_product(b, c));
    constexpr double relative_tolerance = 1e-14;
    if (std::abs(det) <= relative_tolerance * geometry_scale(m_list_edges)) {
        throw std::domain_error("Cannot compute a gradient in a degenerate tetrahedron.");
    }
    return (cross_product(b, c) * du1 + cross_product(c, a) * du2 + cross_product(a, b) * du3) / det;
}

void Tetra::compute_gradient_energy_at_bands() {
    m_gradient_energy_per_band.clear();
    m_gradient_norm_per_band.clear();
    m_nb_bands = m_list_vertices[0]->get_number_bands();
    m_gradient_energy_per_band.reserve(m_nb_bands);
    m_gradient_norm_per_band.reserve(m_nb_bands);
    for (std::size_t band_index = 0; band_index < m_nb_bands; band_index++) {
        const std::array<double, 4> energies_at_vertices = get_band_energies_at_vertices(band_index);
        const vector3 gradient = compute_gradient_at_tetra(energies_at_vertices);
        m_gradient_energy_per_band.push_back(gradient);
        m_gradient_norm_per_band.push_back(gradient.norm());
    }
}

/**
 * @brief Compute the minimum and maximum energies of the bands within the tetrahedron.
 *  *
 */
void Tetra::compute_min_max_energies_at_bands() {
    m_min_energy_per_band.clear();
    m_max_energy_per_band.clear();
    m_nb_bands = m_list_vertices[0]->get_number_bands();
    m_min_energy_per_band.reserve(m_nb_bands);
    m_max_energy_per_band.reserve(m_nb_bands);
    for (std::size_t idx_band = 0; idx_band < m_nb_bands; ++idx_band) {
        const std::array<double, 4> energies = {
            m_list_vertices[0]->get_energy_at_band(idx_band),
            m_list_vertices[1]->get_energy_at_band(idx_band),
            m_list_vertices[2]->get_energy_at_band(idx_band),
            m_list_vertices[3]->get_energy_at_band(idx_band),
        };
        auto minmax   = std::minmax_element(energies.begin(), energies.end());
        m_min_energy_per_band.push_back(*minmax.first);
        m_max_energy_per_band.push_back(*minmax.second);
    }
}

/**
 * @brief Compute the signed volume of the tetrahedron.
 *
 * @return double
 */
double Tetra::compute_signed_volume() const {
    // std::cout << m_list_edges[0] << std::endl;
    // std::cout << m_list_edges[1] << std::endl;
    return (1.0 / 6.0) * scalar_triple_product(m_list_edges[0], m_list_edges[1], m_list_edges[2]);
}

/**
 * @brief Return the values of the energy of the index_band valence band at the 4 vertices of the tetrahedra.
 *
 * @param index_band
 * @return std::vector<double>
 */
std::array<double, 4> Tetra::get_band_energies_at_vertices(std::size_t index_band) const {
    return {m_list_vertices[0]->get_energy_at_band(index_band),
            m_list_vertices[1]->get_energy_at_band(index_band),
            m_list_vertices[2]->get_energy_at_band(index_band),
            m_list_vertices[3]->get_energy_at_band(index_band)};
}

/**
 * @brief Compute the edge vector between two vertices of the tetrahedra.
 * The result is: vtx_1 - vtx_2.
 *
 * @param index_vtx_1
 * @param index_vtx_2
 * @return vector3
 */
vector3 Tetra::compute_edge(std::size_t index_vtx_1, std::size_t index_vtx_2) const {
    if (index_vtx_1 > 3 || index_vtx_2 > 3) {
        throw std::invalid_argument("In Tetra::compute_edge, the index of vertex must be between 0 and 3.");
    }
    return m_list_vertices[index_vtx_1]->get_position() - m_list_vertices[index_vtx_2]->get_position();
}

/**
 * @brief Compute the barycentric coordinate of a given location within the tetrahedra.
 * The returned array of size 4 contains the barycentric coordinates with respect to the vertices in the following order
 * : 0, 1, 2 and 3, respectively.
 *
 * @warning warning message: Do not use this function to check if the location lies in the tetrahedra,
 * The computation relies on the hypothesis that the location do lies in it. Use Tetra::is_location_inside instead.
 *
 * @param location
 * @return std::array<double, 4>
 */
std::array<double, 4> Tetra::compute_barycentric_coordinates(const vector3& location) const {
    const vector3 v_loc1            = location - m_list_vertices[0]->get_position();
    const double  tetra_determinant = 6.0 * m_signed_volume;
    constexpr double relative_tolerance = 1e-14;
    if (std::abs(tetra_determinant) <= relative_tolerance * geometry_scale(m_list_edges)) {
        throw std::domain_error("Cannot compute barycentric coordinates in a degenerate tetrahedron.");
    }

    const double lambda_2 = scalar_triple_product(v_loc1, m_list_edges[1], m_list_edges[2]) / tetra_determinant;
    const double lambda_3 = scalar_triple_product(v_loc1, m_list_edges[2], m_list_edges[0]) / tetra_determinant;
    const double lambda_4 = scalar_triple_product(v_loc1, m_list_edges[0], m_list_edges[1]) / tetra_determinant;
    const double lambda_1 = 1.0 - lambda_2 - lambda_3 - lambda_4;

    return {lambda_1, lambda_2, lambda_3, lambda_4};
}

double Tetra::interpolate_scalar_at_position(const vector3& location, const std::vector<double>& scalar_field) const {
    const auto barycentric_coord = compute_barycentric_coordinates(location);
    return scalar_field[0] * barycentric_coord[0] + scalar_field[1] * barycentric_coord[1] +
           scalar_field[2] * barycentric_coord[2] + scalar_field[3] * barycentric_coord[3];
}

/**
 * @brief Compute the linear interpolation of the energy of the band band_index at the point location.
 *
 * @param location
 * @param band_index
 * @return double
 */
double Tetra::interpolate_energy_at_band(const vector3& location, std::size_t band_index) const {
    const auto                  barycentric_coord    = compute_barycentric_coordinates(location);
    const std::array<double, 4> energies_at_vertices = get_band_energies_at_vertices(band_index);
    return energies_at_vertices[0] * barycentric_coord[0] + energies_at_vertices[1] * barycentric_coord[1] +
           energies_at_vertices[2] * barycentric_coord[2] + energies_at_vertices[3] * barycentric_coord[3];
}

vector3 Tetra::interpolate_gradient_energy_at_band(const vector3& location, std::size_t band_index) const {
    const auto barycentric_coord = compute_barycentric_coordinates(location);

    vector3 gradient_at_location{0.0, 0.0, 0.0};
    for (int i = 0; i < 4; ++i) {
        gradient_at_location += m_list_vertices[i]->get_energy_gradient_at_band(band_index) * barycentric_coord[i];
    }
    // DEBUG
    if (std::isnan(gradient_at_location.x()) || std::isnan(gradient_at_location.y()) ||
        std::isnan(gradient_at_location.z())) {
        std::cerr << "Warning: NaN gradient at location " << location << " for band " << band_index
                  << ". This may indicate a van Hove singularity or insufficient mesh resolution." << std::endl;
        for (int i = 0; i < 4; ++i) {
            std::cerr << "Volume: " << m_signed_volume << std::endl;
            std::cerr << "  Vertex " << i << ": k = " << m_list_vertices[i]->get_position()
                      << ", energy = " << m_list_vertices[i]->get_energy_at_band(band_index)
                      << ", gradient = " << m_list_vertices[i]->get_energy_gradient_at_band(band_index)
                      << " barycentric coord = " << barycentric_coord[i] << std::endl;
        }
    }
    return gradient_at_location;
}

/**
 * @brief Check if a given location lies inside the tetrahedra.
 *
 * @param location
 * @return true
 * @return false
 */
bool Tetra::is_location_inside(const vector3& location) const {
    const vector3    v_loc1            = location - m_list_vertices[0]->get_position();
    const vector3    v_loc2            = location - m_list_vertices[1]->get_position();
    const double     tetra_determinant = 6.0 * m_signed_volume;
    constexpr double relative_tolerance = 1e-14;
    if (std::abs(tetra_determinant) <= relative_tolerance * geometry_scale(m_list_edges)) {
        return false;
    }
    const double     lambda_1 = scalar_triple_product(v_loc2, m_list_edges[4], m_list_edges[3]) / tetra_determinant;
    const double     lambda_2 = scalar_triple_product(v_loc1, m_list_edges[1], m_list_edges[2]) / tetra_determinant;
    const double     lambda_3 = scalar_triple_product(v_loc1, m_list_edges[2], m_list_edges[0]) / tetra_determinant;
    const double     lambda_4 = scalar_triple_product(v_loc1, m_list_edges[0], m_list_edges[1]) / tetra_determinant;
    constexpr double barycentric_tolerance = 1e-12;
    return (lambda_1 >= -barycentric_tolerance && lambda_2 >= -barycentric_tolerance &&
            lambda_3 >= -barycentric_tolerance && lambda_4 >= -barycentric_tolerance);
}

/**
 * @brief Compute the euclidean position from barycentric coordinates.
 *
 * @param barycentric_coordinates
 * @return vector3
 */
vector3 Tetra::compute_euclidean_coordinates(const std::array<double, 4>& barycentric_coordinates) const {
    return (barycentric_coordinates[0] * m_list_vertices[0]->get_position() +
            barycentric_coordinates[1] * m_list_vertices[1]->get_position() +
            barycentric_coordinates[2] * m_list_vertices[2]->get_position() +
            barycentric_coordinates[3] * m_list_vertices[3]->get_position());
}

/**
 * @brief Compute the euclidean position from barycentric coordinates, with a given vertices order,
 * that might be different from the vertices of the tetrahedra.
 *
 * @param barycentric_coordinates
 * @param indices_vertex
 * @return vector3
 */
vector3 Tetra::compute_euclidean_coordinates_with_indices(const std::array<double, 4>& barycentric_coordinates,
                                                          const std::array<int, 4>&    indices_vertex) const {
    return (barycentric_coordinates[0] * m_list_vertices[indices_vertex[0]]->get_position() +
            barycentric_coordinates[1] * m_list_vertices[indices_vertex[1]]->get_position() +
            barycentric_coordinates[2] * m_list_vertices[indices_vertex[2]]->get_position() +
            barycentric_coordinates[3] * m_list_vertices[indices_vertex[3]]->get_position());
}

/**
 * @brief Precompute a list of indices a, b, c, d such as, for the conduction band with index index_band,
 * we have Vtx_a <= Vtx_b <= Vtx_c <= Vtx_d in term of energy.
 *
 * This function is written explicitely instead of using std::sort functions, because the sorting is done
 * with the minimum number of operations for a 4 values sorting. Other solution might be tested later.
 *
 * @param index_band
 * @return std::array<int, 4>
 */
void Tetra::pre_compute_sorted_slots_per_band() {
    m_sorted_slots_per_band.clear();
    m_sorted_slots_per_band.reserve(m_nb_bands);
    for (std::size_t band_index = 0; band_index < m_nb_bands; band_index++) {
        std::array<double, 4> energies_at_vertices = get_band_energies_at_vertices(band_index);
        std::array<int, 4>    sorted_index         = {0, 1, 2, 3};
        if (energies_at_vertices[0] > energies_at_vertices[1]) {
            std::swap(energies_at_vertices[0], energies_at_vertices[1]);
            std::swap(sorted_index[0], sorted_index[1]);
        }
        if (energies_at_vertices[2] > energies_at_vertices[3]) {
            std::swap(energies_at_vertices[2], energies_at_vertices[3]);
            std::swap(sorted_index[2], sorted_index[3]);
        }
        if (energies_at_vertices[0] > energies_at_vertices[2]) {
            std::swap(energies_at_vertices[0], energies_at_vertices[2]);
            std::swap(sorted_index[0], sorted_index[2]);
        }
        if (energies_at_vertices[1] > energies_at_vertices[3]) {
            std::swap(energies_at_vertices[1], energies_at_vertices[3]);
            std::swap(sorted_index[1], sorted_index[3]);
        }
        if (energies_at_vertices[1] > energies_at_vertices[2]) {
            std::swap(energies_at_vertices[1], energies_at_vertices[2]);
            std::swap(sorted_index[1], sorted_index[2]);
        }
        m_sorted_slots_per_band.push_back(sorted_index);
    }
}

/**
 * @brief Compute the iso-energy surface within the tetrahedra for a given energy of a given band.
 * The surface is returned as a list of points (3 when the surface is a triangle, 4 when it is a quadrangle).
 *
 * The case of energy being smaller than the minimum energy of the tetrahedra is not taken into account.
 * Same thing for the case of energy being greater than the maximum energy of the tetrahedra.
 * Those two cases are handle by the caller function. This is done to avoid computing the sorted index which is
 * computationallly intensive. The minimum and maximum energies are stored in the member variables
 * m_min_energy_at_vertices and m_max_energy_at_vertices at the construction of the tetrahedra.
 *
 * This is very important because those 2 trivial cases represent usually more than 95% of the cases.
 *
 *
 * @param iso_energy
 * @param band_index
 * @return std::vector<vector3>
 */
std::vector<vector3> Tetra::compute_band_iso_energy_surface(double iso_energy, std::size_t band_index) const {
    if (band_index >= m_nb_bands) {
        throw std::out_of_range("Band index out of range in tetrahedron iso-surface computation.");
    }

    const auto energies = get_band_energies_at_vertices(band_index);
    const auto minmax   = std::minmax_element(energies.begin(), energies.end());
    const double energy_scale = std::max({std::abs(*minmax.first), std::abs(*minmax.second), 1.0});
    const double energy_tolerance = 1e-12 * energy_scale;
    if (iso_energy < *minmax.first - energy_tolerance || iso_energy > *minmax.second + energy_tolerance) {
        return {};
    }
    if (*minmax.second - *minmax.first <= energy_tolerance) {
        return {};
    }

    constexpr std::array<std::array<std::size_t, 2>, 6> edge_vertices = {
        {{{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}, {{1, 3}}, {{2, 3}}}};
    const double point_tolerance = 1e-12 * std::max(m_bbox.get_diagonal_size(), 1.0);
    std::vector<vector3> intersections;
    intersections.reserve(4);

    for (const auto& edge : edge_vertices) {
        const std::size_t i = edge[0];
        const std::size_t j = edge[1];
        const double di = energies[i] - iso_energy;
        const double dj = energies[j] - iso_energy;
        const bool i_on = std::abs(di) <= energy_tolerance;
        const bool j_on = std::abs(dj) <= energy_tolerance;

        if (i_on) {
            add_unique_point(intersections, m_list_vertices[i]->get_position(), point_tolerance);
        }
        if (j_on) {
            add_unique_point(intersections, m_list_vertices[j]->get_position(), point_tolerance);
        }
        if (!i_on && !j_on && ((di < 0.0) != (dj < 0.0))) {
            const double fraction = (iso_energy - energies[i]) / (energies[j] - energies[i]);
            const vector3 point = (1.0 - fraction) * m_list_vertices[i]->get_position() +
                                  fraction * m_list_vertices[j]->get_position();
            add_unique_point(intersections, point, point_tolerance);
        }
    }

    return intersections;
}

// Area of triangle (A,B,C) in 3D: 0.5 * || (B-A) × (C-A) ||
inline double triangle_area(const vector3& A, const vector3& B, const vector3& C) noexcept {
    const vector3 AB = B - A;
    const vector3 AC = C - A;
    return 0.5 * cross_product(AB, AC).norm();
}

// Returns vertices ordered cyclically in the plane they lie on
inline std::vector<vector3> order_cyclic(const std::vector<vector3>& pts) {
    if (pts.size() < 3) {
        return pts;
    }

    // Compute centroid
    vector3 centroid = std::accumulate(pts.begin(), pts.end(), vector3{0, 0, 0});
    centroid /= static_cast<double>(pts.size());

    // Compute polygon normal from first 3 distinct points
    vector3 n   = cross_product(pts[1] - pts[0], pts[2] - pts[0]);
    double  len = n.norm();
    if (len > 0.0) {
        n /= len;  // normalize
    }

    // Choose an in-plane axis u
    vector3 u    = pts[0] - centroid;
    double  ulen = u.norm();
    if (ulen > 0.0) {
        u /= ulen;
    } else {
        u = vector3{1, 0, 0};  // fallback
    }

    // v = n × u (second in-plane axis)
    vector3 v = cross_product(n, u);

    struct VertexAngle {
        vector3 p;
        double  angle;
    };
    std::vector<VertexAngle> with_angles;
    with_angles.reserve(pts.size());

    for (auto& p : pts) {
        vector3 d     = p - centroid;
        double  x     = dot(d, u);
        double  y     = dot(d, v);
        double  angle = std::atan2(y, x);  // -pi .. pi
        with_angles.push_back({p, angle});
    }

    std::sort(with_angles.begin(), with_angles.end(), [](auto& a, auto& b) { return a.angle < b.angle; });

    std::vector<vector3> ordered;
    ordered.reserve(pts.size());
    for (auto& wa : with_angles) {
        ordered.push_back(wa.p);
    }

    return ordered;
}

/**
 * @brief Compute the area of a polygon defined by its vertices.
 *
 * @param pts
 * @return double
 */
inline double polygon_area(const std::vector<vector3>& pts) {
    auto pts_ordered = order_cyclic(pts);
    if (pts_ordered.size() < 3) {
        return 0.0;
    }

    vector3 sum{0, 0, 0};
    for (size_t i = 0; i < pts_ordered.size(); ++i) {
        const vector3& p0 = pts_ordered[i];
        const vector3& p1 = pts_ordered[(i + 1) % pts_ordered.size()];
        sum += cross_product(p0, p1);
    }
    return 0.5 * sum.norm();
}

double Tetra::compute_tetra_dos_energy_band_reference(double energy_eV, std::size_t band_index) const {
    if (band_index >= m_nb_bands) {
        throw std::out_of_range("Band index out of range in tetrahedron DOS computation.");
    }
    // Early-out if E outside tetra range for this band
    if (energy_eV < m_min_energy_per_band[band_index] || energy_eV > m_max_energy_per_band[band_index]) {
        return 0.0;
    }

    // Intersect isosurface E(k)=energy with tetra edges -> polygon vertices in k (m^-1)
    const vector3 gradient = compute_gradient_at_tetra(get_band_energies_at_vertices(band_index));
    const double gradient_norm = gradient.norm();
    if (!(gradient_norm > 0.0) || !std::isfinite(gradient_norm)) {
        return 0.0;
    }

    std::vector<vector3> iso = compute_band_iso_energy_surface(energy_eV, band_index);
    if (iso.size() < 3) {
        return 0.0;  // no area
    }

    // Surface area in k-space (m^-2)
    const double A = polygon_area(iso);
    if (!(A > 0.0)) {
        return 0.0;
    }

    // Prefactor 1/(2π)^3
    constexpr double pref = 1.0 / (8.0 * uepm::constants::pi * uepm::constants::pi * uepm::constants::pi);

    return pref * A / gradient_norm;  // states / (eV · m^3), without spin degeneracy
}

double Tetra::compute_tetra_dos_energy_band(double energy_eV, std::size_t band_index) const {
    if (band_index >= m_nb_bands) {
        throw std::out_of_range("Band index out of range in tetrahedron DOS computation.");
    }
    if (energy_eV < m_min_energy_per_band[band_index] || energy_eV > m_max_energy_per_band[band_index]) {
        return 0.0;
    }

    IsoEnergyPolygon polygon = compute_band_iso_energy_polygon(energy_eV, band_index);
    return compute_tetra_dos_energy_band(energy_eV, band_index, polygon);
}

double Tetra::compute_tetra_dos_energy_band(double                  energy_eV,
                                            std::size_t             band_index,
                                            const IsoEnergyPolygon& polygon) const {
    if (band_index >= m_nb_bands) {
        throw std::out_of_range("Band index out of range in tetrahedron DOS computation.");
    }
    if (energy_eV < m_min_energy_per_band[band_index] || energy_eV > m_max_energy_per_band[band_index]) {
        return 0.0;
    }
    const double gradient_norm =
        band_index < m_gradient_norm_per_band.size()
            ? m_gradient_norm_per_band[band_index]
            : compute_gradient_at_tetra(get_band_energies_at_vertices(band_index)).norm();
    if (!(gradient_norm > 0.0) || !std::isfinite(gradient_norm)) {
        return 0.0;
    }
    if (!(polygon.area > 0.0)) {
        return 0.0;
    }

    constexpr double pref = 1.0 / (8.0 * uepm::constants::pi * uepm::constants::pi * uepm::constants::pi);
    return pref * polygon.area / gradient_norm;
}

IsoEnergyPolygon Tetra::compute_band_iso_energy_polygon(double iso_energy, std::size_t band_index) const {
    if (band_index >= m_nb_bands) {
        throw std::out_of_range("Band index out of range in tetrahedron iso-energy computation.");
    }
    IsoEnergyPolygon polygon = compute_band_iso_energy_polygon_impl(
        m_list_vertices, m_bbox, get_band_energies_at_vertices(band_index), iso_energy);
    polygon.area = order_iso_polygon_and_compute_area(polygon);
    return polygon;
}

void Tetra::precompute_dos_on_energy_grid_per_band(double energy_step, double energy_max) {
    (void)energy_step;
    (void)energy_max;
    // m_dos_per_band.clear();
    // m_nb_bands = m_list_vertices[0]->get_number_bands();
    // m_dos_per_band.assign(m_nb_bands, UniformDos{});

    // for (std::size_t b = 0; b < m_nb_bands; ++b) {
    //     auto&        T    = m_dos_per_band[b];
    //     const double Emin = m_min_energy_per_band[b];
    //     const double Emax = m_max_energy_per_band[b];

    //     // keep index alignment; mark invalid instead of skipping
    //     if (Emin > energy_max) {
    //         T.valid = false;
    //         continue;
    //     }

    //     // integer number of steps; enforce a small minimum
    //     std::size_t           nb_steps  = static_cast<std::size_t>(std::ceil((Emax - Emin) / energy_step));
    //     constexpr std::size_t min_steps = 5;
    //     if (nb_steps < min_steps) {
    //         nb_steps = min_steps;
    //     }
    //     const double dx = (Emax - Emin) / static_cast<double>(nb_steps);

    //     T.valid  = true;
    //     T.E0     = Emin;
    //     T.Emax   = Emax;
    //     T.inv_dx = 1.0 / dx;
    //     T.N      = static_cast<uint32_t>(nb_steps + 1);

    //     T.D.resize(T.N);
    //     T.D[0]       = 0.0;
    //     T.D[T.N - 1] = 0.0;

    //     for (std::size_t idx_energy = 1; idx_energy < nb_steps; ++idx_energy) {
    //         const double e  = Emin + idx_energy * dx;
    //         T.D[idx_energy] = static_cast<float>(compute_tetra_dos_energy_band(e, b));
    //     }
    // }
}

/**
 * @brief Interpolate the density of states (DOS) at a given energy for a specific band.
 * The DOS is precomputed on a uniform energy grid during the tetrahedron initialization.
 * If the energy is outside the precomputed range, the function returns 0.0.
 * The interpolation is linear between the two nearest grid points.
 *
 * @param energy
 * @param band_index
 * @return double
 */
double Tetra::interpolate_dos_at_energy_per_band(double energy, std::size_t band_index) const {
    // The precomputed cache is not implemented yet. Falling back to the exact
    // tetrahedral expression preserves the physics instead of silently returning zero.
    return compute_tetra_dos_energy_band(energy, band_index);
}

/**
 * @brief Draw a random point on the iso-energy surface within the tetrahedra.
 *
 * @param iso_energy
 * @param band_index
 * @param rng
 * @return vector3
 */
vector3 Tetra::draw_random_uniform_point_at_energy(double iso_energy, std::size_t band_index, std::mt19937& rng) const {
    // std::cout << "Draw random point at energy: " << iso_energy << " in band: " << band_index << std::endl;
    if (iso_energy < m_min_energy_per_band[band_index] || iso_energy > m_max_energy_per_band[band_index]) {
        std::cout << "Band index: " << band_index << std::endl;
        std::cout << "Energie bound: " << m_min_energy_per_band[band_index] << " " << m_max_energy_per_band[band_index]
                  << std::endl;
        std::cout << "iso_energy: " << iso_energy << std::endl;
        throw std::invalid_argument(
            "Energy is not in the band for this tetrahedron. Cannot draw a random point at this energy.");
    }
    IsoEnergyPolygon polygon = compute_band_iso_energy_polygon(iso_energy, band_index);
    return draw_random_uniform_point_at_energy(polygon, rng);
}

vector3 Tetra::draw_random_uniform_point_at_energy(const IsoEnergyPolygon& polygon, std::mt19937& rng) const {
    if (polygon.size == 0) {
        throw std::invalid_argument(
            "Energy is not in the band for this tetrahedron. Cannot draw a random point at this energy.");
    } else if (polygon.size == 3) {
        IsoTriangle triangle(polygon.points[0], polygon.points[1], polygon.points[2], 0.0);
        auto        point = triangle.draw_random_uniform_point_in_triangle(rng);
        return point;
    } else {
        // If the iso-energy shape is a quadrilateral, the point is drawn uniformly in the quadrilateral.
        // To do so, we randomly select on of the triangle, with a probability following the area of the triangle.
        // Then we draw a point in the selected triangle, and return the point.
        if (polygon.size != 4) {
            throw std::runtime_error("A linear tetrahedron iso-energy surface must be a triangle or quadrilateral.");
        }
        IsoTriangle triangle1(polygon.points[0], polygon.points[1], polygon.points[2], 0.0);
        IsoTriangle triangle2(polygon.points[0], polygon.points[2], polygon.points[3], 0.0);
        const double surface_triangle1 = triangle1.get_signed_surface();
        const double surface_triangle2 = triangle2.get_signed_surface();
        std::uniform_real_distribution<double> dist(0.0, surface_triangle1 + surface_triangle2);
        return dist(rng) < surface_triangle1 ? triangle1.draw_random_uniform_point_in_triangle(rng)
                                             : triangle2.draw_random_uniform_point_in_triangle(rng);
    }
}

bool Tetra::is_energy_inside_band(double energy, std::size_t index_band) const {
    return (energy >= m_min_energy_per_band[index_band] && energy <= m_max_energy_per_band[index_band]);
}

bool Tetra::does_intersect_band_energy_range(double e_min, double e_max, std::size_t index_band) const {
    return !(e_max < m_min_energy_per_band[index_band] || e_min > m_max_energy_per_band[index_band]);
}

std::array<double, 8> Tetra::get_tetra_electron_phonon_rates(int band_index) const {
    std::array<double, 8> mean_rates;
    std::fill(mean_rates.begin(), mean_rates.end(), 0.0);
    for (std::size_t i = 0; i < 4; i++) {
        if (band_index >= static_cast<int>(m_list_vertices[i]->get_number_bands())) {
            throw std::invalid_argument("In Tetra::get_tetra_electron_phonon_rates, the band index is out of range.");
        }
        const std::array<double, 8>& rates = m_list_vertices[i]->get_electron_phonon_rates(band_index);
        std::transform(mean_rates.begin(), mean_rates.end(), rates.begin(), mean_rates.begin(), std::plus<double>());
    }
    std::transform(mean_rates.begin(), mean_rates.end(), mean_rates.begin(), [](double val) { return val / 4.0; });
    return mean_rates;
}

std::array<double, 8> Tetra::interpolate_phonon_scattering_rate_at_location(const vector3&     location,
                                                                            const std::size_t& band_index) const {
    if (is_location_inside(location) == false) {
        throw std::invalid_argument(
            "In Tetra::interpolate_phonon_scattering_rate_at_location, the location is not inside the tetrahedra.");
    }
    const auto            barycentric_coord = compute_barycentric_coordinates(location);
    std::array<double, 8> interpolated_rates;
    std::fill(interpolated_rates.begin(), interpolated_rates.end(), 0.0);
    for (std::size_t i = 0; i < 4; i++) {
        const std::array<double, 8>& rates = m_list_vertices[i]->get_electron_phonon_rates(band_index);
        for (std::size_t j = 0; j < 8; j++) {
            interpolated_rates[j] += rates[j] * barycentric_coord[i];
        }
    }
    return interpolated_rates;
}

}  // namespace uepm::mesh_bz
