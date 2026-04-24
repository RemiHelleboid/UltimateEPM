/**
 * @file mesh.hpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief Header file for mesh class implementation.
 * @version 0.1
 * @date 2021-07-10
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <gmsh.h>

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "bbox.hpp"
#include "dataset.hpp"
#include "element.hpp"
#include "element1d.hpp"
#include "element2d.hpp"
#include "element3d.hpp"
#include "function.hpp"
#include "quadrangular_interpolators.hpp"
#include "region.hpp"
#include "tree_node.hpp"
#include "utils_mesh.hpp"
#include "vector.hpp"
#include "vertex.hpp"

namespace uepm {

namespace mesh {

/**
 * @class mesh class
 * @brief Mesh class.
 *
 * The class contains all the regions and geometics objects that defined it.
 * It also contained the diffrents dataset with scalar or vector field.
 *
 */

class writer;
class grid_data;

using sp_element = std::shared_ptr<element>;

using sp_scalar_dataset   = std::shared_ptr<dataset<double>>;
using sp_vector_dataset   = std::shared_ptr<dataset<vector3>>;
using list_scalar_dataset = std::vector<sp_scalar_dataset>;

using sp_scalar_function = std::shared_ptr<function<double>>;
using sp_vector_function = std::shared_ptr<function<vector3>>;

static constexpr std::array<double, 9> no_transformation_matrix = {1u, 0u, 0u, 0u, 1u, 0u, 0u, 0u, 1u};
static constexpr std::array<double, 3> no_translation_vector    = {0u, 0u, 0u};

class mesh {
 private:
    unsigned int m_dimension   = 0;
    unsigned int m_nb_regions  = 0;
    unsigned int m_nb_vertices = 0;

    std::vector<vertex> m_ListVertices;

    std::unique_ptr<tree_node> m_p_search_tree;

    std::array<double, 9> m_transformation_matrix = no_transformation_matrix;
    std::array<double, 3> m_translation_vector    = no_translation_vector;

    std::vector<region_bulk>      m_ListRegionsBulk;
    std::vector<region_interface> m_ListRegionsInterface;
    std::vector<region_contact>   m_ListRegionsContact;

    std::vector<sp_scalar_function> m_list_scalar_functions;
    std::vector<sp_vector_function> m_list_vector_functions;

 public:
    friend class writer;

    mesh() = default;

    //   Mesh
    void                  reset_all();
    unsigned int          get_dimension() const { return m_dimension; }
    std::array<double, 9> get_tranformation_matrix() const { return (m_transformation_matrix); }
    std::array<double, 3> get_translation_vector() const { return (m_translation_vector); }
    void                  set_dimension(unsigned int dim) { m_dimension = dim; }
    void                  set_transformation_matrix(std::array<double, 9> transform_mat) { m_transformation_matrix = transform_mat; }
    void                  set_translation_vector(std::array<double, 3> translat_vector) { m_translation_vector = translat_vector; }
    bbox                  get_bounding_box() const;
    void                  build_search_tree();
    void                  reset_search_tree() { m_p_search_tree = nullptr; }
    bbox                  compute_bounding_box() const;

    mesh &&move() { return std::move(*this); }

    //   Vertices
    [[deprecated("Slow")]] vertex *get_p_vertex_safe(std::size_t index);
    unsigned int                   get_nb_vertices() const { return m_ListVertices.size(); }
    vertex                    *get_p_vertex(std::size_t index) { return index < m_ListVertices.size() ? &m_ListVertices[index] : nullptr; }
    const std::vector<vertex> &get_list_vertices() const { return m_ListVertices; }
    void                       add_vertex(double x, double y);
    void                       add_vertex(double x, double y, double z);
    void                       add_vertex(double x, double y, double z, std::size_t index_vertex);
    void                       set_nb_vertices(unsigned int nb_vertices) { m_nb_vertices = nb_vertices; }
    void                       sort_vertices();
    void                       transform_vertices(std::function<double(double, double)> function);
    void                       reset_all_space_charge_on_vertex();
    void                       update_charge_density_from_vertex_values();

    //   Regions
    unsigned int get_nb_regions() const { return (m_ListRegionsBulk.size() + m_ListRegionsContact.size() + m_ListRegionsInterface.size()); }
    region      *get_p_region(std::size_t index);
    const region      *get_p_region(const std::string &region_name) const;
    region            *get_p_region(const std::string &region_name);
    const region_bulk *get_p_region_bulk(std::size_t index) const;
    region_interface  *get_p_interface_region_from_bulk_indices(std::size_t index_bulk_1, std::size_t index_bulk_2);
    const region      &get_region_of_index(std::size_t search_index) const;
    void               compute_interface_region_betwwen_two_bulks(std::size_t index_bulk_1, std::size_t index_bulk_2);
    void               compute_interface_region_betwwen_two_bulks(std::size_t                     index_bulk_1,
                                                                  std::size_t                     index_bulk_2,
                                                                  const std::vector<std::size_t> &list_index_element_region_2_to_check);

    std::vector<int>                 get_vertices_belonging_to_other_regions(unsigned int idx_region);
    std::vector<region const *>      get_all_p_region() const;
    std::vector<region_bulk const *> get_all_p_bulk_region() const;
    void                             add_bulk_region(const region_bulk &new_region) { m_ListRegionsBulk.push_back(std::move(new_region)); }
    void                             add_contact_region(const region_contact &new_region) { m_ListRegionsContact.push_back(new_region); }
    std::vector<std::size_t>         get_idx_bulk_elements_adjacent_to_contact_region(const std::string &region_name) const;
    void add_interface_region(const region_interface &new_region) { m_ListRegionsInterface.push_back(new_region); }
    void set_nb_regions(unsigned int nb_region) { m_nb_regions = nb_region; }
    void remove_region(std::size_t index_region);
    const std::vector<region_bulk>      &get_list_bulk_region() const { return m_ListRegionsBulk; }
    const std::vector<region_interface> &get_list_interface_region() const { return m_ListRegionsInterface; }
    const std::vector<region_contact>   &get_list_contact_region() const { return m_ListRegionsContact; }

    //   Elements
    std::vector<sp_element> get_list_bulk_element() const;
    std::vector<element *>  get_list_p_bulk_element() const;
    void                    transfer_element_to_other_region(sp_element p_element, region *new_region);
    static void             transfer_element_to_other_region(sp_element p_element, region *origin_region, region *new_region);
    static void transfer_elements_to_other_region(std::vector<std::size_t> list_element_indexes, region *origin_region, region *new_region);

    void compare_bulks_region_with_bounding_box() const;

    //   Datasets
    std::vector<sp_scalar_dataset> get_list_scalar_datasets() const;
    std::vector<sp_vector_dataset> get_list_vector_datasets() const;
    dataset<double>               *get_p_scalar_dataset(const std::string &name, unsigned int region_index) const;
    dataset<double>               *get_p_vector_dataset(const std::string &name, unsigned int region_index) const;
    double                         get_scalar_data_at_vertex(std::size_t idx_vertex, const std::string &dataset_name) const;
    double                         get_scalar_data_at_element(std::size_t idx_element, const std::string &dataset_name) const;
    vector3                        get_vector_data_at_vertex(std::size_t idx_vertex, const std::string &dataset_name) const;
    std::vector<double>            get_all_scalar_dataset_values(const std::string &dataset_name) const;
    std::vector<vector3>           get_all_vector_dataset_values(const std::string &dataset_name) const;
    std::size_t                    get_number_scalar_datasets() const { return get_list_scalar_datasets().size(); }
    std::size_t                    get_number_vector_datasets() const { return get_list_vector_datasets().size(); }
    std::size_t                    get_total_number_dataset() const { return get_number_scalar_datasets() + get_number_vector_datasets(); };
    void                           create_scalar_datasets_from_idx_vertex_and_values(const std::string              &dataset_name,
                                                                                     DataLocationType                data_location_type,
                                                                                     const std::vector<std::size_t> &index_vertices,
                                                                                     const std::vector<double>      &data_values);
    void                           create_vector_datasets_from_idx_vertex_and_values(const std::string              &dataset_name,
                                                                                     DataLocationType                data_location_type,
                                                                                     const std::vector<std::size_t> &index_vertices,
                                                                                     const std::vector<vector3>     &data_values);
    void                           create_scalar_datasets_from_idx_cells_and_values(const std::string              &dataset_name,
                                                                                    DataLocationType                data_location_type,
                                                                                    const std::vector<std::size_t> &index_cells,
                                                                                    const std::vector<double>      &data_values);
    void                           add_scalar_dataset(sp_scalar_dataset new_dataset);
    void                           add_vector_dataset(sp_vector_dataset new_dataset);
    void                           add_scalar_data_to_vertices(const dataset<double> &dtset);
    void                           add_scalar_data_to_elements(const dataset<double> &dtset);
    void                           add_vector_data_to_vertices(const dataset<vector3> &dtset);
    void                           add_vector_data_to_elements(const dataset<vector3> &dtset);
    void                           add_scalar_data_to_all_vertices();
    void                           add_vector_data_to_all_vertices();
    void                           add_doping_concentration_to_vertices(const std::string &doping_fieldname);
    void                           add_electric_field_to_vertices(const std::string &electric_field_fieldname, double factor = 1.0);
    void add_diffusion_gradient_to_vertices(const std::string &e_gradient_fieldname, const std::string &h_gradient_fieldname);
    void add_space_charge_to_vertices(const std::string &space_charge_fieldname);
    void add_charge_density_to_vertices(const std::string &electron_charge_density_fieldname,
                                        const std::string &hole_charge_density_fieldname);
    void recompute_space_charge();
    void update_charge_density_from_elements_values(int nb_iter_windows);
    bool has_datasets() const { return !(get_list_scalar_datasets().empty() && get_list_scalar_datasets().empty()); }
    bool scalar_dataset_exists(const std::string &) const;
    bool vector_dataset_exists(const std::string &) const;
    void re_index_datasets();
    void remove_scalar_datasets(const std::string &name);
    void remove_vector_datasets(const std::string &name);
    void remove_all_scalar_function();
    void remove_all_vector_function();

    //   Functions
    std::vector<std::string>        get_all_functions_names() const;
    std::vector<sp_scalar_function> get_list_scalar_functions() const { return m_list_scalar_functions; }
    std::vector<std::string>        get_scalar_functions_name() const;
    std::vector<std::string>        get_vector_functions_name() const;
    std::vector<sp_vector_function> get_list_vector_functions() const { return m_list_vector_functions; }
    sp_scalar_function              get_sp_scalar_function(const std::string &name) const;
    sp_vector_function              get_sp_vector_function(const std::string &name) const;
    void                            remove_scalar_function(const std::string &name);
    void                            remove_vector_function(const std::string &name);
    bool                            scalar_function_exists(const std::string &name) const;
    bool                            vector_function_exists(const std::string &name) const;
    void                            add_scalar_function_to_vertices(const std::string &name);
    void                            add_vector_function_to_vertices(const std::string &name);
    void                            creates_all_scalar_functions_from_datasets();
    void                            creates_all_vector_functions_from_datasets();
    void                            creates_all_functions_from_datasets();

    void create_scalar_function_from_values_on_vertex(const std::string &name, const std::vector<double> &values);
    void create_vector_function_from_values_on_vertex(const std::string &name, const std::vector<vector3> &values);
    void create_scalar_function_from_values_on_element(const std::string &name, const std::vector<double> &values);
    void create_vector_function_from_values_on_element(const std::string &name, const std::vector<vector3> &values);

    void convert_element_function_to_vertex_function_scalar(const std::string &name_element_function,
                                                            const std::string &name_new_vertex_function);

    void create_null_scalar_function(const std::string &name);
    void create_null_vector_function(const std::string &name);
    void create_scalar_function_from_function(const std::string &new_name, const std::string &origin_function_name);
    void create_vector_function_from_function(const std::string &new_name, const function<vector3> &origin_function);
    void create_test_functions(const std::string &name_test_function);
    void create_test_function_on_elements(const std::string &name_test_function);
    void create_scalar_function_from_grid_data(const std::string &new_name,
                                               const grid_data   &my_grid_data,
                                               const std::string &interpolation_method = "linear");
    void create_scalar_function_from_csv_grid_data(const std::string &new_name,
                                                   const std::string &filename,
                                                   const std::string &interpolation_method = "linear");
    void create_density_function_from_list_positions_gaussian(const std::string          &new_fieldname,
                                                              const std::vector<vector3> &list_positions,
                                                              double                      gaussian_factor = 1.0);
    void create_density_function_from_list_positions_element_method(const std::string          &new_fieldname,
                                                                    const std::vector<vector3> &list_positions,
                                                                    double                      conversion_factor = 1.0);

    void addition_scalar_function_to_function(const std::string &name_function_to_add, const std::string &name_function_to_add_to);

    void multiplication_scalar_function_to_function(const std::string &name_new_function,
                                                    const std::string &name_function_1,
                                                    const std::string &name_function_2);

    void apply_scalar_function_to_function(const std::function<double(double)> &scalar_function,
                                           const std::string                   &name_function_to_apply_to);

    std::pair<std::vector<std::size_t>, std::vector<double>> get_vertices_index_value_of_scalar_function(
        const std::string &function_name) const;
    std::pair<std::vector<std::size_t>, std::vector<double>> get_cells_index_value_of_scalar_function(
        const std::string &function_name) const;
    std::pair<std::vector<std::size_t>, std::vector<vector3>> get_vertices_index_value_of_vector_function(
        const std::string &function_name) const;
    void create_gradient_function(const std::string &scalar_field, const std::string &new_name);

    // Poisson related functions
    void compute_n_charge_on_elements(const std::vector<vector3> &list_positions);
    void compute_p_charge_on_elements(const std::vector<vector3> &list_positions);
    void convert_charge_on_element_into_charge_at_vtx(double factor = 1.0);
    void create_space_charge_function_from_vtx_values(const std::string &new_fieldname);

    // Interpolation
    element                   *find_element_at_location(const vector3 &location) const;
    double                     interpolate_scalar_at_location(const std::string &fieldname, const vector3 &location) const;
    vector3                    interpolate_vector_at_location(const std::string &fieldname, const vector3 &location) const;
    vector3                    interpolate_gradient_at_location(const std::string &fieldname, const vector3 &location) const;
    const region_bulk         *get_p_region_at_location(const vector3 &location);
    std::string                get_region_name_at_location(const vector3 &location);
    std::string                get_material_name_at_location(const vector3 &location);
    std::string                get_material_name_at_element(element *element_location) const;
    std::pair<vector3, double> get_argmax_max_of_function(const std::string &fieldname) const;
    std::pair<vector3, double> get_argmin_min_of_function(const std::string &fieldname) const;
    void                       create_debye_length_function(const std::string &new_fieldname,
                                                            const std::string &donor_fieldname,
                                                            const std::string &acceptor_fieldname);

    // Find elements
    std::vector<element *> find_elements_overlapping_box(const bbox &my_box) const;

    // Intersection
    std::optional<vector3> find_location_line_boundaries_intersection(const vector3 &point_A, const vector3 &point_B) const;
    std::optional<vector3> find_location_line_region_interface_intersection(const vector3 &point_A, const vector3 &point_B) const;

    std::optional<std::pair<sp_element, vector3>> find_element_face_line_boundaries_intersection(const vector3 &point_A,
                                                                                                 const vector3 &point_B) const;

    std::optional<std::pair<sp_element, vector3>> find_element_face_line_region_interface_intersection(const vector3 &point_A,
                                                                                                       const vector3 &point_B) const;

    std::optional<std::pair<sp_element, vector3>> find_line_first_intersection(const vector3 &point_A, const vector3 &point_B) const;

    // Integration
    double integrate_over_mesh(const std::string &fieldname) const;
    double integrate_over_mesh_element_data(const std::string &fieldname) const;
    double integrate_over_region(const std::string &fieldname, const std::string &region_name) const;
    double integrate_over_material(const std::string &fieldname, const std::string &material_name) const;
    double compute_total_mesh_volume() const;
    double compute_region_volume(const std::string &region_name) const;
    double compute_material_volume(const std::string &material_name) const;

    //   I/O
    void print_regions_info() const;
    void print_datasets_info() const;
    void print_functions_info() const;
    void export_vertices_to_csv(const std::string &filename) const;
    void export_vertices_data_to_csv(const std::string &filename, const std::string &dataset_name) const;
    void export_all_vertices_data_to_csv(const std::string &filename) const;
    void export_on_grid(const std::string &filename, double dx, double dy, double dz = 1.0) const;

    void export_x_cut(const std::string &filename, double y_const, double z_const, double dx) const;
    void export_y_cut(const std::string &filename, double x_const, double z_const, double dy) const;
    void export_z_cut(const std::string &filename, double x_const, double y_const, double dz) const;

    std::vector<vector3> generate_mesh_grid(double dx, double dy, double dz) const;
    std::vector<vector3> generate_mesh_grid(int N_x, int N_y, int N_z) const;
    std::vector<vector3> generate_inner_mesh_grid_2d(int N_x, int N_y) const {
        return get_bounding_box().generate_inner_mesh_grid_2d(N_x, N_y);
    }
    std::vector<vector3> generate_inner_mesh_grid_3d(int N_x, int N_y, int N_z) const {
        return get_bounding_box().generate_inner_mesh_grid_3d(N_x, N_y, N_z);
    }

    // Other
    void mundfabisation(const std::string &mundfab_data_name, const double value_silicon, const std::string &silicon_region_name);
};

}  // namespace mesh

}  // namespace uepm