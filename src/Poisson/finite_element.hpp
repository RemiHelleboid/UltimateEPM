/**
 * @file poisson2d.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-09-30
 *
 * @copyright Copyright (c) 2021
 *
 */
#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

#include <plog/Log.h>

#include "materials.hpp"
#include "mesh.hpp"
#include "physical_constants.hpp"

namespace uepm::fem {

typedef Eigen::SparseMatrix<double> EigenSparseMatrix;
typedef Eigen::VectorXd             EigenVector;

enum class fem_status { None = -1, Assembled = -2, SolvedFailed = 1, SolvedSuccess = 0 };
inline std::string get_status_string(fem_status status);

class FiniteElementSystem {
 protected:
    mesh::mesh*                                                    m_p_mesh;
    EigenSparseMatrix                                              m_matrix_lhs;
    Eigen::SparseLU<EigenSparseMatrix, Eigen::COLAMDOrdering<int>> m_solver;
    Eigen::ConjugateGradient<EigenSparseMatrix, Eigen::Lower>      m_iter_solver;
    EigenVector                                                    m_second_member;
    EigenVector                                                    m_solution;
    fem_status                                                     m_status;

 public:
    static const double very_large_value;
    FiniteElementSystem()
        : m_p_mesh{nullptr},
          m_matrix_lhs{EigenSparseMatrix(0, 0)},
          m_second_member{EigenVector(0)},
          m_solution{EigenVector(0)},
          m_status{fem_status::None} {}

    FiniteElementSystem(std::size_t system_size)
        : m_p_mesh{nullptr},
          m_matrix_lhs{EigenSparseMatrix(system_size, system_size)},
          m_second_member{EigenVector(system_size)},
          m_solution{EigenVector(system_size)},
          m_status{fem_status::None} {}

    FiniteElementSystem(mesh::mesh* p_mesh, std::size_t system_size)
        : m_p_mesh{p_mesh},
          m_matrix_lhs{EigenSparseMatrix(system_size, system_size)},
          m_second_member{EigenVector(system_size)},
          m_solution{EigenVector(system_size)},
          m_status{fem_status::None} {}

    virtual ~FiniteElementSystem() = default;
    void reset_system();

    const EigenSparseMatrix& get_lhs_matrix() const { return m_matrix_lhs; }
    const EigenVector&       get_second_member() const { return m_second_member; }
    const EigenVector&       get_solution() const { return m_solution; }
    fem_status               get_status() const { return m_status; }
    std::size_t              get_system_size() const { return m_solution.size(); }

    virtual void compute_stiffness_matrix()                   = 0;
    virtual void compute_mass_matrix()                        = 0;
    virtual void compute_second_member(double constant_value) = 0;
    void         decompose_matrix();
    void         solve_system();

    double solver_L2_error() const;
    void   print_infos();
    void   export_solution_csv(const std::string& filename);
    void   add_second_member_to_mesh_functions(const std::string& function_name);
    void   add_solution_to_mesh_functions(const std::string& function_name, bool add_gradient = true);
    void   export_solution_to_vertices();
};

inline void print_time_different_units(double time_nanoseconds) {
    std::cout << "TIME TO RUN THE SOLVING ITERATIONS           : " << time_nanoseconds << " nanoseconds" << std::endl;
    std::cout << "TIME TO RUN THE SOLVING ITERATIONS           : " << time_nanoseconds / 1.0e3 << " microseconds" << std::endl;
    std::cout << "TIME TO RUN THE SOLVING ITERATIONS           : " << time_nanoseconds / 1.0e6 << " miliseconds" << std::endl;
    std::cout << "TIME TO RUN THE SOLVING ITERATIONS           : " << time_nanoseconds / 1.0e9 << " seconds" << std::endl;
}

}  // namespace uepm::fem