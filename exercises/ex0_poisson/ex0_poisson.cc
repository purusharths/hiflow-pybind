// Copyright (C) 2011-2021 Vincent Heuveline
//
// HiFlow3 is free software: you can redistribute it and/or modify it under the
// terms of the European Union Public Licence (EUPL) v1.2 as published by the
// European Union or (at your option) any later version.
//
// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for
// more details.
//
// You should have received a copy of the European Union Public Licence (EUPL)
// v1.2 along with HiFlow3.  If not, see
// <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

#include "ex0_poisson.h"

static const char *PARAM_FILENAME = "ex0_poisson.xml";
#ifndef MESHES_DATADIR
#define MESHES_DATADIR "./"
#endif
static const char *DATADIR = MESHES_DATADIR;

// Main application class ///////////////////////////////////

class PoissonTutorial {
public:
  PoissonTutorial(const std::string &param_filename,
                  const std::string &path_mesh)
      : path_mesh(path_mesh), comm_(MPI_COMM_WORLD), rank_(-1),
        num_partitions_(-1),
        params_(param_filename, MASTER_RANK, MPI_COMM_WORLD), 
        is_done_(false), refinement_level_(0){
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &num_partitions_);
    
    // Setup Parallel Output / Logging
    if (rank_ == 0) {
      INFO = true;
    }
    else {
      INFO = false;
    }
  }

  // Main algorithm

  void run() {
    
    build_initial_mesh(); // Construct / read in the initial mesh.
    
    Timer timer;
    
    // Initialize space and linear algebra.
    timer.start();
    LOG_INFO ("do", "Prepare System");  
    prepare_system();
    timer.stop();
    LOG_INFO("duration", timer.get_duration()); 
    timer.reset();

    // Compute the stiffness matrix and right-hand side.
    timer.start();
    LOG_INFO ("do", "Assemble System ");
    assemble_system();
    timer.stop();
    LOG_INFO("duration", timer.get_duration()); 
    timer.reset();

    // Solve the linear system.
    timer.start();    
    LOG_INFO ("do", "Solve System ");
    solve_system();
    timer.stop();
    LOG_INFO("duration", timer.get_duration()); 
    timer.reset();
    
    // Compute the error to the exact solution.
    // LOG_INFO ("do", "Compute Error ");
    // timer.start();  
    // compute_error();      
    // timer.stop();
    // LOG_INFO("duration", timer.get_duration()); 
    // timer.reset();

    // Visualize the solution and the errors.
    timer.start();    
    LOG_INFO ("do", "Visualize Solution ");  
    visualize();
    timer.stop();
    LOG_INFO("duration", timer.get_duration()); 
    timer.reset();
  }

  ~PoissonTutorial() {
  }

private:
  // Member functions

  // Read and distribute mesh.
  std::string path_mesh;
  void build_initial_mesh();
  void prepare_system();  // Setup space, linear algebra, and compute Dirichlet values.
  void assemble_system();  // Compute the matrix and rhs.
  void solve_system();  // Compute solution x.
  void compute_error();  // Compute errors compared to exact solution.
  void visualize();   // Visualize the results.
  // void adapt();   // Adapt the space (mesh and/or degree).

  // ----- member variables
  MPI_Comm comm_;   // MPI communicator.
  int rank_, num_partitions_;   // Local process rank and number of processes.
  PropertyTree params_;   // Parameter data read in from file.
  MeshPtr mesh_, mesh_without_ghost_, master_mesh_;   // Local mesh and mesh on master process.
  VectorSpace< DataType, DIM > space_; // Solution space.
  VectorType rhs_, sol_;  // Vectors for solution and load vector.
  MatrixType matrix_;   // System matrix.
  std::vector< DataType > L2_err_, H1_err_; // Vectors for error norms
  std::vector< int > dirichlet_dofs_; // Dof id:s for Dirichlet boundary conditions.
  std::vector< DataType > dirichlet_values_; // Dof values for Dirichlet boundary conditions.

  // Global assembler.
  StandardGlobalAssembler< DataType, DIM > global_asm_;
  
  // bool is_done_; // Flag for stopping adaptive loop.
  // int refinement_level_; // Current refinement level.
  
}; // end class PoissonTutorial

// Program entry point

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  // set default parameter file
  std::string param_filename(PARAM_FILENAME);
  std::string path_mesh;
  // if set take parameter file specified on console
  if (argc > 1) {
    param_filename = std::string(argv[1]);
  }
  // if set take mesh following path specified on console
  if (argc > 2) {
    path_mesh = std::string(argv[2]);
  }
  try {
    // Create log files for INFO and DEBUG output
    std::ofstream info_log("poisson_tutorial_info_log");
    std::ofstream warning_log("poisson_tutorial_warning_log");
    std::ofstream debug_log("poisson_tutorial_debug_log");
    std::ofstream error_log("poisson_tutorial_error_log");
        
    LogKeeper::get_log("info").set_target(&(std::cout));
    LogKeeper::get_log("debug").set_target(&(std::cout));
    LogKeeper::get_log("error").set_target(&(std::cout));
    LogKeeper::get_log("warning").set_target(&(std::cout));

    // Create application object and run it
    PoissonTutorial app(param_filename, path_mesh);
    app.run();

  } catch (std::exception &e) {
    std::cerr << "\nProgram ended with uncaught exception.\n";
    std::cerr << e.what() << "\n";
    return -1;
  }
  
#ifdef WITH_GPERF
  ProfilerStop();
#endif
  MPI_Finalize();
  return 0;
}

//////////////// PoissonTutorial implementation //////////////

void PoissonTutorial::build_initial_mesh() {
  mesh::IMPL impl = mesh::IMPL_DBVIEW;

  // Read in the mesh on the master process. The mesh is chosen according to the
  // DIM of the problem.
  if (rank_ == MASTER_RANK) 
  {
    std::string mesh_name;

    switch (DIM) 
    {
      case 2: 
      {
        mesh_name = params_["Mesh"]["Filename2"].get< std::string >("unit_square.inp");
        break;
      }
      case 3: 
      {
        mesh_name = params_["Mesh"]["Filename3"].get< std::string >("unit_cube.inp");
        break;
      }

      default:
        assert(0);
    }
    std::string mesh_filename;
    if (path_mesh.empty()) 
    {
      mesh_filename = std::string(DATADIR) + mesh_name;
    } 
    else 
    {
      mesh_filename = path_mesh + mesh_name;
    }

    std::vector< MasterSlave > period(0, MasterSlave(0., 0., 0., 0));
    // read the mesh
    master_mesh_ = read_mesh_from_file(mesh_filename, DIM, DIM, 0, impl, period);

    // Refine the mesh until the initial refinement level is reached.
    const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get< int >(3);

    if (initial_ref_lvl > 0) 
    {
      master_mesh_ = master_mesh_->refine_uniform_seq(initial_ref_lvl);
    }
    refinement_level_ = initial_ref_lvl;
  }

  MPI_Bcast(&refinement_level_, 1, MPI_INT, MASTER_RANK, comm_);

  // partition the mesh and distribute the subdomains across all processes
  int uniform_ref_steps;
  mesh_without_ghost_ = partition_and_distribute(master_mesh_, MASTER_RANK, comm_, uniform_ref_steps, impl);

  assert(mesh_without_ghost_ != 0);
  
  // compute ghost cells
  SharedVertexTable shared_verts;
  mesh_ = compute_ghost_cells(*mesh_without_ghost_, comm_, shared_verts, impl);

}

void PoissonTutorial::prepare_system() 
{
  // Assign degrees to each element.
  const int nb_fe_var = 1;
  
  const int fe_degree = params_["FESpace"]["FeDegree"].get< int >(1);
  std::vector< int > fe_params(nb_fe_var, fe_degree);
   
  std::vector< FEType > fe_ansatz (nb_fe_var, FE_TYPE_LAGRANGE);
  std::vector< bool > is_cg (nb_fe_var, true);
  
  // Initialize the VectorSpace object.
  space_.Init(*mesh_, fe_ansatz, is_cg, fe_params, hiflow::doffem::HIFLOW_CLASSIC);
  
  LOG_INFO("nb nonlin trafos", space_.fe_manager().nb_nonlinear_trafos());
  
  // Compute the matrix sparsity structure
  SparsityStructure sparsity;
  compute_sparsity_structure(space_, sparsity);

  // initialize matrix object
  matrix_.Init(comm_, space_.la_couplings());
  matrix_.InitStructure(sparsity);
  matrix_.Zeros();
  
  // initialize vector objects for solution coefficitons and 
  // right hand side
  sol_.Init(comm_, space_.la_couplings());
  rhs_.Init(comm_, space_.la_couplings());
  rhs_.InitStructure();
  sol_.InitStructure();  
  rhs_.Zeros();
  sol_.Zeros();

  // Compute Dirichlet BC dofs and values using known exact solution.
  dirichlet_dofs_.clear();
  dirichlet_values_.clear();

#ifdef SUBEX_D
  // **************************
  // TODO exercise d)

  // **************************
#else
  DirichletConstant bc_zero(0.);
  compute_dirichlet_dofs_and_values(bc_zero, space_, 0, dirichlet_dofs_, dirichlet_values_);
#endif

}

void PoissonTutorial::assemble_system() 
{
  // Assemble matrix and right-hand-side vector.
  LocalPoissonAssembler< ExactSol > local_asm;
  global_asm_.assemble_matrix(space_, local_asm, matrix_);
  global_asm_.assemble_vector(space_, local_asm, rhs_);

  // Correct Dirichlet dofs.
  if (!dirichlet_dofs_.empty()) 
  {
    matrix_.diagonalize_rows(vec2ptr(dirichlet_dofs_), dirichlet_dofs_.size(), 1.0);
    rhs_.SetValues(vec2ptr(dirichlet_dofs_), dirichlet_dofs_.size(), vec2ptr(dirichlet_values_));
    sol_.SetValues(vec2ptr(dirichlet_dofs_), dirichlet_dofs_.size(), vec2ptr(dirichlet_values_));
  }
}

void PoissonTutorial::solve_system() 
{
  LinearSolver< LAD > *solver_;
  LinearSolverFactory< LAD > SolFact;
  solver_ = SolFact.Get(params_["LinearSolver"]["Name"].get< std::string >("CG"))->params(params_["LinearSolver"]);
  solver_->SetupOperator(matrix_);
  solver_->Solve(rhs_, &sol_);
  
  // *************************************
  // TODO ex c): insert your solution here


  // *************************************
  
  delete solver_;
}

void PoissonTutorial::visualize() {
  
  // interpolate exact solution onto FE space 
  ExactSol exact_sol;
  VectorType exact_dofs;
  exact_dofs.CloneFromWithoutContent(sol_);
    
  FeInterNodal<DataType, DIM, ExactSol > fe_inter (space_, &exact_sol, 0);
  fe_inter.interpolate (exact_dofs);
    
  // Setup visualization object.
  int num_sub_intervals = 1;
  CellVisualization< DataType, DIM > visu(space_, num_sub_intervals);

  // collect cell-wise data from mesh object
  const int tdim = mesh_->tdim();
  const int num_cells = mesh_->num_entities(tdim);
  std::vector< DataType > remote_index(num_cells, 0);
  std::vector< DataType > sub_domain(num_cells, 0);
  std::vector< DataType > material_number(num_cells, 0);

  // loop through all cells in the mesh
  for (mesh::EntityIterator it = mesh_->begin(tdim); it != mesh_->end(tdim); ++it) 
  {
    int temp1, temp2;
    const int cell_index = it->index();
    if (DIM > 1) 
    {
      mesh_->get_attribute_value("_remote_index_", tdim, cell_index, &temp1);
      mesh_->get_attribute_value("_sub_domain_", tdim, cell_index, &temp2);
      remote_index.at(cell_index) = temp1;
      sub_domain.at(cell_index) = temp2;
    }
    material_number.at(cell_index) = mesh_->get_material_number(tdim, cell_index);
  }

  // visualize finite element function corresponding to 
  // coefficient vector sol_
  visu.visualize(sol_, 0, "u");

  // visualize finite element function corresponding to 
  // coefficient vector exact_dofs
  visu.visualize(exact_dofs, 0, "exact_u");
  
  // visualize error measures
  visu.visualize_cell_data(L2_err_, "L2");
  visu.visualize_cell_data(H1_err_, "H1");
  
  // visualize some mesh data
  visu.visualize_cell_data(remote_index, "_remote_index_");
  visu.visualize_cell_data(sub_domain, "_sub_domain_");
  visu.visualize_cell_data(material_number, "Material Id");
  
  // write out data
  std::stringstream name;
  name << "ex0_solution" << refinement_level_;
  
  VTKWriter< DataType, DIM> vtk_writer (visu, this->comm_, MASTER_RANK);
  vtk_writer.write(name.str());    
}


// not used for now
void PoissonTutorial::compute_error() 
{
  // prepare sol_ for post processing
  sol_.Update();

  L2_err_.clear();
  H1_err_.clear();

  // Compute square of the L2 error on each element, putting the
  // values into L2_err_.
  L2ErrorIntegrator< ExactSol > L2_int(sol_);
  global_asm_.assemble_scalar(space_, L2_int, L2_err_);

  // Create attribute with L2 error for output.
  AttributePtr L2_err_attr(new DoubleAttribute(L2_err_));
  mesh_->add_attribute("L2 error", DIM, L2_err_attr);
  DataType total_L2_err = std::accumulate(L2_err_.begin(), L2_err_.end(), 0.);
  DataType global_L2_err = 0.;
  MPI_Reduce(&total_L2_err, &global_L2_err, 1, MPI_DOUBLE, MPI_SUM, MASTER_RANK, comm_);
  PLOG_INFO(rank_, "error", "Local L2 error on partition " << rank_ << " = "
                                                   << std::sqrt(total_L2_err));

  // Compute square of the H1 error on each element, putting the
  // values into H1_err_.

  H1ErrorIntegrator< ExactSol > H1_int(sol_);
  global_asm_.assemble_scalar(space_, H1_int, H1_err_);

  // Create attribute with H1 error for output.
  AttributePtr H1_err_attr(new DoubleAttribute(H1_err_));
  mesh_->add_attribute("H1 error", DIM, H1_err_attr);
  DataType total_H1_err = std::accumulate(H1_err_.begin(), H1_err_.end(), 0.);
  DataType global_H1_err = 0.;
  MPI_Reduce(&total_H1_err, &global_H1_err, 1, MPI_DOUBLE, MPI_SUM, MASTER_RANK, comm_);
  PLOG_INFO(rank_, "error", "Local H1 error on partition " << rank_ << " = " << std::sqrt(total_H1_err));

  LOG_INFO("Global L2 error", std::sqrt(global_L2_err));
  LOG_INFO("Global H1 error", std::sqrt(global_H1_err));
}

