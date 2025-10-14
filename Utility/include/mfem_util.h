/*lL.lH.
#
# Copyright (c) 2025, Lawrence Livermore National Security, LLC
# and other MoFA project developers. All Rights reserved.
# See files LICENSE and NOTICE for details. LLNL-CODE-2006961.
#
# This file is part of the MoFA Project. For more information
# and source code availability visit https://github.com/pietrzyk1/MoFA.
#
# SPDX-License-Identifier: BSD-3-Clause
#
lL.lH.*/

#pragma once
#include <iostream>
#include "mfem.hpp"
#include "JSON_IO.h"



using namespace std;    
using namespace mfem;



// ===============================================================
//   While implementing Dirichlett boundary conditions, rows in
//   a system's matrix are zeroed out (and diagonals set to 1). The
//   corresponding columns are also zeroed (to retain symmetry),
//   and the values that were zeroed get multiplied by the
//   solution and subtracted from the RHS. If the BCs depend on
//   time, this process must be repeated each time step.
//
//   This class finds a matrix that can be multiplied by the
//   vdofs of the time-varying BC to find the contributions to
//   the RHS at each time step (more efficient than reconstructing
//   the system matrix). BCs are assumed (for the purposes of MoFA,
//   but not required in general) to be of the form
//   u(x,t) = X(x)T(t), where u and X can have vdim > 1, and
//   vdim = 1 for T.
// ===============================================================
class EnforcedSolsUpdateRHSOperator
{
protected:
    const FiniteElementSpace &fespace;
    const Array<int> attr_marker;
    
    Array<int> marked_ind, unmarked_ind, unmarked_ind_effected;
    
    // Initialize the matrix that finds the contribution to the RHS entries (corresponding to the unmarked vdofs) due to the essential BCs (i.e., marked vdofs).
    // This matrix operators on vectors of marked solution vdofs and produces vectors of size "unmarked_ind".
    DenseMatrix *BC_mat = nullptr;

    bool isOperatorSet = false;
    
    // Private/Protected functions
    void BuildBC_mat(const SparseMatrix &varf_SpMat);
    
public:
    // Class constructors
    EnforcedSolsUpdateRHSOperator(const FiniteElementSpace &fespace_, const Array<int> &attr_marker_) : fespace(fespace_), attr_marker(attr_marker_)
    {
        // Get an array, vdofs_marker, that is -1 for vdofs that are part of the BC, and 0 for vdofs that are not
        Array<int> vdofs_marker;
        fespace.GetEssentialVDofs(attr_marker, vdofs_marker);

        // Use vdofs_marker to get the indices for marked and unmarked vdofs
        for(int i = 0; i < vdofs_marker.Size(); i++)
        {
            if (vdofs_marker[i] == -1) { marked_ind.Append(i); }
            else { unmarked_ind.Append(i); }
        }
    }

    // Function set the values of BC_mat (the operator) from a BilinearForm
    void SetOperator(const BilinearForm &varf);
    void SetOperator(const MixedBilinearForm &varf);

    // Define a function for removing zero rows from a dense matrix
    DenseMatrix RemoveZeroRows(const DenseMatrix &A, Array<int> &kept_rows_inds, const double tol = 0.0);
    
    // Function that takes in the solution vector and b vector (i.e., RHS) and updates b reflect the effects of the temporal BC
    void UpdateLinearFormVector(const Vector &x, Vector &b);

    // Function to return the finite element space in the class
    //const FiniteElementSpace &GetFESpace() const { return fespace; }
    
    // Function to return the vdofs marked as BC vdofs
    //Array<int> GetMarkedVDofs() const { return marked_ind; }

    // Function to return the vdofs not marked as BC vdofs
    //Array<int> GetUnmarkedVDofs() const
    //{
    //    if (!isOperatorSet)
    //    {
    //        cout << "CRITICAL ERROR: mfem_util.h: TimeDependentBC2::GetUnmarkedVDofs: Function called before 'SetOperator'. 'unmarked_ind' is not finalized until 'SetOperator' is called." << endl;
    //        exit(EXIT_FAILURE);
    //    }
    //    return unmarked_ind;
    //}

    // Function to return the vdofs not marked as BC vdofs
    //DenseMatrix GetBC_mat() const { return *BC_mat; }
    
    // Class destructor
    ~EnforcedSolsUpdateRHSOperator()
    {
        delete BC_mat;
    }
};





// ===============================================================
//   This class is the front end of EnforcedSolsUpdateRHSOperator.
//   It has additional functions for integrating better with MFEM
//   functions and objects. Users should define DirichlettBC or
//   VectorDirichlettBC to use EnforcedSolsUpdateRHSOperator.
//   Here, BCs are assumed (for the purposes of MoFA, but not
//   required in general) to be of the form
//   u(x,t) = X(x)T(t), where u, X, and T have vdim = 1.
// ===============================================================
class DirichlettBC : public Coefficient, public EnforcedSolsUpdateRHSOperator
{
private:
    function<void(const Vector&, double&)> BC_func_space;
    bool isSpaceFuncSet = false;

    function<void(const double&, double&)> BC_func_time;
    bool isTimeFuncSet = false;

    // variable 'time' is owned by Coefficient

public:
    // Inherit the constructors
    using Coefficient::Coefficient;
    using EnforcedSolsUpdateRHSOperator::EnforcedSolsUpdateRHSOperator;
    
    // Define the class constructors
    DirichlettBC(const FiniteElementSpace &fespace_, const Array<int> &attr_marker_) : EnforcedSolsUpdateRHSOperator(fespace_, attr_marker_), Coefficient() {}

    // function 'SetTime(double t) { time = t; }' is owned by Coefficient
    
    // Sets the spatially-dependent part/function of the BC to be evaluated (i.e., X(x) in BCs of the form u(x,t) = X(x)T(t))
    void SetSpatiallyDependentFunction(function<void(const Vector&, double&)> X) { BC_func_space = X; isSpaceFuncSet = true; }

    // Sets the temporally-dependent part/function of the BC to be evaluated (i.e., T(t) in BCs of the form u(x,t) = X(x)T(t))
    void SetTemporallyDependentFunction(function<void(const double&, double&)> T) { BC_func_time = T; isTimeFuncSet = true; }


    // Function to project the BC to the corresponding vdofs of a gridfunction
    void ProjectToGridFunction(GridFunction &gf);

    // The handle for evaluating the BC. Used by functions that project the BC to gridfunctions, etc.
    double Eval(ElementTransformation &T, const IntegrationPoint &ip);
    
    // For defining the Eval function in parts: a temporally-dependent part and a spatially-dependent part (MoFA considers space-time separable BCs).
    double EvalTemporallyDependentFunction(const double &t);

    // For defining the Eval function in parts: a temporally-dependent part and a spatially-dependent part (MoFA considers space-time separable BCs).
    double EvalSpatiallyDependentFunction(ElementTransformation &T, const IntegrationPoint &ip);
};





// ===============================================================
//   This class is the front end of EnforcedSolsUpdateRHSOperator.
//   It has additional functions for integrating better with MFEM
//   functions and objects. Users should define DirichlettBC or
//   VectorDirichlettBC to use EnforcedSolsUpdateRHSOperator.
//   Here, BCs are assumed (for the purposes of MoFA, but not
//   required in general) to be of the form
//   u(x,t) = X(x)T(t), where u and X have vdim > 1 and T has
//   vdim = 1.
// ===============================================================
class VectorDirichlettBC : public VectorCoefficient, public EnforcedSolsUpdateRHSOperator
{
private:
    function<void(const Vector&, Vector&)> BC_func_space;
    bool isSpaceFuncSet = false;

    function<void(const double&, double&)> BC_func_time;
    bool isTimeFuncSet = false;

    // variable 'time' is owned by VectorCoefficient

public:
    // Inherit the constructors
    using VectorCoefficient::VectorCoefficient;
    using EnforcedSolsUpdateRHSOperator::EnforcedSolsUpdateRHSOperator;
    
    // Define the class constructors
    VectorDirichlettBC(const FiniteElementSpace &fespace_, const Array<int> &attr_marker_) : EnforcedSolsUpdateRHSOperator(fespace_, attr_marker_), VectorCoefficient(fespace_.GetVDim()) {}

    // function 'SetTime(double t) { time = t; }' is owned by VectorCoefficient
    
    // Sets the spatially-dependent part/function of the BC to be evaluated (i.e., X(x) in BCs of the form u(x,t) = X(x)T(t))
    void SetSpatiallyDependentFunction(function<void(const Vector&, Vector&)> X) { BC_func_space = X; isSpaceFuncSet = true; }

    // Sets the temporally-dependent part/function of the BC to be evaluated (i.e., T(t) in BCs of the form u(x,t) = X(x)T(t))
    void SetTemporallyDependentFunction(function<void(const double&, double&)> T) { BC_func_time = T; isTimeFuncSet = true; }


    // Function to project the BC to the corresponding vdofs of a gridfunction
    void ProjectToGridFunction(GridFunction &gf);

    // The handle for evaluating the BC. Used by functions that project the BC to gridfunctions, etc.
    void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip);
    
    // For defining the Eval function in parts: a temporally-dependent part and a spatially-dependent part (MoFA considers space-time separable BCs).
    double EvalTemporallyDependentFunction(const double &t);

    // For defining the Eval function in parts: a temporally-dependent part and a spatially-dependent part (MoFA considers space-time separable BCs).
    void EvalSpatiallyDependentFunction(Vector &V, ElementTransformation &T, const IntegrationPoint &ip);
};





// ===============================================================
//   This class acts as a Block Preconditioner for Saddle Point
//   Problems (e.g., Stokes flow). Assuming a system of the form
//                        [ A  B^T]
//                        [ B  D  ],
//   where D can be non-zero, this class creates the Schur
//   complement preconditioner for the system. 
// ===============================================================
class SaddlePointBlockPreconditioner : public BlockDiagonalPreconditioner
{
private:
    SparseMatrix S;
    Solver *invA = NULL, *invS = NULL;

public:
    // Inherited constructors
    using BlockDiagonalPreconditioner::BlockDiagonalPreconditioner;
    
    // Class constructors
    SaddlePointBlockPreconditioner(const Array<int> &offsets_) : BlockDiagonalPreconditioner(offsets_)
    {
        assert(offsets_.Size() == 3);
    }

    // Function for creating the saddle point problem block preconditioner.
    // Assume the system is [ A B^T ]
    //                      [ B  0  ]
    void BuildPreconditioner(const SparseMatrix &A, const SparseMatrix &B, const string A_solver = "DSmoother", const string S_solver = "GSSmoother");

    // Function for creating the saddle point problem block preconditioner.
    // Assume the system is [ A B^T ]
    //                      [ B  D  ]
    void BuildPreconditioner(const SparseMatrix &A, const SparseMatrix &B, const SparseMatrix &D, const string A_solver = "DSmoother", const string S_solver = "GSSmoother");

    // Function for creating the Schur complement
    static SparseMatrix CreateSchurComplement(const SparseMatrix &A, const SparseMatrix &B);
    
    // Function for creating the Schur complement
    static SparseMatrix CreateSchurComplement(const SparseMatrix &A, const SparseMatrix &B, const SparseMatrix &D);
    
    // Class destructor
    ~SaddlePointBlockPreconditioner()
    {
        delete invA;
        delete invS;
    }
};





// ===============================================================
//   Declare the class for handling initial conditions.
// ===============================================================
class InitialCondition : public Coefficient
{
private:
    function<void(const Vector&, double&)> IC_func;
    bool isFuncSet = false;

public:
    // Constuctors
    InitialCondition() {}

    // Sets the IC function to be evaluated
    void SetICFunction(function<void(const Vector&, double&)> X) { IC_func = X; isFuncSet = true; }
    
    // The handle for evaluating the IC. Used by functions that project the IC to gridfunctions, etc.
    double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};





// ===============================================================
//   Declare the class for creating the averaging operator from
//   the finite element space of a variable.
// ===============================================================
class AveragingOperator
{
private:
    int N_AR;
    int vdim;
    int DOFs_per_vdim; // This is for the dummy avg space (before making it size N_AR * vdim)
    
    unique_ptr<FiniteElementSpace> fespace = nullptr;
    unique_ptr<FiniteElementSpace> fespace_avg = nullptr;
    unique_ptr<FiniteElementCollection> fec_avg = nullptr;
    unique_ptr<MixedBilinearForm> varf_avg = nullptr;
    
    vector<vector<int>> AR_Elem_inds;
    vector<vector<vector<int>>> AR_Dof_inds;

    SparseMatrix avg_mat;
    Array<double> AR_areas;

    SparseMatrix weighted_avg_mat;
    bool correctARAreas = false;


public:
    // Class constructors
    AveragingOperator(FiniteElementSpace *fespace_) : fespace(fespace_)
    {
        N_AR = fespace_->GetMesh()->attributes.Max(); // Get the number of averaging regions defined by the mesh
        vdim = fespace_->GetVDim(); // Get vdim (i.e., the number of vector components) used to define the finite element space
        
        for (int i = 0; i < N_AR; i++) { AR_Elem_inds.push_back(vector<int>()); } // Initialize AR_Elem_inds as an array of arrays
        
        AR_areas.SetSize(N_AR); // Initialize the size of the AR_areas vector

        // Initialize AR_Dof_inds
        for (int i_vd = 0; i_vd < vdim; i_vd++)
        {
            AR_Dof_inds.push_back(vector<vector<int>>());
            for (int i_AR = 0; i_AR < N_AR; i_AR++) { AR_Dof_inds[i_vd].push_back(vector<int>()); }
        }

        // Carry out the initialization/avg. operator creation procedure
        CreateDummyFESpace(fespace_->GetMesh(), vdim);
        CreateBilinearForm(fespace_);
        GetARIndices(fespace_->GetMesh());
        CreateAvgOperator();
        correctARAreas = true;
    }
    AveragingOperator(FiniteElementSpace *fespace_, FiniteElementSpace *fespace_avg_, bool domain_avg = false) : fespace(fespace_), fespace_avg(fespace_avg_)
    {
        // Get the number of averaging regions that should be defined over the mesh
        if (!domain_avg) { N_AR = fespace_->GetMesh()->attributes.Max(); } 
        else { N_AR = 1; }

        // Get vdim (i.e., the number of vector components) used to define the finite element space
        vdim = fespace_->GetVDim();
        
        // Initialize AR_Elem_inds as an array of arrays
        for (int i = 0; i < N_AR; i++) { AR_Elem_inds.push_back(vector<int>()); } 
        
        AR_areas.SetSize(N_AR); // Initialize the size of the AR_areas vector

        // Initialize AR_Dof_inds
        for (int i_vd = 0; i_vd < vdim; i_vd++)
        {
            AR_Dof_inds.push_back(vector<vector<int>>());
            for (int i_AR = 0; i_AR < N_AR; i_AR++) { AR_Dof_inds[i_vd].push_back(vector<int>()); }
        }

        // Carry out the initialization/avg. operator creation procedure
        DOFs_per_vdim = fespace_avg->FEColl()->GetFE(fespace->GetMesh()->GetElementBaseGeometry(0), 0)->GetDof(); // Get the number of DOFs per vector component of FiniteElementCollection
        CreateBilinearForm(fespace_);
        GetARIndices(fespace->GetMesh());
        PrepareUnweightedAvgOperator();
        correctARAreas = false;
    }

    // Function for creating the dummy avg finite element space that will be used to generate the averaging operator
    void CreateDummyFESpace(Mesh *mesh_, const int vdim_);

    // Function for creating the bilinear form used to generate the averaging operator
    void CreateBilinearForm(FiniteElementSpace *fespace_);

    // Function for getting 1.) the indices of elements inside each AR, and 2.) the indices of dofs (in the dummy space) of each element inside each AR
    void GetARIndices(Mesh *mesh_);

    // Function for creating the averaging operator from the avg space and mixed bilinear form
    void CreateAvgOperator(const double tol = 1e-15);

    // Function for creating (but not finalizing) the unweighted averaging operator from the avg space and mixed bilinear form,
    // as well as computing the AR pore-space areas.
    void PrepareUnweightedAvgOperator(const double tol = 1e-15);

    // Function for creating (and finalizing) the weighted averaging operator from the non-finalized unweighted averaging operator.
    void CreateWeightedAvgOperator(const GridFunction &weight_func);

    // Function that zeroes-out the columns of the avging operator corresponding to DOFs marked by "bdr_attr_is_ess".
    // For the zeroed entries, the function multiplies the "pre-zeroed" entry value by the corresponding value in "sol" and subtracts the product from "rhs".
    void EliminateTrialEssentialBC(const Array<int> &bdr_attr_is_ess, const Vector &sol, Vector &rhs);

    // Function for applying the averaging operator
    void ApplyAvgOperator(const Vector &x, Vector &avg);
    void ApplyAvgOperator(const Vector &x, Vector &avg, const vector<double> &porosities);
    void ApplyWeightedIntegralOperator(const Vector &x, Vector &avg);


    // Static method for computing the AR porosities
    static void ComputePorosities(const vector<double> &pore_areas, const vector<double> &AR_areas, vector<double> &porosities);
    

    // Function to get the number of averaging regions
    int GetNAR() { return N_AR; }

    // Function to get the vector dimension of the operator
    int Getvdim() { return vdim; }

    // Function to get the averaging region areas
    Array<double> GetAR_areas() { if (correctARAreas) {return AR_areas;} else { cout << "AveragingOperator: GetAR_areas(): AR areas should only be obtained from an unweighted AveragingOperator." << endl; exit(1); } }

    // Function to get the averaging operator matrix
    SparseMatrix Getavg_mat() { return avg_mat; }

    // Function meant to mimic MFEM's SparseMatrix return on the (Mixed)BilinearForm class
    SparseMatrix &SpMat() { return avg_mat; }

    // Function to get a reference to AR_Dof_inds
    vector<vector<vector<int>>> &GetAR_Dof_inds() { return AR_Dof_inds; }
};




// POTENTIALLY NOT USED ANY MORE
/*
// ===============================================================
//   Declare the operator for handling time stepping for linear
//   equations
// ===============================================================
class TransportOperator : public TimeDependentOperator
{
    //This operator is M * du/dt + K * u = b

private:
    BlockOperator &M, &K; // M is the mass matrix and K is the stiffness matrix
    BlockVector &b;
    const double dt;

    BlockDiagonalPreconditioner &M_PC; // Preconditioner to mass matrix M
    BlockDiagonalPreconditioner &T_PC; // Preconditioner to the mass matrix M + K*dt (for implicit Euler) 
    
    CGSolver M_inv; // Solver that applies the inverse of symmetric, positive-definite mass matrix M
    GMRESSolver T_inv; // Solver that applied the inverse of a non-symmetic, non-positive-definite matrix T (for implicit Euler)
    
    unique_ptr<BlockOperator> M_Kdt = nullptr;

    
public:
    // Constructor method
    TransportOperator(BlockOperator &M_, BlockOperator &K_, BlockVector &b_, double dt_, BlockDiagonalPreconditioner &M_PC_, BlockDiagonalPreconditioner &T_PC_);
    
    // Static method for summing block operators
    static unique_ptr<BlockOperator> AddBlockOperators(const double A_coeff, const BlockOperator &A, const double B_coeff, const BlockOperator &B);
    
    // Compute du/dt of the ODE system    M * du/dt + K * u = b
    void Mult(const Vector &u, Vector &du_dt) const override;
    
    // Compute du/dt of the ODE system   (M + K*dt) * du/dt + K * u = b
    void ImplicitSolve(const real_t dt, const Vector &u, Vector &du_dt) override;

    // Class destructor
    ~TransportOperator() {}
};
*/




// ===============================================================
//   Declare the operator for handling the time stepping of
//   linear systems.
// ===============================================================
class LinearTimeDependentOperator : public TimeDependentOperator
{
    // For linear systems of the form:   M * du/dt + A * u = b
    // This operator takes in the bilinear forms corresponding to the
    // mass matrix (M) and linear system (A), as well as the linear form (b).

private:
    SparseMatrix &M, &A;
    Vector &b;
    double dt;

    Array<int> block_offsets_;
    
    SparseMatrix *F = NULL;
    BlockOperator *Op = NULL;
    BlockDiagonalPreconditioner *PC = NULL;

    DSmoother *M_inv = NULL;
    CGSolver explicitSolver; // Solver that applies the inverse of symmetric, positive-definite mass matrix M
    
    Solver *F_inv = NULL; 
    GMRESSolver gmresSolver; // Solver that applies the inverse of a non-symmetic, non-positive-definite matrix (for implicit Euler)
    CGSolver cgSolver; // Solver that applies the inverse of a symmetic, positive-definite matrix (for implicit Euler)
    Solver *implicitSolver = NULL;
    
public:
    // Constructor method
    LinearTimeDependentOperator(SparseMatrix &M_, SparseMatrix &A_, Vector &b_, double dt_) : M(M_), A(A_), b(b_), dt(dt_)
    {
        // Initiate the preconditioner
        block_offsets_.SetSize(2);
        block_offsets_[0] = 0;
        block_offsets_[1] = A.Height();
        block_offsets_.PartialSum();
        Op = new BlockOperator(block_offsets_);
        PC = new BlockDiagonalPreconditioner(block_offsets_);
    }
    

    // Function to create the preconditioner for Explicit Euler and set up the solver
    void PrepareExplicitEuler();
    
    // Function to create the preconditioner for Implicit Euler and set up the solver
    void PrepareImplicitEuler(string solver_type = "GMRESSolver", string PC_type = "BlockILU");
    
    // Compute du/dt of the ODE system    M * du/dt + A * u = b
    void Mult(const Vector &u, Vector &du_dt) const override;
    
    // Compute du/dt of the ODE system   (M + A*dt) * du/dt + A * u = b
    void ImplicitSolve(const real_t dt, const Vector &u, Vector &du_dt) override;

    // Class destructor
    ~LinearTimeDependentOperator()
    {
        delete F;
        delete Op;
        delete PC;
        delete F_inv;
        delete M_inv;
        delete implicitSolver;
    }
};

