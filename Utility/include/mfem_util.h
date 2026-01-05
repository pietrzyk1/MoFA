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
        for(int i = 0; i < vdofs_marker.Size(); i++) {
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

#ifdef MPI_BUILD
    class ParEnforcedSolsUpdateRHSOperator
    {
    protected:
        const Array<int> attr_marker;
        
        HypreParMatrix *BC_elim = nullptr;

        bool isOperatorSet = false;
        
    public:
        // Class constructors
        ParEnforcedSolsUpdateRHSOperator(const Array<int> &attr_marker_) : attr_marker(attr_marker_) {}

        // Function set the values of BC_mat (the operator) from a BilinearForm
        void SetOperator(HypreParMatrix &varf_mat, ParBilinearForm &varf);
        
        // Function that takes in the solution vector and b vector (i.e., RHS) and updates b reflect the effects of the temporal BC
        void UpdateLinearFormVector(const Vector &x, Vector &b);

        // Class destructor
        ~ParEnforcedSolsUpdateRHSOperator()
        {
            delete BC_elim;
        }
    };
#endif





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

#ifdef MPI_BUILD
class ParDirichlettBC : public Coefficient, public ParEnforcedSolsUpdateRHSOperator
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
    using ParEnforcedSolsUpdateRHSOperator::ParEnforcedSolsUpdateRHSOperator;
    
    // Define the class constructors
    ParDirichlettBC(const Array<int> &attr_marker_) : ParEnforcedSolsUpdateRHSOperator(attr_marker_), Coefficient() {}

    // function 'SetTime(double t) { time = t; }' is owned by Coefficient
    
    // Sets the spatially-dependent part/function of the BC to be evaluated (i.e., X(x) in BCs of the form u(x,t) = X(x)T(t))
    void SetSpatiallyDependentFunction(function<void(const Vector&, double&)> X) { BC_func_space = X; isSpaceFuncSet = true; }

    // Sets the temporally-dependent part/function of the BC to be evaluated (i.e., T(t) in BCs of the form u(x,t) = X(x)T(t))
    void SetTemporallyDependentFunction(function<void(const double&, double&)> T) { BC_func_time = T; isTimeFuncSet = true; }


    // Function to project the BC to the corresponding vdofs of a gridfunction
    void ProjectToGridFunction(ParGridFunction &gf);

    // The handle for evaluating the BC. Used by functions that project the BC to gridfunctions, etc.
    double Eval(ElementTransformation &T, const IntegrationPoint &ip);
    
    // For defining the Eval function in parts: a temporally-dependent part and a spatially-dependent part (MoFA considers space-time separable BCs).
    double EvalTemporallyDependentFunction(const double &t);

    // For defining the Eval function in parts: a temporally-dependent part and a spatially-dependent part (MoFA considers space-time separable BCs).
    double EvalSpatiallyDependentFunction(ElementTransformation &T, const IntegrationPoint &ip);
};
#endif





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
    string upscaling_theory;

    int N_AR;
    int vdim;
    int DOFs_per_vdim; // This is for the dummy avg space (before making it size N_AR * vdim)
    
    FiniteElementSpace* fespace = nullptr;
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
    AveragingOperator(FiniteElementSpace *fespace_, string upscaling_theory_ = "MoFA") : fespace(fespace_), upscaling_theory(upscaling_theory_)
    {
        // Get the number of averaging regions defined by the mesh
        if (upscaling_theory == "MoFA") { N_AR = fespace_->GetMesh()->attributes.Max(); } 
        else if (upscaling_theory == "Homogenization") { N_AR = 1; }

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
    AveragingOperator(FiniteElementSpace *fespace_, vector<int> AR_tags, string upscaling_theory_ = "MoFA") : fespace(fespace_), upscaling_theory(upscaling_theory_)
    {
        // Get the number of averaging regions defined by the mesh
        if (upscaling_theory == "MoFA") { N_AR = AR_tags.size(); } 
        else if (upscaling_theory == "Homogenization") { N_AR = 1; }

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
        GetARIndices(fespace_->GetMesh(), AR_tags);
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
    void GetARIndices(Mesh *mesh_, vector<int> AR_tags = {-1});

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
    vector<double> GetAR_areas_vector()
    {
        Array<double> AR_areas_Array(this->GetAR_areas());
        vector<double> AR_areas_vector(AR_areas_Array.GetData(), AR_areas_Array.GetData() + AR_areas_Array.Size());
        return AR_areas_vector;
    }

    // Function to get the averaging operator matrix
    SparseMatrix Getavg_mat() { return avg_mat; }

    // Function meant to mimic MFEM's SparseMatrix return on the (Mixed)BilinearForm class
    SparseMatrix &SpMat() { return avg_mat; }

    // Function to get a reference to AR_Dof_inds
    vector<vector<vector<int>>> &GetAR_Dof_inds() { return AR_Dof_inds; }

    // Function for getting a diagonal matrix of the AR pore areas
    SparseMatrix GetAR_areas_SpMat();
};





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
        delete F_inv;
        delete PC;
        delete M_inv;
    }
};

#ifdef MPI_BUILD
    // ===============================================================
    //   Declare the operator for handling the time stepping of
    //   linear systems in MPI parallel.
    // ===============================================================
    class ParLinearTimeDependentOperator : public TimeDependentOperator
    {
        // For linear systems of the form:   M * du/dt + A * u = b
        // This operator takes in the bilinear forms corresponding to the
        // mass matrix (M) and linear system (A), as well as the linear form (b).

    private:
        MPI_Comm comm;

        HypreParMatrix &M, &A;
        Vector &b;
        double dt;

        Array<int> block_offsets_true;
        
        HypreParMatrix *F = NULL;
        BlockOperator *Op = NULL;
        BlockDiagonalPreconditioner *PC = NULL;

        HypreSmoother *M_inv = NULL;
        CGSolver explicitSolver; // Solver that applies the inverse of symmetric, positive-definite mass matrix M
        
        Solver *F_inv = NULL; 
        GMRESSolver gmresSolver; // Solver that applies the inverse of a non-symmetic, non-positive-definite matrix (for implicit Euler)
        CGSolver cgSolver; // Solver that applies the inverse of a symmetic, positive-definite matrix (for implicit Euler)
        Solver *implicitSolver = NULL;
        
    public:
        // Constructor method
        ParLinearTimeDependentOperator(HypreParMatrix &M_, HypreParMatrix &A_, Vector &b_, double dt_, MPI_Comm comm_ ) : M(M_), A(A_), b(b_), dt(dt_), comm(comm_)
        {
            // Initiate the preconditioner
            block_offsets_true.SetSize(2);
            block_offsets_true[0] = 0;
            block_offsets_true[1] = A.Height();
            block_offsets_true.PartialSum();
            Op = new BlockOperator(block_offsets_true);
            PC = new BlockDiagonalPreconditioner(block_offsets_true);
        }
        

        // Function to create the preconditioner for Explicit Euler and set up the solver
        void PrepareExplicitEuler();
        
        // Function to create the preconditioner for Implicit Euler and set up the solver
        void PrepareImplicitEuler(string solver_type = "GMRESSolver", string PC_type = "BlockILU");
        /*
        // Compute du/dt of the ODE system    M * du/dt + A * u = b
        void Mult(const Vector &u, Vector &du_dt) const override;
        */
        // Compute du/dt of the ODE system   (M + A*dt) * du/dt + A * u = b
        void ImplicitSolve(const real_t dt, const Vector &u, Vector &du_dt) override;
        
        // Class destructor
        ~ParLinearTimeDependentOperator()
        {
            delete F;
            delete Op;
            delete F_inv;
            delete PC;
            delete M_inv;
        }
    };
#endif





// Function for projecting an average solution vector onto a GridFunction. We assume the
// GridFunction has a pore-scale finite element space.
void CreateAvgSolGridFunction(const Vector &sol, GridFunction &u, const string avg_area_type);





// ===============================================================
//   Declare a class for saving the results from the closure
//   problems. This includes saving the gridfunctions of the
//   closure variables, as well as the output text file containing
//   closure residuals and the file names of the gridfunctions.
// ===============================================================
class ResultsSaver
{
private:
    int N_residual_sets = 0;
    string solve_mode = "serial";
    
    string sim_key = "default_sim_key";
    string AR_number = "";
    int recursive_iter = -1;
    vector<string> eq_keys;
    vector<string> CFNs; // = {closure_output_file_name_u, closure_output_file_name_p};
    //vector<int> block_IDs_for_saving_residuals; // = {2, 3};
    //vector<int> fespace_vdim; // = {2, 1};

    void VerifyVectorSize(int vecSize) { if (N_residual_sets == 0) { N_residual_sets = vecSize; } else { assert (vecSize == N_residual_sets); } }

    vector<vector<vector<double>>> residuals, alphas, betas, gammas;

    JSONDict residual_dict, time_func_dict, iter_dict;
    vector<JSONDict> eq_dicts, res_dicts, CFN_dicts, alpha_dicts, beta_dicts, gamma_dicts;
    vector<vector<JSONDict>> comp_dicts, alpha_comp_dicts, beta_comp_dicts, gamma_comp_dicts;
    bool initializedJSONDicts = false;

public:
    // Constructor methods
    ResultsSaver() {}
    

    // Functions for setting the various variables
    void SetVariables(string sim_key_, int recursive_iter_, const vector<string> &eq_keys_, const vector<string> &CFNs_, string AR_number_)
    {
        SetSimKey(sim_key_);
        SetRecursiveIter(recursive_iter_);
        SetEqKeys(eq_keys_);
        SetClosureFileNames(CFNs_);
        SetARNumber(AR_number_);
    }

    void SetSolveMode(string solve_mode_) { solve_mode = solve_mode_; }
    void SetSimKey(string sim_key_) { sim_key = sim_key_; }
    void SetRecursiveIter(int recursive_iter_) { recursive_iter = recursive_iter_; }
    void SetEqKeys(const vector<string> &eq_keys_) { VerifyVectorSize(eq_keys_.size()); eq_keys = eq_keys_; }
    void SetClosureFileNames(const vector<string> &CFNs_) { VerifyVectorSize(CFNs_.size()); CFNs = CFNs_; }
    void SetARNumber(string AR_number_) { AR_number = AR_number_; }
    

    // Function for obtaining the closure residuals from the solution block vector
    void ObtainClosureResiduals(const BlockVector &sol_BLK, const vector<int> &block_IDs_for_saving_residuals, const vector<int> &fespace_vdim, int N_AR);
    void ObtainClosureResiduals(const BlockVector &sol_BLK, const vector<int> &block_IDs_for_saving_residuals, const vector<int> &fespace_vdim, int N_AR, const vector<int> AR_inds);

    // Function for initializing the JSONDict vectors at the equation level.
    // This includes eq_dicts, res_dicts, CFN_dicts, and comp_dicts
    void InitializeJSONDicts_EqLevel(const vector<int> &fespace_vdim);

    // Function for loading the residuals and corresponding information/labeling into the various
    // dictionaries, compiling them into a main residual dictionary, and saving that dictionary
    // as a text file
    void SaveResidualDictionary(const string &residual_output_file_path);

    // Function for filling component dictionaries (i.e., either residuals, alphas, betas, gammas, etc.)
    void FillComponentDictionary(vector<vector<JSONDict>> &dicts, vector<vector<vector<double>>> &val);
    // Function for filling parameter dictionaries (i.e., either res_dicts, alpha_dicts, beta_dicts, gamma_dicts, etc.)
    void FillParamDictionary(vector<JSONDict> &dicts, vector<vector<JSONDict>> &val);

    // Function for filling the residual parameters' dictionaries
    void StoreResidualParameters(const string &param_name, const vector<vector<vector<double>>> &param);
    void StoreResidualParameters(const string &param_name, const vector<vector<double>> &param, int N_AR);
    void StoreResidualParameters(const string &param_name, const vector<vector<double>> &param, int N_AR, const vector<int> AR_inds);

    
    // Function for loading the closure residuals/closure file names from a preexisting residual
    // results file into the dictionaries for appending the new closure residuals
    void LoadPreexistingResidualsFile(const string &residual_output_file_path);
};





// ===============================================================
//   Declare a class for manipulating meshes. This includes making
//   a mesh periodic and creating a sub mesh from a parent mesh
//   (for closure problems that have local solutions).
// ===============================================================
class MeshManager
{
private:
    std::shared_ptr<Mesh> parent_mesh, mesh;
    std::shared_ptr<SubMesh> local_mesh;

    vector<Vector> translations = {};

    vector<int> AR_tags_local;
    vector<int> AR_inds_loc2glob;
    int N_AR_local = -1;
    vector<int> central_AR_inds_local;

    bool usingLocalMesh;

    #ifdef MPI_BUILD
        std::shared_ptr<ParMesh> par_mesh;
        int* partitioning;
    #endif

public:
    // Constructor methods
    MeshManager(string mesh_file_path)
        : parent_mesh(std::make_shared<Mesh>(mesh_file_path)) {
            UseParentMesh();
        }
    

    // Function for refining the mesh
    //mesh->UniformRefinement(); // Refine the mesh once uniformly

    // Function for making the mesh periodic (if isPeriodic calls for it). Note, this function creates a new Mesh object that is periodic and deletes the original
    void MakePeriodic(const vector<int> &isPeriodic, const vector<double> &L);

    // Function for creating a "local" mesh around a specified AR index from the "parent" mesh in the class (in variable "mesh")
    void MakeLocalMesh(vector<int> central_AR_inds, int N_neighbor_layers, const vector<int> &AR_tags, const vector<vector<int>> &AR_neighbors);
    
    // Function for clearing and reassigning variable "mesh" to the object pointed at by "parent_mesh"
    void UseParentMesh();

    // Function for clearing and reassigning variable "mesh" to the object pointed at by "local_mesh", but as a Mesh object, not a SubMesh object
    void UseLocalMesh();


    // Functions for providing raw pointers to the "mesh" and "local_mesh" variables
    Mesh* GetParentMesh() { return parent_mesh.get(); }
    Mesh* GetMesh() { return mesh.get(); }
    SubMesh* GetLocalMesh() { return local_mesh.get(); }
    

    // Functions for providing access to other variables in the class
    vector<Vector> get_translations() { return translations; }
    int get_N_AR_local() { return N_AR_local; }
    vector<int> get_AR_kept() { return AR_tags_local; }
    vector<int> get_AR_inds_loc2glob() { return AR_inds_loc2glob; }
    vector<int> get_central_AR_inds_local() { return central_AR_inds_local; }


    // Define functions for MPI-parallel
    #ifdef MPI_BUILD
        // Function for creating a parallel mesh by partitioning the serial mesh (serial mesh could be periodic/local)
        void MakeParallelMesh(MPI_Comm comm);

        // Function for clearing and reassigning variable "mesh" to the object pointed at by "par_mesh", but as a Mesh object, not a ParMesh object
        //void UseParallelMesh();

        // Function for providing a raw pointer to the "mesh" variable
        ParMesh* GetParMesh() { return par_mesh.get(); }
        int* GetPartitioning() { return partitioning; }
    #endif
};





// ===============================================================
//   Declare a class for managing GridFunctions loaded from .gf
//   files. This includes forcing functions, as well as fluid
//   velocity functions for transport simulations. This class
//   can be generlized to manager other GridFunctions as well.
//   This class is used to provide core functionality to classes
//   GridFunctionCoefficientManager and
//   VectorGridFunctionCoefficientManager.
// ===============================================================
class GridFunctionManager
{
protected:
    std::unique_ptr<ifgzstream> gf_ifgz_stream;
    
    FiniteElementSpace *fes;
    
    string fec_name;
    int fec_order;
    int fes_dim;

    std::shared_ptr<GridFunction> gf;
        
    std::unique_ptr<FiniteElementCollection> fec_local;
    std::unique_ptr<FiniteElementSpace> fes_local;
    std::shared_ptr<GridFunction> gf_local;
    
    bool manage_gf_mesh = false;
    Mesh *gf_mesh = nullptr;
    
    #ifdef MPI_BUILD
        ParGridFunction *par_gf;
        FiniteElementSpace *par_fes;
    #endif

public:
    // Constructor methods
    GridFunctionManager(string &gf_file_path, Mesh* gf_mesh_, double gf_scale = 1.0)
    {
        gf_mesh = gf_mesh_;
        LoadGridFunction(gf_file_path, gf_mesh);
        ScaleGridFunction(gf_scale);
    }
    GridFunctionManager(string &gf_file_path, string gf_mesh_file_path, double gf_scale = 1.0)
    {
        // Load in the mesh for the grid function
        gf_mesh = new Mesh(gf_mesh_file_path);
        manage_gf_mesh = true;
        LoadGridFunction(gf_file_path, gf_mesh);
        ScaleGridFunction(gf_scale);
    }
    

    // Function for loading a grid function from a file with path "gf_file_path" onto a mesh with file path "gfc_mesh_file_path"
    void LoadGridFunction(string gf_file_path, Mesh *gf_mesh);

    // Function for loading other grid function information into the class
    void UpdateGridFunctionInfo();

    // Function for scaling the grid function "gf" by "gf_scale"
    void ScaleGridFunction(double gf_scale = 1.0) { *gf *= gf_scale; }

    // Function for transfering the gridfunction "gf" to the local mesh provided
    //void MakeLocalGridFunction(int order, const SubMesh &local_mesh, int dim, string fec_name = "H1_FECollection");
    void MakeLocalGridFunction(SubMesh *local_mesh);

    // Function for making a local finite element collection and finite element space for the grid function coefficient
    void MakeLocalFECAndFES(SubMesh *local_mesh);

    // Function for clearing and reassigning variable "gf" to "gf_local"
    void UseLocalGridFunction();

    
    // Functions for providing access to variables in the class
    GridFunction* GetGridFunction() { return gf.get(); }
    FiniteElementSpace* GetGridFunctionFES() { return fes; }
    Mesh* GetGridFunctionMesh() { return &*gf_mesh; }


    // Define functions for MPI-parallel
    #ifdef MPI_BUILD
        // Function for transfering the gridfunction "gf" to a ParGridFunction
        void MakeParGridFunction(ParMesh* mesh_, int* partitioning_);
    
        // Functions for providing access to variables in the class
        ParGridFunction* GetParGridFunction() { return par_gf; }
        FiniteElementSpace* GetParGridFunctionFES() { return par_fes; }
    #endif


    // Class destructor
    ~GridFunctionManager()
    {
        // Delete the GridFunction mesh, if it is managed by the manager
        if (manage_gf_mesh) { delete gf_mesh; }
        
        // Delete the ParGridFunction
        #ifdef MPI_BUILD
            delete par_gf;
        #endif
    }
};





// ===============================================================
//   Declare a class for managing GridFunctionCoefficients, where
//   the underlying grid function lives in a GridFunctionManager.
// ===============================================================
class GridFunctionCoefficientManager : public GridFunctionManager
{
private:
    GridFunctionCoefficient gfc;

public:
    // Inherit the constructors
    using GridFunctionManager::GridFunctionManager;
    
    // Constructor methods
    GridFunctionCoefficientManager(string &gf_file_path, Mesh* gf_mesh, double gf_scale = 1.0)
        : GridFunctionManager(gf_file_path, gf_mesh, gf_scale) {}
    GridFunctionCoefficientManager(string &gf_file_path, string gf_mesh_file_path, double gf_scale = 1.0)
        : GridFunctionManager(gf_file_path, gf_mesh_file_path, gf_scale) {}


    // Function for making a grid function from "gf" in GridFunctionManager
    void MakeGridFunctionCoefficient() { assert (fes_dim == 1); gfc = GridFunctionCoefficient(gf.get()); }

    
    // Functions for providing access to variables in the class
    GridFunctionCoefficient* GetGridFunctionCoefficient() { return &gfc; }

    
    // Define functions for MPI-parallel
    #ifdef MPI_BUILD
        // Function for making a grid function coefficient from "par_gf" in GridFunctionManager
        void MakeParGridFunctionCoefficient() { assert (fes_dim == 1); gfc = GridFunctionCoefficient(par_gf); }
    #endif
};





// ===============================================================
//   Declare a class for managing VectorGridFunctionCoefficients,
//   where the underlying grid function lives in a
//   GridFunctionManager.
// ===============================================================
class VectorGridFunctionCoefficientManager : public GridFunctionManager
{
private:
    VectorGridFunctionCoefficient vgfc;

public:
    // Inherit the constructors
    using GridFunctionManager::GridFunctionManager;
    
    // Constructor methods
    VectorGridFunctionCoefficientManager(string &gf_file_path, Mesh* gf_mesh, double gf_scale = 1.0)
        : GridFunctionManager(gf_file_path, gf_mesh, gf_scale) {}
    VectorGridFunctionCoefficientManager(string &gf_file_path, string gf_mesh_file_path, double gf_scale = 1.0)
        : GridFunctionManager(gf_file_path, gf_mesh_file_path, gf_scale) {}


    // Function for making a vector grid function coefficient from "gf" in GridFunctionManager
    void MakeVectorGridFunctionCoefficient() { assert (fes_dim == 2 || fes_dim == 3); vgfc = VectorGridFunctionCoefficient(gf.get()); }

    
    // Functions for providing access to variables in the class
    VectorGridFunctionCoefficient* GetVectorGridFunctionCoefficient() { return &vgfc; }


    // Define functions for MPI-parallel
    #ifdef MPI_BUILD
        // Function for making a vector grid function coefficient from "par_gf" in GridFunctionManager
        void MakeParVectorGridFunctionCoefficient() { assert (fes_dim == 2 || fes_dim == 3); vgfc = VectorGridFunctionCoefficient(par_gf); }
    #endif
};

