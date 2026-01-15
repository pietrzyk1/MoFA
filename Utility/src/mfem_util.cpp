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

#include <iostream>
#include <cassert>
#include <chrono>
#include <map>
#include "mfem_util.h"

// ML packages
#include <random>
#include <fstream>




using namespace std;
using namespace mfem;




// ===============================================================
//   Define the functions for the EnforcedSolsUpdateRHSOperator class.
// ===============================================================
// Function to set the values of BC_mat (the matrix) from a BilinearForm
void EnforcedSolsUpdateRHSOperator::SetOperator(const BilinearForm &varf)
{
    // Assert that the trial space of the input bilinear form is the same as the previously provided finite element space
    assert(varf.FESpace() == &fespace);// "mfem_util.cpp: TimeDependentBC::SetOperator: The finite element space of the class and the provided bilinear form must be the same.");
    
    // Get a reference to the sparse matrix of the variational form
    const SparseMatrix &varf_SpMat(varf.SpMat());

    // Build the matrix that obtains the contribution to the RHS from the BC
    BuildBC_mat(varf_SpMat);
    
    // Set isOperatorSet to 'true', now that the operator is set
    isOperatorSet = true;
}
void EnforcedSolsUpdateRHSOperator::SetOperator(const MixedBilinearForm &varf)
{
    // Assert that the trial space of the input (mixed) bilinear form is the same as the previously provided finite element space
    assert(varf.TrialFESpace() == &fespace);// "mfem_util.cpp: TimeDependentBC::SetOperator: The finite element space of the class and the provided bilinear form must be the same.");
    
    // Reset the 'unmarked_ind' array to an array that counts 1 -> N_vdofs of the test space
    unmarked_ind.SetSize(0); for (int i = 0; i < varf.TestFESpace()->GetNDofs(); i++) { unmarked_ind.Append(i); }
    
    // Get a reference to the sparse matrix of the variational form
    const SparseMatrix &varf_SpMat(varf.SpMat());
    
    // Build the matrix that obtains the contribution to the RHS from the BC
    BuildBC_mat(varf_SpMat);
    
    // Set isOperatorSet to 'true', now that the operator is set
    isOperatorSet = true;
}

void EnforcedSolsUpdateRHSOperator::BuildBC_mat(const SparseMatrix &varf_SpMat)
{
    // Initialize the matrix that finds the contribution to the RHS entries (corresponding to the unmarked vdofs) due to the essential BCs (i.e., marked vdofs)
    BC_mat = new DenseMatrix(unmarked_ind.Size(), marked_ind.Size());
    
    // From the sparse matrix of the variational form, get the submatrix of 'potentially-effected vdofs/equations' x 'marked vdofs'
    varf_SpMat.GetSubMatrix(unmarked_ind, marked_ind, *BC_mat);
    
    // Check the submatrix rows for zero-rows. A zero-row means the corresponding 'potentially-effected equation' is not affected. Get an array of indices
    // (w.r.t. unmarked_ind; see next line) for these rows so they can be ignored in the application of BC_mat. Also rebuild BC_mat to exclude these rows.
    *BC_mat = RemoveZeroRows(*BC_mat, unmarked_ind_effected);
    
    // unmarked_ind_effected in the previous line is w.r.t. unmarked_ind; make it w.r.t. BC_mat
    for (int i = 0; i < unmarked_ind_effected.Size(); i++) { unmarked_ind_effected[i] = unmarked_ind[unmarked_ind_effected[i]]; }

    // Bring the bilinear entries (assumed to be on the LHS) to the RHS
    *BC_mat *= -1;
}

// Define a function for removing zero rows from a dense matrix
DenseMatrix EnforcedSolsUpdateRHSOperator::RemoveZeroRows(const DenseMatrix &A, Array<int> &kept_rows_inds, const double tol)
{
    kept_rows_inds.DeleteAll();
    int rows = A.NumRows(), cols = A.NumCols();
    bool is_zero;
    
    // Find rows that are not all zeros
    for (int i = 0; i < rows; i++)
    {
        is_zero = true;
        for (int j = 0; j < cols; j++)
        {
            if (abs(A.Elem(i, j)) > tol)
            {
                is_zero = false;
                break;
            }
        }
        if (!is_zero) { kept_rows_inds.Append(i); }   
    }

    // Create a new DenseMatrix with only the nonzero rows
    int N_nonzero_rows = kept_rows_inds.Size(), row_idx;
    DenseMatrix A_new(N_nonzero_rows, cols);
    for (int i = 0; i < N_nonzero_rows; i++)
    {
        row_idx = kept_rows_inds[i];
        for (int j = 0; j < cols; j++)
        {
            A_new.Elem(i, j) = A.Elem(row_idx, j);
        }
    }
    
    return A_new;
}

// Function that takes in the solution vector and b vector (i.e., RHS) and updates b to reflect the effects of the temporal BC (using BC_mat)
void EnforcedSolsUpdateRHSOperator::UpdateLinearFormVector(const Vector &x, Vector &b)
{
    // Initialize vectors for the marked solution vdofs and unmarked (but affected) RHS entries
    Vector x_marked(marked_ind.Size()), b_unmarked(unmarked_ind_effected.Size()), b_temp(unmarked_ind_effected.Size());

    // Get the marked vdofs from the solution vector; these will be used to get the contributions to the RHS/b vector for the unmarked vdofs
    x.GetSubVector(marked_ind, x_marked);

    // Apply the matrix to the marked vdofs to get the contributions to the RHS/b vector
    BC_mat->Mult(x_marked, b_temp);

    // Get the unmarked vdofs from the RHS/b vector that are affected by the BC vdofs, add the contribution (b_temp) to them, and put them back into the b vector
    b.GetSubVector(unmarked_ind_effected, b_unmarked);
    b_unmarked += b_temp;
    b.SetSubVector(unmarked_ind_effected, b_unmarked);
}


#ifdef MPI_BUILD
    // ===============================================================
    //   Define the functions for the ParEnforcedSolsUpdateRHSOperator class.
    // ===============================================================
    // Function to set the values of BC_mat (the matrix) from a BilinearForm
    void ParEnforcedSolsUpdateRHSOperator::SetOperator(HypreParMatrix &varf_mat, ParBilinearForm &varf)
    {
        // Separate-out the rows and columns (into BC_elim) of varf_mat that correspond to the true vdofs of the BC
        BC_elim = varf.ParallelEliminateEssentialBC(attr_marker, varf_mat);
        *BC_elim *= -1;
                
        // Set isOperatorSet to 'true', now that the operator is set
        isOperatorSet = true;
    }

    // Function that takes in the solution vector and b vector (i.e., RHS) and updates b to reflect the effects of the temporal BC (using BC_mat)
    void ParEnforcedSolsUpdateRHSOperator::UpdateLinearFormVector(const Vector &x, Vector &b)
    {
        Vector b_dummy(b.Size()); b_dummy = 0.0;
        BC_elim->Mult(x, b_dummy);
        b += b_dummy;
    }
#endif





// ===============================================================
//   Define the functions for the DirichlettBC class.
// ===============================================================
// Function to project the BC to the corresponding vdofs of a gridfunction. NOTE: This function requires that DirichlettBC has
// the method double Eval(ElementTransformation &T, const IntegrationPoint &ip)
void DirichlettBC::ProjectToGridFunction(GridFunction &gf)
{
    gf.ProjectBdrCoefficient(*this, attr_marker); // Requires the method Eval(ElementTransformation &T, const IntegrationPoint &ip)
}

// The handle for evaluating the BC. Used by functions that project the BC to gridfunctions, etc.
double DirichlettBC::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
    if (!isSpaceFuncSet || !isTimeFuncSet)
    {
        // In the future, add an option here to evaluate a function for general BCs of the form u(x,t) = F(x,t)
        cout << "CRITICAL ERROR: mfem_util.cpp: DirichlettBC::Eval: Either the space or time function (or both) of the BC have not been set. Set these functions with the 'SetSpatiallyDependentFunction' and 'SetTemporallyDependentFunction' methods." << endl;
        exit(EXIT_FAILURE);
    }
    return EvalTemporallyDependentFunction(time) * EvalSpatiallyDependentFunction(T, ip);
}

// The handle for evaluating the BC. Used by functions that project the BC to gridfunctions, etc.
double DirichlettBC::EvalTemporallyDependentFunction(const double &t)
{
    double sol = 0.0;
    BC_func_time(t, sol);
    return sol;
}

// The handle for evaluating the BC. Used by functions that project the BC to gridfunctions, etc.
double DirichlettBC::EvalSpatiallyDependentFunction(ElementTransformation &T, const IntegrationPoint &ip)
{
    // Get the spatial coordinates in x
    Vector x;
    T.Transform(ip, x); // x can be called like x(0), x(1), x(2) to get x, y, and z coordinates
    
    double sol = 0;
    BC_func_space(x, sol);
    return sol;
}


#ifdef MPI_BUILD
    // Function to project the BC to the corresponding vdofs of a gridfunction. NOTE: This function requires that DirichlettBC has
    // the method double Eval(ElementTransformation &T, const IntegrationPoint &ip)
    void ParDirichlettBC::ProjectToGridFunction(ParGridFunction &gf)
    {
        gf.ProjectBdrCoefficient(*this, attr_marker); // Requires the method Eval(ElementTransformation &T, const IntegrationPoint &ip)
    }

    // The handle for evaluating the BC. Used by functions that project the BC to gridfunctions, etc.
    double ParDirichlettBC::Eval(ElementTransformation &T, const IntegrationPoint &ip)
    {
        if (!isSpaceFuncSet || !isTimeFuncSet)
        {
            // In the future, add an option here to evaluate a function for general BCs of the form u(x,t) = F(x,t)
            cout << "CRITICAL ERROR: mfem_util.cpp: ParDirichlettBC::Eval: Either the space or time function (or both) of the BC have not been set. Set these functions with the 'SetSpatiallyDependentFunction' and 'SetTemporallyDependentFunction' methods." << endl;
            exit(EXIT_FAILURE);
        }
        return EvalTemporallyDependentFunction(time) * EvalSpatiallyDependentFunction(T, ip);
    }

    // The handle for evaluating the BC. Used by functions that project the BC to gridfunctions, etc.
    double ParDirichlettBC::EvalTemporallyDependentFunction(const double &t)
    {
        double sol = 0.0;
        BC_func_time(t, sol);
        return sol;
    }

    // The handle for evaluating the BC. Used by functions that project the BC to gridfunctions, etc.
    double ParDirichlettBC::EvalSpatiallyDependentFunction(ElementTransformation &T, const IntegrationPoint &ip)
    {
        // Get the spatial coordinates in x
        Vector x;
        T.Transform(ip, x); // x can be called like x(0), x(1), x(2) to get x, y, and z coordinates
        
        double sol = 0;
        BC_func_space(x, sol);
        return sol;
    }
#endif




// ===============================================================
//   Define the functions for the VectorDirichlettBC class.
// ===============================================================
// Function to project the BC to the corresponding vdofs of a gridfunction. NOTE: This function requires that VectorDirichlettBC has
// the method void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip)
void VectorDirichlettBC::ProjectToGridFunction(GridFunction &gf)
{
    gf.ProjectBdrCoefficient(*this, attr_marker); // Requires the method Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip)
}

// The handle for evaluating the BC. Used by functions that project the BC to gridfunctions, etc.
void VectorDirichlettBC::Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip)
{
    if (!isSpaceFuncSet || !isTimeFuncSet)
    {
        // In the future, add an option here to evaluate a function for general BCs of the form u(x,t) = F(x,t)
        cout << "CRITICAL ERROR: mfem_util.cpp: VectorDirichlettBC::Eval: Either the space or time function (or both) of the BC have not been set. Set these functions with the 'SetSpatiallyDependentFunction' and 'SetTemporallyDependentFunction' methods." << endl;
        exit(EXIT_FAILURE);
    }

    // NOTE: FOR VECTOR EVALS, YOU HAVE TO SET THE VECTOR SIZE BEFORE USING/ASSIGNING IT A VALUE
    V.SetSize(this->GetVDim());
    
    EvalSpatiallyDependentFunction(V, T, ip);
    V *= EvalTemporallyDependentFunction(time);
}

// The handle for evaluating the BC. Used by functions that project the BC to gridfunctions, etc.
double VectorDirichlettBC::EvalTemporallyDependentFunction(const double &t)
{
    double sol = 0.0;
    BC_func_time(t, sol);
    return sol;
}

// The handle for evaluating the BC. Used by functions that project the BC to gridfunctions, etc.
void VectorDirichlettBC::EvalSpatiallyDependentFunction(Vector &V, ElementTransformation &T, const IntegrationPoint &ip)
{
    // Get the spatial coordinates in x
    Vector x;
    T.Transform(ip, x); // x can be called like x(0), x(1), x(2) to get x, y, and z coordinates
    BC_func_space(x, V);
}





// ===============================================================
//   Define the functions for the SaddlePointBlockPreconditioner class.
// ===============================================================
// Function for creating the saddle point problem block preconditioner.
// Assume the system is [ A B^T ]
//                      [ B  0  ]
void SaddlePointBlockPreconditioner::BuildPreconditioner(const SparseMatrix &A, const SparseMatrix &B, const string A_solver, const string S_solver)
{
    // Approximate the Schur complement of the matrix
    S = CreateSchurComplement(A, B);

    // Use smoother functions to define our approximate inverse matrices. Turn off their iterative modes
    if (A_solver == "BlockILU") { invA = new BlockILU(A); }
    else { invA = new DSmoother(A); }
    if (S_solver == "BlockILU") { invS = new BlockILU(S); }
    else { invS = new GSSmoother(S); }
    
    invA->iterative_mode = false;
    invS->iterative_mode = false;
    
    // Set the block diagonals of the block preconditioner
    this->SetDiagonalBlock(0, invA);
    this->SetDiagonalBlock(1, invS);
}

void SaddlePointBlockPreconditioner::BuildPreconditioner(const SparseMatrix &A, const SparseMatrix &B, const SparseMatrix &D, const string A_solver, const string S_solver)
{
    // Approximate the Schur complement of the matrix
    S = CreateSchurComplement(A, B, D);
    
    // Use smoother functions to define our approximate inverse matrices. Turn off their iterative modes
    if (A_solver == "BlockILU") { invA = new BlockILU(A); }
    else { invA = new DSmoother(A); }
    if (S_solver == "BlockILU") { invS = new BlockILU(S); }
    else { invS = new GSSmoother(S); }

    invA->iterative_mode = false;
    invS->iterative_mode = false;
    
    // Set the block diagonals of the block preconditioner
    this->SetDiagonalBlock(0, invA);
    this->SetDiagonalBlock(1, invS);
}

// Function for generating the Schur complement
SparseMatrix SaddlePointBlockPreconditioner::CreateSchurComplement(const SparseMatrix &A, const SparseMatrix &B)
{
    // Get the diagonal of A to approximate the application of its inverse onto B^T
    Vector A_diag(A.Height()); A.GetDiag(A_diag);
    
    // Approximate the application of A^(-1) on B^T
    unique_ptr<SparseMatrix> A_BT(Transpose(B));
    for (int i = 0; i < A_diag.Size(); i++)
    {
        A_BT->ScaleRow(i, 1./A_diag(i));
    }

    // Multiply by B to complete the approximation of the Schur complement
    return *mfem::Mult(B, *A_BT);
}

// Function for generating the Schur complement
SparseMatrix SaddlePointBlockPreconditioner::CreateSchurComplement(const SparseMatrix &A, const SparseMatrix &B, const SparseMatrix &D)
{
    SparseMatrix sol = CreateSchurComplement(A, B);
    sol *= -1.0;
    sol += D;
    return sol;
}





// ===============================================================
//   Define the functions for the InitialCondition class.
// ===============================================================    
// Define the Eval method for the InitialCondition class
double InitialCondition::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
    // Get the spatial coordinates in x
    Vector x;
    T.Transform(ip, x); // x can be called like x(0), x(1), x(2) to get x, y, and z coordinates
    
    if (isFuncSet)
    {
        double sol = 0;
        IC_func(x, sol);
        return sol;
    }
    else
    {
        bool condition = (x(0) < 0.1 && x(0) > -0.1);
        //bool condition = (x(0) < 0.0);
        if (condition) { return 1; } else { return 0; }

        //bool condition = (x(0) < 0.0);
        //if (condition) { return 0.01; } else { return 0; }
    }    
}





// ===============================================================
//   Define the functions for the AveragingOperator class.
// ===============================================================
// Function for creating the dummy avg finite element space that will be used to generate the averaging operator
void AveragingOperator::CreateDummyFESpace(Mesh *mesh_, const int vdim_)
{
    // Finite element space that will be used to create the average space
    //fec_avg = new L2_FECollection(0, mesh_->Dimension());
    fec_avg = make_unique<L2_FECollection>(0, mesh_->Dimension());
    //fespace_avg = new FiniteElementSpace(mesh_, &*fec_avg, vdim_);
    fespace_avg = make_unique<FiniteElementSpace>(mesh_, &*fec_avg, vdim_);

    // Get the number of DOFs per vector component of FiniteElementCollection
    DOFs_per_vdim = fec_avg->GetFE(mesh_->GetElementBaseGeometry(0), 0)->GetDof();
}

// Function for creating the bilinear form used to generate the averaging operator
void AveragingOperator::CreateBilinearForm(FiniteElementSpace *fespace_)
{
    // Initiate a mixed bilinear form between the concentration space and the dummy average space
    //varf_avg = new MixedBilinearForm(fespace_, &*fespace_avg);
    varf_avg = make_unique<MixedBilinearForm>(fespace_, &*fespace_avg);
    ConstantCoefficient one(1);
    if (fespace_avg->GetVDim() == 1) { varf_avg->AddDomainIntegrator(new MixedScalarMassIntegrator(one)); }
    else { varf_avg->AddDomainIntegrator(new VectorMassIntegrator(one)); }
    varf_avg->Assemble();
    varf_avg->Finalize();
}

// Function for getting 1.) the indices of elements inside each AR, and 2.) the indices of dofs (in the dummy space) of each element inside each AR
void AveragingOperator::GetARIndices(Mesh *mesh_, vector<int> AR_tags)
{
    if (upscaling_theory == "MoFA")
    {
        // Define default AR_tags (i.e., {1 -> N_AR}) if not provided
        if (AR_tags.size() == 1 && AR_tags[0] == -1) {
            AR_tags[0] = 1;
            for (int i = 1; i < N_AR; i++) { AR_tags.push_back( i + 1 ); }
        }
        
        // Get the indices of the elements inside each averaging region
        for (int i_Elem = 0; i_Elem < mesh_->GetNE(); i_Elem++) {
            for (int i_AR = 0; i_AR < N_AR; i_AR++) {
                if (mesh_->GetAttribute(i_Elem) == AR_tags[i_AR]) {
                    AR_Elem_inds[i_AR].push_back(i_Elem);
                    break;
                }
            }
        }
    }
    else if (upscaling_theory == "Homogenization")
    {
        // Get the indices of all elements inside the domain (as the average is taken over the whole unit-cell)
        for (int i_Elem = 0; i_Elem < mesh_->GetNE(); i_Elem++) { AR_Elem_inds[0].push_back(i_Elem); }
    }

    // Get the indices of the DOFS in the dummy avg-space of each element inside each averaging region
    Array<int> Elem_DOFs;
    for (int i_AR = 0; i_AR < N_AR; i_AR++)
    {
        for (int i_Elem = 0; i_Elem < AR_Elem_inds[i_AR].size(); i_Elem++)
        {
            fespace_avg->GetElementVDofs(AR_Elem_inds[i_AR][i_Elem], Elem_DOFs); // GetElementVDofs apparently outputs dofs like [dof_x^1, dof_x^2, dof_y^1, dof_y^2, ...]
            assert (Elem_DOFs.Size() == DOFs_per_vdim * vdim);
            for (int i_vd = 0; i_vd < vdim; i_vd++)
            {
                AR_Dof_inds[i_vd][i_AR].insert(AR_Dof_inds[i_vd][i_AR].end(), Elem_DOFs.begin() + i_vd*DOFs_per_vdim, Elem_DOFs.begin() + (i_vd + 1)*DOFs_per_vdim);
            }
        }
    }
}

// Function for creating the averaging operator from the avg space and mixed bilinear form
void AveragingOperator::CreateAvgOperator(const double tol)
{
    /*
    // Get the sparse matrix of the bilinear form
    SparseMatrix &varf_avg_SpMat(varf_avg->SpMat());
    
    // Sum together the rows of varf_avg_SpMat (i.e., the matrix of the mixed bilinear form) corresponding to the dummy avg-space DOFs
    // in each averaging region (i.e., these are AR_Dof_inds). Insert the sums as rows in avg_mat.
    // Finally, sum together the components of each row and store it in rhs_for_avg; these values will be used to implement non-zero
    // averaging conditions.
    
    SparseMatrix avg_mat_(N_AR * vdim, varf_avg_SpMat.NumCols()); // Initialize the size of the averaging matrix
    Vector row(varf_avg_SpMat.NumCols()), row_sum(varf_avg_SpMat.NumCols());
    
    Array<double> AR_areas_(N_AR);
    Array<int> col_inds, all_col_indices(varf_avg_SpMat.NumCols()); // Create an Array of indices from 0 to the number of columns in cavg_SpMat
    for(int i = 0; i < all_col_indices.Size(); i++) { all_col_indices[i] = i; } // Defines cavg_SpMat
    
    for (int i_vd = 0; i_vd < vdim; i_vd++) // For each vector dimension...
    {
        for (int i_AR = 0; i_AR < N_AR; i_AR++) // ...and each averaging region...
        {
            row_sum = 0.0;
            for (int i = 0; i < AR_Dof_inds[i_vd][i_AR].size(); i++) // ...and each DOF...
            {
                varf_avg_SpMat.GetRow(AR_Dof_inds[i_vd][i_AR][i], col_inds, row); // ...get the row of cavg_SpMat...
                for (int i_sum = 0; i_sum < col_inds.Size(); i_sum++) { row_sum[col_inds[i_sum]] += row[i_sum]; } // ...and add it row_sum to sum the DOF rows.
            }
            for (int i = 0; i < row_sum.Size(); i++) { if (abs(row_sum[i]) < tol) { row_sum[i] = 0; }} // If a component of row_sum is really small, just set it to 0; it is likely numerical error.
            avg_mat_.SetRow(i_AR + N_AR*i_vd, all_col_indices, row_sum); // Insert the row_sum in avg_mat, the matrix discribing the averaging operator
            if (i_vd == 0) { AR_areas_[i_AR] = row_sum.Sum(); } // Keep the sum of the row for implementing averaging conditions (only do this for the first component; assume it will be the same for the others)
        }
    }
    
    avg_mat_.Finalize();
    
    avg_mat = avg_mat_;
    AR_areas = AR_areas_;
    */
    PrepareUnweightedAvgOperator(tol);
    avg_mat.Finalize();
}

// Function for creating the averaging operator from the avg space and mixed bilinear form. Depending on the avg space/mixed bilinear form used,
// the produced avg_mat is either \integral(1 * \basis) or \integral(\basis * \basis) summed over the elements of each AR box.
void AveragingOperator::PrepareUnweightedAvgOperator(const double tol)
{
    // Get the sparse matrix of the bilinear form
    SparseMatrix &varf_avg_SpMat(varf_avg->SpMat());
    
    // Sum together the rows of varf_avg_SpMat (i.e., the matrix of the mixed bilinear form) corresponding to the dummy avg-space DOFs
    // in each averaging region (i.e., these are AR_Dof_inds). Insert the sums as rows in avg_mat.
    // Finally, sum together the components of each row and store it in rhs_for_avg; these values will be used to implement non-zero
    // averaging conditions.
    
    avg_mat = SparseMatrix(N_AR * vdim, varf_avg_SpMat.NumCols());
    //SparseMatrix avg_mat_(N_AR * vdim, varf_avg_SpMat.NumCols()); // Initialize the size of the averaging matrix
    
    Vector row(varf_avg_SpMat.NumCols()), row_sum(varf_avg_SpMat.NumCols());
    
    //Array<double> AR_areas_(N_AR);
    AR_areas = Array<double>(N_AR);
    
    Array<int> col_inds, all_col_indices(varf_avg_SpMat.NumCols()); // Create an Array of indices from 0 to the number of columns in cavg_SpMat
    for(int i = 0; i < all_col_indices.Size(); i++) { all_col_indices[i] = i; } // Defines cavg_SpMat
    
    for (int i_vd = 0; i_vd < vdim; i_vd++) // For each vector dimension...
    {
        for (int i_AR = 0; i_AR < N_AR; i_AR++) // ...and each averaging region...
        {
            row_sum = 0.0;
            for (int i = 0; i < AR_Dof_inds[i_vd][i_AR].size(); i++) // ...and each DOF...
            {
                varf_avg_SpMat.GetRow(AR_Dof_inds[i_vd][i_AR][i], col_inds, row); // ...get the row of cavg_SpMat...
                for (int i_sum = 0; i_sum < col_inds.Size(); i_sum++) { row_sum[col_inds[i_sum]] += row[i_sum]; } // ...and add it row_sum to sum the DOF rows.
            }
            for (int i = 0; i < row_sum.Size(); i++) { if (abs(row_sum[i]) < tol) { row_sum[i] = 0; }} // If a component of row_sum is really small, just set it to 0; it is likely numerical error.
            //avg_mat_.SetRow(i_AR + N_AR*i_vd, all_col_indices, row_sum); // Insert the row_sum in avg_mat, the matrix discribing the averaging operator
            avg_mat.SetRow(i_AR + N_AR*i_vd, all_col_indices, row_sum); // Insert the row_sum in avg_mat, the matrix discribing the averaging operator
            //if (i_vd == 0) { AR_areas_[i_AR] = row_sum.Sum(); } // Keep the sum of the row for implementing averaging conditions (only do this for the first component; assume it will be the same for the others)
            if (i_vd == 0) { AR_areas[i_AR] = row_sum.Sum(); } // Keep the sum of the row for implementing averaging conditions (only do this for the first component; assume it will be the same for the others)
        }
    }
    
    //avg_mat = avg_mat_;
    //AR_areas = AR_areas_;
}

// Function for finalizing the weighted averaging operator (i.e., multiply the rows of avg_mat by the weight function's DOFs in a component-wise fashion).
void AveragingOperator::CreateWeightedAvgOperator(const GridFunction &weight_func)
{
    // Assert that the weight function's finite element space is the same as the dummy fespace used in creating the operator
    assert (weight_func.FESpace() == fespace_avg.get());
    
    Array<int> col_inds;
    Vector row(avg_mat.NumCols()), weighted_row(avg_mat.NumCols());

    // Copy the unweighted average matrix into weighted_avg_mat
    weighted_avg_mat = avg_mat;

    // Go through the rows of weighted_avg_mat and weight each component with the weight function
    for (int i_row = 0; i_row < N_AR * vdim; i_row++)
    {
        weighted_avg_mat.GetRow(i_row, col_inds, row);
        for (int i = 0; i < col_inds.Size(); i++) { weighted_row[col_inds[i]] = weight_func[col_inds[i]] * row[i]; }
        weighted_avg_mat.SetRow(i_row, col_inds, weighted_row);
    }

    // Finalize the weighted average matrix
    weighted_avg_mat.Finalize();
}

// Function that zeroes-out the columns of the avging operator corresponding to DOFs marked by "bdr_attr_is_ess".
// For the zeroed entries, the function multiplies the "pre-zeroed" entry value by the corresponding value in "sol" and subtracts the product from "rhs".
void AveragingOperator::EliminateTrialEssentialBC(const Array<int> &bdr_attr_is_ess, const Vector &sol, Vector &rhs)
{
   Array<int> trial_ess_dofs;
   fespace->GetEssentialVDofs(bdr_attr_is_ess, trial_ess_dofs);
   avg_mat.EliminateCols(trial_ess_dofs, &sol, &rhs);
}

// Function for applying the averaging operator
void AveragingOperator::ApplyAvgOperator(const Vector &x, Vector &avg)
{
    // Define a porosity vector of ones
    vector<double> porosities;
    for (int i = 0; i < N_AR; i++)
    {
        porosities.push_back(1.0);
    }
    
    // Compute the average with a porosity of 1.0
    AveragingOperator::ApplyAvgOperator(x, avg, porosities);
}
void AveragingOperator::ApplyAvgOperator(const Vector &x, Vector &avg, const vector<double> &porosities)
{
    assert(avg.Size() == vdim * N_AR);
    assert(porosities.size() == N_AR);

    // Apply the averaging operator
    avg_mat.Mult(x, avg);
    
    // Divide by the AR pore area and multiply by the AR porosity (default porosity is 1.0)
    int count = 0;
    for (int i_vdim = 0; i_vdim < vdim; i_vdim++)
    {
        for (int i = 0; i < N_AR; i++)
        {
            avg.Elem(count) *= porosities[i] / AR_areas[i];
            count += 1;
        }
    }
}

void AveragingOperator::ApplyWeightedIntegralOperator(const Vector &x, Vector &avg)
{
    assert(avg.Size() == vdim * N_AR);

    // Apply the weighted averaging operator
    weighted_avg_mat.Mult(x, avg);
}

// Static method for computing the AR porosities
void AveragingOperator::ComputePorosities(const vector<double> &pore_areas, const vector<double> &AR_areas, vector<double> &porosities)
{
    assert(pore_areas.size() == AR_areas.size());
    
    porosities.clear();

    for (int i = 0; i < pore_areas.size(); i++)
    {
        porosities.push_back( pore_areas[i] / AR_areas[i] );
    }
}

// Function for getting a diagonal matrix of the AR pore areas
SparseMatrix AveragingOperator::GetAR_areas_SpMat()
{
    SparseMatrix avg_area_mat(vdim * N_AR, vdim * N_AR);
    for (int i_vd = 0; i_vd < vdim; i_vd++) {
        for (int i_ar = 0; i_ar < N_AR; i_ar++) {
            avg_area_mat.Add(i_ar + i_vd*N_AR, i_ar + i_vd*N_AR, AR_areas[i_ar]);
        }
    }
    return avg_area_mat;
}





#ifdef MPI_BUILD
    // ===============================================================
    //   Define the functions for the ParAveragingOperator class.
    // ===============================================================
    // Function for identifying the rows of the final averaging operator that each rank will receive. The output will be
    // a map row_ind -> rank ID
    void ParAveragingOperator::InitializeReceivingRankInfo()
    {
        // Identify which ranks will be receiving the row data of the final averaging operator from all other ranks
        //recv_ranks.resize(N_rows); for (int i = 0; i < N_rows; i++) { recv_ranks[i] = -1; }
        for (int i_res = 0; i_res < N_rows; i_res++) {
            recv_ranks.push_back( -1 );
            for (int i_rank = 0; i_rank < N_ranks; i_rank++) {
                if (i_res >= row_starts[i_rank] && i_res < row_starts[i_rank + 1]) { recv_ranks[i_res] = i_rank; break; } }
            if (recv_ranks[i_res] == -1) { cerr << "ParAveragingOperator::GetReceivingRanks(): CRITICAL ERROR: Could not identify the receiving rank of row " << i_res << " in the averaging matrix." << endl; exit(1); }
        }

        // From recv_ranks, determine how many rows of the final averaging operator are owned by this rank. Also
        // determine which AR and vdim each DOF/row belongs to (this will be used for applying the averaging operator in parallel)
        for (int i_row = 0; i_row < recv_ranks.size(); i_row++) {
            if (recv_ranks[i_row] == rank) {
                N_rows_owned += 1;
                local_row_to_AR_index.push_back( i_row % N_AR );
                //local_row_to_vdim_index.push_back( (i_row - (i_row % N_AR)) / N_AR );
            }
        }
        assert ( N_rows_owned == row_starts[rank + 1] - row_starts[rank] );
    }
    
    // Function for creating the dummy avg finite element space that will be used to generate the averaging operator
    void ParAveragingOperator::CreateDummyFESpace(ParMesh *mesh_, const int vdim_)
    {
        // Finite element space that will be used to create the average space
        fec_avg = make_unique<L2_FECollection>(0, mesh_->Dimension());
        fespace_avg = make_unique<ParFiniteElementSpace>(mesh_, &*fec_avg, vdim_);

        // Get the number of DOFs per vector component of FiniteElementCollection
        DOFs_per_vdim = fec_avg->GetFE(mesh_->GetElementBaseGeometry(0), 0)->GetDof();
    }

    // Function for creating the bilinear form used to generate the averaging operator
    void ParAveragingOperator::CreateBilinearForm(ParFiniteElementSpace *fespace_)
    {
        // Initiate a mixed bilinear form between the concentration space and the dummy average space
        varf_avg = make_unique<ParMixedBilinearForm>(fespace_, &*fespace_avg);
        ConstantCoefficient one(1);
        if (fespace_avg->GetVDim() == 1) { varf_avg->AddDomainIntegrator(new MixedScalarMassIntegrator(one)); }
        else { varf_avg->AddDomainIntegrator(new VectorMassIntegrator(one)); }
        varf_avg->Assemble();
        varf_avg->Finalize();
    }

    // Function for getting 1.) the indices of elements inside each AR, and 2.) the indices of dofs (in the dummy space) of each element inside each AR
    void ParAveragingOperator::GetARIndices(ParMesh *mesh_, vector<int> AR_tags)
    {
        if (upscaling_theory == "MoFA")
        {
            // Define default AR_tags (i.e., {1, 2, ..., N_AR}) if not provided
            if (AR_tags.size() == 1 && AR_tags[0] == -1) { AR_tags[0] = 1; for (int i = 1; i < N_AR; i++) { AR_tags.push_back( i + 1 ); } }
            
            // Get the indices of the elements inside each averaging region
            for (int i_Elem = 0; i_Elem < mesh_->GetNE(); i_Elem++) {
                for (int i_AR = 0; i_AR < N_AR; i_AR++) {
                    if (mesh_->GetAttribute(i_Elem) == AR_tags[i_AR]) {
                        AR_Elem_inds[i_AR].push_back(i_Elem);
                        break;
                    }
                }
            }
        }
        else if (upscaling_theory == "Homogenization")
        {
            // Get the indices of all elements inside the domain (as the average is taken over the whole unit-cell)
            for (int i_Elem = 0; i_Elem < mesh_->GetNE(); i_Elem++) { AR_Elem_inds[0].push_back(i_Elem); }
        }

        // Get the indices of the DOFS in the dummy avg-space of each element inside each averaging region
        Array<int> Elem_DOFs;
        for (int i_AR = 0; i_AR < N_AR; i_AR++)
        {
            for (int i_Elem = 0; i_Elem < AR_Elem_inds[i_AR].size(); i_Elem++)
            {
                fespace_avg->GetElementVDofs(AR_Elem_inds[i_AR][i_Elem], Elem_DOFs); // GetElementVDofs apparently outputs dofs like [dof_x^1, dof_x^2, dof_y^1, dof_y^2, ...]
                assert (Elem_DOFs.Size() == DOFs_per_vdim * vdim);
                for (int i_vd = 0; i_vd < vdim; i_vd++)
                {
                    AR_Dof_inds[i_vd][i_AR].insert(AR_Dof_inds[i_vd][i_AR].end(), Elem_DOFs.begin() + i_vd*DOFs_per_vdim, Elem_DOFs.begin() + (i_vd + 1)*DOFs_per_vdim);
                }
            }
        }
    }

    // Function for creating the averaging operator from the avg space and mixed bilinear form
    void ParAveragingOperator::CreateAvgOperator(const double tol)
    {
        PrepareUnweightedAvgOperator(tol);
    }

    // Function for creating the averaging operator from the avg space and mixed bilinear form. Depending on the avg space/mixed bilinear form used,
    // the produced avg_mat is either \integral(1 * \basis) or \integral(\basis * \basis) summed over the elements of each AR box.
    void ParAveragingOperator::PrepareUnweightedAvgOperator(const double tol)
    {
        // Get the HypreParMatrix of the bilinear form
        HypreParMatrix *varf_avg_HPM = varf_avg->ParallelAssemble();

        // Get the diagonal matrix and corresponding information from the bilinear form HypreParMatrix
        SparseMatrix varf_Dmat; varf_avg_HPM->GetDiag(varf_Dmat);
        int N_cols_owned = varf_Dmat.NumCols();
        
        // Get the off-diagonal matrix and corresponding information from the bilinear form HypreParMatrix
        HYPRE_BigInt *varf_colMap; SparseMatrix varf_ODmat;
        varf_avg_HPM->GetOffd(varf_ODmat, varf_colMap);
        HYPRE_BigInt *varf_col_starts = varf_avg_HPM->GetColStarts(); // NOTE: This function only provides 2 components: the start (inclusive) and end (exclusive) columns of the rank. Do not try to access more than components 0 and 1



        // Initialize the matrices for holding the summed rows (according to AR_Dof_inds) of the diagonal and off-diagonal matrices.
        // Also initialize the vector for holding the AR areas
        SparseMatrix varf_Dmat_summed(N_rows, N_cols_owned);
        SparseMatrix varf_ODmat_summed(N_rows, varf_ODmat.NumCols());
        AR_areas = Array<double>(N_AR);
        
        // Create the summed-row matrices from the diagonal and off-diagonal matrices. Also alter the AR_areas vector
        CreateLocalAvgOperator(varf_Dmat, AR_Dof_inds, varf_Dmat_summed, AR_areas, tol);
        CreateLocalAvgOperator(varf_ODmat, AR_Dof_inds, varf_ODmat_summed, AR_areas, tol);

        // Finalize the summed-row matrices before integrating them into the HypreParMatrix
        varf_Dmat_summed.Finalize();
        varf_ODmat_summed.Finalize();



        // Send AR_areas to receiving_rank, and sum them on receiving_rank. Then, broadcast the result so that all ranks have the same AR_area
        int receiving_rank = 0;
        vector<Array<double>> AR_areas_all;
        if (rank == receiving_rank) {
            // Receive and store the data sent from the other ranks
            for (int i_rank = 0; i_rank < N_ranks; i_rank++) {
                if (i_rank == receiving_rank) { AR_areas_all.push_back( AR_areas ); } // Store the rank's own data
                else { // Receive the other ranks' data
                    AR_areas_all.push_back( Array<double>(N_AR) );
                    MPI_Recv(AR_areas_all[AR_areas_all.size() - 1].GetData(), N_AR, MPI_DOUBLE, i_rank, 0, comm, MPI_STATUS_IGNORE); // Receive the data
                }
            }

            // Sum together the AR_areas for each AR
            AR_areas = 0.0;
            for (int i_AR = 0; i_AR < N_AR; i_AR++) {
                for (int i_rank = 0; i_rank < N_ranks; i_rank++) { AR_areas[i_AR] += AR_areas_all[i_rank][i_AR]; }
            }

            // Broadcast the resulting AR_areas vector
            for (int i_rank = 0; i_rank < N_ranks; i_rank++) {
                if (i_rank != receiving_rank) { MPI_Send(AR_areas.GetData(), N_AR, MPI_DOUBLE, i_rank, 1, comm); }
            }
        }
        else {
            // Send data to receiving_rank
            MPI_Send(AR_areas.GetData(), N_AR, MPI_DOUBLE, receiving_rank, 0, comm);
            // Receive data from receiving_rank
            MPI_Recv(AR_areas.GetData(), N_AR, MPI_DOUBLE, receiving_rank, 1, comm, MPI_STATUS_IGNORE);
        }
        


        // Initialize the structure that will hold the data received from other ranks (as well as this rank's own data)
        vector<vector<pair<HYPRE_BigInt, double>>> recv_global_data;
        
        // For each row in the global averaging operator, pass the column indices and corresponding values to the rank declared to "own" the row
        for (int i_res = 0; i_res < N_rows; i_res++)
        {
            // Get the values/columns of the current row from both the diagonal and off-diagonal matrices. Convert columns indices from both matrices to global column indices
            Array<int> col_inds; vector<double> row_vec;
            {
                Array<int> col_inds_D; Vector row_D; varf_Dmat_summed.GetRow(i_res, col_inds_D, row_D);
                for (int i = 0; i < col_inds_D.Size(); i++) { col_inds.Append( col_inds_D[i] + varf_col_starts[0] ); }
                row_vec.insert(row_vec.end(), row_D.begin(), row_D.end());
            }
            {
                Array<int> col_inds_OD; Vector row_OD; varf_ODmat_summed.GetRow(i_res, col_inds_OD, row_OD);
                for (int i = 0; i < col_inds_OD.Size(); i++) { col_inds.Append( varf_colMap[col_inds_OD[i]] ); }
                row_vec.insert(row_vec.end(), row_OD.begin(), row_OD.end());
            }
            Vector row(row_vec.data(), row_vec.size());
            
            // Create local integer for the size of the data to be sent
            int N_data = col_inds.Size();
            

            // Receive/Send N_data, col_inds, and row data depending on the row ownership defined by recv_ranks
            if (rank == recv_ranks[i_res])
            {
                // Initialize variable for receiving the data from other ranks
                vector<pair<HYPRE_BigInt, double>> dummy_global_data;
                
                // Store the relevant part of the rank's own data
                {
                    vector<pair<HYPRE_BigInt, double>> temp; temp.resize(col_inds.Size());
                    for(int i = 0; i < col_inds.Size(); i++) { temp[i] = {col_inds[i], row[i]}; }
                    dummy_global_data.insert(dummy_global_data.end(), temp.begin(), temp.end());
                }
                
                // Receive and store the data sent from the other ranks
                for (int i_rank = 0; i_rank < N_ranks; i_rank++) {
                    if (i_rank != rank) {
                        // Receive the column data/column data size
                        vector<int> col_inds_recv; int N_data_recv;
                        MPI_Recv(&N_data_recv, 1, MPI_INT, i_rank, 0, comm, MPI_STATUS_IGNORE); // Receive the length of the data being sent
                        col_inds_recv.resize(N_data_recv);
                        MPI_Recv(col_inds_recv.data(), N_data_recv, MPI_INT, i_rank, 1, comm, MPI_STATUS_IGNORE); // Receive the data
                        
                        // Receive the value data (the size is the same as the column data)
                        vector<double> row_data_recv(N_data_recv);
                        MPI_Recv(row_data_recv.data(), N_data_recv, MPI_DOUBLE, i_rank, 2, comm, MPI_STATUS_IGNORE); // Receive the data
                        
                        // Save the column and value data in a vector of pairs for easy sorting later on
                        vector<pair<HYPRE_BigInt, double>> temp; temp.resize(N_data_recv);
                        for(int i = 0; i < N_data_recv; i++) { temp[i] = {col_inds_recv[i], row_data_recv[i]}; }
                        dummy_global_data.insert(dummy_global_data.end(), temp.begin(), temp.end());
                    }
                }

                // Store the received data for the current row i_res as a component in recv_global_data (i.e., each component of recv_global_data stores the data for a different row owned by the same rank)
                recv_global_data.push_back( dummy_global_data );
            }
            else
            {
                // Send the column data/size of column data
                MPI_Send(&N_data, 1, MPI_INT, recv_ranks[i_res], 0, comm); // Send the length of the data being sent
                MPI_Send(col_inds.GetData(), N_data, MPI_INT, recv_ranks[i_res], 1, comm); // Send the data
                // Send the value data
                MPI_Send(row.GetData(), N_data, MPI_DOUBLE, recv_ranks[i_res], 2, comm); // Send the data
            }
        }
        

        
        // Sort the pairs in recv_global_data by the column indices (first component of the pairs). Then, extract the column-value pair data
        // into global_cols and global_vals. To do this, sum values of similar column indices
        vector<vector<int>> global_cols; vector<vector<double>> global_vals;
        for (int i_DOF = 0; i_DOF < N_rows_owned; i_DOF++) {
            // Sort the component of recv_global_data by the column data
            sort(recv_global_data[i_DOF].begin(), recv_global_data[i_DOF].end(),
                [](const pair<HYPRE_BigInt, double> &a, const pair<HYPRE_BigInt, double> &b)
                { return a.first < b.first; }
            );
            
            // Go through the data of the i_DOF component of recv_global_data. Sum values of similar column indices
            vector<int> cols; vector<double> vals;
            for (int i_pair = 0; i_pair < recv_global_data[i_DOF].size(); i_pair++) {
                if (i_pair == 0) {
                    cols.push_back( recv_global_data[i_DOF][i_pair].first );
                    vals.push_back( recv_global_data[i_DOF][i_pair].second );
                    continue;
                }
                if (recv_global_data[i_DOF][i_pair].first == cols[cols.size() - 1]) {
                    vals[vals.size() - 1] += recv_global_data[i_DOF][i_pair].second;
                } else {
                    cols.push_back( recv_global_data[i_DOF][i_pair].first );
                    vals.push_back( recv_global_data[i_DOF][i_pair].second );
                }
            }
            global_cols.push_back( cols );
            global_vals.push_back( vals );
        }



        // The (global) column indices and corresponding values that each rank needs to build their owned row in the final averaging operator
        // is now in global_cols and global_vals. Divide this col/val data by which pairs would reside in the rank's diagonal and off-diagonal
        // matrices.

        // Initialize variables to store the diagonal and off-diagonal matrix col/val data
        vector<vector<double>> Dmat_vals, ODmat_vals;
        vector<Array<HYPRE_BigInt>> Dmat_cols, ODmat_cols;
        
        // For each owned row...
        for (int i_DOF = 0; i_DOF < N_rows_owned; i_DOF++)
        {
            // Initiate temporary variables to store the row data
            vector<double> Dmat_vals_row, ODmat_vals_row;
            Array<HYPRE_BigInt> Dmat_cols_row, ODmat_cols_row;
            
            // ...go through the col/val pairs.
            for (int i_pair = 0; i_pair < global_vals[i_DOF].size(); i_pair++) {
                // If the pair's global column index is within the range of vdofs owned by the rank, append it to the col/val data for the rank's diagonal matrix.
                // Subtract varf_col_starts[0] from the global column index to make it local
                if (global_cols[i_DOF][i_pair] >= varf_col_starts[0] && global_cols[i_DOF][i_pair] < varf_col_starts[1]) {
                    Dmat_cols_row.Append(global_cols[i_DOF][i_pair] - varf_col_starts[0]);
                    Dmat_vals_row.push_back(global_vals[i_DOF][i_pair]);
                }
                // If the pair's global column index is not within the range of vdofs owned by the rank, append it to the col/val data for the rank's off-diagonal matrix
                else {
                    ODmat_cols_row.Append(global_cols[i_DOF][i_pair]);
                    ODmat_vals_row.push_back(global_vals[i_DOF][i_pair]);
                }
            }

            // Store the col/val data for each row
            Dmat_cols.push_back( Dmat_cols_row );
            ODmat_cols.push_back( ODmat_cols_row );
            Dmat_vals.push_back( Dmat_vals_row );
            ODmat_vals.push_back( ODmat_vals_row );
        }
        

        
        // Column maps (global-to-local and local-to-global) need to be created for the off-diagonal matrix.
        vector<HYPRE_BigInt> OD_colMap; // Given a local column index (e.g., 0, 1, ..., ODmat.Width()), it outputs a global column index (e.g., 0, ..., varf_avg_HPM->GetGlobalNumCols())
        map<HYPRE_BigInt, int> OD_colMap_inv; // Given a global column index (e.g., 0, ..., varf_avg_HPM->GetGlobalNumCols()), it outputs a local column index (e.g., 0, 1, ..., ODmat.Width())
        
        // Compile the column indices for the off-diagonal matrix from all owned rows into OD_colMap
        for (int i_DOF = 0; i_DOF < ODmat_cols.size(); i_DOF++) { OD_colMap.insert(OD_colMap.end(), ODmat_cols[i_DOF].begin(), ODmat_cols[i_DOF].end()); }
        
        // Remove the duplicate global columns from OD_colMap
        {
            sort(OD_colMap.begin(), OD_colMap.end(),
                [](const HYPRE_BigInt &a, const HYPRE_BigInt &b)
                { return a < b; }
            );
            auto last = unique(OD_colMap.begin(), OD_colMap.end());
            OD_colMap.erase(last, OD_colMap.end());
        }
        
        // Create the inverse column map
        for (int i = 0; i < OD_colMap.size(); i++) { OD_colMap_inv[OD_colMap[i]] = i; }

        

        // Initiate the rank's diagonal and off-diagonal matrices
        Dmat = SparseMatrix(N_rows_owned, N_cols_owned);
        ODmat = SparseMatrix(N_rows_owned, OD_colMap.size());
        
        // Fill the rank's diagonal and off-diagonal matrices
        for (int i_DOF = 0; i_DOF < N_rows_owned; i_DOF++) {
            // Create an array of the local indices for the off-diagonal matrix
            Array<int> OD_cols_temp;
            for (int i = 0; i < ODmat_cols[i_DOF].Size(); i++) {
                OD_cols_temp.Append( OD_colMap_inv[ODmat_cols[i_DOF][i]] );
            }

            // Set the rows in the diagonal and off-diagonal matrices. Use 
            Vector Dmat_vals_Vec(Dmat_vals[i_DOF].data(), Dmat_vals[i_DOF].size());
            Dmat.SetRow(i_DOF, Dmat_cols[i_DOF], Dmat_vals_Vec);
            Vector ODmat_vals_Vec(ODmat_vals[i_DOF].data(), ODmat_vals[i_DOF].size());
            ODmat.SetRow(i_DOF, OD_cols_temp, ODmat_vals_Vec);
        }
        
        // Finalize the diagonal and off-diagonal matrices
        Dmat.Finalize(); ODmat.Finalize();

        
        
        // Create the HypreParMatrix for the averaging operator
        vector<HYPRE_BigInt> row_starts_local = {row_starts[rank], row_starts[rank + 1]};
        avg_mat = new HypreParMatrix(comm, (HYPRE_BigInt)N_rows, varf_avg_HPM->GetGlobalNumCols(), row_starts_local.data(), varf_avg_HPM->GetColStarts(), &Dmat, &ODmat, OD_colMap.data(), false);
        
        delete varf_avg_HPM; // This HypreParMatix came from ParallelAssemble, which calls for external owner (i.e., the caller deletes the matrix)
    }

    // Function for creating "averaging operators" from the diagonal and off-diagonal HypreParMatrices of the mixed bilinear form
    void ParAveragingOperator::CreateLocalAvgOperator(SparseMatrix &mat, vector<vector<vector<int>>> &row_inds, SparseMatrix &mat_summed_rows, Array<double> &AR_areas_, const double tol)
    {
        // row_inds contains row indices of VDOFs that belong to each vector dimension (vdim) and AR. We sum the rows of mat
        // over the indices provided in each inner-most vector of row_inds to gain one row in the mat_summed_rows, the averaging operator
        int set_row = 0;
        for (int i_vd = 0; i_vd < row_inds.size(); i_vd++) { // For each vector dimension...
            
            // This function will support general row_inds, but the second-level vectors of row_inds should often be N_AR in size (or at least equal in size).
            // A warning will be issued if the second-level vectors are not the same size
            if (i_vd > 0) {
                if (row_inds[i_vd].size() != row_inds[i_vd - 1].size() ) {
                    cout << "Rank " << rank << ": ParAveragingOperator::CreateLocalAvgOperator(): WARNING: the summed row indices provided do not show a consistent number of averaging regions between vector components (vdim). ";
                    cout << "row_inds[" << i_vd << "] = " << row_inds[i_vd].size() << ", row_inds[" << i_vd - 1 << "] = " << row_inds[i_vd - 1].size() << "." << endl; } }
            
            for (int i_AR = 0; i_AR < row_inds[i_vd].size(); i_AR++) { // ...and each averaging region...
                // ...sum over the corresponding rows in mat and insert it into mat_summed_rows
                SumColsGivenRows(mat, row_inds[i_vd][i_AR], mat_summed_rows, set_row, tol);
                
                // Sum the row_sum, as it pertains to the area of the averaging region (only do this for the first vector component (vdim); it should be the same for the other vdim)
                Array<int> col_inds; Vector row; mat_summed_rows.GetRow(set_row, col_inds, row);
                if (i_vd == 0) { AR_areas_[i_AR] += row.Sum(); }
                
                // Iterate the row being set in mat_summed_rows
                set_row += 1;
            }
        }
    }
    
    // Function for summing indicated rows of a matrix mat and inserting the sum into a matrix mat_summed_rows
    void ParAveragingOperator::SumColsGivenRows(SparseMatrix &mat, vector<int> &row_inds, SparseMatrix &mat_summed_rows, int set_ind, const double tol)
    {
        // Initiate a Vector for summing the rows in row_inds
        Vector row_sum(mat.NumCols()); row_sum = 0.0;
        
        // Get the row_inds rows from mat and add them to row_sum
        for (int i_DOF = 0; i_DOF < row_inds.size(); i_DOF++) {
            Array<int> col_inds; Vector row; mat.GetRow(row_inds[i_DOF], col_inds, row); // Get the row from mat
            for (int i_sum = 0; i_sum < col_inds.Size(); i_sum++) { row_sum[col_inds[i_sum]] += row[i_sum]; } // Add the row to row_sum
        }

        // If a component of row_sum is really small, just set it to 0; it is likely numerical error. Simultaneously create an Array of column indices (i.e., {0, ..., N_cols})
        Array<int> all_col_indices(mat.NumCols());
        for (int i = 0; i < row_sum.Size(); i++) { if (abs(row_sum[i]) < tol) { row_sum[i] = 0.0; } all_col_indices[i] = i; }
        
        // Insert the row_sum in mat_summed_rows, the matrix describing the averaging operator (or part of it in the case of MPI-parallel)
        mat_summed_rows.SetRow(set_ind, all_col_indices, row_sum);
    }

    // Static method for computing the AR porosities
    void ParAveragingOperator::ComputePorosities(const vector<double> &pore_areas, const vector<double> &AR_areas, vector<double> &porosities)
    {
        assert(pore_areas.size() == AR_areas.size());
        porosities.clear();
        for (int i = 0; i < pore_areas.size(); i++) {
            porosities.push_back( pore_areas[i] / AR_areas[i] );
        }
    }

    // Function for applying the averaging operator
    void ParAveragingOperator::ApplyAvgOperator(const Vector &x, Vector &avg)
    {
        // Define a porosity vector of ones
        vector<double> porosities;
        for (int i_AR = 0; i_AR < N_AR; i_AR++) { porosities.push_back(1.0); }
        
        // Compute the average with a porosity of 1.0
        ParAveragingOperator::ApplyAvgOperator(x, avg, porosities);
    }
    void ParAveragingOperator::ApplyAvgOperator(const Vector &x, Vector &avg, const vector<double> &porosities)
    {
        assert(porosities.size() == N_AR);

        // Reset the avg output vector
        avg.SetSize(N_rows_owned); avg = 0.0;

        // Apply the averaging operator
        avg_mat->Mult(x, avg);
        
        for (int i_row = 0; i_row < avg.Size(); i_row++) {
            avg.Elem(i_row) *= porosities[local_row_to_AR_index[i_row]] / AR_areas[local_row_to_AR_index[i_row]];
        }
    }

    // Function for piecing together the rank-local vectors that result from applying the averaging operator in parallel.
    // The combine vector is kept on a single rank (by default, this rank is 0)
    void ParAveragingOperator::GetSerialParAveragedVector(double* data, Vector &combined_data, int receiving_rank)
    {
        // Send AR_areas to receiving_rank, and sum them on receiving_rank. Then, broadcast the result so that all ranks have the same AR_area
        if (rank == receiving_rank)
        {
            // Receive and store the data sent from the other ranks
            vector<vector<double>> data_all;
            for (int i_rank = 0; i_rank < N_ranks; i_rank++) {
                if (i_rank == receiving_rank) { // Store the rank's own data
                    vector<double> temp(data, data + N_rows_owned);
                    data_all.push_back( temp );
                }
                else { // Receive the other ranks' data
                    int N_recv_data;
                    MPI_Recv(&N_recv_data, 1, MPI_INT, i_rank, 0, comm, MPI_STATUS_IGNORE);
                    data_all.push_back( vector<double>(N_recv_data) );
                    MPI_Recv(data_all[data_all.size() - 1].data(), N_recv_data, MPI_DOUBLE, i_rank, 1, comm, MPI_STATUS_IGNORE);
                }
            }
            
            // Combine the received data in the order in which rows are owned by rank
            combined_data.SetSize(recv_ranks.size()); combined_data = 0.0;
            vector<int> ind_count(N_ranks); for (int i = 0; i < N_ranks; i++) { ind_count[i] = 0; }
            for (int i_row = 0; i_row < recv_ranks.size(); i_row++) {
                combined_data[i_row] = data_all[ recv_ranks[i_row] ][ ind_count[recv_ranks[i_row]] ];
                ind_count[recv_ranks[i_row]] += 1;
            }
        }
        else
        {
            // Send the size of the data to be received by receiving_rank
            MPI_Send(&N_rows_owned, 1, MPI_INT, receiving_rank, 0, comm);
            // Send data to receiving_rank
            MPI_Send(data, N_rows_owned, MPI_DOUBLE, receiving_rank, 1, comm);
        }
    }    
#endif





// ===============================================================
//   Define the functions for the LinearTimeDependentOperator class.
// ===============================================================
// Define the function that prepares the operator, preconditioner, and solver for Explicitr Euler time stepping
void LinearTimeDependentOperator::PrepareExplicitEuler()
{
    // For reference:   M * du/dt + A * u = b
    // We solve for du/dt.

    // Solve parameters
    int maxIter(5000);
    real_t rtol(1.0e-7);
    real_t atol(1.0e-9);

    // Create the preconditioner (Note: 'dt' is in the preconditioner)
    M_inv = new DSmoother(M);
    M_inv->iterative_mode = false;
    PC->SetDiagonalBlock(0, M_inv);

    // Build the block operator
    Op->SetBlock(0, 0, &M);
    
    // For Explicit Euler
    explicitSolver.iterative_mode = false;
    explicitSolver.SetRelTol(rtol);
    explicitSolver.SetAbsTol(atol);
    explicitSolver.SetMaxIter(maxIter);
    explicitSolver.SetPrintLevel(0);
    explicitSolver.SetOperator(M);
    explicitSolver.SetPreconditioner(*PC);
}

// Define the function that prepares the operator, preconditioner, and solver for Implicit Euler time stepping
void LinearTimeDependentOperator::PrepareImplicitEuler(string solver_type, string PC_type)
{
    // For reference:   (M + A*dt) * du/dt + A * u = b
    // We solve for du/dt.
    
    // Solve parameters
    int maxIter(8000);
    //real_t rtol(1.0e-7);
    //real_t atol(1.0e-9);
    real_t rtol(1.0e-4);
    real_t atol(1.0e-6);
    
    // Create F = M + A dt  (i.e., the inverted matrix)
    F = Add(1.0, M, dt, A);
    F->Finalize();

    // Create the preconditioner (Note: 'dt' is in the preconditioner)
    if (PC_type == "DSmoother")
    {
        F_inv = new DSmoother(*F);
        F_inv->iterative_mode = false;
        PC->SetDiagonalBlock(0, F_inv);
    }
    else if (PC_type == "BlockILU")
    {
        F_inv = new BlockILU(*F);
        F_inv->iterative_mode = false;
        PC->SetDiagonalBlock(0, F_inv);
    }
    else
    {
        cout << "mfem_util.cpp: LinearTimeDependentOperator::PrepareImplicitEuler(): CRITICAL ERROR: PC_type not recognized." << endl;
        exit(1);
    }


    // Build the block operator (Note: 'dt' is in the operator)
    Op->SetBlock(0, 0, F);
    
    // Prepare the solver
    if (solver_type == "CGSolver")
    {
        cgSolver.iterative_mode = true;
        cgSolver.SetRelTol(rtol);
        cgSolver.SetAbsTol(atol);
        cgSolver.SetMaxIter(maxIter);
        cgSolver.SetPrintLevel(0);
        cgSolver.SetOperator(*Op);
        cgSolver.SetPreconditioner(*PC);
        implicitSolver = &cgSolver;
    }
    else if (solver_type == "GMRESSolver")
    {
        gmresSolver.iterative_mode = true;
        gmresSolver.SetRelTol(rtol);
        gmresSolver.SetAbsTol(atol);
        gmresSolver.SetMaxIter(maxIter);
        gmresSolver.SetPrintLevel(0);
        gmresSolver.SetOperator(*Op);
        gmresSolver.SetPreconditioner(*PC);
        implicitSolver = &gmresSolver;
    }
    else
    {
        cout << "mfem_util.cpp: LinearTimeDependentOperator::PrepareImplicitEuler(): CRITICAL ERROR: solver_type not recognized." << endl;
        exit(1);
    }
}

// Define the explicit Mult function for the operator
void LinearTimeDependentOperator::Mult(const Vector &u, Vector &du_dt) const
{
    // For reference:   M * du/dt + A * u = b
    // We solve for du/dt.

    A.Mult(u, du_dt); // Compute A u
    du_dt *= -1; // Move A u to RHS
    du_dt += b; // -A u + b
    explicitSolver.Mult(du_dt, du_dt); // Compute du/dt = (M^{-1}) * (-A u + b)
}

// Define the implicit Mult function for the operator
void LinearTimeDependentOperator::ImplicitSolve(const real_t dt, const Vector &u, Vector &du_dt)
{
    // For reference:   (M + A*dt) * du/dt + A * u = b
    // We solve for du/dt.

    A.Mult(u, du_dt); // Compute  A u
    du_dt *= -1; // Move A u to RHS
    du_dt += b; // -A u + b
    implicitSolver->Mult(du_dt, du_dt); // Compute du/dt = ((M + A*dt)^{-1}) * (-A u + b)
}





// Function for projecting an average solution vector onto a GridFunction. We assume the
// GridFunction has a pore-scale finite element space. Use a DG space
void CreateAvgSolGridFunction(const Vector &sol, GridFunction &u, const string avg_area_type)
{
    // Save relevant GridFunction entities to variables to make manipulation easier
    FiniteElementSpace *fes = u.FESpace();
    
    Mesh *mesh = fes->GetMesh();
    int N_elem = mesh->GetNE();
    
    // Ensure the solution is all zero to start
    u = 0.0;

    // Get the indices of the elements inside each averaging region
    for (int i_Elem = 0; i_Elem < N_elem; i_Elem++)
    {
        int AR_ind;
        if (avg_area_type == "AR") { AR_ind = i_Elem;}
        else { AR_ind = mesh->GetAttribute(i_Elem) - 1;} // Get the AR index of the element

        // Get the vdofs of the element
        Array<int> vdofs;
        fes->GetElementVDofs(i_Elem, vdofs);

        // Set each vdof to the average solution
        for (int i_qp = 0; i_qp < vdofs.Size(); i_qp++) {
            u(vdofs[i_qp]) = sol.Elem(AR_ind);
        }
    }
}





// ===============================================================
//   Define the functions for the ResultsSaver class.
// ===============================================================
// Function for obtaining the closure residuals from the solution block vector
void ResultsSaver::ObtainClosureResiduals(const BlockVector &sol_BLK, const vector<int> &block_IDs_for_saving_residuals, const vector<int> &fespace_vdim, int N_AR)
{
    // Make sure the required variables have been set
    assert (N_residual_sets > 0);
    assert (block_IDs_for_saving_residuals.size() == N_residual_sets);
    
    // Go through the closure residual blocks in the solution and save the residuals
    for (int i_BI = 0; i_BI < N_residual_sets; i_BI++) // i_BI < block_IDs_for_saving_residuals.size()
    {
        residuals.push_back( vector<vector<double>>() );
        if (solve_mode == "serial") { cout << "Closure residuals (placed on RHS): Equation " << i_BI << endl; }
        for (int i_vd = 0; i_vd < fespace_vdim[i_BI]; i_vd++)
        {
            residuals[i_BI].push_back( vector<double>() );
            for (int i_AR = 0; i_AR < N_AR; i_AR++)
            {
                residuals[i_BI][i_vd].push_back( sol_BLK.GetBlock(block_IDs_for_saving_residuals[i_BI]).Elem(i_AR + i_vd * N_AR) );
                if (solve_mode == "serial") { cout << residuals[i_BI][i_vd][i_AR] << endl; }
            }
        }
    }

    // Initialize the JSONDicts with fespace_vdim
    InitializeJSONDicts_EqLevel(fespace_vdim);
}
void ResultsSaver::ObtainClosureResiduals(const Vector &sol_BLK, const vector<int> &fespace_vdim, int N_AR)
{
    // Make sure the required variables have been set
    assert (N_residual_sets > 0);
    assert (N_residual_sets == 1);
    
    // Go through the closure residual blocks in the solution and save the residuals
    for (int i_BI = 0; i_BI < N_residual_sets; i_BI++) // i_BI < block_IDs_for_saving_residuals.size()
    {
        residuals.push_back( vector<vector<double>>() );
        if (solve_mode == "serial") { cout << "Closure residuals (placed on RHS): Equation " << i_BI << endl; }
        for (int i_vd = 0; i_vd < fespace_vdim[i_BI]; i_vd++)
        {
            residuals[i_BI].push_back( vector<double>() );
            for (int i_AR = 0; i_AR < N_AR; i_AR++)
            {
                residuals[i_BI][i_vd].push_back( sol_BLK.Elem(i_AR + i_vd * N_AR) );
                if (solve_mode == "serial") { cout << residuals[i_BI][i_vd][i_AR] << endl; }
            }
        }
    }

    // Initialize the JSONDicts with fespace_vdim
    InitializeJSONDicts_EqLevel(fespace_vdim);
}
void ResultsSaver::ObtainClosureResiduals(const BlockVector &sol_BLK, const vector<int> &block_IDs_for_saving_residuals, const vector<int> &fespace_vdim, int N_AR, const vector<int> AR_inds)
{
    // Make sure the required variables have been set
    assert (N_residual_sets > 0);
    assert (block_IDs_for_saving_residuals.size() == N_residual_sets);
    
    // Go through the closure residual blocks in the solution and save the residuals
    for (int i_BI = 0; i_BI < N_residual_sets; i_BI++) // i_BI < block_IDs_for_saving_residuals.size()
    {
        residuals.push_back( vector<vector<double>>() );
        if (solve_mode == "serial") { cout << "Closure residuals (placed on RHS): Equation " << i_BI << endl; }
        for (int i_vd = 0; i_vd < fespace_vdim[i_BI]; i_vd++)
        {
            residuals[i_BI].push_back( vector<double>() );
            for (int i_AR = 0; i_AR < N_AR; i_AR++)
            {
                bool isZero = true;
                int i_AR_loc;
                for (int i_loc = 0; i_loc < AR_inds.size(); i_loc++) { if (i_AR == AR_inds[i_loc]) { isZero = false; i_AR_loc = i_loc; break; } }
                if (isZero) { residuals[i_BI][i_vd].push_back( 0.0 ); }
                else { residuals[i_BI][i_vd].push_back( sol_BLK.GetBlock(block_IDs_for_saving_residuals[i_BI]).Elem(i_AR_loc + i_vd * AR_inds.size()) ); }
                if (solve_mode == "serial") { cout << residuals[i_BI][i_vd][i_AR] << endl; }
            }
        }
    }

    // Initialize the JSONDicts with fespace_vdim
    InitializeJSONDicts_EqLevel(fespace_vdim);
}

// Function for initializing the JSONDict vectors at the equation level.
// This includes eq_dicts, res_dicts, CFN_dicts, and comp_dicts
void ResultsSaver::InitializeJSONDicts_EqLevel(const vector<int> &fespace_vdim)
{
    // Initialize the dictionary vectors
    for (int i_RS = 0; i_RS < N_residual_sets; i_RS++) { //i < eq_keys.size()
        eq_dicts.push_back(JSONDict());
        CFN_dicts.push_back(JSONDict());
        res_dicts.push_back(JSONDict());
        alpha_dicts.push_back(JSONDict());
        beta_dicts.push_back(JSONDict());
        gamma_dicts.push_back(JSONDict());
        
        comp_dicts.push_back(vector<JSONDict>());
        alpha_comp_dicts.push_back(vector<JSONDict>());
        beta_comp_dicts.push_back(vector<JSONDict>());
        gamma_comp_dicts.push_back(vector<JSONDict>());

        for (int i_vd = 0; i_vd < fespace_vdim[i_RS]; i_vd++) {
            comp_dicts[i_RS].push_back(JSONDict());
            alpha_comp_dicts[i_RS].push_back(JSONDict());
            beta_comp_dicts[i_RS].push_back(JSONDict());
            gamma_comp_dicts[i_RS].push_back(JSONDict());
        }
    }
    initializedJSONDicts = true;
}

// Function for filling the residual parameters' dictionaries
void ResultsSaver::StoreResidualParameters(const string &param_name, const vector<vector<vector<double>>> &param)
{
    // Make sure the required variables have been set
    assert (N_residual_sets > 0);
    assert (param.size() == N_residual_sets);
    assert (param_name == "gamma");
    gammas = param;
}
void ResultsSaver::StoreResidualParameters(const string &param_name, const vector<vector<double>> &param, int N_AR)
{
    // Make sure the required variables have been set
    assert (N_residual_sets > 0);
    assert (param.size() == N_residual_sets);
    assert (param_name == "alpha" || param_name == "beta");

    vector<vector<vector<double>>> *ref_vec = nullptr;
    if (param_name == "alpha") { ref_vec = &alphas; }
    else if (param_name == "beta") { ref_vec = &betas; }    

    for (int i_BI = 0; i_BI < N_residual_sets; i_BI++) {
        ref_vec->push_back( vector<vector<double>>() );
        for (int i_vd = 0; i_vd < param[i_BI].size(); i_vd++) {
            (*ref_vec)[i_BI].push_back( vector<double>() );
            for (int i_AR = 0; i_AR < N_AR; i_AR++) {
                (*ref_vec)[i_BI][i_vd].push_back( param[i_BI][i_AR] );
            }
        }
    }
}
void ResultsSaver::StoreResidualParameters(const string &param_name, const vector<vector<double>> &param, int N_AR, const vector<int> AR_inds)
{
    // Make sure the required variables have been set
    assert (N_residual_sets > 0);
    assert (param.size() == N_residual_sets);
    assert (param_name == "alpha" || param_name == "beta");

    vector<vector<vector<double>>> *ref_vec = nullptr;
    if (param_name == "alpha") { ref_vec = &alphas; }
    else if (param_name == "beta") { ref_vec = &betas; }    

    for (int i_BI = 0; i_BI < N_residual_sets; i_BI++) {
        ref_vec->push_back( vector<vector<double>>() );
        for (int i_vd = 0; i_vd < param[i_BI].size(); i_vd++) {
            (*ref_vec)[i_BI].push_back( vector<double>() );
            for (int i_AR = 0; i_AR < N_AR; i_AR++) {
                bool isZero = true;
                for (int i_loc = 0; i_loc < AR_inds.size(); i_loc++) { if (i_AR == AR_inds[i_loc]) { isZero = false; break; } }
                if (isZero) { (*ref_vec)[i_BI][i_vd].push_back( 0.0 ); }
                else { (*ref_vec)[i_BI][i_vd].push_back( param[i_BI][i_AR] ); }
            }
        }
    }
}

// Function for filling component dictionaries (i.e., either residuals, alphas, betas, gammas, etc.)
void ResultsSaver::FillComponentDictionary(vector<vector<JSONDict>> &dicts, vector<vector<vector<double>>> &val)
{
    for (int i_eq = 0; i_eq < N_residual_sets; i_eq++) {
        for (int i_comp = 0; i_comp < dicts[i_eq].size(); i_comp++) {
            dicts[i_eq][i_comp][AR_number] = val[i_eq][i_comp];
        }
    }
}

// Function for filling parameter dictionaries (i.e., either res_dicts, alpha_dicts, beta_dicts, gamma_dicts, etc.)
void ResultsSaver::FillParamDictionary(vector<JSONDict> &dicts, vector<vector<JSONDict>> &val)
{
    for (int i_eq = 0; i_eq < N_residual_sets; i_eq++) {
        for (int i_comp = 0; i_comp < val[i_eq].size(); i_comp++) {
            dicts[i_eq]["component " + to_string(i_comp)] = &val[i_eq][i_comp];
        }
    }
}

// Function for loading the residuals and corresponding information/labeling into the various
// dictionaries, compiling them into a main residual dictionary, and saving that dictionary
// as a text file
void ResultsSaver::SaveResidualDictionary(const string &residual_output_file_path)
{
    // Make sure the required variables have been set
    assert (N_residual_sets > 0);
    assert (eq_keys.size() > 0);
    assert (recursive_iter > -1);
    assert (AR_number != "");
    assert (initializedJSONDicts);

    // See if at least one of the parameters was filled. If so, save them in the JSON dict
    bool saveParams = !(alphas.size() == 0 && betas.size() == 0 && gammas.size() == 0);

    // If the residual file exists, load it into the dictionaries
    ifstream file(residual_output_file_path);
    if (file) { LoadPreexistingResidualsFile(residual_output_file_path); }

    // Load up the residual component dictionaries
    FillComponentDictionary(comp_dicts, residuals);
    if (saveParams) {
        FillComponentDictionary(alpha_comp_dicts, alphas);
        FillComponentDictionary(beta_comp_dicts, betas);
        FillComponentDictionary(gamma_comp_dicts, gammas);
    }
    
    // Create the new residual and file name dictionaries
    FillParamDictionary(res_dicts, comp_dicts);
    if (saveParams) {
        FillParamDictionary(alpha_dicts, alpha_comp_dicts);
        FillParamDictionary(beta_dicts, beta_comp_dicts);
        FillParamDictionary(gamma_dicts, gamma_comp_dicts);
    }
    
    for (int i_eq = 0; i_eq < N_residual_sets; i_eq++) { CFN_dicts[i_eq][AR_number] = CFNs[i_eq]; }
    
    // Save the new dictionaries into the loaded dictionaries
    for (int i_eq = 0; i_eq < N_residual_sets; i_eq++) {
        eq_dicts[i_eq]["closure file name"] = &CFN_dicts[i_eq];
        eq_dicts[i_eq]["residuals"] = &res_dicts[i_eq];
        if (saveParams) {
            eq_dicts[i_eq]["alphas"] = &alpha_dicts[i_eq];
            eq_dicts[i_eq]["betas"] = &beta_dicts[i_eq];
            eq_dicts[i_eq]["gammas"] = &gamma_dicts[i_eq];
        }
        iter_dict[eq_keys[i_eq]] = &eq_dicts[i_eq];
    }
    time_func_dict["iter_" + to_string(recursive_iter)] = &iter_dict;
    residual_dict[sim_key] = &time_func_dict;
    residual_dict.saveToFile(residual_output_file_path);
}

// Function for loading the closure residuals/closure file names from a preexisting residual
// results file into the dictionaries for appending the new closure residuals
void ResultsSaver::LoadPreexistingResidualsFile(const string &residual_output_file_path)
{
    residual_dict.loadFromFile(residual_output_file_path);
    residual_dict.getValue(sim_key, time_func_dict);
    time_func_dict.getValue("iter_" + to_string(recursive_iter), iter_dict);
    for (int i_eq = 0; i_eq < N_residual_sets; i_eq++) // i_eq < eq_keys.size()
    {
        iter_dict.getValue(eq_keys[i_eq], eq_dicts[i_eq]);
        
        eq_dicts[i_eq].getValue("closure file name", CFN_dicts[i_eq]);
        
        eq_dicts[i_eq].getValue("residuals", res_dicts[i_eq]);
        eq_dicts[i_eq].getValue("alphas", alpha_dicts[i_eq]);
        eq_dicts[i_eq].getValue("betas", beta_dicts[i_eq]);
        eq_dicts[i_eq].getValue("gammas", gamma_dicts[i_eq]);
        for (int i_comp = 0; i_comp < comp_dicts[i_eq].size(); i_comp++)
        {
            res_dicts[i_eq].getValue("component " + to_string(i_comp), comp_dicts[i_eq][i_comp]);
            if (!alpha_dicts[i_eq].empty()) { alpha_dicts[i_eq].getValue("component " + to_string(i_comp), alpha_comp_dicts[i_eq][i_comp]); }
            if (!beta_dicts[i_eq].empty()) { beta_dicts[i_eq].getValue("component " + to_string(i_comp), beta_comp_dicts[i_eq][i_comp]); }
            if (!gamma_dicts[i_eq].empty()) { gamma_dicts[i_eq].getValue("component " + to_string(i_comp), gamma_comp_dicts[i_eq][i_comp]); }
        }
    }
}





// ===============================================================
//   Define the functions for the MeshManager class.
// ===============================================================
// Function for making the mesh periodic (if isPeriodic calls for it). Note, this function creates a new Mesh object that is periodic and deletes the original
void MeshManager::MakePeriodic(const vector<int> &isPeriodic, const vector<double> &L)
{
    // Create a periodic mesh if required
    for (int i = 0; i < parent_mesh->Dimension(); i++) {
        if (isPeriodic[i]) {
            Vector translation(3);
            for (int j = 0; j < translation.Size(); j++) { translation[j] = 0.0; }
            translation[i] = L[i];
            translations.push_back( translation );
        }
    }
    if (translations.size() > 0) {
        Mesh periodic_mesh = Mesh::MakePeriodic(*parent_mesh, parent_mesh->CreatePeriodicVertexMapping(translations));
        parent_mesh.reset();
        parent_mesh = std::make_shared<Mesh>(periodic_mesh);
        if (!usingLocalMesh) { UseParentMesh(); }
        // This one liner might be able to replace the first 3 lines above in the if statment
        //parent_mesh = std::make_shared<Mesh>(Mesh::MakePeriodic(*parent_mesh, parent_mesh->CreatePeriodicVertexMapping(translations)));
    }
}

// Function for creating a "local" mesh around a specified AR index from the parent mesh in the class
void MeshManager::MakeLocalMesh(vector<int> central_AR_inds, int N_neighbor_layers, const vector<int> &AR_tags, const vector<vector<int>> &AR_neighbors)
{
    // Initiate variables
    for (int i_aA = 0; i_aA < central_AR_inds.size(); i_aA++) {
        AR_tags_local.push_back( AR_tags[central_AR_inds[i_aA]] );
    }
    central_AR_inds_local = AR_tags_local; // central_AR_inds_local is initially filled with the corresponding tags, and then replaced later with the indices
    vector<int> AR_inds = central_AR_inds;
    AR_inds_loc2glob = AR_inds;
    
    for (int i_nl = 0; i_nl < N_neighbor_layers; i_nl++) {
        // Get the AR tags of the neighbors to the ARs in AR_inds
        vector<int> new_AR_tags;
        for (int i_ari = 0; i_ari < AR_inds.size(); i_ari++) {
            new_AR_tags.insert(new_AR_tags.end(), AR_neighbors[AR_inds[i_ari]].begin(), AR_neighbors[AR_inds[i_ari]].end());
        }

        // Delete the duplicates in new_AR_tags
        std::sort(new_AR_tags.begin(), new_AR_tags.end());
        auto last = std::unique(new_AR_tags.begin(), new_AR_tags.end());
        new_AR_tags.erase(last, new_AR_tags.end());

        // Store new_AR_tags (i.e., the neighbor AR tags) in AR_kept, and clear AR_inds
        AR_tags_local.insert(AR_tags_local.end(), new_AR_tags.begin(), new_AR_tags.end());
        AR_inds.clear();
        
        // Get the indices of the AR tags in new_AR_tags for the next iteration
        for (int i_art = 0; i_art < new_AR_tags.size(); i_art++) {
            bool wasFound = false;
            for (int i_ar = 0; i_ar < AR_tags.size(); i_ar++) {
                if (new_AR_tags[i_art] == AR_tags[i_ar]) {
                    AR_inds.push_back( i_ar );
                    wasFound = true;
                    break;
                }
            }
            if (!wasFound) { cerr << "MeshManager::MakeLocalMesh(): CRITICAL ERROR: Could not find the AR index from AR tag " << new_AR_tags[i_art] << "." << endl; exit(1); }
        }

        // Store AR_inds in AR_inds_loc2glob
        AR_inds_loc2glob.insert(AR_inds_loc2glob.end(), AR_inds.begin(), AR_inds.end());
    }
    
    // Delete the duplicates in AR_kept
    std::sort(AR_tags_local.begin(), AR_tags_local.end());
    auto last = std::unique(AR_tags_local.begin(), AR_tags_local.end());
    AR_tags_local.erase(last, AR_tags_local.end());
    
    // Delete the duplicates in AR_inds_loc2glob
    std::sort(AR_inds_loc2glob.begin(), AR_inds_loc2glob.end());
    auto last2 = std::unique(AR_inds_loc2glob.begin(), AR_inds_loc2glob.end());
    AR_inds_loc2glob.erase(last2, AR_inds_loc2glob.end());
    
    // Record the number of AR in the local mesh
    N_AR_local = AR_tags_local.size();
    // Get the indices of the central AR in AR_tags_local. These are are new 'active_AR' in the solver codes
    for (int i_cal = 0; i_cal < central_AR_inds_local.size(); i_cal++) {
        bool wasfound = false;
        for (int i_AR = 0; i_AR < AR_tags_local.size(); i_AR++) {
            if (AR_tags_local[i_AR] == central_AR_inds_local[i_cal]) { central_AR_inds_local[i_cal] = i_AR; wasfound = true; break; }
        }
        if (!wasfound) { cerr << "mfem_util.cpp: MeshManager::MakeLocalMesh(): Could not find a tag in 'central_AR_inds_local' in 'AR_tags_local', but should be able to. Tag was: " << central_AR_inds_local[i_cal] << endl; exit(1); }
    }

    // Create the submesh with the kept ARs
    Array<int> keep(AR_tags_local.data(), AR_tags_local.size());
    local_mesh.reset();
    local_mesh = std::make_shared<SubMesh>(SubMesh::CreateFromDomain(*mesh, keep));
    
    //return SubMesh::CreateFromDomain(*mesh, keep);
    
    //SubMesh smesh = SubMesh::CreateFromDomain(*mesh, keep);
    //mesh = new Mesh(smesh);
    //mesh = new Mesh(smesh);
    //SubMesh *mm = &smesh;
}

// Function for clearing and reassigning variable "mesh" to the object pointed at by "parent_mesh"
void MeshManager::UseParentMesh()
{   
    mesh.reset();
    mesh = parent_mesh;
    usingLocalMesh = false;
}

// Function for clearing and reassigning variable "mesh" to the object pointed at by "local_mesh", but as a Mesh object, not a SubMesh object
void MeshManager::UseLocalMesh()
{
    mesh.reset();
    mesh = std::make_shared<Mesh>(Mesh(*local_mesh.get()));
    usingLocalMesh = true;
}


// Define functions for MPI-parallel
#ifdef MPI_BUILD
    // Function for creating a parallel mesh by partitioning the serial mesh (serial mesh could be periodic/local)
    void MeshManager::MakeParallelMesh(MPI_Comm comm)
    {
        //par_mesh = std::make_shared<ParMesh>(ParMesh(comm, *mesh));
        partitioning = mesh->GeneratePartitioning(Mpi::WorldSize());
        par_mesh = std::make_shared<ParMesh>(ParMesh(comm, *mesh, partitioning));
    }

    // Function for clearing and reassigning variable "mesh" to the object pointed at by "par_mesh", but as a Mesh object, not a ParMesh object
    //void MeshManager::UseParallelMesh()
    //{
    //    mesh.reset();
    //    mesh = std::make_shared<ParMesh>(ParMesh(*par_mesh.get()));
    //    // if (usingLocalMesh); // clear serial meshes?
    //}
#endif





// ===============================================================
//   Define the functions for the GridFunctionManager class.
// ===============================================================
// Function for loading a grid function from a file with path "gf_file_path" onto a mesh with file path "gfc_mesh_file_path"
void GridFunctionManager::LoadGridFunction(string gf_file_path, Mesh *gf_mesh)
{
    // Create an ifgzstream pointer to read in the file containing the grid function
    gf_ifgz_stream = std::make_unique<ifgzstream>(gf_file_path);
    if (!(*gf_ifgz_stream)) {
        cerr << "mfem_utils.cpp: GridFunctionCoefficientManager::LoadGridFunction(): CRITICAL ERROR: Unable to open grid function file. File path is: " << gf_file_path << endl;
        exit(1);
    }

    // Define the grid function using the mesh and the stream
    gf = std::make_shared<GridFunction>(gf_mesh, *gf_ifgz_stream);
    
    // For debugging
    //gf->Save("fluid_velocity_verification_plot.gf");
    //gf_mesh->Save("fluid_velocity_verification_mesh.mesh");
    //exit(1);

    // Save a non-const pointer to the grid function's finite element space (this is needed for some MFEM functions, like LinearForm)
    fes = gf->FESpace();

    // Load other grid function information into the class for easier manipulation later
    UpdateGridFunctionInfo();
}

// Function for loading other grid function information into the class
void GridFunctionManager::UpdateGridFunctionInfo()
{
    // Determine the name of the grid function's finite element collection 
    const FiniteElementCollection *fec = gf->FESpace()->FEColl();
    if (dynamic_cast<const H1_FECollection*>(fec)) { fec_name = "H1"; }
    //else if (dynamic_cast<const DG_FECollection*>(fec)) { fec_name = "DG"; }
    //else if (dynamic_cast<const ND_FECollection*>(fec)) { fec_name = "ND"; }
    //else if (dynamic_cast<const RT_FECollection*>(fec)) { fec_name = "RT"; }
    else {
        cerr << "mfem_util.cpp: GridFunctionCoefficientManager::UpdateGridFunctionInfo(): CRITICAL ERROR: Unknown finite element collection (fec) was provided. If the fec is valid, it should be added as an option in this function. Provided fec: " << fec->Name() << endl;
        exit(1);
    }

    // Determine the order of the grid function's finite element collection
    fec_order = gf->FESpace()->FEColl()->GetOrder();

    // Determine the dimension (vdim) of the grid function's finite element space
    fes_dim = gf->FESpace()->GetVDim();
}

// Function for transfering the gridfunction "gf" to the local mesh provided
void GridFunctionManager::MakeLocalGridFunction(SubMesh *local_mesh)
{
    // Make the local finite element collection and finite element space
    MakeLocalFECAndFES(local_mesh);

    // Make the local grid function using the local finite element space
    gf_local = std::make_shared<GridFunction>(fes_local.get());
    
    // Make the transfer map for the local mesh, and map "gf" to "gf_local"
    TransferMap map(*gf, *gf_local);
    map.Transfer(*gf, *gf_local);
}

// Function for making a local finite element collection and finite element space for the grid function coefficient
void GridFunctionManager::MakeLocalFECAndFES(SubMesh *local_mesh)
{
    // Create the finite element collection based on the class's info about the loaded grid function
    if (fec_name == "H1") { fec_local = std::make_unique<H1_FECollection>(fec_order, local_mesh->Dimension()); }
    else {
        cerr << "mfem_util.cpp: GridFunctionCoefficientManager::MakeLocalFECAndFES(): CRITICAL ERROR: Unknown finite element collection (fec_name) was provided. If 'fec_name' provides a valid FEC, it should be added as an option in this function. Provided fec_name: " << fec_name << endl;
        exit(1);
    }

    // Create the finite element space
    fes_local = std::make_unique<FiniteElementSpace>(local_mesh, fec_local.get(), fes_dim);
}

// Function for clearing and reassigning variable "gf" to "gf_local"
void GridFunctionManager::UseLocalGridFunction()
{
    gf.reset();
    gf = gf_local;
}


// Define functions for MPI-parallel
#ifdef MPI_BUILD
    // Function for transfering the gridfunction "gf" to a ParGridFunction
    void GridFunctionManager::MakeParGridFunction(ParMesh* mesh_, int* partitioning_)
    {
        par_gf = new ParGridFunction(mesh_, gf.get(), partitioning_);
        par_fes = par_gf->FESpace();
    }
#endif




#ifdef MPI_BUILD
    // ===============================================================
    //   Define the functions for the ParLinearTimeDependentOperator class.
    // ===============================================================
    // Define the function that prepares the operator, preconditioner, and solver for Explicitr Euler time stepping
    void ParLinearTimeDependentOperator::PrepareExplicitEuler()
    {
        // For reference:   M * du/dt + A * u = b
        // We solve for du/dt.

        // Solve parameters
        int maxIter(5000);
        real_t rtol(1.0e-7);
        real_t atol(1.0e-9);

        // Create the preconditioner (Note: 'dt' is in the preconditioner)
        M_inv = new HypreSmoother(M, 0); // 0 for DSmoother
        M_inv->iterative_mode = false;
        PC->SetDiagonalBlock(0, M_inv);

        // Build the block operator
        Op->SetBlock(0, 0, &M);
        
        // For Explicit Euler
        explicitSolver = CGSolver(comm);
        explicitSolver.iterative_mode = false;
        explicitSolver.SetRelTol(rtol);
        explicitSolver.SetAbsTol(atol);
        explicitSolver.SetMaxIter(maxIter);
        explicitSolver.SetPrintLevel(0);
        explicitSolver.SetOperator(M);
        explicitSolver.SetPreconditioner(*PC);
    }
    
    // Define the function that prepares the operator, preconditioner, and solver for Implicit Euler time stepping
    void ParLinearTimeDependentOperator::PrepareImplicitEuler(string solver_type, string PC_type)
    {
        // For reference:   (M + A*dt) * du/dt + A * u = b
        // We solve for du/dt.
        
        // Solve parameters
        int maxIter(8000);
        //real_t rtol(1.0e-7);
        //real_t atol(1.0e-9);
        real_t rtol(1.0e-4);
        real_t atol(1.0e-6);
        
        // Create F = M + A dt  (i.e., the inverted matrix)
        F = Add(1.0, M, dt, A);

        // Create the preconditioner (Note: 'dt' is in the preconditioner)
        if (PC_type == "DSmoother")
        {
            F_inv = new HypreSmoother(*F, 0); // 0 for DSmoother
            F_inv->iterative_mode = false;
            PC->SetDiagonalBlock(0, F_inv);
        }
        else if (PC_type == "BlockILU")
        {
            F_inv = new HypreILU();
            //F_inv->SetPrintLevel(0);
            F_inv->SetOperator(*F);
            PC->SetDiagonalBlock(0, F_inv);
        }
        else
        {
            cout << "mfem_util.cpp: ParLinearTimeDependentOperator::PrepareImplicitEuler(): CRITICAL ERROR: PC_type not recognized." << endl;
            exit(1);
        }


        // Build the block operator (Note: 'dt' is in the operator)
        Op->SetBlock(0, 0, F);
        
        // Prepare the solver
        if (solver_type == "CGSolver")
        {
            cgSolver = CGSolver(comm);
            cgSolver.iterative_mode = true;
            cgSolver.SetRelTol(rtol);
            cgSolver.SetAbsTol(atol);
            cgSolver.SetMaxIter(maxIter);
            cgSolver.SetPrintLevel(0);
            cgSolver.SetOperator(*Op);
            cgSolver.SetPreconditioner(*PC);
            implicitSolver = &cgSolver;
        }
        else if (solver_type == "GMRESSolver")
        {
            gmresSolver = GMRESSolver(comm);
            gmresSolver.iterative_mode = true;
            gmresSolver.SetRelTol(rtol);
            gmresSolver.SetAbsTol(atol);
            gmresSolver.SetMaxIter(maxIter);
            gmresSolver.SetPrintLevel(0);
            gmresSolver.SetOperator(*Op);
            gmresSolver.SetPreconditioner(*PC);
            implicitSolver = &gmresSolver;
        }
        else
        {
            cout << "mfem_util.cpp: ParLinearTimeDependentOperator::PrepareImplicitEuler(): CRITICAL ERROR: solver_type not recognized." << endl;
            exit(1);
        }
    }
    /*
    // Define the explicit Mult function for the operator
    void LinearTimeDependentOperator::Mult(const Vector &u, Vector &du_dt) const
    {
        // For reference:   M * du/dt + A * u = b
        // We solve for du/dt.

        A.Mult(u, du_dt); // Compute A u
        du_dt *= -1; // Move A u to RHS
        du_dt += b; // -A u + b
        explicitSolver.Mult(du_dt, du_dt); // Compute du/dt = (M^{-1}) * (-A u + b)
    }
    */
    // Define the implicit Mult function for the operator
    void ParLinearTimeDependentOperator::ImplicitSolve(const real_t dt, const Vector &u, Vector &du_dt)
    {
        // For reference:   (M + A*dt) * du/dt + A * u = b
        // We solve for du/dt.
        A.Mult(u, du_dt); // Compute  A u
        du_dt *= -1; // Move A u to RHS
        du_dt += b; // -A u + b
        implicitSolver->Mult(du_dt, du_dt); // Compute du/dt = ((M + A*dt)^{-1}) * (-A u + b)
    }
    
#endif
