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
void AveragingOperator::GetARIndices(Mesh *mesh_)
{
    // Get the indices of the elements inside each averaging region
    for (int i_Elem = 0; i_Elem < mesh_->GetNE(); i_Elem++)
    {
        for (int i_AR = 0; i_AR < N_AR; i_AR++)
        {
            if (mesh_->GetAttribute(i_Elem) == i_AR + 1)
            {
                AR_Elem_inds[i_AR].push_back(i_Elem);
                break;
            }
        }
    }

    // Get the indices of the DOFS in the dummy avg-space of each element inside each averaging region
    Array<int> Elem_DOFs;
    for (int i_AR = 0; i_AR < N_AR; i_AR++)
    {
        for (int i_Elem = 0; i_Elem < AR_Elem_inds[i_AR].size(); i_Elem++)
        {
            //fespace_avg->GetElementDofs(AR_Elem_inds[i_AR][i_Elem], Elem_DOFs);
            fespace_avg->GetElementVDofs(AR_Elem_inds[i_AR][i_Elem], Elem_DOFs); // GetElementVDofs apparently outputs dofs like [dof_x^1, dof_x^2, dof_y^1, dof_y^2, ...]
            assert (Elem_DOFs.Size() == DOFs_per_vdim * vdim);
            for (int i_vd = 0; i_vd < vdim; i_vd++)
            {
                AR_Dof_inds[i_vd][i_AR].insert(AR_Dof_inds[i_vd][i_AR].end(), Elem_DOFs.begin() + i_vd*DOFs_per_vdim, Elem_DOFs.begin() + (i_vd + 1)*DOFs_per_vdim);
            }
            //AR_Dof_inds[i_AR].insert(AR_Dof_inds[i_AR].end(), Elem_DOFs.begin(), Elem_DOFs.end());
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
    
    avg_mat = avg_mat_;
    AR_areas = AR_areas_;
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




/*
// ===============================================================
//   Define the functions for the TransportOperator class.
// ===============================================================
// Define the class constructor
TransportOperator::TransportOperator(BlockOperator &M_, BlockOperator &K_, BlockVector &b_, double dt_, BlockDiagonalPreconditioner &M_PC_, BlockDiagonalPreconditioner &T_PC_)
                 : TimeDependentOperator(K_.Width(), K_.Height(), 0.0), M(M_), K(K_), b(b_), dt(dt_), M_PC(M_PC_), T_PC(T_PC_)
{
    int maxIter(5000);
    real_t rtol(1.0e-7);
    real_t atol(0.0);

    // For Explicit Euler
    M_inv.iterative_mode = false;
    M_inv.SetRelTol(rtol);
    M_inv.SetAbsTol(atol);
    M_inv.SetMaxIter(maxIter);
    M_inv.SetPrintLevel(0);
    M_inv.SetOperator(M);
    if (&M_PC) { M_inv.SetPreconditioner(M_PC); }
    
    // For Implicit Euler
    T_inv.iterative_mode = true;
    T_inv.SetRelTol(rtol);
    T_inv.SetAbsTol(atol);
    T_inv.SetMaxIter(maxIter);
    T_inv.SetPrintLevel(0);
    M_Kdt = AddBlockOperators(1.0, M, dt, K);
    T_inv.SetOperator(*M_Kdt);
    if (&T_PC) { T_inv.SetPreconditioner(T_PC); }
};

// Define the Mult function for the operator
void TransportOperator::Mult(const Vector &u, Vector &du_dt) const
{
    // For reference:   M * du/dt + K * u = b
    // We solve for du/dt.

    K.Mult(u, du_dt); // Compute K u
    du_dt *= -1; // Move K u to RHS
    du_dt += b; // -K u + b
    M_inv.Mult(du_dt, du_dt); // Compute du/dt = (M^{-1}) * (-K u + b)
}

// Define the implicit Mult function for the operator
void TransportOperator::ImplicitSolve(const real_t dt, const Vector &u, Vector &du_dt)
{
    // For reference:   M * du/dt + K * (u + dt * du/dt) = b
    // We solve for du/dt.

    Vector RHS = u;
    K.Mult(u, RHS); // Compute  K u
    RHS *= -1; // Move K u to RHS
    RHS += b; // RHS = -K u + b
    T_inv.Mult(RHS, du_dt); // Compute du/dt = ((M + K*dt)^{-1}) * (-K u + b)
}

// Primarily written by ChatGPT.
unique_ptr<BlockOperator> TransportOperator::AddBlockOperators(const double A_coeff, const BlockOperator &A, const double B_coeff, const BlockOperator &B)
{
    MFEM_VERIFY(A.NumRowBlocks() == B.NumRowBlocks(), "Row block size mismatch");
    MFEM_VERIFY(A.NumColBlocks() == B.NumColBlocks(), "Col block size mismatch");


    auto row_offsets = A.RowOffsets();
    auto col_offsets = A.ColOffsets(); 
    
    unique_ptr<BlockOperator> C(new BlockOperator(row_offsets, col_offsets));
    

    for (int i = 0; i < A.NumRowBlocks(); ++i)
    {
        for (int j = 0; j < A.NumColBlocks(); ++j)
        {
            const Operator &Ai = A.GetBlock(i, j);
            const Operator &Bi = B.GetBlock(i, j);
            
            // If both blocks exist
            if (&Ai && &Bi)
            {
                const auto &As = dynamic_cast<const SparseMatrix *>(&Ai);
                const auto &Bs = dynamic_cast<const SparseMatrix *>(&Bi);
                MFEM_VERIFY(&As && &Bs, "Non-sparse matrix blocks not supported");
                MFEM_VERIFY(&As != nullptr && &Bs != nullptr, "One of the added matrices is null");
                
                SparseMatrix *Cij = Add(A_coeff, *As, B_coeff, *Bs);
                C->SetBlock(i, j, Cij);
            }
            else if (&Ai)
            {
                const auto &As = dynamic_cast<const SparseMatrix *>(&Ai);
                MFEM_VERIFY(&As, "Non-sparse block in A");
                
                SparseMatrix *Cij = new SparseMatrix(*As);
                *Cij *= A_coeff;
                C->SetBlock(i, j, Cij); // copy
            }
            else if (&Bi)
            {
                const auto &Bs = dynamic_cast<const SparseMatrix *>(&Bi);
                MFEM_VERIFY(&Bs, "Non-sparse block in B");
                
                SparseMatrix *Cij = new SparseMatrix(*Bs);
                *Cij *= B_coeff; 
                C->SetBlock(i, j, Cij); // copy
            }
            else
            {
                C->SetBlock(i, j, nullptr); // no block
            }
        }
    }
 
    return C;
}
*/




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


