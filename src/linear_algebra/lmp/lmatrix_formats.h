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

/// @author Dimitar Lukarski

#ifndef __LMATRIX_FORMATS_H
#define __LMATRIX_FORMATS_H

#include <vector>
#include <cassert>
#include <iostream>
#include <stdlib.h>
#include <typeinfo>
#include <cstring>

namespace hiflow {
namespace la {

/// @brief COO Matrix structure
/// @author Dimitar Lukarski

template < typename ValueType, typename IndexType = int > 
struct COO_lMatrixType {
  IndexType *row;       // row index
  IndexType *col;       // col index
  ValueType *val; // values
};

template < typename ValueType, typename IndexType = int > 
struct COO_lMatrixTypeConst {
  const IndexType *row;       // row index
  const IndexType *col;       // col index
  const ValueType *val; // values
};

/// @brief CSR Matrix structure
/// @author Dimitar Lukarski

template < typename ValueType, typename IndexType = int > 
struct CSR_lMatrixType {
  IndexType *row;       // row pointer
  IndexType *col;       // col index
  ValueType *val; // values
};

template < typename ValueType, typename IndexType = int > 
struct CSR_lMatrixTypeConst {
  const IndexType *row;       // row pointer
  const IndexType *col;       // col index
  const ValueType *val; // values
};

template < typename ValueType, typename IndexType = int > 
struct CSC_lMatrixType {
  IndexType *row;       // row pointer
  IndexType *col;       // col index
  ValueType *val; // values
};

template < typename ValueType, typename IndexType = int > 
struct CSC_lMatrixTypeConst {
  const IndexType *row;       // row pointer
  const IndexType *col;       // col index
  const ValueType *val; // values
};

/// @brief MCSR Matrix structure
/// @author Dimitar Lukarski

template < typename ValueType > struct MCSR_lMatrixType {
  int *row;        // row pointer
  int *col;        // col index
  ValueType *val;  // values
  ValueType *diag; // diagonal values
};

/// @brief DIA Matrix structure
/// @author Dimitar Lukarski

template < typename ValueType > struct DIA_lMatrixType {
  // extra matrix info
  int num_diagonals; // number of diagonals

  // data
  int *offsets;    // the offsets wrt the diagonal
  ValueType *data; // values
};

/// @brief ELL Matrix structure
/// @author Dimitar Lukarski

template < typename ValueType > struct ELL_lMatrixType {
  // extra matrix info
  int num_cols; // max number of column per row

  // data
  int *index;      //  column index set
  ValueType *data; // values
};

/// @brief DENSE Matrix structure
/// @author Niels Wegh

template < typename ValueType > struct DENSE_lMatrixType {
  ValueType *val; // values
};


/// @brief convert csr structure to csc structure
/// @author Philipp Gerstner

template < typename ValueType, typename IndexType >
void csr_2_csc (const CSR_lMatrixTypeConst<ValueType, IndexType>& input,
                IndexType nrows, IndexType ncols, IndexType nnz,
                bool allocate_output,
                CSC_lMatrixType<ValueType, IndexType>& output)
{
  assert (input.row != nullptr);
  assert (input.col != nullptr);
  
  if (allocate_output)
  {
    if (output.col != nullptr)
      delete[] output.col;
    if (output.row != nullptr)
      delete[] output.row;
    if (output.val != nullptr)
      delete[] output.val;
    
    output.col = new IndexType[ncols+1];
    output.row = new IndexType[nnz];
    output.val = new ValueType[nnz];
  }
  
  assert (output.row != nullptr);
  assert (output.col != nullptr);
  
  memset(output.col, 0, (ncols + 1) * sizeof(IndexType));
  memset(output.row, 0, nnz * sizeof(IndexType));
  
  if (output.val != nullptr)
  {
    memset(output.val, 0, nnz * sizeof(ValueType));
  }
  
  // fill column pointer
  std::vector<IndexType> nnz_per_col (ncols, 0);
  for (IndexType c=0; c!= nnz; ++c)
  {
    ++(nnz_per_col[input.col[c]]);
  }
  for (IndexType j=1; j!= ncols+1; ++j)
  {
    output.col[j] = output.col[j-1] + nnz_per_col[j];
  }
  
  // fill row ptr and nonzero values
  // loop through csr structure: 
  nnz_per_col.clear();
  nnz_per_col.resize(ncols, 0);
  if ( (output.val != nullptr) && (input.val != nullptr) )
  {
    assert (input.val != nullptr);

    for (IndexType i=0; i!= nrows; ++i)
    {
      // i = row
      const IndexType first_nz = input.row[i];
      const IndexType last_nz = input.row[i+1];
  
      for (IndexType c=first_nz; c!=last_nz; ++c)
      {
        // j = col, v = A(i,j)
        const IndexType j = input.col[c];
        const ValueType v = input.val[c];
        const IndexType pos = output.col[j]+nnz_per_col[j]; 
      
        output.row[pos] = i;
        output.val[pos] = v;
  
        ++(nnz_per_col[j]);
      }
    } 
  }
  else
  {
    for (IndexType i=0; i!= nrows; ++i)
    {
      // i = row
      const IndexType first_nz = input.row[i];
      const IndexType last_nz = input.row[i+1];
  
      for (IndexType c=first_nz; c!=last_nz; ++c)
      {
        // j = col, v = A(i,j)
        const IndexType j = input.col[c];
        const IndexType pos = output.col[j]+nnz_per_col[j]; 
      
        output.row[pos] = i; 
        ++(nnz_per_col[j]);
      }
    } 
  }
}

template < typename ValueType, typename IndexType >
void csrval_2_cscval (const CSR_lMatrixTypeConst<ValueType, IndexType>& input,
                      IndexType nrows, IndexType ncols, IndexType nnz,
                      CSC_lMatrixType<ValueType, IndexType>& output)
{
  assert (input.row != nullptr);
  assert (input.col != nullptr);
  assert (input.val != nullptr);

  assert (output.row != nullptr);
  assert (output.col != nullptr);
  
  if (output.val != nullptr)
  {
    memset(output.val, 0, nnz * sizeof(ValueType));
  }
  else
  {
    output.val = new ValueType[nnz];
  }
  
  // fill row ptr and nonzero values
  // loop through csr structure: 
  std::vector<IndexType> nnz_per_col(ncols, 0);

  for (IndexType i=0; i!= nrows; ++i)
  {
    // i = row
    const IndexType first_nz = input.row[i];
    const IndexType last_nz = input.row[i+1];
  
    for (IndexType c=first_nz; c!=last_nz; ++c)
    {
      // j = col, v = A(i,j)
      const IndexType j = input.col[c];
      const ValueType v = input.val[c];
      const IndexType pos = output.col[j]+nnz_per_col[j]; 

      output.val[pos] = v;
  
      ++(nnz_per_col[j]);
    }
  } 
}

template < typename ValueType, typename IndexType >
void csr_2_coo (const CSR_lMatrixTypeConst<ValueType, IndexType>& input,
                IndexType nrows, IndexType ncols, IndexType nnz,
                bool allocate_output,
                COO_lMatrixType<ValueType, IndexType>& output)
{
  assert (input.row != nullptr);
  assert (input.col != nullptr);

  if (allocate_output)
  {
    if (output.col != nullptr)
      delete[] output.col;
    if (output.row != nullptr)
      delete[] output.row;
    if (output.val != nullptr)
      delete[] output.val;
    
    output.col = new IndexType[nnz];
    output.row = new IndexType[nnz];
    output.val = new ValueType[nnz];
  }

  assert (output.row != nullptr);
  assert (output.col != nullptr);
  
  for (IndexType i = 0; i != nrows; ++i) 
  {
    // i = row
    const IndexType first_nz = input.row[i];
    const IndexType last_nz = input.row[i+1];
    
    for (IndexType c=first_nz; c!=last_nz; ++c)
    {
      output.row[c] = i;
    }
  }

  if ((output.val != nullptr) && (input.val != nullptr))
  {
    assert (input.val != nullptr);
    for (IndexType c = 0; c != nnz; ++c) 
    {
      output.col[c] = input.col[c];
      output.val[c] = input.val[c];
    }
  }
  else
  {
    for (IndexType c = 0; c != nnz; ++c) 
    {
      output.col[c] = input.col[c];
    }
  }
}

template < typename ValueType, typename IndexType >
void coo_2_csr (const COO_lMatrixTypeConst<ValueType, IndexType>& input,
                IndexType nrows, IndexType ncols, IndexType nnz,
                bool allocate_output, 
                CSR_lMatrixType<ValueType, IndexType>& output)
{
  assert (input.row != nullptr);
  assert (input.col != nullptr);

  if (allocate_output)
  {
    if (output.col != nullptr)
      delete[] output.col;
    if (output.row != nullptr)
      delete[] output.row;
    if (output.val != nullptr)
      delete[] output.val;
    
    output.row = new IndexType[nrows+1];
    output.col = new IndexType[nnz];
    output.val = new ValueType[nnz];
  }

  assert (output.row != nullptr);
  assert (output.col != nullptr);
  
  memset(output.row, 0, (nrows + 1) * sizeof(IndexType));
  memset(output.col, 0, nnz * sizeof(IndexType));
  
  if (output.val != nullptr)
  {
    memset(output.val, 0, nnz * sizeof(ValueType));
  }
  
  // build the row offsets
  for (IndexType i = 0; i != nnz; ++i)
  {
    ++(output.row[input.row[i]]);
  }
  output.row[nrows] = nnz;

  IndexType *index_row = new IndexType[nrows]; // auxiliary index for accessing the CSR cols, vals

  // accumulate the row offsets
  IndexType acum = 0;
  for (IndexType i = 0; i != nrows; ++i) 
  {
    IndexType current_row = output.row[i];
    output.row[i] = acum;
    index_row[i] = acum; // = output.row[i]
    acum = acum + current_row;
  }

  // for init structure only - no data
  if (output.val == nullptr || input.val == nullptr) 
  {
    for (IndexType i = 0; i != nnz; ++i) 
    {
      // inverse mapping
      output.col[index_row[input.row[i]]] = input.col[i];

      // keep the track of the current row index
      ++(index_row[input.row[i]]); // = index_row[ rows[i] ] + 1 ;
    }
  } 
  else 
  {
    for (IndexType i = 0; i != nnz; ++i) 
    {
      const IndexType temp = index_row[input.row[i]];
      // inverse mapping
      output.col[temp] = input.col[i];
      output.val[temp] = input.val[i];

      // keep the track of the current row index
      ++(index_row[input.row[i]]); // = index_row[ rows[i] ] + 1 ;
    }
  }

  delete[] index_row;

  // Sorting the Cols (per Row)
  ValueType tv;
  IndexType ti;

  if (output.val != nullptr)
  {
    for (IndexType k = 0, k_e = nrows; k != k_e; ++k) 
    {
      for (IndexType i = output.row[k], i_e = output.row[k + 1]; i != i_e; ++i) 
      {
        IndexType current_min_index = i;
        IndexType current_min = output.col[i];
        for (IndexType j = i + 1; j < i_e; ++j) 
        {
          if (output.col[j] < current_min) 
          {
            current_min_index = j;
            current_min = output.col[j];
          }
        }
        ti = output.col[i];
        output.col[i] = output.col[current_min_index];
        output.col[current_min_index] = ti;

        tv = output.val[i];
        output.val[i] = output.val[current_min_index];
        output.val[current_min_index] = tv;
      }
    }
  }
  else
  {
    for (IndexType k = 0, k_e = nrows; k != k_e; ++k) 
    {
      for (IndexType i = output.row[k], i_e = output.row[k + 1]; i != i_e; ++i) 
      {
        IndexType current_min_index = i;
        IndexType current_min = output.col[i];
        for (IndexType j = i + 1; j < i_e; ++j) 
        {
          if (output.col[j] < current_min) 
          {
            current_min_index = j;
            current_min = output.col[j];
          }
        }
        ti = output.col[i];
        output.col[i] = output.col[current_min_index];
        output.col[current_min_index] = ti;
      }
    }
  }
}

// check if column indices are in scending order
template < typename ValueType, typename IndexType >
bool csr_check_ascending_order (const CSR_lMatrixTypeConst<ValueType, IndexType>& input,
                                IndexType nrows, IndexType nnz)
{
  assert (input.row != nullptr);
  assert (input.col != nullptr);
  
  for (IndexType i=0; i!= nrows; ++i)
  {
    // i = row
    const IndexType first_nz = input.row[i];
    const IndexType last_nz = input.row[i+1];
    IndexType prev_j = -1;
    
    for (IndexType c=first_nz; c!=last_nz; ++c)
    {
      // j = col, v = A(i,j)
      const IndexType j = input.col[c];
      
      if (j <= prev_j)
      {
        return false;
      }
      prev_j = j;
    }
  }
  return true; 
}

} // namespace la
} // namespace hiflow

#endif
