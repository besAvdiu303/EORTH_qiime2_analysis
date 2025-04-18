// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup op_sum
//! @{


class op_sum
  : public traits_op_xvec
  {
  public:
  
  // dense matrices
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op< T1,                 op_sum >& in);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op< eOp<T1,eop_square>, op_sum >& in);

  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op< eOp<T1,eop_pow   >, op_sum >& in);

  template<typename eT>
  inline static void apply_mat_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim);
  
  template<typename eT>
  inline static void apply_mat_square_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim);

  template<typename T1>
  inline static void apply_proxy_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword dim);
  
  
  // cubes
  
  template<typename T1>
  inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1, op_sum>& in);
  
  template<typename eT>
  inline static void apply_cube_noalias(Cube<eT>& out, const Cube<eT>& X, const uword dim);
  
  template<typename T1>
  inline static void apply_proxy_noalias(Cube<typename T1::elem_type>& out, const ProxyCube<T1>& P, const uword dim);
  };


//! @}
