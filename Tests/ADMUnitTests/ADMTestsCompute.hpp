/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMTESTSCOMPUTE_HPP_
#define ADMTESTSCOMPUTE_HPP_

#include "ADMFixedBGVars.hpp"
#include "KerrSchildFixedBG.hpp"
#include "VarsTools.hpp"
#include "UserVariables.hpp"


//template <class data_t> using Vars = CCZ4Vars::VarsWithGauge<data_t>;
template <class deriv_t> class ADMTestsCompute
{
  protected:
    const deriv_t m_deriv;
    const double m_dx; 
    //KerrSchildFixedBG::params_t m_bg_params;
    

  public:
    ADMTestsCompute(double dx) : m_deriv(dx), m_dx(dx) {}
    
    //template <class data_t> struct Vars// using Vars = ADMFixedBGVars::Vars<data_t>;
    
    //template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;
    template <class data_t> struct Vars
    {
      
    // ADM vars needed in matter only rhs (ok for Proca and SF)
    
    Tensor<2, data_t> gamma;
    Tensor<2, data_t> K_tensor;
    data_t K;
    data_t lapse;
    Tensor<1, data_t> shift;
    Tensor<2, Tensor<1, data_t>> d1_gamma;
    Tensor<1, data_t> d1_lapse;
    Tensor<2, data_t> d1_shift;

    //Extra variables needed for tensor field case
    Tensor<2, data_t> gamma_UU;
    Tensor<1, data_t> d1_K;    
    Tensor<3, data_t> d1_K_tensor; 
    Tensor<2, data_t> d2_lapse;
    Tensor<3, data_t> d2_shift;
    Tensor<4, data_t> d2_gamma;
    Tensor<3, data_t> d1_gamma_UU;
    Tensor<4, data_t> d1_chris_phys; 
    Tensor<4, data_t> riemann_phys_ULLL;
    Tensor<2, data_t> ricci_phys;

    Tensor<1, data_t> drdx;   
    Tensor<2, data_t> d2rdx2;
    
    
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                        // VarsTools
            //define_enum_mapping(mapping_function, GRInterval<c_gamma_UU11,c_gamma_UU33>(), gamma_UU);
            //define_enum_mapping(mapping_function, c_gamma_UU11, gamma_UU[0][0]);
            //define_enum_mapping(mapping_function, c_gamma_UU11, gamma_UU[0][0]);
            //define_enum_mapping(mapping_function, c_d1_gamma_UU111, gamma_UU[0][0]);
            define_symmetric_enum_mapping(mapping_function, GRInterval<c_d1_gamma_UU11, c_d1_gamma_UU33>(), gamma_UU);
            define_symmetric_enum_mapping(mapping_function, GRInterval<c_d2_gamma11, c_d2_gamma33>(), gamma);
            define_symmetric_enum_mapping(mapping_function, GRInterval<c_d1_K_tensor11, c_d1_K_tensor33>(), K_tensor);
            

            define_enum_mapping(mapping_function, c_d2_lapse, lapse);

            /*
            define_enum_mapping(mapping_function, c_d1_shift1, shift[0]);
            define_enum_mapping(mapping_function, c_d1_shift2, shift[1]);
            define_enum_mapping(mapping_function, c_d1_shift3, shift[2]);
            */

            define_enum_mapping(mapping_function, GRInterval<c_d2_shift1, c_d2_shift3>(), shift);

        }
        
    };
    

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
   
        const auto out_d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto out_d2 = m_deriv.template diff2<Vars>(current_cell);
   
        current_cell.store_vars(out_d1.gamma_UU[0][0][0], c_d1_gamma_UU11);
        current_cell.store_vars(out_d1.gamma_UU[0][1][0], c_d1_gamma_UU12);
        current_cell.store_vars(out_d1.gamma_UU[0][2][0], c_d1_gamma_UU13);
        current_cell.store_vars(out_d1.gamma_UU[1][1][0], c_d1_gamma_UU22);
        current_cell.store_vars(out_d1.gamma_UU[1][2][0], c_d1_gamma_UU23);
        current_cell.store_vars(out_d1.gamma_UU[2][2][0], c_d1_gamma_UU33);

        current_cell.store_vars(out_d2.gamma[0][0][0][1], c_d2_gamma11);
        current_cell.store_vars(out_d2.gamma[0][1][0][1], c_d2_gamma12);
        current_cell.store_vars(out_d2.gamma[0][2][0][1], c_d2_gamma13);
        current_cell.store_vars(out_d2.gamma[1][1][0][1], c_d2_gamma22);
        current_cell.store_vars(out_d2.gamma[1][2][0][1], c_d2_gamma23);
        current_cell.store_vars(out_d2.gamma[2][2][0][1], c_d2_gamma33);

        current_cell.store_vars(out_d1.K_tensor[0][0][2], c_d1_K_tensor11);
        current_cell.store_vars(out_d1.K_tensor[0][1][2], c_d1_K_tensor12);
        current_cell.store_vars(out_d1.K_tensor[0][2][2], c_d1_K_tensor13);
        current_cell.store_vars(out_d1.K_tensor[1][1][2], c_d1_K_tensor22);
        current_cell.store_vars(out_d1.K_tensor[1][2][2], c_d1_K_tensor23);
        current_cell.store_vars(out_d1.K_tensor[2][2][2], c_d1_K_tensor33);

        current_cell.store_vars(out_d2.lapse[0][0], c_d2_lapse);

        current_cell.store_vars(out_d1.shift[0][0], c_d1_shift1);
        current_cell.store_vars(out_d1.shift[1][0], c_d1_shift2);
        current_cell.store_vars(out_d1.shift[2][0], c_d1_shift3);        

        current_cell.store_vars(out_d2.shift[0][0][0], c_d2_shift1);
        current_cell.store_vars(out_d2.shift[1][0][0], c_d2_shift2);
        current_cell.store_vars(out_d2.shift[2][0][0], c_d2_shift3);
 
        
 
    }
};

#endif /* ADMTESTSCOMPUTE_HPP_ */
