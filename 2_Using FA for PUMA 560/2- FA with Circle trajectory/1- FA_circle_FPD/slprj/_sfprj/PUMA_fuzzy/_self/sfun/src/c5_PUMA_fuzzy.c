/* Include files */

#include <stddef.h>
#include "blas.h"
#include "PUMA_fuzzy_sfun.h"
#include "c5_PUMA_fuzzy.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "PUMA_fuzzy_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c5_debug_family_names[52] = { "e", "pe", "vLR", "vSR", "vZ",
  "vSL", "vLL", "LAN", "SMN", "Ze", "SMP", "LAP", "FAL", "SLL", "Zpe", "SLR",
  "FAR", "muLAN", "muSMN", "muZe", "muSMP", "muLAP", "muFAL", "muSLL", "muZpe",
  "muSLR", "muFAR", "outLL1", "outSL1", "outZ1", "outSR1", "outLR1", "outLL",
  "outSL", "outZ", "outSR", "outLR", "woutLL", "woutSL", "woutZ", "woutSR",
  "woutLR", "i", "sumout", "sumwout", "nargin", "nargout", "teta", "desiredteta",
  "dteta", "desireddteta", "T" };

/* Function Declarations */
static void initialize_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance);
static void initialize_params_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance);
static void enable_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct *chartInstance);
static void disable_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct *chartInstance);
static void c5_update_debugger_state_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct *
  chartInstance);
static void set_sim_state_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance, const mxArray *c5_st);
static void finalize_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct *chartInstance);
static void sf_gateway_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance);
static void c5_chartstep_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance);
static void initSimStructsc5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c5_machineNumber, uint32_T
  c5_chartNumber, uint32_T c5_instanceNumber);
static const mxArray *c5_sf_marshallOut(void *chartInstanceVoid, void *c5_inData);
static real_T c5_emlrt_marshallIn(SFc5_PUMA_fuzzyInstanceStruct *chartInstance,
  const mxArray *c5_T, const char_T *c5_identifier);
static real_T c5_b_emlrt_marshallIn(SFc5_PUMA_fuzzyInstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static const mxArray *c5_b_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData);
static void c5_c_emlrt_marshallIn(SFc5_PUMA_fuzzyInstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, real_T c5_y[7]);
static void c5_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static void c5_info_helper(const mxArray **c5_info);
static const mxArray *c5_emlrt_marshallOut(const char * c5_u);
static const mxArray *c5_b_emlrt_marshallOut(const uint32_T c5_u);
static void c5_eml_scalar_eg(SFc5_PUMA_fuzzyInstanceStruct *chartInstance);
static const mxArray *c5_c_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData);
static int32_T c5_d_emlrt_marshallIn(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static uint8_T c5_e_emlrt_marshallIn(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance, const mxArray *c5_b_is_active_c5_PUMA_fuzzy, const char_T
  *c5_identifier);
static uint8_T c5_f_emlrt_marshallIn(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void init_dsm_address_info(SFc5_PUMA_fuzzyInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance)
{
  chartInstance->c5_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c5_is_active_c5_PUMA_fuzzy = 0U;
}

static void initialize_params_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c5_update_debugger_state_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct *
  chartInstance)
{
  const mxArray *c5_st;
  const mxArray *c5_y = NULL;
  real_T c5_hoistedGlobal;
  real_T c5_u;
  const mxArray *c5_b_y = NULL;
  uint8_T c5_b_hoistedGlobal;
  uint8_T c5_b_u;
  const mxArray *c5_c_y = NULL;
  real_T *c5_T;
  c5_T = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c5_st = NULL;
  c5_st = NULL;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_createcellmatrix(2, 1), false);
  c5_hoistedGlobal = *c5_T;
  c5_u = c5_hoistedGlobal;
  c5_b_y = NULL;
  sf_mex_assign(&c5_b_y, sf_mex_create("y", &c5_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c5_y, 0, c5_b_y);
  c5_b_hoistedGlobal = chartInstance->c5_is_active_c5_PUMA_fuzzy;
  c5_b_u = c5_b_hoistedGlobal;
  c5_c_y = NULL;
  sf_mex_assign(&c5_c_y, sf_mex_create("y", &c5_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c5_y, 1, c5_c_y);
  sf_mex_assign(&c5_st, c5_y, false);
  return c5_st;
}

static void set_sim_state_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance, const mxArray *c5_st)
{
  const mxArray *c5_u;
  real_T *c5_T;
  c5_T = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c5_doneDoubleBufferReInit = true;
  c5_u = sf_mex_dup(c5_st);
  *c5_T = c5_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 0)),
    "T");
  chartInstance->c5_is_active_c5_PUMA_fuzzy = c5_e_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 1)),
     "is_active_c5_PUMA_fuzzy");
  sf_mex_destroy(&c5_u);
  c5_update_debugger_state_c5_PUMA_fuzzy(chartInstance);
  sf_mex_destroy(&c5_st);
}

static void finalize_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance)
{
  real_T *c5_T;
  real_T *c5_teta;
  real_T *c5_desiredteta;
  real_T *c5_dteta;
  real_T *c5_desireddteta;
  c5_desireddteta = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c5_dteta = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c5_desiredteta = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c5_teta = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c5_T = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 4U, chartInstance->c5_sfEvent);
  chartInstance->c5_sfEvent = CALL_EVENT;
  c5_chartstep_c5_PUMA_fuzzy(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_PUMA_fuzzyMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  _SFD_DATA_RANGE_CHECK(*c5_T, 0U);
  _SFD_DATA_RANGE_CHECK(*c5_teta, 1U);
  _SFD_DATA_RANGE_CHECK(*c5_desiredteta, 2U);
  _SFD_DATA_RANGE_CHECK(*c5_dteta, 3U);
  _SFD_DATA_RANGE_CHECK(*c5_desireddteta, 4U);
}

static void c5_chartstep_c5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance)
{
  real_T c5_hoistedGlobal;
  real_T c5_b_hoistedGlobal;
  real_T c5_c_hoistedGlobal;
  real_T c5_d_hoistedGlobal;
  real_T c5_teta;
  real_T c5_desiredteta;
  real_T c5_dteta;
  real_T c5_desireddteta;
  uint32_T c5_debug_family_var_map[52];
  real_T c5_e;
  real_T c5_pe;
  real_T c5_vLR;
  real_T c5_vSR;
  real_T c5_vZ;
  real_T c5_vSL;
  real_T c5_vLL;
  real_T c5_LAN;
  real_T c5_SMN;
  real_T c5_Ze;
  real_T c5_SMP;
  real_T c5_LAP;
  real_T c5_FAL;
  real_T c5_SLL;
  real_T c5_Zpe;
  real_T c5_SLR;
  real_T c5_FAR;
  real_T c5_muLAN;
  real_T c5_muSMN;
  real_T c5_muZe;
  real_T c5_muSMP;
  real_T c5_muLAP;
  real_T c5_muFAL;
  real_T c5_muSLL;
  real_T c5_muZpe;
  real_T c5_muSLR;
  real_T c5_muFAR;
  real_T c5_outLL1[7];
  real_T c5_outSL1[7];
  real_T c5_outZ1[7];
  real_T c5_outSR1[7];
  real_T c5_outLR1[7];
  real_T c5_outLL;
  real_T c5_outSL;
  real_T c5_outZ;
  real_T c5_outSR;
  real_T c5_outLR;
  real_T c5_woutLL;
  real_T c5_woutSL;
  real_T c5_woutZ;
  real_T c5_woutSR;
  real_T c5_woutLR;
  real_T c5_i;
  real_T c5_sumout;
  real_T c5_sumwout;
  real_T c5_nargin = 4.0;
  real_T c5_nargout = 1.0;
  real_T c5_T;
  int32_T c5_i0;
  int32_T c5_i1;
  int32_T c5_i2;
  int32_T c5_i3;
  int32_T c5_i4;
  real_T c5_A;
  real_T c5_x;
  real_T c5_b_x;
  real_T c5_c_x;
  real_T c5_y;
  real_T c5_varargin_2;
  real_T c5_varargin_3;
  real_T c5_b_y;
  real_T c5_c_y;
  real_T c5_yk;
  real_T c5_d_y;
  real_T c5_b_A;
  real_T c5_d_x;
  real_T c5_e_x;
  real_T c5_f_x;
  real_T c5_e_y;
  real_T c5_c_A;
  real_T c5_g_x;
  real_T c5_h_x;
  real_T c5_i_x;
  real_T c5_f_y;
  real_T c5_varargin_1;
  real_T c5_b_varargin_2;
  real_T c5_c_varargin_2;
  real_T c5_b_varargin_3;
  real_T c5_j_x;
  real_T c5_g_y;
  real_T c5_k_x;
  real_T c5_h_y;
  real_T c5_xk;
  real_T c5_b_yk;
  real_T c5_l_x;
  real_T c5_i_y;
  real_T c5_minval;
  real_T c5_d_varargin_2;
  real_T c5_c_varargin_3;
  real_T c5_j_y;
  real_T c5_k_y;
  real_T c5_c_yk;
  real_T c5_l_y;
  real_T c5_d_A;
  real_T c5_m_x;
  real_T c5_n_x;
  real_T c5_o_x;
  real_T c5_m_y;
  real_T c5_e_A;
  real_T c5_p_x;
  real_T c5_q_x;
  real_T c5_r_x;
  real_T c5_n_y;
  real_T c5_b_varargin_1;
  real_T c5_e_varargin_2;
  real_T c5_f_varargin_2;
  real_T c5_d_varargin_3;
  real_T c5_s_x;
  real_T c5_o_y;
  real_T c5_t_x;
  real_T c5_p_y;
  real_T c5_b_xk;
  real_T c5_d_yk;
  real_T c5_u_x;
  real_T c5_q_y;
  real_T c5_b_minval;
  real_T c5_g_varargin_2;
  real_T c5_e_varargin_3;
  real_T c5_r_y;
  real_T c5_s_y;
  real_T c5_e_yk;
  real_T c5_t_y;
  real_T c5_f_A;
  real_T c5_v_x;
  real_T c5_w_x;
  real_T c5_x_x;
  real_T c5_u_y;
  real_T c5_g_A;
  real_T c5_y_x;
  real_T c5_ab_x;
  real_T c5_bb_x;
  real_T c5_v_y;
  real_T c5_c_varargin_1;
  real_T c5_h_varargin_2;
  real_T c5_i_varargin_2;
  real_T c5_f_varargin_3;
  real_T c5_cb_x;
  real_T c5_w_y;
  real_T c5_db_x;
  real_T c5_x_y;
  real_T c5_c_xk;
  real_T c5_f_yk;
  real_T c5_eb_x;
  real_T c5_y_y;
  real_T c5_c_minval;
  real_T c5_j_varargin_2;
  real_T c5_g_varargin_3;
  real_T c5_ab_y;
  real_T c5_bb_y;
  real_T c5_g_yk;
  real_T c5_cb_y;
  real_T c5_h_A;
  real_T c5_fb_x;
  real_T c5_gb_x;
  real_T c5_hb_x;
  real_T c5_db_y;
  real_T c5_k_varargin_2;
  real_T c5_h_varargin_3;
  real_T c5_eb_y;
  real_T c5_fb_y;
  real_T c5_h_yk;
  real_T c5_gb_y;
  real_T c5_i_A;
  real_T c5_ib_x;
  real_T c5_jb_x;
  real_T c5_kb_x;
  real_T c5_hb_y;
  real_T c5_l_varargin_2;
  real_T c5_i_varargin_3;
  real_T c5_ib_y;
  real_T c5_jb_y;
  real_T c5_i_yk;
  real_T c5_kb_y;
  real_T c5_j_A;
  real_T c5_lb_x;
  real_T c5_mb_x;
  real_T c5_nb_x;
  real_T c5_lb_y;
  real_T c5_k_A;
  real_T c5_ob_x;
  real_T c5_pb_x;
  real_T c5_qb_x;
  real_T c5_mb_y;
  real_T c5_d_varargin_1;
  real_T c5_m_varargin_2;
  real_T c5_n_varargin_2;
  real_T c5_j_varargin_3;
  real_T c5_rb_x;
  real_T c5_nb_y;
  real_T c5_sb_x;
  real_T c5_ob_y;
  real_T c5_d_xk;
  real_T c5_j_yk;
  real_T c5_tb_x;
  real_T c5_pb_y;
  real_T c5_d_minval;
  real_T c5_o_varargin_2;
  real_T c5_k_varargin_3;
  real_T c5_qb_y;
  real_T c5_rb_y;
  real_T c5_k_yk;
  real_T c5_sb_y;
  real_T c5_l_A;
  real_T c5_ub_x;
  real_T c5_vb_x;
  real_T c5_wb_x;
  real_T c5_tb_y;
  real_T c5_m_A;
  real_T c5_xb_x;
  real_T c5_yb_x;
  real_T c5_ac_x;
  real_T c5_ub_y;
  real_T c5_e_varargin_1;
  real_T c5_p_varargin_2;
  real_T c5_q_varargin_2;
  real_T c5_l_varargin_3;
  real_T c5_bc_x;
  real_T c5_vb_y;
  real_T c5_cc_x;
  real_T c5_wb_y;
  real_T c5_e_xk;
  real_T c5_l_yk;
  real_T c5_dc_x;
  real_T c5_xb_y;
  real_T c5_e_minval;
  real_T c5_r_varargin_2;
  real_T c5_m_varargin_3;
  real_T c5_yb_y;
  real_T c5_ac_y;
  real_T c5_m_yk;
  real_T c5_bc_y;
  real_T c5_n_A;
  real_T c5_ec_x;
  real_T c5_fc_x;
  real_T c5_gc_x;
  real_T c5_cc_y;
  real_T c5_o_A;
  real_T c5_hc_x;
  real_T c5_ic_x;
  real_T c5_jc_x;
  real_T c5_dc_y;
  real_T c5_f_varargin_1;
  real_T c5_s_varargin_2;
  real_T c5_t_varargin_2;
  real_T c5_n_varargin_3;
  real_T c5_kc_x;
  real_T c5_ec_y;
  real_T c5_lc_x;
  real_T c5_fc_y;
  real_T c5_f_xk;
  real_T c5_n_yk;
  real_T c5_mc_x;
  real_T c5_gc_y;
  real_T c5_f_minval;
  real_T c5_u_varargin_2;
  real_T c5_o_varargin_3;
  real_T c5_hc_y;
  real_T c5_ic_y;
  real_T c5_o_yk;
  real_T c5_jc_y;
  real_T c5_p_A;
  real_T c5_nc_x;
  real_T c5_oc_x;
  real_T c5_pc_x;
  real_T c5_kc_y;
  real_T c5_v_varargin_2;
  real_T c5_p_varargin_3;
  real_T c5_lc_y;
  real_T c5_mc_y;
  real_T c5_p_yk;
  real_T c5_nc_y;
  real_T c5_g_varargin_1;
  real_T c5_w_varargin_2;
  real_T c5_x_varargin_2;
  real_T c5_q_varargin_3;
  real_T c5_qc_x;
  real_T c5_oc_y;
  real_T c5_rc_x;
  real_T c5_pc_y;
  real_T c5_g_xk;
  real_T c5_q_yk;
  real_T c5_sc_x;
  real_T c5_qc_y;
  real_T c5_g_minval;
  real_T c5_h_varargin_1;
  real_T c5_y_varargin_2;
  real_T c5_ab_varargin_2;
  real_T c5_r_varargin_3;
  real_T c5_tc_x;
  real_T c5_rc_y;
  real_T c5_uc_x;
  real_T c5_sc_y;
  real_T c5_h_xk;
  real_T c5_r_yk;
  real_T c5_vc_x;
  real_T c5_tc_y;
  real_T c5_h_minval;
  real_T c5_i_varargin_1;
  real_T c5_bb_varargin_2;
  real_T c5_cb_varargin_2;
  real_T c5_s_varargin_3;
  real_T c5_wc_x;
  real_T c5_uc_y;
  real_T c5_xc_x;
  real_T c5_vc_y;
  real_T c5_i_xk;
  real_T c5_s_yk;
  real_T c5_yc_x;
  real_T c5_wc_y;
  real_T c5_i_minval;
  real_T c5_j_varargin_1;
  real_T c5_db_varargin_2;
  real_T c5_eb_varargin_2;
  real_T c5_t_varargin_3;
  real_T c5_ad_x;
  real_T c5_xc_y;
  real_T c5_bd_x;
  real_T c5_yc_y;
  real_T c5_j_xk;
  real_T c5_t_yk;
  real_T c5_cd_x;
  real_T c5_ad_y;
  real_T c5_j_minval;
  int32_T c5_b_i;
  real_T c5_k_varargin_1;
  real_T c5_fb_varargin_2;
  real_T c5_gb_varargin_2;
  real_T c5_u_varargin_3;
  real_T c5_dd_x;
  real_T c5_bd_y;
  real_T c5_ed_x;
  real_T c5_cd_y;
  real_T c5_k_xk;
  real_T c5_u_yk;
  real_T c5_fd_x;
  real_T c5_dd_y;
  real_T c5_k_minval;
  real_T c5_l_varargin_1;
  real_T c5_hb_varargin_2;
  real_T c5_ib_varargin_2;
  real_T c5_v_varargin_3;
  real_T c5_gd_x;
  real_T c5_ed_y;
  real_T c5_hd_x;
  real_T c5_fd_y;
  real_T c5_l_xk;
  real_T c5_v_yk;
  real_T c5_id_x;
  real_T c5_gd_y;
  real_T c5_l_minval;
  real_T c5_m_varargin_1;
  real_T c5_jb_varargin_2;
  real_T c5_kb_varargin_2;
  real_T c5_w_varargin_3;
  real_T c5_jd_x;
  real_T c5_hd_y;
  real_T c5_kd_x;
  real_T c5_id_y;
  real_T c5_m_xk;
  real_T c5_w_yk;
  real_T c5_ld_x;
  real_T c5_jd_y;
  real_T c5_m_minval;
  real_T c5_n_varargin_1;
  real_T c5_lb_varargin_2;
  real_T c5_mb_varargin_2;
  real_T c5_x_varargin_3;
  real_T c5_md_x;
  real_T c5_kd_y;
  real_T c5_nd_x;
  real_T c5_ld_y;
  real_T c5_n_xk;
  real_T c5_x_yk;
  real_T c5_od_x;
  real_T c5_md_y;
  real_T c5_n_minval;
  real_T c5_o_varargin_1;
  real_T c5_nb_varargin_2;
  real_T c5_ob_varargin_2;
  real_T c5_y_varargin_3;
  real_T c5_pd_x;
  real_T c5_nd_y;
  real_T c5_qd_x;
  real_T c5_od_y;
  real_T c5_o_xk;
  real_T c5_y_yk;
  real_T c5_rd_x;
  real_T c5_pd_y;
  real_T c5_o_minval;
  real_T c5_p_varargin_1;
  real_T c5_pb_varargin_2;
  real_T c5_qb_varargin_2;
  real_T c5_ab_varargin_3;
  real_T c5_sd_x;
  real_T c5_qd_y;
  real_T c5_td_x;
  real_T c5_rd_y;
  real_T c5_p_xk;
  real_T c5_ab_yk;
  real_T c5_ud_x;
  real_T c5_sd_y;
  real_T c5_p_minval;
  int32_T c5_c_i;
  real_T c5_q_varargin_1;
  real_T c5_rb_varargin_2;
  real_T c5_sb_varargin_2;
  real_T c5_bb_varargin_3;
  real_T c5_vd_x;
  real_T c5_td_y;
  real_T c5_wd_x;
  real_T c5_ud_y;
  real_T c5_q_xk;
  real_T c5_bb_yk;
  real_T c5_xd_x;
  real_T c5_vd_y;
  real_T c5_q_minval;
  real_T c5_r_varargin_1;
  real_T c5_tb_varargin_2;
  real_T c5_ub_varargin_2;
  real_T c5_cb_varargin_3;
  real_T c5_yd_x;
  real_T c5_wd_y;
  real_T c5_ae_x;
  real_T c5_xd_y;
  real_T c5_r_xk;
  real_T c5_cb_yk;
  real_T c5_be_x;
  real_T c5_yd_y;
  real_T c5_r_minval;
  real_T c5_s_varargin_1;
  real_T c5_vb_varargin_2;
  real_T c5_wb_varargin_2;
  real_T c5_db_varargin_3;
  real_T c5_ce_x;
  real_T c5_ae_y;
  real_T c5_de_x;
  real_T c5_be_y;
  real_T c5_s_xk;
  real_T c5_db_yk;
  real_T c5_ee_x;
  real_T c5_ce_y;
  real_T c5_s_minval;
  real_T c5_t_varargin_1;
  real_T c5_xb_varargin_2;
  real_T c5_yb_varargin_2;
  real_T c5_eb_varargin_3;
  real_T c5_fe_x;
  real_T c5_de_y;
  real_T c5_ge_x;
  real_T c5_ee_y;
  real_T c5_t_xk;
  real_T c5_eb_yk;
  real_T c5_he_x;
  real_T c5_fe_y;
  real_T c5_t_minval;
  real_T c5_u_varargin_1;
  real_T c5_ac_varargin_2;
  real_T c5_bc_varargin_2;
  real_T c5_fb_varargin_3;
  real_T c5_ie_x;
  real_T c5_ge_y;
  real_T c5_je_x;
  real_T c5_he_y;
  real_T c5_u_xk;
  real_T c5_fb_yk;
  real_T c5_ke_x;
  real_T c5_ie_y;
  real_T c5_u_minval;
  int32_T c5_d_i;
  real_T c5_v_varargin_1;
  real_T c5_cc_varargin_2;
  real_T c5_dc_varargin_2;
  real_T c5_gb_varargin_3;
  real_T c5_le_x;
  real_T c5_je_y;
  real_T c5_me_x;
  real_T c5_ke_y;
  real_T c5_v_xk;
  real_T c5_gb_yk;
  real_T c5_ne_x;
  real_T c5_le_y;
  real_T c5_v_minval;
  real_T c5_w_varargin_1;
  real_T c5_ec_varargin_2;
  real_T c5_fc_varargin_2;
  real_T c5_hb_varargin_3;
  real_T c5_oe_x;
  real_T c5_me_y;
  real_T c5_pe_x;
  real_T c5_ne_y;
  real_T c5_w_xk;
  real_T c5_hb_yk;
  real_T c5_qe_x;
  real_T c5_oe_y;
  real_T c5_w_minval;
  real_T c5_x_varargin_1;
  real_T c5_gc_varargin_2;
  real_T c5_hc_varargin_2;
  real_T c5_ib_varargin_3;
  real_T c5_re_x;
  real_T c5_pe_y;
  real_T c5_se_x;
  real_T c5_qe_y;
  real_T c5_x_xk;
  real_T c5_ib_yk;
  real_T c5_te_x;
  real_T c5_re_y;
  real_T c5_x_minval;
  real_T c5_y_varargin_1;
  real_T c5_ic_varargin_2;
  real_T c5_jc_varargin_2;
  real_T c5_jb_varargin_3;
  real_T c5_ue_x;
  real_T c5_se_y;
  real_T c5_ve_x;
  real_T c5_te_y;
  real_T c5_y_xk;
  real_T c5_jb_yk;
  real_T c5_we_x;
  real_T c5_ue_y;
  real_T c5_y_minval;
  real_T c5_ab_varargin_1;
  real_T c5_kc_varargin_2;
  real_T c5_lc_varargin_2;
  real_T c5_kb_varargin_3;
  real_T c5_xe_x;
  real_T c5_ve_y;
  real_T c5_ye_x;
  real_T c5_we_y;
  real_T c5_ab_xk;
  real_T c5_kb_yk;
  real_T c5_af_x;
  real_T c5_xe_y;
  real_T c5_ab_minval;
  real_T c5_bb_varargin_1;
  real_T c5_mc_varargin_2;
  real_T c5_nc_varargin_2;
  real_T c5_lb_varargin_3;
  real_T c5_bf_x;
  real_T c5_ye_y;
  real_T c5_cf_x;
  real_T c5_af_y;
  real_T c5_bb_xk;
  real_T c5_lb_yk;
  real_T c5_df_x;
  real_T c5_bf_y;
  real_T c5_bb_minval;
  real_T c5_cb_varargin_1;
  real_T c5_oc_varargin_2;
  real_T c5_pc_varargin_2;
  real_T c5_mb_varargin_3;
  real_T c5_ef_x;
  real_T c5_cf_y;
  real_T c5_ff_x;
  real_T c5_df_y;
  real_T c5_cb_xk;
  real_T c5_mb_yk;
  real_T c5_gf_x;
  real_T c5_ef_y;
  real_T c5_cb_minval;
  int32_T c5_e_i;
  real_T c5_db_varargin_1;
  real_T c5_qc_varargin_2;
  real_T c5_rc_varargin_2;
  real_T c5_nb_varargin_3;
  real_T c5_hf_x;
  real_T c5_ff_y;
  real_T c5_if_x;
  real_T c5_gf_y;
  real_T c5_db_xk;
  real_T c5_nb_yk;
  real_T c5_jf_x;
  real_T c5_hf_y;
  real_T c5_db_minval;
  real_T c5_eb_varargin_1;
  real_T c5_sc_varargin_2;
  real_T c5_tc_varargin_2;
  real_T c5_ob_varargin_3;
  real_T c5_kf_x;
  real_T c5_if_y;
  real_T c5_lf_x;
  real_T c5_jf_y;
  real_T c5_eb_xk;
  real_T c5_ob_yk;
  real_T c5_mf_x;
  real_T c5_kf_y;
  real_T c5_eb_minval;
  real_T c5_fb_varargin_1;
  real_T c5_uc_varargin_2;
  real_T c5_vc_varargin_2;
  real_T c5_pb_varargin_3;
  real_T c5_nf_x;
  real_T c5_lf_y;
  real_T c5_of_x;
  real_T c5_mf_y;
  real_T c5_fb_xk;
  real_T c5_pb_yk;
  real_T c5_pf_x;
  real_T c5_nf_y;
  real_T c5_fb_minval;
  int32_T c5_f_i;
  real_T c5_q_A;
  real_T c5_B;
  real_T c5_qf_x;
  real_T c5_of_y;
  real_T c5_rf_x;
  real_T c5_pf_y;
  real_T c5_sf_x;
  real_T c5_qf_y;
  real_T *c5_b_desireddteta;
  real_T *c5_b_dteta;
  real_T *c5_b_desiredteta;
  real_T *c5_b_teta;
  real_T *c5_b_T;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  boolean_T guard5 = false;
  boolean_T guard6 = false;
  boolean_T guard7 = false;
  boolean_T guard8 = false;
  boolean_T guard9 = false;
  boolean_T guard10 = false;
  boolean_T guard11 = false;
  boolean_T guard12 = false;
  boolean_T guard13 = false;
  boolean_T guard14 = false;
  boolean_T guard15 = false;
  boolean_T guard16 = false;
  boolean_T guard17 = false;
  boolean_T guard18 = false;
  boolean_T guard19 = false;
  boolean_T guard20 = false;
  boolean_T guard21 = false;
  boolean_T guard22 = false;
  boolean_T guard23 = false;
  boolean_T guard24 = false;
  boolean_T guard25 = false;
  boolean_T guard26 = false;
  boolean_T guard27 = false;
  boolean_T guard28 = false;
  boolean_T guard29 = false;
  boolean_T guard30 = false;
  boolean_T guard31 = false;
  boolean_T guard32 = false;
  boolean_T guard33 = false;
  boolean_T guard34 = false;
  boolean_T guard35 = false;
  boolean_T guard36 = false;
  boolean_T guard37 = false;
  boolean_T guard38 = false;
  boolean_T guard39 = false;
  boolean_T guard40 = false;
  boolean_T guard41 = false;
  boolean_T guard42 = false;
  boolean_T guard43 = false;
  boolean_T guard44 = false;
  boolean_T guard45 = false;
  boolean_T guard46 = false;
  boolean_T guard47 = false;
  boolean_T guard48 = false;
  boolean_T guard49 = false;
  boolean_T guard50 = false;
  boolean_T guard51 = false;
  boolean_T guard52 = false;
  boolean_T guard53 = false;
  c5_b_desireddteta = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c5_b_dteta = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c5_b_desiredteta = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c5_b_teta = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c5_b_T = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 4U, chartInstance->c5_sfEvent);
  c5_hoistedGlobal = *c5_b_teta;
  c5_b_hoistedGlobal = *c5_b_desiredteta;
  c5_c_hoistedGlobal = *c5_b_dteta;
  c5_d_hoistedGlobal = *c5_b_desireddteta;
  c5_teta = c5_hoistedGlobal;
  c5_desiredteta = c5_b_hoistedGlobal;
  c5_dteta = c5_c_hoistedGlobal;
  c5_desireddteta = c5_d_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 52U, 52U, c5_debug_family_names,
    c5_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_e, 0U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_pe, 1U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c5_vLR, 2U, c5_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c5_vSR, 3U, c5_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c5_vZ, 4U, c5_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c5_vSL, 5U, c5_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c5_vLL, 6U, c5_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_LAN, 7U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_SMN, 8U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_Ze, 9U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_SMP, 10U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_LAP, 11U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_FAL, 12U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_SLL, 13U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_Zpe, 14U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_SLR, 15U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_FAR, 16U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_muLAN, 17U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_muSMN, 18U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_muZe, 19U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_muSMP, 20U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_muLAP, 21U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_muFAL, 22U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_muSLL, 23U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_muZpe, 24U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_muSLR, 25U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_muFAR, 26U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_outLL1, 27U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_outSL1, 28U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_outZ1, 29U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_outSR1, 30U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_outLR1, 31U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_outLL, 32U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_outSL, 33U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_outZ, 34U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_outSR, 35U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_outLR, 36U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_woutLL, 37U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_woutSL, 38U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_woutZ, 39U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_woutSR, 40U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_woutLR, 41U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_i, 42U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_sumout, 43U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_sumwout, 44U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_nargin, 45U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_nargout, 46U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c5_teta, 47U, c5_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c5_desiredteta, 48U, c5_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c5_dteta, 49U, c5_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c5_desireddteta, 50U, c5_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_T, 51U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 6);
  c5_e = c5_desiredteta - c5_teta;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 7);
  c5_pe = c5_desireddteta - c5_dteta;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 16);
  c5_vLR = -140.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 17);
  c5_vSR = -70.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 18);
  c5_vZ = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 19);
  c5_vSL = 70.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 20);
  c5_vLL = 140.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 23);
  c5_LAN = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 24);
  c5_SMN = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 25);
  c5_Ze = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 26);
  c5_SMP = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 27);
  c5_LAP = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 29);
  c5_FAL = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 30);
  c5_SLL = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 31);
  c5_Zpe = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 32);
  c5_SLR = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 33);
  c5_FAR = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 35);
  c5_muLAN = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 36);
  c5_muSMN = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 37);
  c5_muZe = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 38);
  c5_muSMP = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 39);
  c5_muLAP = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 41);
  c5_muFAL = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 42);
  c5_muSLL = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 43);
  c5_muZpe = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 44);
  c5_muSLR = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 45);
  c5_muFAR = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 47);
  for (c5_i0 = 0; c5_i0 < 7; c5_i0++) {
    c5_outLL1[c5_i0] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 48);
  for (c5_i1 = 0; c5_i1 < 7; c5_i1++) {
    c5_outSL1[c5_i1] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 49);
  for (c5_i2 = 0; c5_i2 < 7; c5_i2++) {
    c5_outZ1[c5_i2] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 50);
  for (c5_i3 = 0; c5_i3 < 7; c5_i3++) {
    c5_outSR1[c5_i3] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 51);
  for (c5_i4 = 0; c5_i4 < 7; c5_i4++) {
    c5_outLR1[c5_i4] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 53);
  c5_outLL = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 54);
  c5_outSL = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 55);
  c5_outZ = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 56);
  c5_outSR = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 57);
  c5_outLR = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 59);
  c5_woutLL = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 60);
  c5_woutSL = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 61);
  c5_woutZ = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 62);
  c5_woutSR = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 63);
  c5_woutLR = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 68);
  c5_A = -1.5707963267948966 - c5_e;
  c5_x = c5_A;
  c5_b_x = c5_x;
  c5_c_x = c5_b_x;
  c5_y = c5_c_x / 1.5707963267948966;
  c5_varargin_2 = c5_y;
  c5_varargin_3 = c5_varargin_2;
  c5_b_y = c5_varargin_3;
  c5_c_y = c5_b_y;
  c5_eml_scalar_eg(chartInstance);
  c5_yk = c5_c_y;
  c5_d_y = c5_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_muLAN = muDoubleScalarMax(0.0, c5_d_y);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 69);
  c5_b_A = c5_e - -3.1415926535897931;
  c5_d_x = c5_b_A;
  c5_e_x = c5_d_x;
  c5_f_x = c5_e_x;
  c5_e_y = c5_f_x / 1.5707963267948966;
  c5_c_A = -c5_e;
  c5_g_x = c5_c_A;
  c5_h_x = c5_g_x;
  c5_i_x = c5_h_x;
  c5_f_y = c5_i_x / 1.5707963267948966;
  c5_varargin_1 = c5_e_y;
  c5_b_varargin_2 = c5_f_y;
  c5_c_varargin_2 = c5_varargin_1;
  c5_b_varargin_3 = c5_b_varargin_2;
  c5_j_x = c5_c_varargin_2;
  c5_g_y = c5_b_varargin_3;
  c5_k_x = c5_j_x;
  c5_h_y = c5_g_y;
  c5_eml_scalar_eg(chartInstance);
  c5_xk = c5_k_x;
  c5_b_yk = c5_h_y;
  c5_l_x = c5_xk;
  c5_i_y = c5_b_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_minval = muDoubleScalarMin(c5_l_x, c5_i_y);
  c5_d_varargin_2 = c5_minval;
  c5_c_varargin_3 = c5_d_varargin_2;
  c5_j_y = c5_c_varargin_3;
  c5_k_y = c5_j_y;
  c5_eml_scalar_eg(chartInstance);
  c5_c_yk = c5_k_y;
  c5_l_y = c5_c_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_muSMN = muDoubleScalarMax(0.0, c5_l_y);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 70);
  c5_d_A = c5_e - -1.5707963267948966;
  c5_m_x = c5_d_A;
  c5_n_x = c5_m_x;
  c5_o_x = c5_n_x;
  c5_m_y = c5_o_x / 1.5707963267948966;
  c5_e_A = 1.5707963267948966 - c5_e;
  c5_p_x = c5_e_A;
  c5_q_x = c5_p_x;
  c5_r_x = c5_q_x;
  c5_n_y = c5_r_x / 1.5707963267948966;
  c5_b_varargin_1 = c5_m_y;
  c5_e_varargin_2 = c5_n_y;
  c5_f_varargin_2 = c5_b_varargin_1;
  c5_d_varargin_3 = c5_e_varargin_2;
  c5_s_x = c5_f_varargin_2;
  c5_o_y = c5_d_varargin_3;
  c5_t_x = c5_s_x;
  c5_p_y = c5_o_y;
  c5_eml_scalar_eg(chartInstance);
  c5_b_xk = c5_t_x;
  c5_d_yk = c5_p_y;
  c5_u_x = c5_b_xk;
  c5_q_y = c5_d_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_b_minval = muDoubleScalarMin(c5_u_x, c5_q_y);
  c5_g_varargin_2 = c5_b_minval;
  c5_e_varargin_3 = c5_g_varargin_2;
  c5_r_y = c5_e_varargin_3;
  c5_s_y = c5_r_y;
  c5_eml_scalar_eg(chartInstance);
  c5_e_yk = c5_s_y;
  c5_t_y = c5_e_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_muZe = muDoubleScalarMax(0.0, c5_t_y);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 71);
  c5_f_A = c5_e;
  c5_v_x = c5_f_A;
  c5_w_x = c5_v_x;
  c5_x_x = c5_w_x;
  c5_u_y = c5_x_x / 1.5707963267948966;
  c5_g_A = 3.1415926535897931 - c5_e;
  c5_y_x = c5_g_A;
  c5_ab_x = c5_y_x;
  c5_bb_x = c5_ab_x;
  c5_v_y = c5_bb_x / 1.5707963267948966;
  c5_c_varargin_1 = c5_u_y;
  c5_h_varargin_2 = c5_v_y;
  c5_i_varargin_2 = c5_c_varargin_1;
  c5_f_varargin_3 = c5_h_varargin_2;
  c5_cb_x = c5_i_varargin_2;
  c5_w_y = c5_f_varargin_3;
  c5_db_x = c5_cb_x;
  c5_x_y = c5_w_y;
  c5_eml_scalar_eg(chartInstance);
  c5_c_xk = c5_db_x;
  c5_f_yk = c5_x_y;
  c5_eb_x = c5_c_xk;
  c5_y_y = c5_f_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_c_minval = muDoubleScalarMin(c5_eb_x, c5_y_y);
  c5_j_varargin_2 = c5_c_minval;
  c5_g_varargin_3 = c5_j_varargin_2;
  c5_ab_y = c5_g_varargin_3;
  c5_bb_y = c5_ab_y;
  c5_eml_scalar_eg(chartInstance);
  c5_g_yk = c5_bb_y;
  c5_cb_y = c5_g_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_muSMP = muDoubleScalarMax(0.0, c5_cb_y);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 72);
  c5_h_A = c5_e - 1.5707963267948966;
  c5_fb_x = c5_h_A;
  c5_gb_x = c5_fb_x;
  c5_hb_x = c5_gb_x;
  c5_db_y = c5_hb_x / 1.5707963267948966;
  c5_k_varargin_2 = c5_db_y;
  c5_h_varargin_3 = c5_k_varargin_2;
  c5_eb_y = c5_h_varargin_3;
  c5_fb_y = c5_eb_y;
  c5_eml_scalar_eg(chartInstance);
  c5_h_yk = c5_fb_y;
  c5_gb_y = c5_h_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_muLAP = muDoubleScalarMax(0.0, c5_gb_y);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 74);
  if (CV_EML_IF(0, 1, 0, c5_e < -3.1415926535897931)) {
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 75);
    c5_muLAN = 1.0;
  } else {
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 76);
    if (CV_EML_IF(0, 1, 1, c5_e > 3.1415926535897931)) {
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 77);
      c5_muLAP = 1.0;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 84);
  c5_i_A = -0.78539816339744828 - c5_pe;
  c5_ib_x = c5_i_A;
  c5_jb_x = c5_ib_x;
  c5_kb_x = c5_jb_x;
  c5_hb_y = c5_kb_x / 0.78539816339744828;
  c5_l_varargin_2 = c5_hb_y;
  c5_i_varargin_3 = c5_l_varargin_2;
  c5_ib_y = c5_i_varargin_3;
  c5_jb_y = c5_ib_y;
  c5_eml_scalar_eg(chartInstance);
  c5_i_yk = c5_jb_y;
  c5_kb_y = c5_i_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_muFAL = muDoubleScalarMax(0.0, c5_kb_y);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 85);
  c5_j_A = c5_pe - -1.5707963267948966;
  c5_lb_x = c5_j_A;
  c5_mb_x = c5_lb_x;
  c5_nb_x = c5_mb_x;
  c5_lb_y = c5_nb_x / 0.78539816339744828;
  c5_k_A = 0.0 - c5_pe;
  c5_ob_x = c5_k_A;
  c5_pb_x = c5_ob_x;
  c5_qb_x = c5_pb_x;
  c5_mb_y = c5_qb_x / 0.78539816339744828;
  c5_d_varargin_1 = c5_lb_y;
  c5_m_varargin_2 = c5_mb_y;
  c5_n_varargin_2 = c5_d_varargin_1;
  c5_j_varargin_3 = c5_m_varargin_2;
  c5_rb_x = c5_n_varargin_2;
  c5_nb_y = c5_j_varargin_3;
  c5_sb_x = c5_rb_x;
  c5_ob_y = c5_nb_y;
  c5_eml_scalar_eg(chartInstance);
  c5_d_xk = c5_sb_x;
  c5_j_yk = c5_ob_y;
  c5_tb_x = c5_d_xk;
  c5_pb_y = c5_j_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_d_minval = muDoubleScalarMin(c5_tb_x, c5_pb_y);
  c5_o_varargin_2 = c5_d_minval;
  c5_k_varargin_3 = c5_o_varargin_2;
  c5_qb_y = c5_k_varargin_3;
  c5_rb_y = c5_qb_y;
  c5_eml_scalar_eg(chartInstance);
  c5_k_yk = c5_rb_y;
  c5_sb_y = c5_k_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_muSLL = muDoubleScalarMax(0.0, c5_sb_y);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 86);
  c5_l_A = c5_pe - -0.78539816339744828;
  c5_ub_x = c5_l_A;
  c5_vb_x = c5_ub_x;
  c5_wb_x = c5_vb_x;
  c5_tb_y = c5_wb_x / 0.78539816339744828;
  c5_m_A = 0.78539816339744828 - c5_pe;
  c5_xb_x = c5_m_A;
  c5_yb_x = c5_xb_x;
  c5_ac_x = c5_yb_x;
  c5_ub_y = c5_ac_x / 0.78539816339744828;
  c5_e_varargin_1 = c5_tb_y;
  c5_p_varargin_2 = c5_ub_y;
  c5_q_varargin_2 = c5_e_varargin_1;
  c5_l_varargin_3 = c5_p_varargin_2;
  c5_bc_x = c5_q_varargin_2;
  c5_vb_y = c5_l_varargin_3;
  c5_cc_x = c5_bc_x;
  c5_wb_y = c5_vb_y;
  c5_eml_scalar_eg(chartInstance);
  c5_e_xk = c5_cc_x;
  c5_l_yk = c5_wb_y;
  c5_dc_x = c5_e_xk;
  c5_xb_y = c5_l_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_e_minval = muDoubleScalarMin(c5_dc_x, c5_xb_y);
  c5_r_varargin_2 = c5_e_minval;
  c5_m_varargin_3 = c5_r_varargin_2;
  c5_yb_y = c5_m_varargin_3;
  c5_ac_y = c5_yb_y;
  c5_eml_scalar_eg(chartInstance);
  c5_m_yk = c5_ac_y;
  c5_bc_y = c5_m_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_muZpe = muDoubleScalarMax(0.0, c5_bc_y);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 87);
  c5_n_A = c5_pe;
  c5_ec_x = c5_n_A;
  c5_fc_x = c5_ec_x;
  c5_gc_x = c5_fc_x;
  c5_cc_y = c5_gc_x / 0.78539816339744828;
  c5_o_A = 1.5707963267948966 - c5_pe;
  c5_hc_x = c5_o_A;
  c5_ic_x = c5_hc_x;
  c5_jc_x = c5_ic_x;
  c5_dc_y = c5_jc_x / 0.78539816339744828;
  c5_f_varargin_1 = c5_cc_y;
  c5_s_varargin_2 = c5_dc_y;
  c5_t_varargin_2 = c5_f_varargin_1;
  c5_n_varargin_3 = c5_s_varargin_2;
  c5_kc_x = c5_t_varargin_2;
  c5_ec_y = c5_n_varargin_3;
  c5_lc_x = c5_kc_x;
  c5_fc_y = c5_ec_y;
  c5_eml_scalar_eg(chartInstance);
  c5_f_xk = c5_lc_x;
  c5_n_yk = c5_fc_y;
  c5_mc_x = c5_f_xk;
  c5_gc_y = c5_n_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_f_minval = muDoubleScalarMin(c5_mc_x, c5_gc_y);
  c5_u_varargin_2 = c5_f_minval;
  c5_o_varargin_3 = c5_u_varargin_2;
  c5_hc_y = c5_o_varargin_3;
  c5_ic_y = c5_hc_y;
  c5_eml_scalar_eg(chartInstance);
  c5_o_yk = c5_ic_y;
  c5_jc_y = c5_o_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_muSLR = muDoubleScalarMax(0.0, c5_jc_y);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 88);
  c5_p_A = c5_pe - 0.78539816339744828;
  c5_nc_x = c5_p_A;
  c5_oc_x = c5_nc_x;
  c5_pc_x = c5_oc_x;
  c5_kc_y = c5_pc_x / 0.78539816339744828;
  c5_v_varargin_2 = c5_kc_y;
  c5_p_varargin_3 = c5_v_varargin_2;
  c5_lc_y = c5_p_varargin_3;
  c5_mc_y = c5_lc_y;
  c5_eml_scalar_eg(chartInstance);
  c5_p_yk = c5_mc_y;
  c5_nc_y = c5_p_yk;
  c5_eml_scalar_eg(chartInstance);
  c5_muFAR = muDoubleScalarMax(0.0, c5_nc_y);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 90);
  if (CV_EML_IF(0, 1, 2, c5_pe < -1.5707963267948966)) {
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 91);
    c5_muFAL = 1.0;
  } else {
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 92);
    if (CV_EML_IF(0, 1, 3, c5_pe > 1.5707963267948966)) {
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 93);
      c5_muFAR = 1.0;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 102);
  if (CV_EML_IF(0, 1, 4, c5_e < -3.1415926535897931)) {
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 104);
    c5_LAN = 1.0;
  } else {
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 106);
    guard50 = false;
    if (CV_EML_COND(0, 1, 0, c5_e > -3.1415926535897931)) {
      if (CV_EML_COND(0, 1, 1, c5_e < -1.5707963267948966)) {
        CV_EML_MCDC(0, 1, 0, true);
        CV_EML_IF(0, 1, 5, true);
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 108);
        c5_LAN = 1.0;
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 109);
        c5_SMN = 1.0;
      } else {
        guard50 = true;
      }
    } else {
      guard50 = true;
    }

    if (guard50 == true) {
      CV_EML_MCDC(0, 1, 0, false);
      CV_EML_IF(0, 1, 5, false);
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 111);
      guard51 = false;
      if (CV_EML_COND(0, 1, 2, c5_e > -1.5707963267948966)) {
        if (CV_EML_COND(0, 1, 3, c5_e < 0.0)) {
          CV_EML_MCDC(0, 1, 1, true);
          CV_EML_IF(0, 1, 6, true);
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 113);
          c5_SMN = 1.0;
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 114);
          c5_Ze = 1.0;
        } else {
          guard51 = true;
        }
      } else {
        guard51 = true;
      }

      if (guard51 == true) {
        CV_EML_MCDC(0, 1, 1, false);
        CV_EML_IF(0, 1, 6, false);
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 116);
        guard52 = false;
        if (CV_EML_COND(0, 1, 4, c5_e > 0.0)) {
          if (CV_EML_COND(0, 1, 5, c5_e < 1.5707963267948966)) {
            CV_EML_MCDC(0, 1, 2, true);
            CV_EML_IF(0, 1, 7, true);
            _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 118);
            c5_Ze = 1.0;
            _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 119);
            c5_SMP = 1.0;
          } else {
            guard52 = true;
          }
        } else {
          guard52 = true;
        }

        if (guard52 == true) {
          CV_EML_MCDC(0, 1, 2, false);
          CV_EML_IF(0, 1, 7, false);
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 121);
          guard53 = false;
          if (CV_EML_COND(0, 1, 6, c5_e > 1.5707963267948966)) {
            if (CV_EML_COND(0, 1, 7, c5_e < 3.1415926535897931)) {
              CV_EML_MCDC(0, 1, 3, true);
              CV_EML_IF(0, 1, 8, true);
              _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 123);
              c5_SMP = 1.0;
              _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 124);
              c5_LAP = 1.0;
            } else {
              guard53 = true;
            }
          } else {
            guard53 = true;
          }

          if (guard53 == true) {
            CV_EML_MCDC(0, 1, 3, false);
            CV_EML_IF(0, 1, 8, false);
            _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 126);
            if (CV_EML_IF(0, 1, 9, c5_e > 3.1415926535897931)) {
              _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 128U);
              c5_LAP = 1.0;
            }
          }
        }
      }
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 135U);
  if (CV_EML_IF(0, 1, 10, c5_pe < -1.5707963267948966)) {
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 137U);
    c5_FAL = 1.0;
  } else {
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 139U);
    guard46 = false;
    if (CV_EML_COND(0, 1, 8, c5_pe > -1.5707963267948966)) {
      if (CV_EML_COND(0, 1, 9, c5_pe < -0.78539816339744828)) {
        CV_EML_MCDC(0, 1, 4, true);
        CV_EML_IF(0, 1, 11, true);
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 141U);
        c5_FAL = 1.0;
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 142U);
        c5_SLL = 1.0;
      } else {
        guard46 = true;
      }
    } else {
      guard46 = true;
    }

    if (guard46 == true) {
      CV_EML_MCDC(0, 1, 4, false);
      CV_EML_IF(0, 1, 11, false);
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 144U);
      guard47 = false;
      if (CV_EML_COND(0, 1, 10, c5_pe > -0.78539816339744828)) {
        if (CV_EML_COND(0, 1, 11, c5_pe < 0.0)) {
          CV_EML_MCDC(0, 1, 5, true);
          CV_EML_IF(0, 1, 12, true);
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 146U);
          c5_SLL = 1.0;
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 147U);
          c5_Zpe = 1.0;
        } else {
          guard47 = true;
        }
      } else {
        guard47 = true;
      }

      if (guard47 == true) {
        CV_EML_MCDC(0, 1, 5, false);
        CV_EML_IF(0, 1, 12, false);
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 149U);
        guard48 = false;
        if (CV_EML_COND(0, 1, 12, c5_pe > 0.0)) {
          if (CV_EML_COND(0, 1, 13, c5_pe < 0.78539816339744828)) {
            CV_EML_MCDC(0, 1, 6, true);
            CV_EML_IF(0, 1, 13, true);
            _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 151U);
            c5_Zpe = 1.0;
            _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 152U);
            c5_SLR = 1.0;
          } else {
            guard48 = true;
          }
        } else {
          guard48 = true;
        }

        if (guard48 == true) {
          CV_EML_MCDC(0, 1, 6, false);
          CV_EML_IF(0, 1, 13, false);
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 154U);
          guard49 = false;
          if (CV_EML_COND(0, 1, 14, c5_pe > 0.78539816339744828)) {
            if (CV_EML_COND(0, 1, 15, c5_pe < 1.5707963267948966)) {
              CV_EML_MCDC(0, 1, 7, true);
              CV_EML_IF(0, 1, 14, true);
              _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 156U);
              c5_SLR = 1.0;
              _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 157U);
              c5_FAR = 1.0;
            } else {
              guard49 = true;
            }
          } else {
            guard49 = true;
          }

          if (guard49 == true) {
            CV_EML_MCDC(0, 1, 7, false);
            CV_EML_IF(0, 1, 14, false);
            _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 159U);
            if (CV_EML_IF(0, 1, 15, c5_pe > 1.5707963267948966)) {
              _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 161U);
              c5_FAR = 1.0;
            }
          }
        }
      }
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 169U);
  guard39 = false;
  guard40 = false;
  guard41 = false;
  guard42 = false;
  guard43 = false;
  guard44 = false;
  guard45 = false;
  if (CV_EML_COND(0, 1, 16, c5_LAN == 1.0)) {
    if (CV_EML_COND(0, 1, 17, c5_FAL == 1.0)) {
      guard44 = true;
    } else {
      guard45 = true;
    }
  } else {
    guard45 = true;
  }

  if (guard45 == true) {
    if (CV_EML_COND(0, 1, 18, c5_LAN == 1.0)) {
      if (CV_EML_COND(0, 1, 19, c5_SLL == 1.0)) {
        guard44 = true;
      } else {
        guard43 = true;
      }
    } else {
      guard43 = true;
    }
  }

  if (guard44 == true) {
    guard42 = true;
  }

  if (guard43 == true) {
    if (CV_EML_COND(0, 1, 20, c5_LAN == 1.0)) {
      if (CV_EML_COND(0, 1, 21, c5_Zpe == 1.0)) {
        guard42 = true;
      } else {
        guard41 = true;
      }
    } else {
      guard41 = true;
    }
  }

  if (guard42 == true) {
    guard40 = true;
  }

  if (guard41 == true) {
    if (CV_EML_COND(0, 1, 22, c5_SMN == 1.0)) {
      if (CV_EML_COND(0, 1, 23, c5_FAL == 1.0)) {
        guard40 = true;
      } else {
        guard39 = true;
      }
    } else {
      guard39 = true;
    }
  }

  if (guard40 == true) {
    CV_EML_MCDC(0, 1, 8, true);
    CV_EML_IF(0, 1, 16, true);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 173U);
    c5_g_varargin_1 = c5_muLAN;
    c5_w_varargin_2 = c5_muFAL;
    c5_x_varargin_2 = c5_g_varargin_1;
    c5_q_varargin_3 = c5_w_varargin_2;
    c5_qc_x = c5_x_varargin_2;
    c5_oc_y = c5_q_varargin_3;
    c5_rc_x = c5_qc_x;
    c5_pc_y = c5_oc_y;
    c5_eml_scalar_eg(chartInstance);
    c5_g_xk = c5_rc_x;
    c5_q_yk = c5_pc_y;
    c5_sc_x = c5_g_xk;
    c5_qc_y = c5_q_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_g_minval = muDoubleScalarMin(c5_sc_x, c5_qc_y);
    c5_outLR1[0] = c5_g_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 174U);
    c5_h_varargin_1 = c5_muLAN;
    c5_y_varargin_2 = c5_muSLL;
    c5_ab_varargin_2 = c5_h_varargin_1;
    c5_r_varargin_3 = c5_y_varargin_2;
    c5_tc_x = c5_ab_varargin_2;
    c5_rc_y = c5_r_varargin_3;
    c5_uc_x = c5_tc_x;
    c5_sc_y = c5_rc_y;
    c5_eml_scalar_eg(chartInstance);
    c5_h_xk = c5_uc_x;
    c5_r_yk = c5_sc_y;
    c5_vc_x = c5_h_xk;
    c5_tc_y = c5_r_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_h_minval = muDoubleScalarMin(c5_vc_x, c5_tc_y);
    c5_outLR1[1] = c5_h_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 175U);
    c5_i_varargin_1 = c5_muLAN;
    c5_bb_varargin_2 = c5_muZpe;
    c5_cb_varargin_2 = c5_i_varargin_1;
    c5_s_varargin_3 = c5_bb_varargin_2;
    c5_wc_x = c5_cb_varargin_2;
    c5_uc_y = c5_s_varargin_3;
    c5_xc_x = c5_wc_x;
    c5_vc_y = c5_uc_y;
    c5_eml_scalar_eg(chartInstance);
    c5_i_xk = c5_xc_x;
    c5_s_yk = c5_vc_y;
    c5_yc_x = c5_i_xk;
    c5_wc_y = c5_s_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_i_minval = muDoubleScalarMin(c5_yc_x, c5_wc_y);
    c5_outLR1[2] = c5_i_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 176U);
    c5_j_varargin_1 = c5_muSMN;
    c5_db_varargin_2 = c5_muFAL;
    c5_eb_varargin_2 = c5_j_varargin_1;
    c5_t_varargin_3 = c5_db_varargin_2;
    c5_ad_x = c5_eb_varargin_2;
    c5_xc_y = c5_t_varargin_3;
    c5_bd_x = c5_ad_x;
    c5_yc_y = c5_xc_y;
    c5_eml_scalar_eg(chartInstance);
    c5_j_xk = c5_bd_x;
    c5_t_yk = c5_yc_y;
    c5_cd_x = c5_j_xk;
    c5_ad_y = c5_t_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_j_minval = muDoubleScalarMin(c5_cd_x, c5_ad_y);
    c5_outLR1[3] = c5_j_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 178U);
    c5_i = 1.0;
    c5_b_i = 0;
    while (c5_b_i < 4) {
      c5_i = 1.0 + (real_T)c5_b_i;
      CV_EML_FOR(0, 1, 0, 1);
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 179U);
      if (CV_EML_IF(0, 1, 17, c5_outLR1[_SFD_EML_ARRAY_BOUNDS_CHECK("outLR1",
            (int32_T)_SFD_INTEGER_CHECK("i", c5_i), 1, 7, 1, 0) - 1] != 0.0)) {
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 180U);
        c5_outLR = c5_outLR1[_SFD_EML_ARRAY_BOUNDS_CHECK("outLR1", (int32_T)
          _SFD_INTEGER_CHECK("i", c5_i), 1, 7, 1, 0) - 1];
      }

      c5_b_i++;
      _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
    }

    CV_EML_FOR(0, 1, 0, 0);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 185U);
    c5_woutLR = c5_outLR * -140.0;
  }

  if (guard39 == true) {
    CV_EML_MCDC(0, 1, 8, false);
    CV_EML_IF(0, 1, 16, false);
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 191U);
  guard28 = false;
  guard29 = false;
  guard30 = false;
  guard31 = false;
  guard32 = false;
  guard33 = false;
  guard34 = false;
  guard35 = false;
  guard36 = false;
  guard37 = false;
  guard38 = false;
  if (CV_EML_COND(0, 1, 24, c5_LAN == 1.0)) {
    if (CV_EML_COND(0, 1, 25, c5_SLR == 1.0)) {
      guard37 = true;
    } else {
      guard38 = true;
    }
  } else {
    guard38 = true;
  }

  if (guard38 == true) {
    if (CV_EML_COND(0, 1, 26, c5_LAN == 1.0)) {
      if (CV_EML_COND(0, 1, 27, c5_FAR == 1.0)) {
        guard37 = true;
      } else {
        guard36 = true;
      }
    } else {
      guard36 = true;
    }
  }

  if (guard37 == true) {
    guard35 = true;
  }

  if (guard36 == true) {
    if (CV_EML_COND(0, 1, 28, c5_SMN == 1.0)) {
      if (CV_EML_COND(0, 1, 29, c5_SLL == 1.0)) {
        guard35 = true;
      } else {
        guard34 = true;
      }
    } else {
      guard34 = true;
    }
  }

  if (guard35 == true) {
    guard33 = true;
  }

  if (guard34 == true) {
    if (CV_EML_COND(0, 1, 30, c5_SMN == 1.0)) {
      if (CV_EML_COND(0, 1, 31, c5_Zpe == 1.0)) {
        guard33 = true;
      } else {
        guard32 = true;
      }
    } else {
      guard32 = true;
    }
  }

  if (guard33 == true) {
    guard31 = true;
  }

  if (guard32 == true) {
    if (CV_EML_COND(0, 1, 32, c5_Ze == 1.0)) {
      if (CV_EML_COND(0, 1, 33, c5_FAL == 1.0)) {
        guard31 = true;
      } else {
        guard30 = true;
      }
    } else {
      guard30 = true;
    }
  }

  if (guard31 == true) {
    guard29 = true;
  }

  if (guard30 == true) {
    if (CV_EML_COND(0, 1, 34, c5_Ze == 1.0)) {
      if (CV_EML_COND(0, 1, 35, c5_SLL == 1.0)) {
        guard29 = true;
      } else {
        guard28 = true;
      }
    } else {
      guard28 = true;
    }
  }

  if (guard29 == true) {
    CV_EML_MCDC(0, 1, 9, true);
    CV_EML_IF(0, 1, 18, true);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 195U);
    c5_k_varargin_1 = c5_muLAN;
    c5_fb_varargin_2 = c5_muSLR;
    c5_gb_varargin_2 = c5_k_varargin_1;
    c5_u_varargin_3 = c5_fb_varargin_2;
    c5_dd_x = c5_gb_varargin_2;
    c5_bd_y = c5_u_varargin_3;
    c5_ed_x = c5_dd_x;
    c5_cd_y = c5_bd_y;
    c5_eml_scalar_eg(chartInstance);
    c5_k_xk = c5_ed_x;
    c5_u_yk = c5_cd_y;
    c5_fd_x = c5_k_xk;
    c5_dd_y = c5_u_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_k_minval = muDoubleScalarMin(c5_fd_x, c5_dd_y);
    c5_outSR1[0] = c5_k_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 196U);
    c5_l_varargin_1 = c5_muLAN;
    c5_hb_varargin_2 = c5_muFAR;
    c5_ib_varargin_2 = c5_l_varargin_1;
    c5_v_varargin_3 = c5_hb_varargin_2;
    c5_gd_x = c5_ib_varargin_2;
    c5_ed_y = c5_v_varargin_3;
    c5_hd_x = c5_gd_x;
    c5_fd_y = c5_ed_y;
    c5_eml_scalar_eg(chartInstance);
    c5_l_xk = c5_hd_x;
    c5_v_yk = c5_fd_y;
    c5_id_x = c5_l_xk;
    c5_gd_y = c5_v_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_l_minval = muDoubleScalarMin(c5_id_x, c5_gd_y);
    c5_outSR1[1] = c5_l_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 197U);
    c5_m_varargin_1 = c5_muSMN;
    c5_jb_varargin_2 = c5_muSLL;
    c5_kb_varargin_2 = c5_m_varargin_1;
    c5_w_varargin_3 = c5_jb_varargin_2;
    c5_jd_x = c5_kb_varargin_2;
    c5_hd_y = c5_w_varargin_3;
    c5_kd_x = c5_jd_x;
    c5_id_y = c5_hd_y;
    c5_eml_scalar_eg(chartInstance);
    c5_m_xk = c5_kd_x;
    c5_w_yk = c5_id_y;
    c5_ld_x = c5_m_xk;
    c5_jd_y = c5_w_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_m_minval = muDoubleScalarMin(c5_ld_x, c5_jd_y);
    c5_outSR1[2] = c5_m_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 198U);
    c5_n_varargin_1 = c5_muSMN;
    c5_lb_varargin_2 = c5_muZpe;
    c5_mb_varargin_2 = c5_n_varargin_1;
    c5_x_varargin_3 = c5_lb_varargin_2;
    c5_md_x = c5_mb_varargin_2;
    c5_kd_y = c5_x_varargin_3;
    c5_nd_x = c5_md_x;
    c5_ld_y = c5_kd_y;
    c5_eml_scalar_eg(chartInstance);
    c5_n_xk = c5_nd_x;
    c5_x_yk = c5_ld_y;
    c5_od_x = c5_n_xk;
    c5_md_y = c5_x_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_n_minval = muDoubleScalarMin(c5_od_x, c5_md_y);
    c5_outSR1[3] = c5_n_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 199U);
    c5_o_varargin_1 = c5_muZe;
    c5_nb_varargin_2 = c5_muFAL;
    c5_ob_varargin_2 = c5_o_varargin_1;
    c5_y_varargin_3 = c5_nb_varargin_2;
    c5_pd_x = c5_ob_varargin_2;
    c5_nd_y = c5_y_varargin_3;
    c5_qd_x = c5_pd_x;
    c5_od_y = c5_nd_y;
    c5_eml_scalar_eg(chartInstance);
    c5_o_xk = c5_qd_x;
    c5_y_yk = c5_od_y;
    c5_rd_x = c5_o_xk;
    c5_pd_y = c5_y_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_o_minval = muDoubleScalarMin(c5_rd_x, c5_pd_y);
    c5_outSR1[4] = c5_o_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 200U);
    c5_p_varargin_1 = c5_muZe;
    c5_pb_varargin_2 = c5_muSLL;
    c5_qb_varargin_2 = c5_p_varargin_1;
    c5_ab_varargin_3 = c5_pb_varargin_2;
    c5_sd_x = c5_qb_varargin_2;
    c5_qd_y = c5_ab_varargin_3;
    c5_td_x = c5_sd_x;
    c5_rd_y = c5_qd_y;
    c5_eml_scalar_eg(chartInstance);
    c5_p_xk = c5_td_x;
    c5_ab_yk = c5_rd_y;
    c5_ud_x = c5_p_xk;
    c5_sd_y = c5_ab_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_p_minval = muDoubleScalarMin(c5_ud_x, c5_sd_y);
    c5_outSR1[5] = c5_p_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 204U);
    c5_i = 1.0;
    c5_c_i = 0;
    while (c5_c_i < 6) {
      c5_i = 1.0 + (real_T)c5_c_i;
      CV_EML_FOR(0, 1, 1, 1);
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 205U);
      if (CV_EML_IF(0, 1, 19, c5_outSR1[_SFD_EML_ARRAY_BOUNDS_CHECK("outSR1",
            (int32_T)_SFD_INTEGER_CHECK("i", c5_i), 1, 7, 1, 0) - 1] != 0.0)) {
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 206U);
        c5_outSR = c5_outSR1[_SFD_EML_ARRAY_BOUNDS_CHECK("outSR1", (int32_T)
          _SFD_INTEGER_CHECK("i", c5_i), 1, 7, 1, 0) - 1];
      }

      c5_c_i++;
      _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
    }

    CV_EML_FOR(0, 1, 1, 0);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 211U);
    c5_woutSR = c5_outSR * -70.0;
  }

  if (guard28 == true) {
    CV_EML_MCDC(0, 1, 9, false);
    CV_EML_IF(0, 1, 18, false);
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 217U);
  guard19 = false;
  guard20 = false;
  guard21 = false;
  guard22 = false;
  guard23 = false;
  guard24 = false;
  guard25 = false;
  guard26 = false;
  guard27 = false;
  if (CV_EML_COND(0, 1, 36, c5_SMN == 1.0)) {
    if (CV_EML_COND(0, 1, 37, c5_SLR == 1.0)) {
      guard26 = true;
    } else {
      guard27 = true;
    }
  } else {
    guard27 = true;
  }

  if (guard27 == true) {
    if (CV_EML_COND(0, 1, 38, c5_SMN == 1.0)) {
      if (CV_EML_COND(0, 1, 39, c5_FAR == 1.0)) {
        guard26 = true;
      } else {
        guard25 = true;
      }
    } else {
      guard25 = true;
    }
  }

  if (guard26 == true) {
    guard24 = true;
  }

  if (guard25 == true) {
    if (CV_EML_COND(0, 1, 40, c5_Ze == 1.0)) {
      if (CV_EML_COND(0, 1, 41, c5_Zpe == 1.0)) {
        guard24 = true;
      } else {
        guard23 = true;
      }
    } else {
      guard23 = true;
    }
  }

  if (guard24 == true) {
    guard22 = true;
  }

  if (guard23 == true) {
    if (CV_EML_COND(0, 1, 42, c5_SMP == 1.0)) {
      if (CV_EML_COND(0, 1, 43, c5_FAL == 1.0)) {
        guard22 = true;
      } else {
        guard21 = true;
      }
    } else {
      guard21 = true;
    }
  }

  if (guard22 == true) {
    guard20 = true;
  }

  if (guard21 == true) {
    if (CV_EML_COND(0, 1, 44, c5_SMP == 1.0)) {
      if (CV_EML_COND(0, 1, 45, c5_SLL == 1.0)) {
        guard20 = true;
      } else {
        guard19 = true;
      }
    } else {
      guard19 = true;
    }
  }

  if (guard20 == true) {
    CV_EML_MCDC(0, 1, 10, true);
    CV_EML_IF(0, 1, 20, true);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 221U);
    c5_q_varargin_1 = c5_muSMN;
    c5_rb_varargin_2 = c5_muSLR;
    c5_sb_varargin_2 = c5_q_varargin_1;
    c5_bb_varargin_3 = c5_rb_varargin_2;
    c5_vd_x = c5_sb_varargin_2;
    c5_td_y = c5_bb_varargin_3;
    c5_wd_x = c5_vd_x;
    c5_ud_y = c5_td_y;
    c5_eml_scalar_eg(chartInstance);
    c5_q_xk = c5_wd_x;
    c5_bb_yk = c5_ud_y;
    c5_xd_x = c5_q_xk;
    c5_vd_y = c5_bb_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_q_minval = muDoubleScalarMin(c5_xd_x, c5_vd_y);
    c5_outZ1[0] = c5_q_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 222U);
    c5_r_varargin_1 = c5_muSMN;
    c5_tb_varargin_2 = c5_muFAR;
    c5_ub_varargin_2 = c5_r_varargin_1;
    c5_cb_varargin_3 = c5_tb_varargin_2;
    c5_yd_x = c5_ub_varargin_2;
    c5_wd_y = c5_cb_varargin_3;
    c5_ae_x = c5_yd_x;
    c5_xd_y = c5_wd_y;
    c5_eml_scalar_eg(chartInstance);
    c5_r_xk = c5_ae_x;
    c5_cb_yk = c5_xd_y;
    c5_be_x = c5_r_xk;
    c5_yd_y = c5_cb_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_r_minval = muDoubleScalarMin(c5_be_x, c5_yd_y);
    c5_outZ1[1] = c5_r_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 223U);
    c5_s_varargin_1 = c5_muZe;
    c5_vb_varargin_2 = c5_muZpe;
    c5_wb_varargin_2 = c5_s_varargin_1;
    c5_db_varargin_3 = c5_vb_varargin_2;
    c5_ce_x = c5_wb_varargin_2;
    c5_ae_y = c5_db_varargin_3;
    c5_de_x = c5_ce_x;
    c5_be_y = c5_ae_y;
    c5_eml_scalar_eg(chartInstance);
    c5_s_xk = c5_de_x;
    c5_db_yk = c5_be_y;
    c5_ee_x = c5_s_xk;
    c5_ce_y = c5_db_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_s_minval = muDoubleScalarMin(c5_ee_x, c5_ce_y);
    c5_outZ1[2] = c5_s_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 224U);
    c5_t_varargin_1 = c5_muSMP;
    c5_xb_varargin_2 = c5_muFAL;
    c5_yb_varargin_2 = c5_t_varargin_1;
    c5_eb_varargin_3 = c5_xb_varargin_2;
    c5_fe_x = c5_yb_varargin_2;
    c5_de_y = c5_eb_varargin_3;
    c5_ge_x = c5_fe_x;
    c5_ee_y = c5_de_y;
    c5_eml_scalar_eg(chartInstance);
    c5_t_xk = c5_ge_x;
    c5_eb_yk = c5_ee_y;
    c5_he_x = c5_t_xk;
    c5_fe_y = c5_eb_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_t_minval = muDoubleScalarMin(c5_he_x, c5_fe_y);
    c5_outZ1[3] = c5_t_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 225U);
    c5_u_varargin_1 = c5_muSMP;
    c5_ac_varargin_2 = c5_muSLL;
    c5_bc_varargin_2 = c5_u_varargin_1;
    c5_fb_varargin_3 = c5_ac_varargin_2;
    c5_ie_x = c5_bc_varargin_2;
    c5_ge_y = c5_fb_varargin_3;
    c5_je_x = c5_ie_x;
    c5_he_y = c5_ge_y;
    c5_eml_scalar_eg(chartInstance);
    c5_u_xk = c5_je_x;
    c5_fb_yk = c5_he_y;
    c5_ke_x = c5_u_xk;
    c5_ie_y = c5_fb_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_u_minval = muDoubleScalarMin(c5_ke_x, c5_ie_y);
    c5_outZ1[4] = c5_u_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 227U);
    c5_i = 1.0;
    c5_d_i = 0;
    while (c5_d_i < 5) {
      c5_i = 1.0 + (real_T)c5_d_i;
      CV_EML_FOR(0, 1, 2, 1);
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 228U);
      if (CV_EML_IF(0, 1, 21, c5_outZ1[_SFD_EML_ARRAY_BOUNDS_CHECK("outZ1",
            (int32_T)_SFD_INTEGER_CHECK("i", c5_i), 1, 7, 1, 0) - 1] != 0.0)) {
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 229U);
        c5_outZ = c5_outZ1[_SFD_EML_ARRAY_BOUNDS_CHECK("outZ1", (int32_T)
          _SFD_INTEGER_CHECK("i", c5_i), 1, 7, 1, 0) - 1];
      }

      c5_d_i++;
      _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
    }

    CV_EML_FOR(0, 1, 2, 0);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 234U);
    c5_woutZ = c5_outZ * 0.0;
  }

  if (guard19 == true) {
    CV_EML_MCDC(0, 1, 10, false);
    CV_EML_IF(0, 1, 20, false);
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 241U);
  guard6 = false;
  guard7 = false;
  guard8 = false;
  guard9 = false;
  guard10 = false;
  guard11 = false;
  guard12 = false;
  guard13 = false;
  guard14 = false;
  guard15 = false;
  guard16 = false;
  guard17 = false;
  guard18 = false;
  if (CV_EML_COND(0, 1, 46, c5_Ze == 1.0)) {
    if (CV_EML_COND(0, 1, 47, c5_SLR == 1.0)) {
      guard17 = true;
    } else {
      guard18 = true;
    }
  } else {
    guard18 = true;
  }

  if (guard18 == true) {
    if (CV_EML_COND(0, 1, 48, c5_Ze == 1.0)) {
      if (CV_EML_COND(0, 1, 49, c5_FAR == 1.0)) {
        guard17 = true;
      } else {
        guard16 = true;
      }
    } else {
      guard16 = true;
    }
  }

  if (guard17 == true) {
    guard15 = true;
  }

  if (guard16 == true) {
    if (CV_EML_COND(0, 1, 50, c5_SMP == 1.0)) {
      if (CV_EML_COND(0, 1, 51, c5_Zpe == 1.0)) {
        guard15 = true;
      } else {
        guard14 = true;
      }
    } else {
      guard14 = true;
    }
  }

  if (guard15 == true) {
    guard13 = true;
  }

  if (guard14 == true) {
    if (CV_EML_COND(0, 1, 52, c5_SMP == 1.0)) {
      if (CV_EML_COND(0, 1, 53, c5_SLR == 1.0)) {
        guard13 = true;
      } else {
        guard12 = true;
      }
    } else {
      guard12 = true;
    }
  }

  if (guard13 == true) {
    guard11 = true;
  }

  if (guard12 == true) {
    if (CV_EML_COND(0, 1, 54, c5_SMP == 1.0)) {
      if (CV_EML_COND(0, 1, 55, c5_FAR == 1.0)) {
        guard11 = true;
      } else {
        guard10 = true;
      }
    } else {
      guard10 = true;
    }
  }

  if (guard11 == true) {
    guard9 = true;
  }

  if (guard10 == true) {
    if (CV_EML_COND(0, 1, 56, c5_LAP == 1.0)) {
      if (CV_EML_COND(0, 1, 57, c5_FAL == 1.0)) {
        guard9 = true;
      } else {
        guard8 = true;
      }
    } else {
      guard8 = true;
    }
  }

  if (guard9 == true) {
    guard7 = true;
  }

  if (guard8 == true) {
    if (CV_EML_COND(0, 1, 58, c5_LAP == 1.0)) {
      if (CV_EML_COND(0, 1, 59, c5_SLL == 1.0)) {
        guard7 = true;
      } else {
        guard6 = true;
      }
    } else {
      guard6 = true;
    }
  }

  if (guard7 == true) {
    CV_EML_MCDC(0, 1, 11, true);
    CV_EML_IF(0, 1, 22, true);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 245U);
    c5_v_varargin_1 = c5_muZe;
    c5_cc_varargin_2 = c5_muSLR;
    c5_dc_varargin_2 = c5_v_varargin_1;
    c5_gb_varargin_3 = c5_cc_varargin_2;
    c5_le_x = c5_dc_varargin_2;
    c5_je_y = c5_gb_varargin_3;
    c5_me_x = c5_le_x;
    c5_ke_y = c5_je_y;
    c5_eml_scalar_eg(chartInstance);
    c5_v_xk = c5_me_x;
    c5_gb_yk = c5_ke_y;
    c5_ne_x = c5_v_xk;
    c5_le_y = c5_gb_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_v_minval = muDoubleScalarMin(c5_ne_x, c5_le_y);
    c5_outSL1[0] = c5_v_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 246U);
    c5_w_varargin_1 = c5_muZe;
    c5_ec_varargin_2 = c5_muFAR;
    c5_fc_varargin_2 = c5_w_varargin_1;
    c5_hb_varargin_3 = c5_ec_varargin_2;
    c5_oe_x = c5_fc_varargin_2;
    c5_me_y = c5_hb_varargin_3;
    c5_pe_x = c5_oe_x;
    c5_ne_y = c5_me_y;
    c5_eml_scalar_eg(chartInstance);
    c5_w_xk = c5_pe_x;
    c5_hb_yk = c5_ne_y;
    c5_qe_x = c5_w_xk;
    c5_oe_y = c5_hb_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_w_minval = muDoubleScalarMin(c5_qe_x, c5_oe_y);
    c5_outSL1[1] = c5_w_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 247U);
    c5_x_varargin_1 = c5_muSMP;
    c5_gc_varargin_2 = c5_muZpe;
    c5_hc_varargin_2 = c5_x_varargin_1;
    c5_ib_varargin_3 = c5_gc_varargin_2;
    c5_re_x = c5_hc_varargin_2;
    c5_pe_y = c5_ib_varargin_3;
    c5_se_x = c5_re_x;
    c5_qe_y = c5_pe_y;
    c5_eml_scalar_eg(chartInstance);
    c5_x_xk = c5_se_x;
    c5_ib_yk = c5_qe_y;
    c5_te_x = c5_x_xk;
    c5_re_y = c5_ib_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_x_minval = muDoubleScalarMin(c5_te_x, c5_re_y);
    c5_outSL1[2] = c5_x_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 248U);
    c5_y_varargin_1 = c5_muSMP;
    c5_ic_varargin_2 = c5_muSLR;
    c5_jc_varargin_2 = c5_y_varargin_1;
    c5_jb_varargin_3 = c5_ic_varargin_2;
    c5_ue_x = c5_jc_varargin_2;
    c5_se_y = c5_jb_varargin_3;
    c5_ve_x = c5_ue_x;
    c5_te_y = c5_se_y;
    c5_eml_scalar_eg(chartInstance);
    c5_y_xk = c5_ve_x;
    c5_jb_yk = c5_te_y;
    c5_we_x = c5_y_xk;
    c5_ue_y = c5_jb_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_y_minval = muDoubleScalarMin(c5_we_x, c5_ue_y);
    c5_outSL1[3] = c5_y_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 249U);
    c5_ab_varargin_1 = c5_muSMP;
    c5_kc_varargin_2 = c5_muFAR;
    c5_lc_varargin_2 = c5_ab_varargin_1;
    c5_kb_varargin_3 = c5_kc_varargin_2;
    c5_xe_x = c5_lc_varargin_2;
    c5_ve_y = c5_kb_varargin_3;
    c5_ye_x = c5_xe_x;
    c5_we_y = c5_ve_y;
    c5_eml_scalar_eg(chartInstance);
    c5_ab_xk = c5_ye_x;
    c5_kb_yk = c5_we_y;
    c5_af_x = c5_ab_xk;
    c5_xe_y = c5_kb_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_ab_minval = muDoubleScalarMin(c5_af_x, c5_xe_y);
    c5_outSL1[4] = c5_ab_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 250U);
    c5_bb_varargin_1 = c5_muLAP;
    c5_mc_varargin_2 = c5_muFAL;
    c5_nc_varargin_2 = c5_bb_varargin_1;
    c5_lb_varargin_3 = c5_mc_varargin_2;
    c5_bf_x = c5_nc_varargin_2;
    c5_ye_y = c5_lb_varargin_3;
    c5_cf_x = c5_bf_x;
    c5_af_y = c5_ye_y;
    c5_eml_scalar_eg(chartInstance);
    c5_bb_xk = c5_cf_x;
    c5_lb_yk = c5_af_y;
    c5_df_x = c5_bb_xk;
    c5_bf_y = c5_lb_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_bb_minval = muDoubleScalarMin(c5_df_x, c5_bf_y);
    c5_outSL1[5] = c5_bb_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 251U);
    c5_cb_varargin_1 = c5_muLAP;
    c5_oc_varargin_2 = c5_muSLL;
    c5_pc_varargin_2 = c5_cb_varargin_1;
    c5_mb_varargin_3 = c5_oc_varargin_2;
    c5_ef_x = c5_pc_varargin_2;
    c5_cf_y = c5_mb_varargin_3;
    c5_ff_x = c5_ef_x;
    c5_df_y = c5_cf_y;
    c5_eml_scalar_eg(chartInstance);
    c5_cb_xk = c5_ff_x;
    c5_mb_yk = c5_df_y;
    c5_gf_x = c5_cb_xk;
    c5_ef_y = c5_mb_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_cb_minval = muDoubleScalarMin(c5_gf_x, c5_ef_y);
    c5_outSL1[6] = c5_cb_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 253U);
    c5_i = 1.0;
    c5_e_i = 0;
    while (c5_e_i < 7) {
      c5_i = 1.0 + (real_T)c5_e_i;
      CV_EML_FOR(0, 1, 3, 1);
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 254U);
      if (CV_EML_IF(0, 1, 23, c5_outSL1[_SFD_EML_ARRAY_BOUNDS_CHECK("outSL1",
            (int32_T)_SFD_INTEGER_CHECK("i", c5_i), 1, 7, 1, 0) - 1] != 0.0)) {
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, MAX_uint8_T);
        c5_outSL = c5_outSL1[_SFD_EML_ARRAY_BOUNDS_CHECK("outSL1", (int32_T)
          _SFD_INTEGER_CHECK("i", c5_i), 1, 7, 1, 0) - 1];
      }

      c5_e_i++;
      _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
    }

    CV_EML_FOR(0, 1, 3, 0);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 260);
    c5_woutSL = c5_outSL * 70.0;
  }

  if (guard6 == true) {
    CV_EML_MCDC(0, 1, 11, false);
    CV_EML_IF(0, 1, 22, false);
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 266);
  guard1 = false;
  guard2 = false;
  guard3 = false;
  guard4 = false;
  guard5 = false;
  if (CV_EML_COND(0, 1, 60, c5_LAP == 1.0)) {
    if (CV_EML_COND(0, 1, 61, c5_Zpe == 1.0)) {
      guard4 = true;
    } else {
      guard5 = true;
    }
  } else {
    guard5 = true;
  }

  if (guard5 == true) {
    if (CV_EML_COND(0, 1, 62, c5_LAP == 1.0)) {
      if (CV_EML_COND(0, 1, 63, c5_SLR == 1.0)) {
        guard4 = true;
      } else {
        guard3 = true;
      }
    } else {
      guard3 = true;
    }
  }

  if (guard4 == true) {
    guard2 = true;
  }

  if (guard3 == true) {
    if (CV_EML_COND(0, 1, 64, c5_LAP == 1.0)) {
      if (CV_EML_COND(0, 1, 65, c5_FAR == 1.0)) {
        guard2 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard2 == true) {
    CV_EML_MCDC(0, 1, 12, true);
    CV_EML_IF(0, 1, 24, true);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 270);
    c5_db_varargin_1 = c5_muLAP;
    c5_qc_varargin_2 = c5_muZpe;
    c5_rc_varargin_2 = c5_db_varargin_1;
    c5_nb_varargin_3 = c5_qc_varargin_2;
    c5_hf_x = c5_rc_varargin_2;
    c5_ff_y = c5_nb_varargin_3;
    c5_if_x = c5_hf_x;
    c5_gf_y = c5_ff_y;
    c5_eml_scalar_eg(chartInstance);
    c5_db_xk = c5_if_x;
    c5_nb_yk = c5_gf_y;
    c5_jf_x = c5_db_xk;
    c5_hf_y = c5_nb_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_db_minval = muDoubleScalarMin(c5_jf_x, c5_hf_y);
    c5_outLL1[0] = c5_db_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 271);
    c5_eb_varargin_1 = c5_muLAP;
    c5_sc_varargin_2 = c5_muSLR;
    c5_tc_varargin_2 = c5_eb_varargin_1;
    c5_ob_varargin_3 = c5_sc_varargin_2;
    c5_kf_x = c5_tc_varargin_2;
    c5_if_y = c5_ob_varargin_3;
    c5_lf_x = c5_kf_x;
    c5_jf_y = c5_if_y;
    c5_eml_scalar_eg(chartInstance);
    c5_eb_xk = c5_lf_x;
    c5_ob_yk = c5_jf_y;
    c5_mf_x = c5_eb_xk;
    c5_kf_y = c5_ob_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_eb_minval = muDoubleScalarMin(c5_mf_x, c5_kf_y);
    c5_outLL1[1] = c5_eb_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 272);
    c5_fb_varargin_1 = c5_muLAP;
    c5_uc_varargin_2 = c5_muFAR;
    c5_vc_varargin_2 = c5_fb_varargin_1;
    c5_pb_varargin_3 = c5_uc_varargin_2;
    c5_nf_x = c5_vc_varargin_2;
    c5_lf_y = c5_pb_varargin_3;
    c5_of_x = c5_nf_x;
    c5_mf_y = c5_lf_y;
    c5_eml_scalar_eg(chartInstance);
    c5_fb_xk = c5_of_x;
    c5_pb_yk = c5_mf_y;
    c5_pf_x = c5_fb_xk;
    c5_nf_y = c5_pb_yk;
    c5_eml_scalar_eg(chartInstance);
    c5_fb_minval = muDoubleScalarMin(c5_pf_x, c5_nf_y);
    c5_outLL1[2] = c5_fb_minval;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 274);
    c5_i = 1.0;
    c5_f_i = 0;
    while (c5_f_i < 3) {
      c5_i = 1.0 + (real_T)c5_f_i;
      CV_EML_FOR(0, 1, 4, 1);
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 275);
      if (CV_EML_IF(0, 1, 25, c5_outLL1[_SFD_EML_ARRAY_BOUNDS_CHECK("outLL1",
            (int32_T)_SFD_INTEGER_CHECK("i", c5_i), 1, 7, 1, 0) - 1] != 0.0)) {
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 276);
        c5_outLL = c5_outLL1[_SFD_EML_ARRAY_BOUNDS_CHECK("outLL1", (int32_T)
          _SFD_INTEGER_CHECK("i", c5_i), 1, 7, 1, 0) - 1];
      }

      c5_f_i++;
      _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
    }

    CV_EML_FOR(0, 1, 4, 0);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 281);
    c5_woutLL = c5_outLL * 140.0;
  }

  if (guard1 == true) {
    CV_EML_MCDC(0, 1, 12, false);
    CV_EML_IF(0, 1, 24, false);
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 290);
  c5_sumout = (((c5_outLL + c5_outSL) + c5_outZ) + c5_outSR) + c5_outLR;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 291);
  c5_sumwout = (((c5_woutLL + c5_woutSL) + c5_woutZ) + c5_woutSR) + c5_woutLR;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 292);
  c5_q_A = c5_sumwout;
  c5_B = c5_sumout;
  c5_qf_x = c5_q_A;
  c5_of_y = c5_B;
  c5_rf_x = c5_qf_x;
  c5_pf_y = c5_of_y;
  c5_sf_x = c5_rf_x;
  c5_qf_y = c5_pf_y;
  c5_T = c5_sf_x / c5_qf_y;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, -292);
  _SFD_SYMBOL_SCOPE_POP();
  *c5_b_T = c5_T;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 4U, chartInstance->c5_sfEvent);
}

static void initSimStructsc5_PUMA_fuzzy(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c5_machineNumber, uint32_T
  c5_chartNumber, uint32_T c5_instanceNumber)
{
  (void)c5_machineNumber;
  (void)c5_chartNumber;
  (void)c5_instanceNumber;
}

static const mxArray *c5_sf_marshallOut(void *chartInstanceVoid, void *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  real_T c5_u;
  const mxArray *c5_y = NULL;
  SFc5_PUMA_fuzzyInstanceStruct *chartInstance;
  chartInstance = (SFc5_PUMA_fuzzyInstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  c5_u = *(real_T *)c5_inData;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", &c5_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, false);
  return c5_mxArrayOutData;
}

static real_T c5_emlrt_marshallIn(SFc5_PUMA_fuzzyInstanceStruct *chartInstance,
  const mxArray *c5_T, const char_T *c5_identifier)
{
  real_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_T), &c5_thisId);
  sf_mex_destroy(&c5_T);
  return c5_y;
}

static real_T c5_b_emlrt_marshallIn(SFc5_PUMA_fuzzyInstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  real_T c5_y;
  real_T c5_d0;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_d0, 1, 0, 0U, 0, 0U, 0);
  c5_y = c5_d0;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_T;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  real_T c5_y;
  SFc5_PUMA_fuzzyInstanceStruct *chartInstance;
  chartInstance = (SFc5_PUMA_fuzzyInstanceStruct *)chartInstanceVoid;
  c5_T = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_T), &c5_thisId);
  sf_mex_destroy(&c5_T);
  *(real_T *)c5_outData = c5_y;
  sf_mex_destroy(&c5_mxArrayInData);
}

static const mxArray *c5_b_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  int32_T c5_i5;
  real_T c5_b_inData[7];
  int32_T c5_i6;
  real_T c5_u[7];
  const mxArray *c5_y = NULL;
  SFc5_PUMA_fuzzyInstanceStruct *chartInstance;
  chartInstance = (SFc5_PUMA_fuzzyInstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  for (c5_i5 = 0; c5_i5 < 7; c5_i5++) {
    c5_b_inData[c5_i5] = (*(real_T (*)[7])c5_inData)[c5_i5];
  }

  for (c5_i6 = 0; c5_i6 < 7; c5_i6++) {
    c5_u[c5_i6] = c5_b_inData[c5_i6];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 0, 0U, 1U, 0U, 2, 1, 7), false);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, false);
  return c5_mxArrayOutData;
}

static void c5_c_emlrt_marshallIn(SFc5_PUMA_fuzzyInstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, real_T c5_y[7])
{
  real_T c5_dv0[7];
  int32_T c5_i7;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv0, 1, 0, 0U, 1, 0U, 2, 1, 7);
  for (c5_i7 = 0; c5_i7 < 7; c5_i7++) {
    c5_y[c5_i7] = c5_dv0[c5_i7];
  }

  sf_mex_destroy(&c5_u);
}

static void c5_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_outLR1;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  real_T c5_y[7];
  int32_T c5_i8;
  SFc5_PUMA_fuzzyInstanceStruct *chartInstance;
  chartInstance = (SFc5_PUMA_fuzzyInstanceStruct *)chartInstanceVoid;
  c5_outLR1 = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_outLR1), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_outLR1);
  for (c5_i8 = 0; c5_i8 < 7; c5_i8++) {
    (*(real_T (*)[7])c5_outData)[c5_i8] = c5_y[c5_i8];
  }

  sf_mex_destroy(&c5_mxArrayInData);
}

const mxArray *sf_c5_PUMA_fuzzy_get_eml_resolved_functions_info(void)
{
  const mxArray *c5_nameCaptureInfo = NULL;
  c5_nameCaptureInfo = NULL;
  sf_mex_assign(&c5_nameCaptureInfo, sf_mex_createstruct("structure", 2, 18, 1),
                false);
  c5_info_helper(&c5_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c5_nameCaptureInfo);
  return c5_nameCaptureInfo;
}

static void c5_info_helper(const mxArray **c5_info)
{
  const mxArray *c5_rhs0 = NULL;
  const mxArray *c5_lhs0 = NULL;
  const mxArray *c5_rhs1 = NULL;
  const mxArray *c5_lhs1 = NULL;
  const mxArray *c5_rhs2 = NULL;
  const mxArray *c5_lhs2 = NULL;
  const mxArray *c5_rhs3 = NULL;
  const mxArray *c5_lhs3 = NULL;
  const mxArray *c5_rhs4 = NULL;
  const mxArray *c5_lhs4 = NULL;
  const mxArray *c5_rhs5 = NULL;
  const mxArray *c5_lhs5 = NULL;
  const mxArray *c5_rhs6 = NULL;
  const mxArray *c5_lhs6 = NULL;
  const mxArray *c5_rhs7 = NULL;
  const mxArray *c5_lhs7 = NULL;
  const mxArray *c5_rhs8 = NULL;
  const mxArray *c5_lhs8 = NULL;
  const mxArray *c5_rhs9 = NULL;
  const mxArray *c5_lhs9 = NULL;
  const mxArray *c5_rhs10 = NULL;
  const mxArray *c5_lhs10 = NULL;
  const mxArray *c5_rhs11 = NULL;
  const mxArray *c5_lhs11 = NULL;
  const mxArray *c5_rhs12 = NULL;
  const mxArray *c5_lhs12 = NULL;
  const mxArray *c5_rhs13 = NULL;
  const mxArray *c5_lhs13 = NULL;
  const mxArray *c5_rhs14 = NULL;
  const mxArray *c5_lhs14 = NULL;
  const mxArray *c5_rhs15 = NULL;
  const mxArray *c5_lhs15 = NULL;
  const mxArray *c5_rhs16 = NULL;
  const mxArray *c5_lhs16 = NULL;
  const mxArray *c5_rhs17 = NULL;
  const mxArray *c5_lhs17 = NULL;
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("mrdivide"), "name", "name", 0);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c5_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 1);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 1);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c5_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 2);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("rdivide"), "name", "name", 2);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c5_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c5_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 4);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c5_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_div"), "name", "name", 5);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c5_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 6);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c5_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 7);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("max"), "name", "name", 7);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1311262516U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c5_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 8);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1378303184U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c5_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 9);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 9);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c5_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 10);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c5_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 11);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 11);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c5_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 12);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c5_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 13);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 13);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c5_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 14);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 14);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c5_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 15);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 15);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c5_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 16);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("min"), "name", "name", 16);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 16);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1311262518U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c5_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 17);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1378303184U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c5_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs17), "lhs", "lhs",
                  17);
  sf_mex_destroy(&c5_rhs0);
  sf_mex_destroy(&c5_lhs0);
  sf_mex_destroy(&c5_rhs1);
  sf_mex_destroy(&c5_lhs1);
  sf_mex_destroy(&c5_rhs2);
  sf_mex_destroy(&c5_lhs2);
  sf_mex_destroy(&c5_rhs3);
  sf_mex_destroy(&c5_lhs3);
  sf_mex_destroy(&c5_rhs4);
  sf_mex_destroy(&c5_lhs4);
  sf_mex_destroy(&c5_rhs5);
  sf_mex_destroy(&c5_lhs5);
  sf_mex_destroy(&c5_rhs6);
  sf_mex_destroy(&c5_lhs6);
  sf_mex_destroy(&c5_rhs7);
  sf_mex_destroy(&c5_lhs7);
  sf_mex_destroy(&c5_rhs8);
  sf_mex_destroy(&c5_lhs8);
  sf_mex_destroy(&c5_rhs9);
  sf_mex_destroy(&c5_lhs9);
  sf_mex_destroy(&c5_rhs10);
  sf_mex_destroy(&c5_lhs10);
  sf_mex_destroy(&c5_rhs11);
  sf_mex_destroy(&c5_lhs11);
  sf_mex_destroy(&c5_rhs12);
  sf_mex_destroy(&c5_lhs12);
  sf_mex_destroy(&c5_rhs13);
  sf_mex_destroy(&c5_lhs13);
  sf_mex_destroy(&c5_rhs14);
  sf_mex_destroy(&c5_lhs14);
  sf_mex_destroy(&c5_rhs15);
  sf_mex_destroy(&c5_lhs15);
  sf_mex_destroy(&c5_rhs16);
  sf_mex_destroy(&c5_lhs16);
  sf_mex_destroy(&c5_rhs17);
  sf_mex_destroy(&c5_lhs17);
}

static const mxArray *c5_emlrt_marshallOut(const char * c5_u)
{
  const mxArray *c5_y = NULL;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c5_u)), false);
  return c5_y;
}

static const mxArray *c5_b_emlrt_marshallOut(const uint32_T c5_u)
{
  const mxArray *c5_y = NULL;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", &c5_u, 7, 0U, 0U, 0U, 0), false);
  return c5_y;
}

static void c5_eml_scalar_eg(SFc5_PUMA_fuzzyInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *c5_c_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  int32_T c5_u;
  const mxArray *c5_y = NULL;
  SFc5_PUMA_fuzzyInstanceStruct *chartInstance;
  chartInstance = (SFc5_PUMA_fuzzyInstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  c5_u = *(int32_T *)c5_inData;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", &c5_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, false);
  return c5_mxArrayOutData;
}

static int32_T c5_d_emlrt_marshallIn(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  int32_T c5_y;
  int32_T c5_i9;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_i9, 1, 6, 0U, 0, 0U, 0);
  c5_y = c5_i9;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_b_sfEvent;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  int32_T c5_y;
  SFc5_PUMA_fuzzyInstanceStruct *chartInstance;
  chartInstance = (SFc5_PUMA_fuzzyInstanceStruct *)chartInstanceVoid;
  c5_b_sfEvent = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_sfEvent),
    &c5_thisId);
  sf_mex_destroy(&c5_b_sfEvent);
  *(int32_T *)c5_outData = c5_y;
  sf_mex_destroy(&c5_mxArrayInData);
}

static uint8_T c5_e_emlrt_marshallIn(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance, const mxArray *c5_b_is_active_c5_PUMA_fuzzy, const char_T
  *c5_identifier)
{
  uint8_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c5_b_is_active_c5_PUMA_fuzzy), &c5_thisId);
  sf_mex_destroy(&c5_b_is_active_c5_PUMA_fuzzy);
  return c5_y;
}

static uint8_T c5_f_emlrt_marshallIn(SFc5_PUMA_fuzzyInstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  uint8_T c5_y;
  uint8_T c5_u0;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_u0, 1, 3, 0U, 0, 0U, 0);
  c5_y = c5_u0;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void init_dsm_address_info(SFc5_PUMA_fuzzyInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c5_PUMA_fuzzy_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(4065800326U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3740360445U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2271397531U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3485960127U);
}

mxArray *sf_c5_PUMA_fuzzy_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("IldpWiAaXSlrDqOr8XalaH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c5_PUMA_fuzzy_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c5_PUMA_fuzzy_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c5_PUMA_fuzzy(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[4],T\"T\",},{M[8],M[0],T\"is_active_c5_PUMA_fuzzy\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c5_PUMA_fuzzy_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc5_PUMA_fuzzyInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc5_PUMA_fuzzyInstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _PUMA_fuzzyMachineNumber_,
           5,
           1,
           1,
           0,
           5,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize ist own list of scripts */
        init_script_number_translation(_PUMA_fuzzyMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_PUMA_fuzzyMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _PUMA_fuzzyMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,2,0,1,"T");
          _SFD_SET_DATA_PROPS(1,1,1,0,"teta");
          _SFD_SET_DATA_PROPS(2,1,1,0,"desiredteta");
          _SFD_SET_DATA_PROPS(3,1,1,0,"dteta");
          _SFD_SET_DATA_PROPS(4,1,1,0,"desireddteta");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,26,0,0,0,5,0,66,13);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,7663);
        _SFD_CV_INIT_EML_IF(0,1,0,1548,1560,1584,1599);
        _SFD_CV_INIT_EML_IF(0,1,1,1584,1599,-1,1599);
        _SFD_CV_INIT_EML_IF(0,1,2,2180,2195,2219,2237);
        _SFD_CV_INIT_EML_IF(0,1,3,2219,2237,-1,2237);
        _SFD_CV_INIT_EML_IF(0,1,4,2416,2428,2514,2543);
        _SFD_CV_INIT_EML_IF(0,1,5,2514,2543,2630,3058);
        _SFD_CV_INIT_EML_IF(0,1,6,2630,2657,2742,3058);
        _SFD_CV_INIT_EML_IF(0,1,7,2742,2768,2856,3058);
        _SFD_CV_INIT_EML_IF(0,1,8,2856,2883,2976,3058);
        _SFD_CV_INIT_EML_IF(0,1,9,2976,2991,-1,2991);
        _SFD_CV_INIT_EML_IF(0,1,10,3126,3141,3202,3234);
        _SFD_CV_INIT_EML_IF(0,1,11,3202,3234,3316,3787);
        _SFD_CV_INIT_EML_IF(0,1,12,3316,3345,3440,3787);
        _SFD_CV_INIT_EML_IF(0,1,13,3440,3468,3567,3787);
        _SFD_CV_INIT_EML_IF(0,1,14,3567,3598,3697,3787);
        _SFD_CV_INIT_EML_IF(0,1,15,3697,3715,-1,3715);
        _SFD_CV_INIT_EML_IF(0,1,16,3861,3964,-1,4462);
        _SFD_CV_INIT_EML_IF(0,1,17,4304,4321,-1,4368);
        _SFD_CV_INIT_EML_IF(0,1,18,4482,4635,-1,5334);
        _SFD_CV_INIT_EML_IF(0,1,19,5175,5192,-1,5239);
        _SFD_CV_INIT_EML_IF(0,1,20,5355,5483,-1,6015);
        _SFD_CV_INIT_EML_IF(0,1,21,5853,5869,-1,5914);
        _SFD_CV_INIT_EML_IF(0,1,22,6036,6215,-1,6825);
        _SFD_CV_INIT_EML_IF(0,1,23,6666,6683,-1,6730);
        _SFD_CV_INIT_EML_IF(0,1,24,6847,6924,-1,7385);
        _SFD_CV_INIT_EML_IF(0,1,25,7226,7243,-1,7290);
        _SFD_CV_INIT_EML_FOR(0,1,0,4280,4292,4379);
        _SFD_CV_INIT_EML_FOR(0,1,1,5151,5163,5251);
        _SFD_CV_INIT_EML_FOR(0,1,2,5829,5841,5926);
        _SFD_CV_INIT_EML_FOR(0,1,3,6642,6654,6742);
        _SFD_CV_INIT_EML_FOR(0,1,4,7202,7214,7302);

        {
          static int condStart[] = { 2522, 2533 };

          static int condEnd[] = { 2529, 2542 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,0,2522,2542,2,0,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 2638, 2651 };

          static int condEnd[] = { 2647, 2656 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,1,2638,2656,2,2,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 2750, 2759 };

          static int condEnd[] = { 2755, 2767 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,2,2750,2767,2,4,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 2864, 2876 };

          static int condEnd[] = { 2872, 2882 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,3,2864,2882,2,6,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 3209, 3223 };

          static int condEnd[] = { 3219, 3233 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,4,3209,3233,2,8,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 3324, 3338 };

          static int condEnd[] = { 3334, 3344 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,5,3324,3344,2,10,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 3448, 3458 };

          static int condEnd[] = { 3454, 3467 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,6,3448,3467,2,12,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 3575, 3588 };

          static int condEnd[] = { 3584, 3597 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,7,3575,3597,2,14,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 3865, 3877, 3891, 3903, 3917, 3929, 3943,
            3955 };

          static int condEnd[] = { 3873, 3885, 3899, 3911, 3925, 3937, 3951,
            3963 };

          static int pfixExpr[] = { 0, 1, -3, 2, 3, -3, -2, 4, 5, -3, -2, 6, 7,
            -3, -2 };

          _SFD_CV_INIT_EML_MCDC(0,1,8,3864,3964,8,16,&(condStart[0]),&(condEnd[0]),
                                15,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 4486, 4498, 4512, 4524, 4538, 4550, 4564,
            4576, 4590, 4601, 4615, 4626 };

          static int condEnd[] = { 4494, 4506, 4520, 4532, 4546, 4558, 4572,
            4584, 4597, 4609, 4622, 4634 };

          static int pfixExpr[] = { 0, 1, -3, 2, 3, -3, -2, 4, 5, -3, -2, 6, 7,
            -3, -2, 8, 9, -3, -2, 10, 11, -3, -2 };

          _SFD_CV_INIT_EML_MCDC(0,1,9,4485,4635,12,24,&(condStart[0]),&(condEnd
            [0]),23,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 5359, 5371, 5385, 5397, 5411, 5422, 5436,
            5448, 5462, 5474 };

          static int condEnd[] = { 5367, 5379, 5393, 5405, 5418, 5430, 5444,
            5456, 5470, 5482 };

          static int pfixExpr[] = { 0, 1, -3, 2, 3, -3, -2, 4, 5, -3, -2, 6, 7,
            -3, -2, 8, 9, -3, -2 };

          _SFD_CV_INIT_EML_MCDC(0,1,10,5358,5483,10,36,&(condStart[0]),
                                &(condEnd[0]),19,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 6040, 6051, 6065, 6076, 6090, 6102, 6116,
            6128, 6142, 6154, 6168, 6180, 6194, 6206 };

          static int condEnd[] = { 6047, 6059, 6072, 6084, 6098, 6110, 6124,
            6136, 6150, 6162, 6176, 6188, 6202, 6214 };

          static int pfixExpr[] = { 0, 1, -3, 2, 3, -3, -2, 4, 5, -3, -2, 6, 7,
            -3, -2, 8, 9, -3, -2, 10, 11, -3, -2, 12, 13, -3, -2 };

          _SFD_CV_INIT_EML_MCDC(0,1,11,6039,6215,14,46,&(condStart[0]),
                                &(condEnd[0]),27,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 6851, 6863, 6877, 6889, 6903, 6915 };

          static int condEnd[] = { 6859, 6871, 6885, 6897, 6911, 6923 };

          static int pfixExpr[] = { 0, 1, -3, 2, 3, -3, -2, 4, 5, -3, -2 };

          _SFD_CV_INIT_EML_MCDC(0,1,12,6850,6924,6,60,&(condStart[0]),&(condEnd
            [0]),11,&(pfixExpr[0]));
        }

        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c5_sf_marshallOut,(MexInFcnForType)c5_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c5_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c5_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c5_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c5_sf_marshallOut,(MexInFcnForType)NULL);

        {
          real_T *c5_T;
          real_T *c5_teta;
          real_T *c5_desiredteta;
          real_T *c5_dteta;
          real_T *c5_desireddteta;
          c5_desireddteta = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
          c5_dteta = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c5_desiredteta = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c5_teta = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          c5_T = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
          _SFD_SET_DATA_VALUE_PTR(0U, c5_T);
          _SFD_SET_DATA_VALUE_PTR(1U, c5_teta);
          _SFD_SET_DATA_VALUE_PTR(2U, c5_desiredteta);
          _SFD_SET_DATA_VALUE_PTR(3U, c5_dteta);
          _SFD_SET_DATA_VALUE_PTR(4U, c5_desireddteta);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _PUMA_fuzzyMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "O6UEVGpRWbGFDcPWl4iRKF";
}

static void sf_opaque_initialize_c5_PUMA_fuzzy(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc5_PUMA_fuzzyInstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c5_PUMA_fuzzy((SFc5_PUMA_fuzzyInstanceStruct*)
    chartInstanceVar);
  initialize_c5_PUMA_fuzzy((SFc5_PUMA_fuzzyInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c5_PUMA_fuzzy(void *chartInstanceVar)
{
  enable_c5_PUMA_fuzzy((SFc5_PUMA_fuzzyInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c5_PUMA_fuzzy(void *chartInstanceVar)
{
  disable_c5_PUMA_fuzzy((SFc5_PUMA_fuzzyInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c5_PUMA_fuzzy(void *chartInstanceVar)
{
  sf_gateway_c5_PUMA_fuzzy((SFc5_PUMA_fuzzyInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c5_PUMA_fuzzy(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c5_PUMA_fuzzy
    ((SFc5_PUMA_fuzzyInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c5_PUMA_fuzzy();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c5_PUMA_fuzzy(SimStruct* S, const mxArray *
  st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c5_PUMA_fuzzy();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c5_PUMA_fuzzy((SFc5_PUMA_fuzzyInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c5_PUMA_fuzzy(SimStruct* S)
{
  return sf_internal_get_sim_state_c5_PUMA_fuzzy(S);
}

static void sf_opaque_set_sim_state_c5_PUMA_fuzzy(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c5_PUMA_fuzzy(S, st);
}

static void sf_opaque_terminate_c5_PUMA_fuzzy(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc5_PUMA_fuzzyInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_PUMA_fuzzy_optimization_info();
    }

    finalize_c5_PUMA_fuzzy((SFc5_PUMA_fuzzyInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc5_PUMA_fuzzy((SFc5_PUMA_fuzzyInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c5_PUMA_fuzzy(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    initialize_params_c5_PUMA_fuzzy((SFc5_PUMA_fuzzyInstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c5_PUMA_fuzzy(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_PUMA_fuzzy_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,5);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,5,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,5,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,5);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,5,4);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,5,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 4; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,5);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3231126885U));
  ssSetChecksum1(S,(725293514U));
  ssSetChecksum2(S,(1199111913U));
  ssSetChecksum3(S,(3087117645U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c5_PUMA_fuzzy(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c5_PUMA_fuzzy(SimStruct *S)
{
  SFc5_PUMA_fuzzyInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc5_PUMA_fuzzyInstanceStruct *)utMalloc(sizeof
    (SFc5_PUMA_fuzzyInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc5_PUMA_fuzzyInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c5_PUMA_fuzzy;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c5_PUMA_fuzzy;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c5_PUMA_fuzzy;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c5_PUMA_fuzzy;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c5_PUMA_fuzzy;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c5_PUMA_fuzzy;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c5_PUMA_fuzzy;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c5_PUMA_fuzzy;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c5_PUMA_fuzzy;
  chartInstance->chartInfo.mdlStart = mdlStart_c5_PUMA_fuzzy;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c5_PUMA_fuzzy;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  crtInfo->instanceInfo = (&(chartInstance->chartInfo));
  crtInfo->isJITEnabled = false;
  ssSetUserData(S,(void *)(crtInfo));  /* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c5_PUMA_fuzzy_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c5_PUMA_fuzzy(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c5_PUMA_fuzzy(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c5_PUMA_fuzzy(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c5_PUMA_fuzzy_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
