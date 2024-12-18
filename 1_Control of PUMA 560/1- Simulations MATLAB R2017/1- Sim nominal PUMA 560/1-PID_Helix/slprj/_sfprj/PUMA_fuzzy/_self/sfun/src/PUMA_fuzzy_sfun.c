/* Include files */

#include "PUMA_fuzzy_sfun.h"
#include "PUMA_fuzzy_sfun_debug_macros.h"
#include "c1_PUMA_fuzzy.h"
#include "c2_PUMA_fuzzy.h"
#include "c3_PUMA_fuzzy.h"
#include "c4_PUMA_fuzzy.h"
#include "c5_PUMA_fuzzy.h"
#include "c6_PUMA_fuzzy.h"
#include "c7_PUMA_fuzzy.h"
#include "c8_PUMA_fuzzy.h"
#include "c9_PUMA_fuzzy.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
uint32_T _PUMA_fuzzyMachineNumber_;

/* Function Declarations */

/* Function Definitions */
void PUMA_fuzzy_initializer(void)
{
}

void PUMA_fuzzy_terminator(void)
{
}

/* SFunction Glue Code */
unsigned int sf_PUMA_fuzzy_method_dispatcher(SimStruct *simstructPtr, unsigned
  int chartFileNumber, const char* specsCksum, int_T method, void *data)
{
  if (chartFileNumber==1) {
    c1_PUMA_fuzzy_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==2) {
    c2_PUMA_fuzzy_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==3) {
    c3_PUMA_fuzzy_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==4) {
    c4_PUMA_fuzzy_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==5) {
    c5_PUMA_fuzzy_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==6) {
    c6_PUMA_fuzzy_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==7) {
    c7_PUMA_fuzzy_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==8) {
    c8_PUMA_fuzzy_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==9) {
    c9_PUMA_fuzzy_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  return 0;
}

unsigned int sf_PUMA_fuzzy_process_check_sum_call( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[20];
  if (nrhs<1 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the checksum */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"sf_get_check_sum"))
    return 0;
  plhs[0] = mxCreateDoubleMatrix( 1,4,mxREAL);
  if (nrhs>1 && mxIsChar(prhs[1])) {
    mxGetString(prhs[1], commandName,sizeof(commandName)/sizeof(char));
    commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
    if (!strcmp(commandName,"machine")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(501110421U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1517896906U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2715824327U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(105509480U);
    } else if (!strcmp(commandName,"exportedFcn")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0U);
    } else if (!strcmp(commandName,"makefile")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(674415114U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2426248082U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1000325487U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3683618400U);
    } else if (nrhs==3 && !strcmp(commandName,"chart")) {
      unsigned int chartFileNumber;
      chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
      switch (chartFileNumber) {
       case 1:
        {
          extern void sf_c1_PUMA_fuzzy_get_check_sum(mxArray *plhs[]);
          sf_c1_PUMA_fuzzy_get_check_sum(plhs);
          break;
        }

       case 2:
        {
          extern void sf_c2_PUMA_fuzzy_get_check_sum(mxArray *plhs[]);
          sf_c2_PUMA_fuzzy_get_check_sum(plhs);
          break;
        }

       case 3:
        {
          extern void sf_c3_PUMA_fuzzy_get_check_sum(mxArray *plhs[]);
          sf_c3_PUMA_fuzzy_get_check_sum(plhs);
          break;
        }

       case 4:
        {
          extern void sf_c4_PUMA_fuzzy_get_check_sum(mxArray *plhs[]);
          sf_c4_PUMA_fuzzy_get_check_sum(plhs);
          break;
        }

       case 5:
        {
          extern void sf_c5_PUMA_fuzzy_get_check_sum(mxArray *plhs[]);
          sf_c5_PUMA_fuzzy_get_check_sum(plhs);
          break;
        }

       case 6:
        {
          extern void sf_c6_PUMA_fuzzy_get_check_sum(mxArray *plhs[]);
          sf_c6_PUMA_fuzzy_get_check_sum(plhs);
          break;
        }

       case 7:
        {
          extern void sf_c7_PUMA_fuzzy_get_check_sum(mxArray *plhs[]);
          sf_c7_PUMA_fuzzy_get_check_sum(plhs);
          break;
        }

       case 8:
        {
          extern void sf_c8_PUMA_fuzzy_get_check_sum(mxArray *plhs[]);
          sf_c8_PUMA_fuzzy_get_check_sum(plhs);
          break;
        }

       case 9:
        {
          extern void sf_c9_PUMA_fuzzy_get_check_sum(mxArray *plhs[]);
          sf_c9_PUMA_fuzzy_get_check_sum(plhs);
          break;
        }

       default:
        ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0.0);
      }
    } else if (!strcmp(commandName,"target")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3031367619U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(4001028638U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3978939492U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(838979348U);
    } else {
      return 0;
    }
  } else {
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3807072667U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3872474013U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2817460162U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1923074366U);
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_PUMA_fuzzy_autoinheritance_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[32];
  char aiChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the autoinheritance_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_autoinheritance_info"))
    return 0;
  mxGetString(prhs[2], aiChksum,sizeof(aiChksum)/sizeof(char));
  aiChksum[(sizeof(aiChksum)/sizeof(char)-1)] = '\0';

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(aiChksum, "CM92fGLI0YW4eEB0eOFczG") == 0) {
          extern mxArray *sf_c1_PUMA_fuzzy_get_autoinheritance_info(void);
          plhs[0] = sf_c1_PUMA_fuzzy_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 2:
      {
        if (strcmp(aiChksum, "JhhfrIDgI3GLImbeYr6vY") == 0) {
          extern mxArray *sf_c2_PUMA_fuzzy_get_autoinheritance_info(void);
          plhs[0] = sf_c2_PUMA_fuzzy_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 3:
      {
        if (strcmp(aiChksum, "IldpWiAaXSlrDqOr8XalaH") == 0) {
          extern mxArray *sf_c3_PUMA_fuzzy_get_autoinheritance_info(void);
          plhs[0] = sf_c3_PUMA_fuzzy_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 4:
      {
        if (strcmp(aiChksum, "IldpWiAaXSlrDqOr8XalaH") == 0) {
          extern mxArray *sf_c4_PUMA_fuzzy_get_autoinheritance_info(void);
          plhs[0] = sf_c4_PUMA_fuzzy_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 5:
      {
        if (strcmp(aiChksum, "IldpWiAaXSlrDqOr8XalaH") == 0) {
          extern mxArray *sf_c5_PUMA_fuzzy_get_autoinheritance_info(void);
          plhs[0] = sf_c5_PUMA_fuzzy_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 6:
      {
        if (strcmp(aiChksum, "12kFc1ERD0HZlaUHP7SUP") == 0) {
          extern mxArray *sf_c6_PUMA_fuzzy_get_autoinheritance_info(void);
          plhs[0] = sf_c6_PUMA_fuzzy_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 7:
      {
        if (strcmp(aiChksum, "NkUVpTvpcIfnQ3dGHoUIiB") == 0) {
          extern mxArray *sf_c7_PUMA_fuzzy_get_autoinheritance_info(void);
          plhs[0] = sf_c7_PUMA_fuzzy_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 8:
      {
        if (strcmp(aiChksum, "NkUVpTvpcIfnQ3dGHoUIiB") == 0) {
          extern mxArray *sf_c8_PUMA_fuzzy_get_autoinheritance_info(void);
          plhs[0] = sf_c8_PUMA_fuzzy_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 9:
      {
        if (strcmp(aiChksum, "NkUVpTvpcIfnQ3dGHoUIiB") == 0) {
          extern mxArray *sf_c9_PUMA_fuzzy_get_autoinheritance_info(void);
          plhs[0] = sf_c9_PUMA_fuzzy_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_PUMA_fuzzy_get_eml_resolved_functions_info( int nlhs, mxArray *
  plhs[], int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[64];
  if (nrhs<2 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the get_eml_resolved_functions_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_eml_resolved_functions_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        extern const mxArray *sf_c1_PUMA_fuzzy_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c1_PUMA_fuzzy_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 2:
      {
        extern const mxArray *sf_c2_PUMA_fuzzy_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c2_PUMA_fuzzy_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 3:
      {
        extern const mxArray *sf_c3_PUMA_fuzzy_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c3_PUMA_fuzzy_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 4:
      {
        extern const mxArray *sf_c4_PUMA_fuzzy_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c4_PUMA_fuzzy_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 5:
      {
        extern const mxArray *sf_c5_PUMA_fuzzy_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c5_PUMA_fuzzy_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 6:
      {
        extern const mxArray *sf_c6_PUMA_fuzzy_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c6_PUMA_fuzzy_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 7:
      {
        extern const mxArray *sf_c7_PUMA_fuzzy_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c7_PUMA_fuzzy_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 8:
      {
        extern const mxArray *sf_c8_PUMA_fuzzy_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c8_PUMA_fuzzy_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 9:
      {
        extern const mxArray *sf_c9_PUMA_fuzzy_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c9_PUMA_fuzzy_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_PUMA_fuzzy_third_party_uses_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the third_party_uses_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_third_party_uses_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "aAm9ijpDHDAauwqf1ezeEH") == 0) {
          extern mxArray *sf_c1_PUMA_fuzzy_third_party_uses_info(void);
          plhs[0] = sf_c1_PUMA_fuzzy_third_party_uses_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "ooUbnnZvOjLciQTjlVYnMC") == 0) {
          extern mxArray *sf_c2_PUMA_fuzzy_third_party_uses_info(void);
          plhs[0] = sf_c2_PUMA_fuzzy_third_party_uses_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "O6UEVGpRWbGFDcPWl4iRKF") == 0) {
          extern mxArray *sf_c3_PUMA_fuzzy_third_party_uses_info(void);
          plhs[0] = sf_c3_PUMA_fuzzy_third_party_uses_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "O6UEVGpRWbGFDcPWl4iRKF") == 0) {
          extern mxArray *sf_c4_PUMA_fuzzy_third_party_uses_info(void);
          plhs[0] = sf_c4_PUMA_fuzzy_third_party_uses_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "O6UEVGpRWbGFDcPWl4iRKF") == 0) {
          extern mxArray *sf_c5_PUMA_fuzzy_third_party_uses_info(void);
          plhs[0] = sf_c5_PUMA_fuzzy_third_party_uses_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "zogkxKwma5JxPXRH1HG7YB") == 0) {
          extern mxArray *sf_c6_PUMA_fuzzy_third_party_uses_info(void);
          plhs[0] = sf_c6_PUMA_fuzzy_third_party_uses_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "fWf8KhLhbwTrgZwrH8S6GF") == 0) {
          extern mxArray *sf_c7_PUMA_fuzzy_third_party_uses_info(void);
          plhs[0] = sf_c7_PUMA_fuzzy_third_party_uses_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "fWf8KhLhbwTrgZwrH8S6GF") == 0) {
          extern mxArray *sf_c8_PUMA_fuzzy_third_party_uses_info(void);
          plhs[0] = sf_c8_PUMA_fuzzy_third_party_uses_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "fWf8KhLhbwTrgZwrH8S6GF") == 0) {
          extern mxArray *sf_c9_PUMA_fuzzy_third_party_uses_info(void);
          plhs[0] = sf_c9_PUMA_fuzzy_third_party_uses_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

unsigned int sf_PUMA_fuzzy_updateBuildInfo_args_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the updateBuildInfo_args_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_updateBuildInfo_args_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "aAm9ijpDHDAauwqf1ezeEH") == 0) {
          extern mxArray *sf_c1_PUMA_fuzzy_updateBuildInfo_args_info(void);
          plhs[0] = sf_c1_PUMA_fuzzy_updateBuildInfo_args_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "ooUbnnZvOjLciQTjlVYnMC") == 0) {
          extern mxArray *sf_c2_PUMA_fuzzy_updateBuildInfo_args_info(void);
          plhs[0] = sf_c2_PUMA_fuzzy_updateBuildInfo_args_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "O6UEVGpRWbGFDcPWl4iRKF") == 0) {
          extern mxArray *sf_c3_PUMA_fuzzy_updateBuildInfo_args_info(void);
          plhs[0] = sf_c3_PUMA_fuzzy_updateBuildInfo_args_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "O6UEVGpRWbGFDcPWl4iRKF") == 0) {
          extern mxArray *sf_c4_PUMA_fuzzy_updateBuildInfo_args_info(void);
          plhs[0] = sf_c4_PUMA_fuzzy_updateBuildInfo_args_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "O6UEVGpRWbGFDcPWl4iRKF") == 0) {
          extern mxArray *sf_c5_PUMA_fuzzy_updateBuildInfo_args_info(void);
          plhs[0] = sf_c5_PUMA_fuzzy_updateBuildInfo_args_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "zogkxKwma5JxPXRH1HG7YB") == 0) {
          extern mxArray *sf_c6_PUMA_fuzzy_updateBuildInfo_args_info(void);
          plhs[0] = sf_c6_PUMA_fuzzy_updateBuildInfo_args_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "fWf8KhLhbwTrgZwrH8S6GF") == 0) {
          extern mxArray *sf_c7_PUMA_fuzzy_updateBuildInfo_args_info(void);
          plhs[0] = sf_c7_PUMA_fuzzy_updateBuildInfo_args_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "fWf8KhLhbwTrgZwrH8S6GF") == 0) {
          extern mxArray *sf_c8_PUMA_fuzzy_updateBuildInfo_args_info(void);
          plhs[0] = sf_c8_PUMA_fuzzy_updateBuildInfo_args_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "fWf8KhLhbwTrgZwrH8S6GF") == 0) {
          extern mxArray *sf_c9_PUMA_fuzzy_updateBuildInfo_args_info(void);
          plhs[0] = sf_c9_PUMA_fuzzy_updateBuildInfo_args_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

void PUMA_fuzzy_debug_initialize(struct SfDebugInstanceStruct* debugInstance)
{
  _PUMA_fuzzyMachineNumber_ = sf_debug_initialize_machine(debugInstance,
    "PUMA_fuzzy","sfun",0,9,0,0,0);
  sf_debug_set_machine_event_thresholds(debugInstance,_PUMA_fuzzyMachineNumber_,
    0,0);
  sf_debug_set_machine_data_thresholds(debugInstance,_PUMA_fuzzyMachineNumber_,0);
}

void PUMA_fuzzy_register_exported_symbols(SimStruct* S)
{
}

static mxArray* sRtwOptimizationInfoStruct= NULL;
mxArray* load_PUMA_fuzzy_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct==NULL) {
    sRtwOptimizationInfoStruct = sf_load_rtw_optimization_info("PUMA_fuzzy",
      "PUMA_fuzzy");
    mexMakeArrayPersistent(sRtwOptimizationInfoStruct);
  }

  return(sRtwOptimizationInfoStruct);
}

void unload_PUMA_fuzzy_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct!=NULL) {
    mxDestroyArray(sRtwOptimizationInfoStruct);
    sRtwOptimizationInfoStruct = NULL;
  }
}
