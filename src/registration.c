#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void odin_model_2sp_initmod_desolve(void *);
extern void odin_model_2sp_output_dde(void *);
extern void odin_model_2sp_rhs_dde(void *);
extern void odin_model_2sp_rhs_desolve(void *);

/* .Call calls */
extern SEXP odin_model_2sp_contents(void *);
extern SEXP odin_model_2sp_create(void *);
extern SEXP odin_model_2sp_initial_conditions(void *, void *);
extern SEXP odin_model_2sp_metadata(void *);
extern SEXP odin_model_2sp_rhs_r(void *, void *, void *);
extern SEXP odin_model_2sp_set_initial(void *, void *, void *, void *);
extern SEXP odin_model_2sp_set_user(void *, void *);

static const R_CMethodDef CEntries[] = {
    {"odin_model_2sp_initmod_desolve", (DL_FUNC) &odin_model_2sp_initmod_desolve, 1},
    {"odin_model_2sp_output_dde",      (DL_FUNC) &odin_model_2sp_output_dde,      1},
    {"odin_model_2sp_rhs_dde",         (DL_FUNC) &odin_model_2sp_rhs_dde,         1},
    {"odin_model_2sp_rhs_desolve",     (DL_FUNC) &odin_model_2sp_rhs_desolve,     1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"odin_model_2sp_contents",           (DL_FUNC) &odin_model_2sp_contents,           1},
    {"odin_model_2sp_create",             (DL_FUNC) &odin_model_2sp_create,             1},
    {"odin_model_2sp_initial_conditions", (DL_FUNC) &odin_model_2sp_initial_conditions, 2},
    {"odin_model_2sp_metadata",           (DL_FUNC) &odin_model_2sp_metadata,           1},
    {"odin_model_2sp_rhs_r",              (DL_FUNC) &odin_model_2sp_rhs_r,              3},
    {"odin_model_2sp_set_initial",        (DL_FUNC) &odin_model_2sp_set_initial,        4},
    {"odin_model_2sp_set_user",           (DL_FUNC) &odin_model_2sp_set_user,           2},
    {NULL, NULL, 0}
};

void R_init_ICDMM(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
