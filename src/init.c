#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void EntropySd(void *, void *, void *);
extern void EntropySharp(void *, void *, void *);
extern void GenSimpSd(void *, void *, void *, void *);
extern void GenSimpSharp(void *, void *, void *, void *);
extern void KlPlugin(void *, void *, void *, void *);
extern void KlSd(void *, void *, void *, void *);
extern void KlSharp(void *, void *, void *, void *);
extern void MISd(void *, void *, void *, void *);
extern void RenyiEqEntropySharp(void *, void *, void *, void *);
extern void RenyiEqSd(void *, void *, void *, void *);
extern void SymSd(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"EntropySd",           (DL_FUNC) &EntropySd,           3},
    {"EntropySharp",        (DL_FUNC) &EntropySharp,        3},
    {"GenSimpSd",           (DL_FUNC) &GenSimpSd,           4},
    {"GenSimpSharp",        (DL_FUNC) &GenSimpSharp,        4},
    {"KlPlugin",            (DL_FUNC) &KlPlugin,            4},
    {"KlSd",                (DL_FUNC) &KlSd,                4},
    {"KlSharp",             (DL_FUNC) &KlSharp,             4},
    {"MISd",                (DL_FUNC) &MISd,                4},
    {"RenyiEqEntropySharp", (DL_FUNC) &RenyiEqEntropySharp, 4},
    {"RenyiEqSd",           (DL_FUNC) &RenyiEqSd,           4},
    {"SymSd",               (DL_FUNC) &SymSd,               4},
    {NULL, NULL, 0}
};

void R_init_EntropyEstimation(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
