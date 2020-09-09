#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "rsBayes.h"

static const R_CallMethodDef CallEntries[] = {
    {"phenoMod0", (DL_FUNC) &phenoMod0, 17},
    {"phenoMod0Pred", (DL_FUNC) &phenoMod0, 11},
    {NULL, NULL, 0}
};

void R_init_rsBayes(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

