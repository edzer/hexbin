#include <R.h>
#include <R_ext/Rdynload.h>

extern void F77_NAME(hbin  )(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(herode)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hsm   )(void *, void *, void *, void *, void *, void *, void *);


static const R_FortranMethodDef FortranEntries[] = {
	{"hbin", (DL_FUNC) &F77_NAME(hbin), 13},
	{"herode", (DL_FUNC) &F77_NAME(herode), 10},
	{"hsm", (DL_FUNC) &F77_NAME(hsm), 7},
	{NULL, NULL, 0}
};

void R_init_hexbin(DllInfo *info) {
	R_registerRoutines(info, NULL, NULL, FortranEntries, NULL);
	R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);
}
