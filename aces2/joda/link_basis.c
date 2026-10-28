/*
 * When the job directory has neither GENBAS nor ZMAT.BAS, link GENBAS from
 * the basis set library into it (and ECPDATA too, when the job directory has
 * none), so that every program run after joda finds the files as before.
 * The library is the first of these directories that holds a GENBAS:
 *
 *   $ACES_BASIS_DIR
 *   ACES_BASIS_DIR, the ACES_PROGRAM/basis of the tree joda was built from
 *   <directory of this executable>/../share/aces/basis
 *
 * A GENBAS, ZMAT.BAS or ECPDATA already in the job directory always wins.
 * If no library is found nothing happens, and joda reports the missing
 * GENBAS as it always has.
 */

#include <stdio.h>	/* for printf */
#include <stdlib.h>	/* for getenv */
#include <string.h>	/* for strrchr */
#include <unistd.h>	/* for access, readlink, symlink */
#include <limits.h>	/* for PATH_MAX */

#include "f77_name.h"	/* for F77_NAME */

#ifndef ACES_BASIS_DIR
#define ACES_BASIS_DIR ""
#endif

static int
has_genbas(const char * szDir)
{
    char szFile[PATH_MAX];
    if (!szDir || !*szDir) return 0;
    snprintf(szFile,PATH_MAX,"%s/GENBAS",szDir);
    return access(szFile,R_OK) == 0;
}

static void
link_file(const char * szDir, const char * szName)
{
    char szFile[PATH_MAX];
    snprintf(szFile,PATH_MAX,"%s/%s",szDir,szName);
    if (access(szName,F_OK) == 0 || access(szFile,R_OK)) return;
    if (symlink(szFile,szName))
        printf("@LINK_BASIS: could not link %s to %s\n",szName,szFile);
    else
        printf("@LINK_BASIS: no %s in the job directory; using %s\n",
               szName,szFile);
}

void
F77_NAME(link_basis,LINK_BASIS)
(void)
{
    char szExe[PATH_MAX], szShare[PATH_MAX];
    const char * szDir = NULL;
    ssize_t n;
    char * p;

    if (access("GENBAS",F_OK) == 0 || access("ZMAT.BAS",F_OK) == 0) return;

    if (has_genbas(getenv("ACES_BASIS_DIR")))
        szDir = getenv("ACES_BASIS_DIR");
    else if (has_genbas(ACES_BASIS_DIR))
        szDir = ACES_BASIS_DIR;
    else
    {
        n = readlink("/proc/self/exe",szExe,PATH_MAX-1);
        if (n > 0)
        {
            szExe[n] = '\0';
            p = strrchr(szExe,'/');
            if (p) *p = '\0';
            snprintf(szShare,PATH_MAX,"%s/../share/aces/basis",szExe);
            if (has_genbas(szShare)) szDir = szShare;
        }
    }
    if (!szDir) return;

    link_file(szDir,"GENBAS");
    link_file(szDir,"ECPDATA");
}
