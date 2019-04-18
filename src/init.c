/** 
 * @file  init.c
 * @brief init file for all C files
 *
 * @author Christophe Dutang
 *
 * Copyright (C) 2018, Christophe Dutang. 
 * All rights reserved.
 *  
 */
 
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "header.h"

 
 //table of registration routines accessed with .C() all defined in header.h
 static const R_CMethodDef CEntries[] = {
   {"do_ts2tsd", (DL_FUNC) &do_ts2tsd, 8}, /*called in util-ts2tsd.R*/
   {"do_gpdlik", (DL_FUNC) &do_gpdlik, 6}, /*called in util-confGPD-prof.R, util-GPD-(g)mle.R*/
   {"do_pplik", (DL_FUNC) &do_pplik, 8}, /*called in uvpot-fitPP.R*/
   {"do_samlmu", (DL_FUNC) &do_samlmu, 4}, /*called in util-lmom.R*/
   {"do_clust", (DL_FUNC) &do_clust, 6}, /*called in graph-clust.R*/
   {"do_cgpdlik", (DL_FUNC) &do_cgpdlik, 8}, 
   {"do_gpdbvlog", (DL_FUNC) &do_gpdbvlog, 13}, /*called in bvpot-fitBvGPD.R*/
   {"do_gpdbvalog", (DL_FUNC) &do_gpdbvalog, 15}, /*called in bvpot-fitBvGPD.R*/
   {"do_gpdbvnlog", (DL_FUNC) &do_gpdbvnlog, 13}, /*called in bvpot-fitBvGPD.R*/
   {"do_gpdbvanlog", (DL_FUNC) &do_gpdbvanlog, 15}, /*called in bvpot-fitBvGPD.R*/
   {"do_gpdbvmix", (DL_FUNC) &do_gpdbvmix, 13}, /*called in bvpot-fitBvGPD.R*/
   {"do_gpdbvamix", (DL_FUNC) &do_gpdbvamix, 14}, /*called in bvpot-fitBvGPD.R*/
   {"do_rbvlog_shi", (DL_FUNC) &do_rbvlog_shi, 3}, /*called in bvGPD-dpqr-fun.R*/
   {"do_rbvalog_shi", (DL_FUNC) &do_rbvalog_shi, 4}, /*called in bvGPD-dpqr-fun.R*/
   {"do_rbvlog", (DL_FUNC) &do_rbvlog, 3}, /*called in bvGPD-dpqr-fun.R, mcpot-simmc.R*/
   {"do_rbvalog", (DL_FUNC) &do_rbvalog, 4}, /*called in bvGPD-dpqr-fun.R, mcpot-simmc.R*/
   {"do_rbvnlog", (DL_FUNC) &do_rbvnlog, 3}, /*called in bvGPD-dpqr-fun.R, mcpot-simmc.R*/
   {"do_rbvanlog", (DL_FUNC) &do_rbvanlog, 4}, /*called in bvGPD-dpqr-fun.R, mcpot-simmc.R*/
   {"do_rbvmix", (DL_FUNC) &do_rbvmix, 3}, /*called in bvGPD-dpqr-fun.R, mcpot-simmc.R*/
   {"do_rbvamix", (DL_FUNC) &do_rbvamix, 4}, /*called in bvGPD-dpqr-fun.R, mcpot-simmc.R*/
   {"do_gpdmclog", (DL_FUNC) &do_gpdmclog, 13}, /*called in mcpot-fitMcGPD.R*/
   {"do_gpdmcalog", (DL_FUNC) &do_gpdmcalog, 15}, /*called in mcpot-fitMcGPD.R*/
   {"do_gpdmcnlog", (DL_FUNC) &do_gpdmcnlog, 13}, /*called in mcpot-fitMcGPD.R*/
   {"do_gpdmcanlog", (DL_FUNC) &do_gpdmcanlog, 15}, /*called in mcpot-fitMcGPD.R*/
   {"do_gpdmcmix", (DL_FUNC) &do_gpdmcmix, 13}, /*called in mcpot-fitMcGPD.R*/
   {"do_gpdmcamix", (DL_FUNC) &do_gpdmcamix, 14}, /*called in mcpot-fitMcGPD.R*/
   /* never defined by M. Ribatet in C code
    * {"do_gpdmcamixtest", (DL_FUNC) &do_gpdmcamixtest, 14}, called in mcpot-fitMcGPD.R*/
   {NULL, NULL, 0}
 };
 
 
 
 
 
 //table of registered routines
 void R_init_POT(DllInfo *info)
 {
   //register method accessed with .C(), no Call(), no Fortran(), no external()
   R_registerRoutines(info, CEntries, NULL, NULL, NULL); 
   //avoid dynamic lookup for routines
   R_useDynamicSymbols(info, FALSE);
   //force the use of symbols and not character strings
   R_forceSymbols(info, TRUE);
 }
 
