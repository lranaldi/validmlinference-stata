*! _build.do -- compile lvmli.mlib from .mata source files
*! Run from ValidMLInference-stata/src/mata/ :  do _build.do

version 17

qui cap mata: mata mlib drop lvmli
qui cap erase lvmli.mlib

do vmli_core.mata

mata: mata mlib create lvmli, replace
mata: mata mlib add lvmli vmli_*() _vmli_*()
mata: mata mlib index

di as result "Built lvmli.mlib in `c(pwd)'"
