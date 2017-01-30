#!/bin/bash
dir=`pwd`
str="HOMEDIR := $dir"
sed -i "/HOMEDIR := /c $str" environment

str="CHARACTER(*), PARAMETER :: resultsDir=\"$1/\""
sed -i "/resultsDir=/c $str" tests/testDAMVW.f95
