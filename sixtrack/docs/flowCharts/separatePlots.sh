#!/bin/bash

texFileName="createFlowCharts"

Files=`grep input ${texFileName}.tex | grep -v '%' | cut -d\{ -f2 | cut -d\} -f1`
Files=( ${Files} )
for (( ii=0; ii<${#Files[@]} ; ii++ )) ; do
    let jj=$ii+1
    pdftops ${texFileName}.pdf - | psselect -p${jj} | ps2pdf14 - ${Files[$ii]}.pdf
done
