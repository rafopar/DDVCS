#/usr/bin/sh

echo "Kuku compiling latex"

echo "* Run = $1"

#cmd="pdflatex -jobname=AnalysisPlots_Run_$1_tOvrThr_$2_tCross_$3_tLRMatch_$4.pdf \"\def\RUN{$1} \def\T_OVER_THR{$2} \def\CRS_DELTA_T{$3} \def\LR_MATCH_DELTA_T{$4} \input{XYHodoAnalysisPlots.tex}\""
pdflatex -jobname=AnaGEMCPlots_Run_$1 "\def\RUN{$1} \input{AnaGEMC.tex}"

