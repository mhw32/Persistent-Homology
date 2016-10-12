#!/usr/bin/env bash

echo "Generating Figure 1"
Rscript gen_fig1.r 
echo "Generating Figure 2"
python gen_fig2.py
Rscript gen_fig2.r
echo "Generating Figure 3"
python gen_fig3.py
echo "Generating Figure 5"
python gen_fig5.py
Rscript gen_fig5.r
echo "Generating Figure 7"
Rscript gen_fig7.r
echo "Generating Figure 8"
python gen_fig8.py
echo "Generating Figure 9"
python gen_fig9.py
echo "Generating Figure 10"
Rscript gen_fig10.r
echo "Generating Figure 11"
Rscript gen_fig11.r
echo "Generating Figure 12"
python gen_fig_12.py
echo "Generating Figure 13"
python gen_fig_13.py
Rscript gen_fig_13.r
echo "Generating Figure A1"
python gen_fig_appendix1.py
echo "Generating Figure A2"
python gen_fig_appendix2.py

echo "Saving files to output/"
rm output/*.pdf
mv *.pdf output
