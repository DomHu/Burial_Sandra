#! /bin/bash

#first we convert all ps to eps
ps2eps *.ps

##now we get a list of all eps files and call eps as many times as required
FILESeps="*.eps"
for f in $FILESeps; do epstopdf $f; done


##and now we crop the pdfs in case they are not cropped already

FILES="*.pdf"
for f in $FILES; do pdfcrop $f $f; done


##and now we get rid of the eps files and ps files
mkdir eps-files
mv *.eps eps-files/.
rm *.ps

#for f in $FILESeps; do rm $f; done
