#!/bin/sh
#
# Usage: eps2png.sh name.eps
#  Writes name.png and name-crop.png to cwd (WARNING overwrites)
# Using ImageMagick
#
# ImageMagick convert no longer seems to work on either cluster (linux)
# or cygwin (on new Lenovo PC) so use gs directly instead
#
# Note 'crop' e.g. pdfcrop is called 'trim' in ImageMagick

if [ $# -ne 1 ]; then
   echo "Usage: $0 name.eps" >&2
   exit 1
fi
name=$1


namebase=`basename ${name} .eps`

gs -dQUIET -dSAFER -dBATCH -dNOPAUSE -dNOPROMPT -dMaxBitmap=500000000 -dAlignToPixels=0 -dGridFitTT=2 '-sDEVICE=pngalpha' -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -r300 -dEPSCrop -sOutputFile=${namebase}.png -f ${name}

convert -trim +repage ${namebase}.png ${namebase}-crop.png
