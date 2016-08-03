#!/bin/bash
FILEBASE=$1_
FILEEXT1=vx.su
FILEEXT2=vy.su
FILEEXT3=vz.su
FILEEXT4=p.su
FILEEXT5=div.su
FILEEXT6=curl.su
FILESTAR=.*

FOLDER_IN=su/
FOLDER_OUT=su/

#usage : bash$: ./sucat.sh test
#all following clauses differ only by the chosen file extension $FILEEXT#
#file names are constructed by adding FILEBASE+FILEEXT+FILESTAR

#merge of *.vx seismogram files
#first number of files with wildcard are counted
NUMOFFFILES=$(ls $FOLDER_IN$FILEBASE$FILEEXT1$FILESTAR 2> /dev/null | wc -l)
#if number is larger than 0, merge will be performed into $FOLDER
if [ "$NUMOFFFILES" != "0" ]
then
   echo "File $FOLDER_IN$FILEBASE$FILEEXT1$FILESTAR exists"
   cat $FOLDER_IN$FILEBASE$FILEEXT1$FILESTAR > $FOLDER_OUT$FILEBASE$FILEEXT1
else
   echo "File $FOLDER_IN$FILEBASE$FILEEXT1$FILESTAR does not exists"
fi

#merge of *.vy seismogram files
NUMOFFFILES=$(ls $FOLDER_IN$FILEBASE$FILEEXT2$FILESTAR 2> /dev/null | wc -l)
if [ "$NUMOFFFILES" != "0" ]
then
   echo "File $FOLDER_IN$FILEBASE$FILEEXT2$FILESTAR exists"
   cat $FOLDER_IN$FILEBASE$FILEEXT2$FILESTAR > $FOLDER_OUT$FILEBASE$FILEEXT2
else
   echo "File $FOLDER_IN$FILEBASE$FILEEXT2$FILESTAR does not exists"
fi

#merge of *.vz seismogram files
NUMOFFFILES=$(ls $FOLDER_IN$FILEBASE$FILEEXT3$FILESTAR 2> /dev/null | wc -l)
if [ "$NUMOFFFILES" != "0" ]
then
   echo "File $FOLDER_IN$FILEBASE$FILEEXT3$FILESTAR exists"
   cat $FOLDER_IN$FILEBASE$FILEEXT3$FILESTAR > $FOLDER_OUT$FILEBASE$FILEEXT3
else
   echo "File $FOLDER_IN$FILEBASE$FILEEXT3$FILESTAR does not exists"
fi


#merge of *.p seismogram files
NUMOFFFILES=$(ls $FOLDER_IN$FILEBASE$FILEEXT4$FILESTAR 2> /dev/null | wc -l)
if [ "$NUMOFFFILES" != "0" ]
then
   echo "File $FOLDER_IN$FILEBASE$FILEEXT4$FILESTAR exists"
   cat $FOLDER_IN$FILEBASE$FILEEXT4$FILESTAR > $FOLDER_OUT$FILEBASE$FILEEXT4
else
   echo "File $FOLDER_IN$FILEBASE$FILEEXT4$FILESTAR does not exists"
fi

#merge of *.div seismogram files
NUMOFFFILES=$(ls $FOLDER_IN$FILEBASE$FILEEXT5$FILESTAR 2> /dev/null | wc -l)
if [ "$NUMOFFFILES" != "0" ]
then
   echo "File $FOLDER_IN$FILEBASE$FILEEXT5$FILESTAR exists"
   cat $FOLDER_IN$FILEBASE$FILEEXT5$FILESTAR > $FOLDER_OUT$FILEBASE$FILEEXT5
else
   echo "File $FOLDER_IN$FILEBASE$FILEEXT3$FILESTAR does not exists"
fi

#merge of *.rot seismogram files
NUMOFFFILES=$(ls $FOLDER_IN$FILEBASE$FILEEXT6$FILESTAR 2> /dev/null | wc -l)
if [ "$NUMOFFFILES" != "0" ]
then
   echo "File $FOLDER_IN$FILEBASE$FILEEXT6$FILESTAR exists"
   cat $FOLDER_IN$FILEBASE$FILEEXT6$FILESTAR > $FOLDER_OUT$FILEBASE$FILEEXT6
else
   echo "File $FOLDER_IN$FILEBASE$FILEEXT6$FILESTAR does not exists"
fi
