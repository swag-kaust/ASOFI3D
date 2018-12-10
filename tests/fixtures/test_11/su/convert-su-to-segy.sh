#!/usr/bin/env bash
# Convert seismograms from the old Seismic Unix (.su) format
# to the SEG-Y format.
# The source of the seismograms is SAVA code.

fatal () {
    # Print an error message to the standard error output and exit.
    echo ERROR: "$1" > /dev/stderr
    exit 1
}

if [ $# -ge 1 ]; then
    filename_su_old=$1
else
    echo 'USAGE: convert_to_segy_format.sh filename_in_old_su_format'
    echo 'For example: convert_to_segy_format.sh fd3d_ortho0001_p.su'
    exit 2
fi

filename=${filename_su_old%.*}
filename_su=${filename}.su.new
filename_sgy=${filename}.sgy

# Convert from the old Seismic Unix format to the new one.
#suoldtonew < "$filename_su_old" > "${filename_su}"
cp "$filename_su_old" "${filename_su}"
if [ $? -ne 0 ]; then
    fatal "Conversion from old Seismic Unix format to the new one failed"
fi

segyhdrs < "${filename_su}" hfile=header bfile=binary

segywrite < "${filename_su}" hfile=header bfile=binary tape="${filename_sgy}"
if [ $? -ne 0 ]; then
    fatal "Conversion from new Seismic Unix format to SEG-Y failed"
fi
