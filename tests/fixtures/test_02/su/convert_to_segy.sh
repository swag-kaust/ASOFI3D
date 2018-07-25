#!/usr/bin/env bash
# Convert seismograms from the old Seismic Unix (.su) format
# to the SEG-Y format.
# The source of the seismograms is the `overnightbuilt` branch
# of the SOFI3D code.

print_error () {
    # Print an error message to the standard error output.
    echo ERROR: "$1" > /dev/stderr
}

if [ $# -ge 1 ]; then
    filename=$1
else
    echo 'USAGE: convert_to_segy_format.sh filename_prefix'
    echo 'For example: convert_to_segy_format.sh test_p'
    exit 1
fi

# Convert from the old Seismic Unix format to the new one.
suoldtonew < "${filename}.su" > "${filename}_new.su"
if [ $? -ne 0 ]; then
    print_error "Conversion from old Seismic Unix format to the new one failed"
fi

segyhdrs < "${filename}_new.su"

segywrite < "${filename}_new.su" tape="${filename}.sgy"
if [ $? -ne 0 ]; then
    print_error "Conversion from new Seismic Unix format to SEG-Y failed"
fi
