# Functions related to testing.

setup() {
    # Setup environment for the execution of a test.
    rm -rf tmp/
    mkdir tmp
    mkdir tmp/in_and_out
    mkdir tmp/log
    mkdir tmp/model
    mkdir tmp/snap
    mkdir tmp/sources
    mkdir tmp/receiver
    mkdir tmp/su
}

error() {
    echo "${TEST_ID}:" "$@" > /dev/stderr
    exit 1
}

compile_code () {
    # Compile code.
    cd src

    # First invoke `make` in question mode to check if
    # `sofi3D` target is up to date.
    make sofi3D --question

    if [ $? -eq 0 ]
    then
        echo "Executable 'sofi3D' is up do date."
    else
        echo "Recompiling 'sofi3D'. See tmp/make.log for details"
        make sofi3D > ../tmp/make.log

        if [ "$?" -ne 0 ]; then
            cd ..
            echo "${TEST_ID}: FAIL" > /dev/stderr
            exit 1
        fi
    fi

    cd ..
}

animate_progress () {
    # Animate long-running tasks visually for usability.
    # USAGE: animate_progress task_id message
    # where
    #     task_id    process id
    #     message    message to show before the animated indicator
    #
    # Examples:
    #     $ animate_progress $! "Task is running"
    task_id="$1"
    message="$2"
    pause=0.2

    if [ "${CI}" ]
    then
        wait $task_id
    else
        while kill -0 "${task_id}" > /dev/null 2>&1
        do
            printf "\r"
            printf "%s |" "${message}"
            sleep $pause
            printf "\r"
            printf "%s /" "${message}"
            sleep $pause
            printf "\r"
            printf "%s -" "${message}"
            sleep $pause
            printf "\r"
            printf "%s \\" "${message}"
            sleep $pause
        done
        printf "\r"
        printf "%s\n" "${message}"
    fi
}

convert_segy_to_rsf() {
    # Convert SEG-Y format to RSF format.
    # USAGE: convert_segy_to_rsf path/filename.sgy
    #
    # Examples:
    #     $ convert_segy_to_rsf tmp/p.sgy
    # will produce RSF file tmp/p.rsf.

    file_sgy="$1"
    dir=$(dirname "$file_sgy")
    filename_base=$(basename -s .sgy "${file_sgy}")

    file_rsf="${dir}/${filename_base}.rsf"
    file_rsf_tfile="${dir}/${filename_base}.tfile.rsf"

    sfsegyread tape="$file_sgy" tfile="$file_rsf_tfile" > "$file_rsf"
}
