# Functions related to testing.

setup () {
    # Setup environment for the execution of a test.
    rm -rf tmp/
    mkdir tmp
    mkdir tmp/in_and_out
    mkdir tmp/log
    mkdir tmp/model
    mkdir tmp/snap
    mkdir tmp/sources
    mkdir tmp/su
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
    task_id=$1
    message=$2

    while kill -0 $task_id > /dev/null 2>&1 ; do
        echo -en "\r"
        echo -en "${message} |"
        sleep 0.5;
        echo -en "\r"
        echo -en "${message} /"
        sleep 0.5;
        echo -en "\r"
        echo -en "${message} -"
        sleep 0.5;
        echo -en "\r"
        echo -en "${message} \\"
        sleep 0.5;
    done
    echo -en "\r"
    echo "${message}"
}
