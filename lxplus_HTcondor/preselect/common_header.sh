function escapeDoublequotes {
    echo "$1" | sed 's/\\/\\\\/g' | sed 's/"/\\"/g'
}

function wrapDryRun {
    if (( dryRun )); then
        echo "$@"
    else
        "$@"
    fi
}

declare -a queuedArgsArray=()
# shellcheck disable=SC2034 # nonflag parameters isn't required
declare -a nonflagArgsArray=()

function shiftArgs {
    local -i nShift=1
    if [[ "$#" -gt 0 ]]; then
        nShift="$1"
    fi
    queuedArgsArray=( "${queuedArgsArray[@]:$nShift}" )
}

function manageShorthands {
    local strRaw="${queuedArgsArray[0]:1}"
    shiftArgs
    local -a subArgs=()
    local -i i=0;
    # Parameterized expansions of for loop
    # doesn't work with 'set -u' (throw unbound variable error)
    # {1..2} doesn't work with variables as the upper bound
    # So use 'seq' from 'coreutils'
    for i in $(seq 0 $((${#strRaw} - 1))); do
        subArgs+=( "-${strRaw:$i:1}" )
    done
    queuedArgsArray=( "${subArgs[@]}" "${queuedArgsArray[@]}" )
}

queuedArgsArray=( "$@" )
