#!/usr/bin/env fish

function rip_suffix
    set strIn $argv[1]
    if test (count (echo "$strIn" | grep "^\\s*..(/[^.]*)?\$")) -gt 0
        echo $strIn
    else if test (count (string match "?*.*" $strIn)) -gt 0
        echo $strIn | sed -r "s/^(.*)\\.[^.]*\$/\\1/"
    else
        echo $strIn
    end
end

function get_suffix
    set strIn $argv[1]
    if test (count (echo "$strIn" | grep "^\\s*..(/[^.]*)?\$")) -gt 0
        return
    end
    if test (count (string match "?*.*" $strIn)) -gt 0
        echo $strIn | sed -r "s/^.*(\\.[^.]*)\$/\\1/"
    end
end

set debug false

set pathInfile $argv[1]
if test $debug = true; echo pathInfile: $pathInfile; end
if test $debug = true; echo rip_suffix \$pathInfile: (rip_suffix $pathInfile); end
if test $debug = true; echo get_suffix \$pathInfile: (get_suffix $pathInfile); end
if test (count $argv) -ge 2
    set pathOutfile $argv[2]
else
    set pathOutfile ""(rip_suffix "$pathInfile")"_withComma"(get_suffix "$pathInfile")
end
if test $debug = true; echo pathOutfile: $pathOutfile; end

set listLines (grep "\\S" "$pathInfile")
if test "$pathInfile" = "-"
    echo -e (string join ",\n" $listLines)
else
    echo -e (string join ",\n" $listLines) > "$pathOutfile"
end
