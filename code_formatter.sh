#!/bin/bash

#######################################

required_command="astyle"
code_directory="src/"

#######################################

usage() {
	echo
    echo -e "\tUsage: $(basename $0) [files]"
	echo
    echo -e "\tIf files are not specified, $(basename $0) formats all ".C" and ".H" files"
	echo -e "\tin source directory; otherwise, it formats all given files."
	echo
	echo -e "\tRequired command: $required_command"
	echo
	exit 0
}


[[ $1 == "-h" ]] && usage

# Test for required program
for comm in $required_command; do
	command -v $comm >/dev/null 2>&1 || {
		echo "I require $comm but it's not installed. Aborting." >&2;
		exit 1
	}
done

# Set the files to format
[[ $# != 0 ]] && src_files=$@ || src_files="--recursive $code_directory*.C,*.H"


# Here the important part: astyle formats the files.
astyle --style=bsd\
       --indent=spaces=4\
       --indent-classes\
       --indent-switches\
       --indent-col1-comments\
       --break-blocks\
       --pad-oper\
       --pad-comma\
       --pad-header\
       --delete-empty-lines\
       --align-pointer=type\
       --align-reference=type\
       --add-braces\
       --convert-tabs\
       --close-templates\
       --max-code-length=80\
       --mode=c\
       $src_files
