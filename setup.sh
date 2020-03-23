if [[ -v BASH_SOURCE ]]; then
	THIS_FOLDER=`dirname $BASH_SOURCE[0]`
else
	THIS_FOLDER=`dirname $0`
fi
THIS_FOLDER=`readlink -f ${THIS_FOLDER}`

export LD_LIBRARY_PATH=${THIS_FOLDER}/lib:${LD_LIBRARY_PATH}
