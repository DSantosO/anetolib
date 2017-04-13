
function check_cmake {
	if [[ $1 != "0" ]];then
		echo
		echo
		echo
		echo "*********************************************"
		echo "*********************************************"
		echo "***   cmake building failed!!! "
		echo "***   Please look carefully its output and check the README file for possible solutions."
		echo "*********************************************"
		exit
	fi
}

function check_make_doc {
	if [[ $1 != "0" ]];then
		echo
		echo
		echo
		echo "*********************************************"
		echo "*********************************************"
		echo "***   make doc failed!!! "
		echo "***   Please check if doxygen is installed."
		echo "***   If this is not the problem look carefully its output and check the README file for possible solutions."
		echo "*********************************************"
		exit
	fi
}

function check_latex {
	if [[ $1 != "0" ]];then
		echo
		echo
		echo
		echo "*********************************************"
		echo "*********************************************"
		echo "***   LaTeX compilation failed!!! "
		echo "***   Please check if you have all the LaTeX packages needed"
		echo "***   If this is not the problem look carefully its output and check the README file for possible solutions."
		echo "*********************************************"
		exit
	fi
}




OPTION=$1
ANETO_PATH="$(pwd)"


	echo "*********************************************"
	echo "***   ANETO library Documentation script"

	mkdir build
	cd build
	cmake ../.

	OUT=$?
	check_cmake ${OUT}

	make doc

	OUT=$?
	check_make_doc ${OUT}

	DIRECTORY="$(pwd)"
	echo
	echo
	echo "*********************************************"
	echo "***   Documentation created succesfully at ${ANETO_PATH}/doc"
	echo "***   HTML files can be checked at ${ANETO_PATH}/doc/html/index.html"

	cd ../doc/latex/
	make
	OUT=$?
	check_latex ${OUT}

	echo "*********************************************"
	echo "***   ANETO library Documentation script"
	echo "***   Documentation created succesfully at ${ANETO_PATH}/doc"
	echo "***   HTML files can be checked at ${ANETO_PATH}/doc/html/index.html"
	echo "***   LaTeX pdf file can be found at ${ANETO_PATH}/doc/latex/refman.pdf"
	echo "*********************************************"

