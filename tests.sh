
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




OPTION=$1


	echo "*********************************************"
	echo "***   ANETO library Tests script"

if [ "${OPTION}" = "build" ]; then
	echo "***   Building the tests"
	echo

	mkdir build
	cd build
	cmake ../.

	OUT=$?
	check_cmake ${OUT}

	make tests

	DIRECTORY="$(pwd)"

	echo
	echo
	echo "*********************************************"
	echo "***   Test executable created at ${DIRECTORY}/tests"
	echo "***   For running the tests execute   '$0 run'"


elif [ "${OPTION}" = "run" ]; then
	if [[ -a "build/tests" ]];then
		echo "***   Running the tests"
		echo "***   For a complete testing check the options with './build/tests -h' "
		cd build
		./tests
	else
		echo "***   The executable does not exist."
		echo "***   Build it with '$0 build'"
	fi

elif [ "${OPTION}" = "clean" ]; then
	echo "***   Cleaning all the compilation files"
	rm -r build
else
	echo "***   Script help"
	echo "***"
	echo "***   For building the tests execute               '$0 build'"
	echo "***   For running the tests execute                '$0 run'"
	echo "***   For cleaning the compilation files execute   '$0 clean'"
fi


echo "*********************************************"





