# Try to find the EIGEN library
# See http://http://eigen.tuxfamily.org
#
# @author Lluis Gesa <gesa@ieec.cat>

include(CheckIncludeFileCXX)

CHECK_INCLUDE_FILE_CXX("eigen3/unsupported/Eigen/FFT" HAVE_EIGEN)

IF(NOT HAVE_EIGEN)
  UNSET(HAVE_EIGEN CACHE)
  message( FATAL_ERROR "\n\nEIGEN library not found, please download from http://eigen.tuxfamily.org/ \n\n

Follow the instructions to change install directories like \n
cmake . -DINCLUDE_INSTALL_DIR=/usr/local/include -DCMAKE_INSTALL_PREFIX=  \n" )
ENDIF()