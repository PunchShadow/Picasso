# Looking for ColPack Library
#set(ColPack_DIR /opt/hsl
#						CACHE PATH "/opt/hsl")
set(ColPack_DIR /home/sferdou/opt/ColPack-bin 
						CACHE PATH "home/sferdou/opt/ColPack-bin")

find_library(ColPack_LIBRARIES libColPack.a
						 PATHS ${ColPack_DIR}/lib64/archive REQUIRED)

set(ColPack_INCDIR ${ColPack_DIR}/include)
set(ColPack_LIBDIR ${ColPack_DIR}/lib64/archive)
set(ColPack_LIBS -lColPack)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
	"ColPack"
	FOUND_VAR ColPack_FOUND
	REQUIRED_VARS ColPack_INCDIR ColPack_LIBDIR ColPack_LIBS
    FAIL_MESSAGE "Couldn't find ColPack"
	)
