// Define library behaviour for MSVC
#ifdef _WIN32
	#ifdef Gripper_EXPORTS
		#define EXPORT __declspec(dllexport)
        #define IMPORT
	#else
		#define EXPORT __declspec(dllimport)
        #define IMPORT extern
	#endif
#else
	#define EXPORT
    #define IMPORT
#endif