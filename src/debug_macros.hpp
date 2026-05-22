#ifndef DEBUG_MACRO_HPP
#define DEBUG_MACRO_HPP

#ifndef __DEBUG
#define __DEBUG 1 
#endif

#if __DEBUG 
#define DEBUG(s) {\
	std::cerr << "[" << __FILE__ <<  "][" << __FUNCTION__ << "][" << __LINE__ << "]: " << s << std::endl;\
}
#else 
#define DEBUG(s)
#endif


#endif // !DEBUG_MACRO_HPP
