// ===========================================================================
//	carb.r			  			PowerPlant 2.2.2		 �2000 Metrowerks Inc.
// ===========================================================================
//
//	Mac OS X looks for a 'carb' resource with ID 0 to determine whether to
//	launch a program under Classic or Carbon
//
//	Add this file to Carbon targets, but NOT to Classic targets

type 'carb' {
	longint = 0;				// Four zero bytes
};


resource 'carb'(0) {
};
