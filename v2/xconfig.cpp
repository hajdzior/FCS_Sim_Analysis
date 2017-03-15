#include "config.h"

int main()
{
	ConfigFile cfg("config.cfg");

	bool exists = cfg.keyExists("blinking");
	std::cout << "blinking key: " << std::boolalpha << exists << "\n";

	std::string someValue = cfg.getValueOfKey<std::string>("mykey", "Unknown");
	std::cout << "value of key mykey: " << someValue << "\n";

	std::string blinking = cfg.getValueOfKey<std::string>("blinking");
	std::cout << "value of key blinking: " << blinking << "\n";

	double D = cfg.getValueOfKey<double>("D");
	std::cout << "value of key D: " << D << "\n\n";


	unsigned int n_particles = cfg.getValueOfKey<unsigned int>("n_particles");
    std::cout << "n_particles= " << n_particles << "\n";   
    
    std::cout << "blinking= " << D << "\n";    
    std::cout << "D= " << D << "\n";

	double t_max = cfg.getValueOfKey<double>("t_max");
	std::cout << "t_max= " << t_max << "\n";

	double w_xy = cfg.getValueOfKey<double>("w_xy");
	std::cout << "w_xy= " << w_xy << "\n";

	double w_z = cfg.getValueOfKey<double>("w_z");
	std::cout << "w_z= " << w_z << "\n";

	return 0;
}
