#include "structures.hpp"

void write_sorted_set(std::vector<point> &buff,
					  unsigned long int num_points,
					  int box_length) {

	for (unsigned long int i = 0; i < num_points; ++i) {

		buff.emplace_back(i / (box_length * box_length),
					      (i / box_length) % box_length,
						  i % box_length);

	}

	std::sort(buff.begin(), buff.end(), [](const point &p1, const point &p2) {
											return p1.peanokey < p2.peanokey;
										});

	if (verbose >= 3) {

		std::cout << "\nProposed point set : " << std::endl;

		for (unsigned long int i = 0; i < num_points; ++i) {

			std::cout << buff[i] << std::endl;
		}
	}

	return;
}
