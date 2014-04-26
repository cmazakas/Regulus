/*

    Welcome to project Regulus! This is a poor man's version of Volker Springel's
    code Arepo and is quite the work in progress. But what is Arepo?
    Arepo is a physics simulator done using an unstructured mesh employing a 
    Vornoi tessellation to solve flux equations. This is the ultimate goal
    of Regulus as well, only a free and open-source implementation which we all
    know means it's obviously a pale imitation. Contributors are always welcome,
    however. Or at least, I haven't been able to find a download link for Arepo...
    This is subject to change.

    As of now, the code is a Delaunay triangulator and, unfortunately, currently
    lacks the ability to receive user-defined point sets. Luckily, this is one
    of the areas of future development.

    For Regulus, I would  currently love to add :

    1. Voronoi Mesh Generation
    2. User-defined Point Sets
    3. Arbitrary Precision Support

    Also, there are various bug fixes and optimizations that still need to
    happen. Namely, when Regulus is run on more than one thread, the point
    set is split across buffers and the largest problem is that each
    individual mesh formed is disparate in the sense that they are not linked
    with each other. This is most likely to be the first bug fix to be
    addressed as it would be nice to have one giant mesh instead of four 
    smaller ones. So for a complete mesh, run with "-np 1" only

    Things you will need to install Regulus :

    1. Eigen 3.0 and higher
    2. GMP
    3. MPFR
    4. C++11-compliant compiler

    I also have no idea if it'll Make on non-Unix based systems. But that's
    open-source for you, now isn't it?

*/

#include "structures.hpp"

double get_wall_time(void);

/*
    main loop
*/

int main(int argc, char **argv) {

    /*
        We choose to handle simple command line input
    */

    int num_threads(0), box_length(0);

    try {

        if (argc == 5) {

            if (strcmp(argv[1], "-np") == 0 && 
                strcmp(argv[3], "-bl") == 0) {

                sscanf(argv[2], "%d", &num_threads);
                sscanf(argv[4], "%d", &box_length);
            } else
                throw std::runtime_error("Incorrect arugments");
        } else
            throw std::runtime_error("Incorrect number of arguments");

    } catch(std::runtime_error &error) {

        std::cout << error.what() << std::endl;
          std::cout << "Use the form :: ./regulus -np X -bl Y \n";
        return -1;
    }

    /*
        ... And we then build a point set to triangulate...
    */

    unsigned long int num_points = std::pow(box_length, 3);
    std::vector<mesh> meshes;

    {

        std::vector<point> tmp_point_buff;
        tmp_point_buff.reserve(num_points);

        write_sorted_set(tmp_point_buff, num_points, box_length);

        auto m = num_points % num_threads;
        auto ppp = (num_points - m) / num_threads; 

        for (int i = 0; i < num_threads; ++i) {

            unsigned long int num_per_tetra = 0;

            if (i != num_threads - 1)
                num_per_tetra = ppp;
            else
                num_per_tetra = ppp + m;

            meshes.emplace_back(num_per_tetra);

            /*
                Allocate root points
            */

            auto tl = 3 * (box_length - 1);

            meshes[i].p_buffer.emplace_back(0, 0, 0);
            meshes[i].p_buffer.emplace_back(tl, 0, 0);
            meshes[i].p_buffer.emplace_back(0, tl, 0);
            meshes[i].p_buffer.emplace_back(0, 0, tl);

            meshes[i].p_position += 4;

            /*
                Allocate root tetrahedron
            */

            new(meshes[i].t_position) tetra(&meshes[i].p_buffer[0],
                                            &meshes[i].p_buffer[1],
                                            &meshes[i].p_buffer[2],
                                            &meshes[i].p_buffer[3]);

            meshes[i].m_position = meshes[i].t_position;
            ++meshes[i].t_position;

            /*
                Copy partitioned points
            */

            for (unsigned long int j = 0; j < num_per_tetra; ++j) {

                auto &cpy = tmp_point_buff[i*ppp + j];
                meshes[i].p_buffer.emplace_back(cpy.x, cpy.y, cpy.z);
            }

            if (verbose >= 3) {

                std::cout << "\nthread : " << i << std::endl;

                for (auto it = meshes[i].p_buffer.begin(); it < meshes[i].p_buffer.end(); ++it)
                    std::cout << *it << std::endl;
            }

            std::cout << std::endl;
        }
    }    

    /*
        This is the crux of regulus
        
        We begin the triangulation
    */

    double start_time = get_wall_time();

    std::vector<std::thread> threads;

    for (int i = 0; i < num_threads; ++i) {

        threads.emplace_back([&meshes, i](void)->void {

            std::vector<std::pair<tetra*, unsigned int>> leaves;
            leaves.reserve(512);

            std::vector<tetra*> pile;
            pile.reserve(2048);

            point *p = meshes[i].p_buffer.data();

            for (unsigned long int j = 5; j < meshes[i].p_buffer.size(); ++j) {

                if (verbose >= 1)
                    std::cout << "\niterative index : " << j - 4 << std::endl;

                walk(leaves, p + j, meshes[i].m_position, &*(meshes[i].t_buffer.begin()));
                fracture(meshes[i], leaves, pile, p + j);
                delaunay(meshes[i], pile);

                leaves.clear();
                pile.clear();
            }
        });
    }

    for (auto &t : threads)
        t.join();

    unsigned long int tot_num_tetra = 0;

    for (auto &m : meshes)
        tot_num_tetra += m.num_tetra;

    double end_time = get_wall_time();

    std::cout << "\nTotal number of tetrahedra triangulated : " << tot_num_tetra << std::endl;
    std::cout<< "\nTriangulated " << num_points << " points in " << end_time - start_time << " seconds" << std::endl;
    std::cout << num_points / ((end_time - start_time) * num_threads) << " points per second per thread" << std::endl;

    return 0;
}

/*
    Simple function found on the internet to return wall time (Unix only)
*/

double get_wall_time(void){

    struct timeval time;

    if (gettimeofday(&time, nullptr)){

        return 0;
    }
    return (double )time.tv_sec + (double )time.tv_usec * .000001;
}
