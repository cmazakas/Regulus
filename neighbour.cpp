#include "structures.hpp"

/*

    This is a simple, brute-force approach to linking two tetrahedra.

*/

int link_tetra(tetra *t1, tetra *t2) {

    assert(t1 && t2 && t1 != t2);

    int sum1 = 0;
    int sum2 = 0;
    int matches = 0;

    for (unsigned int i = 0; i < 4; ++i) {
        for (unsigned int j = 0; j < 4; ++j) {
            if (t1->p[i]->peanokey == t2->p[j]->peanokey) {
                sum1 += i;
                sum2 += j;
                ++matches;
            }
        }
    }

    if (matches >= 4)
        throw std::runtime_error("Identical tetrahedra created in the mesh\n");

    if (matches == 3) {

        switch (sum1) {
        
            case 6 : // 0 - 321
            
                switch (sum2) {
                
                    case 6 : // 321
                    
                        t1->ngb[0] = t2;
                        t2->ngb[0] = t1;
                        
                        t1->nid[0] = 0;
                        t2->nid[0] = 0;
                        
                        break;
                        
                    case 5 : // 023
                    
                        t1->ngb[0] = t2;
                        t2->ngb[1] = t1;
                        
                        t1->nid[0] = 1;
                        t2->nid[1] = 0;
                        
                        break;
                        
                    case 4 : // 031
                    
                        t1->ngb[0] = t2;
                        t2->ngb[2] = t1;
                        
                        t1->nid[0] = 2;
                        t2->nid[2] = 0;
                    
                        break;
                        
                    case 3 : // 012
                    
                        t1->ngb[0] = t2;
                        t2->ngb[3] = t1;
                        
                        t1->nid[0] = 3;
                        t2->nid[3] = 0;
                
                        break;
            }
                
            break;
            
            case 5 : // 1 - 023
            
                switch (sum2) {
            
                    case 6 : // 321
                    
                        t1->ngb[1] = t2;
                        t2->ngb[0] = t1;
                        
                        t1->nid[1] = 0;
                        t2->nid[0] = 1;
                        
                        break;
                        
                    case 5 : // 023
                    
                        t1->ngb[1] = t2;
                        t2->ngb[1] = t1;
                        
                        t1->nid[1] = 1;
                        t2->nid[1] = 1;
                        
                        break;
                        
                    case 4 : // 031
                    
                        t1->ngb[1] = t2;
                        t2->ngb[2] = t1;
                        
                        t1->nid[1] = 2;
                        t2->nid[2] = 1;
                    
                        break;
                        
                    case 3 : // 012
                    
                        t1->ngb[1] = t2;
                        t2->ngb[3] = t1;
                        
                        t1->nid[1] = 3;
                        t2->nid[3] = 1;
                    
                        break;
                }
            
                break;
                
            case 4 : // 2 - 031
            
                switch (sum2) {
        
                    case 6 : // 321
                    
                        t1->ngb[2] = t2;
                        t2->ngb[0] = t1;
                        
                        t1->nid[2] = 0;
                        t2->nid[0] = 2;
                        
                        break;
                        
                    case 5 : // 023
                    
                        t1->ngb[2] = t2;
                        t2->ngb[1] = t1;
                        
                        t1->nid[2] = 1;
                        t2->nid[1] = 2;
                    
                        break;
                        
                    case 4 : // 031
                    
                        t1->ngb[2] = t2;
                        t2->ngb[2] = t1;
                        
                        t1->nid[2] = 2;
                        t2->nid[2] = 2;
                
                        break;
                        
                    case 3 : // 012
                    
                        t1->ngb[2] = t2;
                        t2->ngb[3] = t1;
                        
                        t1->nid[2] = 3;
                        t2->nid[3] = 2;
                    
                        break;
                }
            
                break;
                
            case 3 : // 3 - 012
        
                switch (sum2) {
        
                    case 6 : // 321
                    
                        t1->ngb[3] = t2;
                        t2->ngb[0] = t1;
                        
                        t1->nid[3] = 0;
                        t2->nid[0] = 3;
                        
                        break;
                        
                    case 5 : // 023
                
                        t1->ngb[3] = t2;
                        t2->ngb[1] = t1;
                        
                        t1->nid[3] = 1;
                        t2->nid[1] = 3;
                        
                        break;
                        
                    case 4 : // 031
                
                        t1->ngb[3] = t2;
                        t2->ngb[2] = t1;
                        
                        t1->nid[3] = 2;
                        t2->nid[2] = 3;
                
                        break;
                        
                    case 3 : // 012
                    
                        t1->ngb[3] = t2;
                        t2->ngb[3] = t1;
                        
                        t1->nid[3] = 3;
                        t2->nid[3] = 3;
                    
                        break;
                }
                
                break;        
        }
            
        if (verbose >= 1) {
        
            std::cout << "link established between : " <<
                      "(" << t1 << ", " << t2 << ")" << std::endl;
        }
                    
        return 0;
    }

    return 1;
}
