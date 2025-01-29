// mesh funcs and mesh types
#include "mesh.h"
#include "vector.h"
// CRS lib and CG solver
#include "CRSMat_types.h"
#include "CRSmatfuncs.h"
#include "CGSolver.h"
#include "CGSolver_types.h"
// General funcs
#include "common.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
// Operator
#include "gradient.h"

int main (void){
    int npoin, nelem, *elems;
    double *ptxyz;
    char path [] = {"temp/a06161.1.flds.zfem"};
    CHECK_ERROR(read_zfem(path,&npoin,&nelem,&ptxyz,&elems));
    if (ptxyz == NULL || elems == NULL)
    {
        fprintf(stderr,"Memory allocation (elems/ptxyz) failed.\n");
        return 1;
    }
    // created required data structure for mesh 
    int Nredge = 3;
    int *esurp,*esurp_ptr,*esure,*open;
    save_esurp(npoin,nelem,elems,&esurp,&esurp_ptr,Nredge);
    if (ptxyz == NULL || elems == NULL)
    {
        fprintf(stderr,"Memory allocation (esurp/esurp_ptr) failed.\n");
        return 1;
    }
    save_esure(nelem,elems,esurp_ptr,esurp,&esure,&open,Nredge);
    if (open == NULL || esure == NULL)
    {
        fprintf(stderr,"Memory allocation (esure/open) failed.\n");
        return 1;
    }
    int opencount = 0;
    for (int i=0;i<nelem;i++){
        if (open[i]==0) opencount++;
    }
    (opencount==0) ? printf("! this is open mesh.\n") : printf("* this is close mesh.\n");
    // calc norm of ele
    double *normele;
    CHECK_ERROR(save_normele(nelem,elems,ptxyz,&normele));
    // flip the normal vector to be outward:
    for (int ele = 0; ele < (3 * nelem); ele++)
        normele[ele] = -1 * normele[ele];
    // calculate the centeroid of each element
    double *cen;
    CHECK_ERROR(save_centri3(nelem,elems,ptxyz,&cen));
    if (cen == NULL){
        fprintf(stderr,"Memory allocation (cen) failed.\n");
        return 1;        
    }
    // calculate the center of each edge for each element
    double *cenedge;
    CHECK_ERROR(save_cenedgetri3(nelem,elems,ptxyz,&cenedge));
    // calculate the normal vectors of edges
    double *normedge;
    CHECK_ERROR(save_normedge(nelem,ptxyz,elems,normele,&normedge)); 
    if (normedge == NULL){
        fprintf(stderr,"Memory allocation (normedge) failed.\n");
        return 1;        
    }
    // find a boundary condition
    int *cell_stat = (int *)calloc((size_t)nelem,sizeof(int));
    char bcpath [] = {"/dagon1/achitsaz/poisson/temp/a06161.1.flds.zfem.labels"};
    FILE *fptr = fopen(bcpath,"r");
    if (fptr == NULL){
        fprintf(stderr,"BC file does not open properly\n");
        return 1;
    }
    int buffer = 50, nscan;
    char line[buffer];
    char *str;
    for (int ele =0;ele<nelem;ele++){
        str=edit_endline_character(line,buffer,fptr);
        nscan = sscanf(str, "%d", &cell_stat[ele]);
        if (nscan!=1){
            fprintf(stderr,"! there is error in the line %d of BC file.\n",ele+1);
            return 1;
        }
    }
    fclose(fptr);
    // creat mesh struct 
    mesh *M1 = (mesh *)malloc(sizeof(mesh));
    if (M1)
    {
        *M1 = (mesh){0}; // Set all integer and pointer fields to 0 or NULL
    }
    if (M1 == NULL)
    {
        fprintf(stderr, "Memory allocation failed for M1 pointer\n");
        exit(EXIT_FAILURE);
    }
    M1->elems =elems;M1->npoin =npoin;M1->nelem=nelem;M1->esure=esure;M1->esurp=esurp;M1->esurp_ptr=esurp_ptr;M1->nredge=Nredge;M1->open=open;M1->ptxyz=ptxyz;
    #ifdef DEBUG
        // coordinate 
            M1->numExtraPoints = nelem+3*nelem;
            double *extra_ptxyz = calloc((size_t)M1->numExtraPoints * 3, sizeof(double));
            // points for normal of elements
            for (int i = 0; i < (nelem* 3); i++)
                extra_ptxyz[i] = cen[i];
            // // point for edge p1->p2 
            // for (int i = 0; i < (nelem* 3); i++)
            //     extra_ptxyz[3*nelem+i] = cen[i];
            // // point for edge p2->p3    
            // for (int i = 0; i < (nelem* 3); i++)
            //     extra_ptxyz[6*nelem+i] = cen[i];
            // // point for edge p3->p1    
            // for (int i = 0; i < (nelem* 3); i++)
            //     extra_ptxyz[9*nelem+i] = cen[i]; 
            
            // point for edge p1->p2 
            for (int ele = 0; ele < nelem; ele++){
                for (int i=0;i<3;i++)
                extra_ptxyz[3*nelem+3*ele+i] = cenedge[9*ele+i];
            }
                
            // point for edge p2->p3    
            for (int ele = 0; ele < nelem; ele++){
                for (int i=0;i<3;i++)
                extra_ptxyz[6*nelem+3*ele+i] = cenedge[9*ele+3+i];
            }
            // point for edge p3->p1    
            for (int ele = 0; ele < nelem; ele++){
                for (int i=0;i<3;i++)
                extra_ptxyz[9*nelem+3*ele+i] = cenedge[9*ele+6+i];
            }

            M1->extra_ptxyz = extra_ptxyz;

            double *new_normele = calloc(((size_t)M1->npoin + (size_t)M1->numExtraPoints) * 3, sizeof(double));
            // normal vector of each element
            for (int i = 0; i < nelem; i++)
            {
                for (int j=0;j<3;j++)
                new_normele[3*M1->npoin + 3*i+j] = normele[3*i+j];
            }
            // normal vector of edge1 p1 -> p2
            for (int i = 0; i < nelem; i++)
            {
                for (int j=0;j<3;j++)
                new_normele[3*M1->npoin +3*nelem+ 3*i+j] = normedge[9*i+j];
            }
            // normal vector of edge1 p2 -> p3
            for (int i = 0; i < nelem; i++)
            {
                for (int j=0;j<3;j++)
                new_normele[3*M1->npoin +6*nelem+ 3*i+j] = normedge[9*i+3+j];
            }
            // normal vector of edge1 p3 -> p1
            for (int i = 0; i < nelem; i++)
            {
                for (int j=0;j<3;j++)
                new_normele[3*M1->npoin +9*nelem+ 3*i+j] = normedge[9*i+6+j];
            }
           FunctionWithArgs prtelefield[] =
            {
                {"open", 1, nelem, open, SCA_int_VTK},
                {"BC", 1, nelem, cell_stat, SCA_int_VTK},
            };
        size_t countele = sizeof(prtelefield) / sizeof(prtelefield[0]);
        FunctionWithArgs prtpntfield[] = {
            {"normal_vec", 3, (M1->npoin + M1->numExtraPoints), new_normele, VEC_double_VTK}
        };
        size_t countpnt = 1;
        CHECK_ERROR(SaveVTK("./", "checkmesh", 0, M1, tri3funcVTK, prtelefield, countele, prtpntfield, countpnt));
        free (new_normele);
        free (extra_ptxyz);
    #endif
    // define sparse matrix of coefficient
    int *row_ptr = (int *)calloc((size_t)(nelem + 1), sizeof(int));
    int max_nnz = nelem * (Nredge+1); // At most non-zero entries per row (element+neighbours)
    double *val = (double *)malloc((size_t)max_nnz * sizeof(double));
    int *col_ind = (int *)malloc((size_t)max_nnz * sizeof(int));
    int nnz = 0;
    double *RHS = (double *)calloc((size_t)nelem, sizeof(double));
    if (row_ptr == NULL || val == NULL || col_ind == NULL || RHS == NULL)
    {
        printf("Memory allocation failed.\n");
        return 1;
    }

    // calculate the symmetrical Positive Definite Coefficient martix and Right-Hand-Sided vector
    int order_in_cells [] = {0,1,1,2,2,0};
    for (int ele = 0; ele < nelem; ele++)
    {
        nnz++;
        int IDele = nnz - 1;
        // Known cells 
        if (cell_stat[ele]>0){
            // equation for this element is u=Tb
            val[IDele] = 1;
            col_ind[IDele] = ele;
            RHS[ele] = (cell_stat[ele] == 1) ? 1000 : 0; 
        }else{
        // Unknown cells 
            for (int nei = 0; nei < Nredge; nei++)
            {
                int neighbor = esure[Nredge * ele + nei]; 
                int lp1 = elems [Nredge * ele + order_in_cells[2*nei]]-1; 
                int lp2 = elems [Nredge * ele + order_in_cells[2*nei+1]]-1; 
                double dA = 0;
                dA = SQUARE((ptxyz [3 * lp1 ] - ptxyz [3 * lp2 ] ));          // x-coordinate  
                dA += SQUARE((ptxyz [3 * lp1 + 1] - ptxyz [3 * lp2 +1 ] ));   // y-coordinate
                dA += SQUARE((ptxyz [3 * lp1 + 2] - ptxyz [3 * lp2 +2 ] ));   // z-coordinate
                dA = sqrt (dA); 
                double dl = 0;
                dl = SQUARE(cen [Nredge * ele] - cen [Nredge *neighbor]);     // x-coordinate
                dl += SQUARE(cen [Nredge * ele + 1] - cen [Nredge *neighbor + 1]);     // y-coordinate
                dl += SQUARE(cen [Nredge * ele + 2] - cen [Nredge *neighbor + 2]);     // z-coordinate
                dl = sqrt (dl);
                double coef = dA/dl;
                // Internal flux  
                // contribution for diagonal element of coefficient matrix 
                val[IDele] += coef;
                col_ind[IDele] = ele;
                // contribution for off-diagonal element of coefficient matrix 
                if (cell_stat[neighbor] >0 ){
                    // between known and unknown cells
                    RHS [ele] += (cell_stat[neighbor] == 1) ? 1000 : 0;
                }else{
                    // between 2known cells 
                    val[nnz] = -1 * coef;
                    col_ind[nnz] = neighbor;
                    nnz++; // cause creat new element in the coefficient matrix
                }
            }
        }
        row_ptr[ele + 1] = nnz;
    }
    row_ptr[nelem] = nnz;

    // SOLVER SECTION
    CRSMatrix A;
    A.n = nelem;
    A.nnz = nnz;
    A.row_ptr = row_ptr;
    A.col_index = col_ind;
    A.values = val;

    // unknown vector
    double *u = (double *)calloc((size_t)A.n, sizeof(double));

    // Solve using CG
    clock_t start_time,end_time;
    double cpu_time_used;
    SolverConfig config = {100000, 1e-8, true};
    solver_set_config(config);
    // precond_conjugate_gradient(&A, RHS, u);
    start_time = clock();
    conjugate_gradient(&A, RHS, u);
    end_time = clock();
    cpu_time_used = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("* CG Solver execution time : %.2f seconds\n", cpu_time_used);

    // calculate the gradient scaler on each element
    double *norm_grad;
    CHECK_ERROR(gradient_ele_tri3(nelem,elems,ptxyz,esure,normedge,u,&norm_grad));
    if (norm_grad == NULL ){
        fprintf(stderr,"! ERROR: grad array is empty\n");
        return -1;
    }
    for (int ele=0;ele<nelem;ele++){
        double sum = 0;
        sum += SQUARE (norm_grad[3*ele]);
        sum += SQUARE (norm_grad[3*ele+1]);
        sum += SQUARE (norm_grad[3*ele+2]);
        sum = sqrt (sum);
        if (sum != 0){
        norm_grad [3*ele] /=  sum;
        norm_grad [3*ele+1] /=  sum;
        norm_grad [3*ele+2] /=  sum;
        }
    }



    FunctionWithArgs prtelefield2[] =
        {
            {"open", 1, nelem, open, SCA_int_VTK},
            {"BC", 1, nelem, cell_stat, SCA_int_VTK},
            {"u", 1, nelem, u, SCA_double_VTK},
            {"norm_grad_u", 3, nelem, norm_grad, VEC_double_VTK},
        };
    size_t countele2 = sizeof(prtelefield2) / sizeof(prtelefield2[0]);
    FunctionWithArgs prtpntfield2[] = {0};
    size_t countpnt2 = 0;
    CHECK_ERROR(SaveVTK("./", "checkmesh", 1, M1, tri3funcVTK, prtelefield2, countele2, prtpntfield2, countpnt2));





    // free dynamics arraies
    free (elems);
    free (ptxyz);
    free (open);
    free (esure);
    free (esurp);
    free (esurp_ptr);
    free (row_ptr);
    free (val);
    free (col_ind);
    free (RHS);
    free (cen);
    free (M1);
    free (cell_stat);
    free (u);
    free (cenedge);
    free (normele);
    free (normedge);
    free (norm_grad);
    return 0; // success signal
}