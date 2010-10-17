/* Copyright (C) 2010 Romain Dubessy */
#ifndef OCTTREE_H
#define OCTTREE_H
#define TEST_CONST 1.0
/*! \brief This class implements an octree structure.
 *
 * In the resulting tree, each node has exactly eight children.
 * This structure is used to partition the tree dimensional space by recursively
 * subdviding it into octants.
 * This allow for efficient representation of a bunch of particules by grouping
 * them into clusters.
 *
 * Therefore it is usefull to implement a Barnes-Hut algorithm to
 * approximatively compute pair-wise interaction in coulombian systems, with 
 * complexity O(N*log(N)) (instead of the O(N^2) complexity of the full
 * problem).
 */
#ifdef VECTORIZE
typedef float v4sf __attribute__((vector_size(16)));
#endif //VECTORIZE
class octTree {
    public:
        /*! \brief Default constructor for an empty octTree. */
        octTree(void);
        /*! \brief Construct an octTree from an array of charges. */
        octTree(int, double *, double *);
        /*! \brief Destructor to keep low memory usage. */
        ~octTree(void);
        /*! \brief Method to insert a charge in the tree. */
        void insert(double *, double *, int);
        /*! \brief Initialize the pointers for ordered tree walking. */
        void flatten(void);
        /*! \brief Sets the local nodes values, as barycenter. */
        void finalize(void);
        /*!\brief Prints an ascii version of the tree on the screen. */
        void display(void);
        /*! \brief Method to compute the energy. */
        int energy_grad(double *, double &);
        /*! \brief Method to compute the stepsize. */
        int guess(double &, double &);
#ifdef VECTORIZE
        /*! \brief Overloaded new operator for memory alignement requirements. 
         * 
         * This trick ensures that the addresses of class octtree objects are 16
         * bytes aligned in memory.
         * This is required by some of the vector operators.
         */
        void *operator new(size_t size) {
            void *p;
            posix_memalign(&p,16,size);
            return p;
        };
        /*! \brief Overloaded new array operator for memory alignement
         * requirements.
         */
        void *operator new[](size_t size) { return operator new(size); };
        /*! \brief Overloaded insert method. */
        void insert(v4sf&, v4sf&, int);
#endif //VECTORIZE

    private:
#ifdef VECTORIZE
        /*!\brief Node's charges barycenter coordinates. */
        v4sf mean __attribute__((aligned(16)));
        v4sf gmean __attribute__((aligned(16)));
        v4sf gcmean __attribute__((aligned(16)));
        v4sf gcmeanb __attribute__((aligned(16)));
        float center[3];  /*!<\brief Node's cell center coordinates. */
        float size;       /*!<\brief  Node's cell size. */
        float g2mean;
        float delta;      /*!<\brief Node's center-barycenter distance. */
#else
        double mean[3];   /*!<\brief Node's charges barycenter coordinates. */
        double gmean[3];  /*!<\brief Node's gradient barycenter. */
        /*!\brief Auxiliary node's gradient two point averages. */
        double gcmean[6];
        double center[3]; /*!<\brief Node's cell center coordinates. */
        double size;      /*!<\brief Node's cell size. */
        double g2mean;    /*!<\brief Auxiliary node's gradient average. */
        double delta;     /*!<\brief Node's center-barycenter distance. */
#endif //VECTORIZE
        octTree *next;    /*!<\brief Pointer for tree walking. */
        octTree *skip;    /*!<\brief Pointer for tree walking. */
        octTree *child;   /*!<\brief Array of children. */
        int n_q;          /*!<\brief Number of charges at the node. */
        int indice;       /*!<\brief ID of charge stored in a leaf. */
        int level;        /*!<\brief Level of the node in the tree structure. */
#ifdef VECTORIZE
} __attribute__((aligned(16)));
#else
};
#endif //VECTORIZE
#endif //OCTTREE_H
/* octtree.h */
