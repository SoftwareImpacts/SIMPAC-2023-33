/// @file main.cpp
/// @brief Main function to configure the user's simulation settings


#include "SimulationFunctions.h"   /* complete simulation functions and all headers */


int main(int argc, char *argv[])
{
    // Initialize MPI environment
    MPI_Init (&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    // Prepare MPI for Master-only threading
    //int provided;
    //MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    int rank = 0;
    MPI_Comm_rank(comm,&rank);
    double ti=MPI_Wtime(); // Overall start time

    /** Determine the output directory.                      
     * A "SimResults" folder will be created if non-existent 
     * with a subdirectory named in the identifier format    
     * "yy-mm-dd_hh-MM-ss" that contains the csv files     */
    constexpr auto outputDirectory = "/path/to/directory/";

    if(rank==0 && !filesystem::exists(outputDirectory)) {
        cerr<<"\nOutput directory nonexistent.\n";
        MPI_Abort(comm,1);
    }

    
    //------------ BEGIN OF CONFIGURATION ------------//
    
    ///////////////////// -- 1D -- ///////////////////////
    /** A 1D simulation with specified */
    /*
    //// Specify your settings here ////
    constexpr array <sunrealtype,2> CVodeTolerances={1.0e-16,1.0e-16}; /// - relative and absolute tolerances of the CVode solver
    constexpr int StencilOrder=13;                                     /// - accuracy order of the stencils in the range 1-13
    constexpr sunrealtype physical_sidelength=300e-6;                  /// - physical length of the lattice in meters
    constexpr sunindextype latticepoints=6e3;                          /// - number of lattice points
    constexpr bool periodic=true;                                      /// - periodic or vanishing boundary values
    int processOrder=3;                                                /// - included processes of the weak-field expansion, see README.md
    constexpr sunrealtype simulationTime=100.0e-6l;                    /// -  physical total simulation time
    constexpr int numberOfSteps=100;                                   /// -  discrete time steps
    constexpr int outputStep=100;                                      /// -  output step multiples
    constexpr char outputStyle='c';                                    /// - output in csv (c) or binary (b)

    /// Add electromagnetic waves.
    planewave plane1;                /// A plane wave with
    plane1.k = {1e5,0,0};            /// - wavevector (normalized to \f$ 1/\lambda \f$)
    plane1.p = {0,0,0.1};            /// - amplitude/polarization
    plane1.phi = {0,0,0};            /// - phase shift
    planewave plane2;                /// Another plane wave with
    plane2.k = {-1e6,0,0};           /// - wavevector (normalized to \f$ 1/\lambda \f$)
    plane2.p = {0,0,0.5};            /// - amplitude/polarization
    plane2.phi = {0,0,0};            /// - phase shift
    // Do not comment out this vector, even if no plane wave is used. But if, emplace used plane waves.
    vector<planewave> planewaves;
    //planewaves.emplace_back(plane1);
    //planewaves.emplace_back(plane2);
    
    gaussian1D gauss1;               /// A Gaussian wave with
    gauss1.k = {1.0e6,0,0};          /// - wavevector (normalized to \f$ 1/\lambda \f$)
    gauss1.p = {0,0,0.1};            /// - polarization/amplitude
    gauss1.x0 = {100e-6,0,0};        /// - shift from origin
    gauss1.phig = 5e-6;              /// - width
    gauss1.phi = {0,0,0};            /// - phase shift
    gaussian1D gauss2;               /// Another Gaussian with
    gauss2.k = {-0.2e6,0,0};         /// - wavevector (normalized to \f$ 1/\lambda \f$)
    gauss2.p = {0,0,0.5};            /// - polarization/amplitude
    gauss2.x0 = {200e-6,0,0};        /// - shift from origin
    gauss2.phig = 15e-6;             /// - width
    gauss2.phi = {0,0,0};            /// - phase shift
    // Do not comment out this vector, even if no Gaussian wave is used. But if, emplace used Gaussian waves.
    vector<gaussian1D> Gaussians1D;
    Gaussians1D.emplace_back(gauss1);
    Gaussians1D.emplace_back(gauss2);
  
    //// Do not change this below ////
    int *interactions = &processOrder;
    Sim1D(CVodeTolerances,StencilOrder,physical_sidelength,latticepoints,
            periodic,interactions,simulationTime,numberOfSteps,
            outputDirectory,outputStep,outputStyle,
            planewaves,Gaussians1D);
    */
    ////////////////////////////////////////////////////


    ///////////////////// -- 2D -- ///////////////////////
    /** A 2D simulation with specified */
    /*
    //// Specify your settings here ////
    constexpr array<sunrealtype,2> CVodeTolerances={1.0e-12,1.0e-12};  /// - relative and absolute tolerances of the CVode solver
    constexpr int StencilOrder=13;                                     /// - accuracy order of the stencils in the range 1-13
    constexpr array<sunrealtype,2> physical_sidelengths={80e-6,80e-6}; /// - physical length of the lattice in the given dimensions in meters
    constexpr array<sunindextype,2> latticepoints_per_dim={800,800};   /// - number of lattice points per dimension
    constexpr array<int,2> patches_per_dim={2,2};                      /// - slicing of discrete dimensions into patches
    constexpr bool periodic=true;                                      /// - periodic or vanishing boundary values
    int processOrder=3;                                                /// - included processes of the weak-field expansion, see README.md
    constexpr sunrealtype simulationTime=40e-6l;                       /// - physical total simulation time
    constexpr int numberOfSteps=100;                                   /// - discrete time steps
    constexpr int outputStep=100;                                      /// - output step multiples
    constexpr char outputStyle='c';                                    /// - output in csv (c) or binary (b)

    /// Add electromagnetic waves.
    planewave plane1;                  /// A plane wave with 
    plane1.k = {1e5,0,0};              /// - wavevector (normalized to \f$ 1/\lambda \f$)
    plane1.p = {0,0,0.1};              /// - amplitude/polarization
    plane1.phi = {0,0,0};              /// - phase shift
    planewave plane2;                  /// Another plane wave with
    plane2.k = {-1e6,0,0};             /// - wavevector
    plane2.p = {0,0,0.5};              /// - amplitude/polarization
    plane2.phi = {0,0,0};              /// - phase shift
    // Do not comment out this vector, even if no plane wave is used. But if, emplace used plane waves.
    vector<planewave> planewaves;
    //planewaves.emplace_back(plane1);
    //planewaves.emplace_back(plane2);

    gaussian2D gauss1;                 /// A Gaussian wave with
    gauss1.x0 = {40e-6,40e-6};         /// - center it approaches
    gauss1.axis = {1,0};               /// - normalized direction _from_ which the wave approaches the center
    gauss1.amp = 0.5;                  /// - amplitude
    gauss1.phip = 2*atan(0);           /// - polarization rotation from TE-mode (z-axis)
    gauss1.w0 = 2.3e-6;                /// - taille
    gauss1.zr = 16.619e-6;             /// - Rayleigh length
    /// the wavelength is determined by the relation \f$ \lambda = \pi*w_0^2/z_R \f$
    gauss1.ph0 = 2e-5;                 /// - beam center
    gauss1.phA = 0.45e-5;              /// - beam length
    gaussian2D gauss2;                 /// Another Gaussian wave with
    gauss2.x0 = {40e-6,40e-6};         /// - center it approaches
    gauss2.axis = {-0.7071,0.7071};    /// - normalized direction from which the wave approaches the center
    gauss2.amp = 0.5;                  /// - amplitude
    gauss2.phip = 2*atan(0);           /// - polarization rotation fom TE-mode (z-axis)
    gauss2.w0 = 2.3e-6;                /// - taille
    gauss2.zr = 16.619e-6;             /// - Rayleigh length
    gauss2.ph0 = 2e-5;                 /// - beam center
    gauss2.phA = 0.45e-5;              /// - beam length
    // Do not comment out this vector, even if no Gaussian wave is used. But if, emplace used Gaussian waves.
    vector<gaussian2D> Gaussians2D;
    Gaussians2D.emplace_back(gauss1);
    Gaussians2D.emplace_back(gauss2);

    //// Do not change this below ////
    static_assert(latticepoints_per_dim[0]%patches_per_dim[0]==0 &&
            latticepoints_per_dim[1]%patches_per_dim[1]==0,
            "The number of lattice points in each dimension must be "
            "divisible by the number of patches in that direction.");
    int * interactions = &processOrder;
    Sim2D(CVodeTolerances,StencilOrder,physical_sidelengths,
            latticepoints_per_dim,patches_per_dim,periodic,interactions,
            simulationTime,numberOfSteps,outputDirectory,outputStep,
            outputStyle,planewaves,Gaussians2D);
    */
    ////////////////////////////////////////////////////


    ///////////////////// -- 3D -- ///////////////////////
    /** A 3D simulation with specified */
    /*
    //// Specify your settings here ////
    constexpr array<sunrealtype,2> CVodeTolerances={1.0e-12,1.0e-12};        /// - relative and absolute tolerances of the CVode solver
    constexpr int StencilOrder=13;                                           /// - accuracy order of the stencils in the range 1-13
    constexpr array<sunrealtype,3> physical_sidelengths={80e-6,80e-6,20e-6}; /// - physical dimensions in meters
    constexpr array<sunindextype,3> latticepoints_per_dim={1200,1200,300};   /// - number of lattice points in any dimension
    constexpr array<int,3> patches_per_dim= {12,12,3};                       /// - slicing of discrete dimensions into patches
    constexpr bool periodic=true;                                            /// - perodic or non-periodic boundaries
    int processOrder=3;                                                      /// - processes of the weak-field expansion, see README.md
    constexpr sunrealtype simulationTime=40e-6;                              /// - physical total simulation time
    constexpr int numberOfSteps=40;                                          /// - discrete time steps
    constexpr int outputStep=10;                                             /// - output step multiples
    constexpr char outputStyle='b';                                          /// - output in csv (c) or binary (b)
    /// Add electromagnetic waves.
    planewave plane1;                   /// A plane wave with
    plane1.k = {1e5,0,0};               /// - wavevector (normalized to \f$ 1/\lambda \f$)
    plane1.p = {0,0,0.1};               /// - amplitude/polarization
    plane1.phi = {0,0,0};               /// - phase shift
    planewave plane2;                   /// Another plane wave with
    plane2.k = {-1e6,0,0};              /// - wavevector (normalized to \f$ 1/\lambda \f$)
    plane2.p = {0,0,0.5};               /// - amplitude/polarization
    plane2.phi = {0,0,0};               /// - phase shift
    // Do not comment out this vector, even if no plane wave is used. But if, emplace used plane waves.
    vector<planewave> planewaves;
    //planewaves.emplace_back(plane1);
    //planewaves.emplace_back(plane2);

    gaussian3D gauss1;                  /// A Gaussian wave with
    gauss1.x0 = {40e-6,40e-6,10e-6};    /// - center it approaches
    gauss1.axis = {1,0,0};              /// - normalized direction _from_ which the wave approaches the center
    gauss1.amp = 0.05;                  /// - amplitude
    gauss1.phip = 2*atan(0);            /// - polarization rotation from TE-mode (z-axis)
    gauss1.w0 = 3.5e-6;                 /// - taille
    gauss1.zr = 19.242e-6;              /// - Rayleigh length
    /// the wavelength is determined by the relation \f$ \lambda = \pi*w_0^2/z_R \f$
    gauss1.ph0 = 2e-5;                  /// - beam center
    gauss1.phA = 0.45e-5;               /// - beam length
    gaussian3D gauss2;                  /// Another Gaussian wave with
    gauss2.x0 = {40e-6,40e-6,10e-6};    /// - center it approaches
    gauss2.axis = {0,1,0};              /// - normalized direction from which the wave approaches the center
    gauss2.amp = 0.05;                  /// - amplitude
    gauss2.phip = 2*atan(0);            /// - polarization rotation from TE-mode (z-axis)
    gauss2.w0 = 3.5e-6;                 /// - taille
    gauss2.zr = 19.242e-6;              /// - Rayleigh length
    gauss2.ph0 = 2e-5;                  /// - beam center
    gauss2.phA = 0.45e-5;               /// - beam length
    // Do not comment out this vector, even if no Gaussian wave is used. But if, emplace used Gaussian waves.
    vector<gaussian3D> Gaussians3D;
    Gaussians3D.emplace_back(gauss1);
    Gaussians3D.emplace_back(gauss2);
    
    //// Do not change this below ////
    static_assert(latticepoints_per_dim[0]%patches_per_dim[0]==0 &&
            latticepoints_per_dim[1]%patches_per_dim[1]==0 &&
            latticepoints_per_dim[2]%patches_per_dim[2]==0,
            "The number of lattice points in each dimension must be "
            "divisible by the number of patches in that direction.");
    int *interactions = &processOrder;
    Sim3D(CVodeTolerances,StencilOrder,physical_sidelengths,
            latticepoints_per_dim,patches_per_dim,periodic,interactions,
            simulationTime,numberOfSteps,outputDirectory,outputStep,
            outputStyle,planewaves,Gaussians3D);
    */
    ////////////////////////////////////////////////////

    //------------- END OF CONFIGURATION -------------//

    double tf=MPI_Wtime(); // Overall finish time
    if(rank==0) {cout<<endl; timer(ti,tf);} // Print the elapsed time

    // Finalize MPI environment
    MPI_Finalize();

    return 0;
}
