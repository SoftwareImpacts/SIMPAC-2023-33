/// @file main.cpp
/// @brief Main function to configure the user's simulation settings


#include "SimulationFunctions.h"   /* complete simulation functions and all headers */

using namespace std;

int main(int argc, char *argv[])
{
    /** Determine the output directory.                      
     * A "SimResults" folder will be created if non-existent 
     * with a subdirectory named in the identifier format    
     * "yy-mm-dd_hh-MM-ss" that contains the csv files      */
    string outputDirectory = argv[1];  // command line

    if(!filesystem::exists(outputDirectory)) {
        cerr<<"\nOutput directory nonexistent.\n";
        exit(1);
    }

    // Initialize MPI environment
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    
    //------------ BEGIN OF CONFIGURATION ------------//
    
    if (stoi(argv[2]) == 1) {
    ///////////////////// -- 1D -- ///////////////////////
    /** A 1D simulation with specified */
    array <sunrealtype,2> CVodeTolerances={stod(argv[3]),stod(argv[4])}; /// - relative and absolute tolerances of the CVode solver
    int StencilOrder=stoi(argv[5]);  /// - accuracy order of the stencils in the range 1-13
    sunrealtype physical_sidelength=stod(argv[6]);  /// - physical length of the lattice in meters
    sunindextype latticepoints=stoi(argv[9]);  /// - number of lattice points
    constexpr bool periodic=true;  /// - periodic or vanishing boundary values
    int processOrder=stoi(argv[15]);  /// - included processes of the weak-field expansion, see README.md
    sunrealtype simulationTime=stod(argv[16]);  /// -  physical total simulation time
    int numberOfSteps=stoi(argv[17]);  /// -  discrete time steps
    int outputStep=stoi(argv[18]);  /// -  output step multiples
    char outputStyle=*(argv[19]); /// - output in csv (c) or binary (b) format

    /// Add electromagnetic waves.
    vector<planewave> planewaves;
    vector<gaussian1D> Gaussians1D;
    if (stoi(argv[20])) {  // if use of plane waves
        planewave plane1;  /// A plane wave with
        plane1.k = {stod(argv[22]),stod(argv[23]),stod(argv[24])};  /// - wavevector (normalized to \f$ 1/\lambda \f$)
        plane1.p = {stod(argv[25]),stod(argv[26]),stod(argv[27])};  /// - amplitude/polarization
        plane1.phi = {stod(argv[28]),stod(argv[29]),stod(argv[30])};  /// - phase shift
        planewave plane2;  /// Another plane wave with
        plane2.k = {stod(argv[31]),stod(argv[32]),stod(argv[33])};  /// - wavevector (normalized to \f$ 1/\lambda \f$)
        plane2.p = {stod(argv[34]),stod(argv[35]),stod(argv[36])};  /// - amplitude/polarization
        plane2.phi = {stod(argv[37]),stod(argv[38]),stod(argv[39])};  /// - phase shift
        if (stoi(argv[21])>2) {
            cerr<<"You can only use up to two plane waves." << endl;
            exit(1);
        }
        if (stoi(argv[21])>0) planewaves.emplace_back(plane1);
        if (stoi(argv[21])==2) planewaves.emplace_back(plane2);
    }
    if (stoi(argv[40])) {  // if use of Gaussians
        gaussian1D gauss1;  /// A Gaussian wave with
        gauss1.k = {stod(argv[42]),stod(argv[43]),stod(argv[44])};  /// - wavevector (normalized to \f$ 1/\lambda \f$)
        gauss1.p = {stod(argv[45]),stod(argv[46]),stod(argv[47])};  /// - polarization/amplitude
        gauss1.x0 = {stod(argv[48]),stod(argv[49]),stod(argv[50])};  /// - shift from origin
        gauss1.phig = stod(argv[51]);  /// - width
        gauss1.phi = {stod(argv[52]),stod(argv[53]),stod(argv[54])};  /// - phase shift
        gaussian1D gauss2;  /// Another Gaussian with
        gauss2.k = {stod(argv[55]),stod(argv[56]),stod(argv[57])};  /// - wavevector (normalized to \f$ 1/\lambda \f$)
        gauss2.p = {stod(argv[58]),stod(argv[59]),stod(argv[60])};  /// - polarization/amplitude
        gauss2.x0 = {stod(argv[61]),stod(argv[62]),stod(argv[63])};  /// - shift from origin
        gauss2.phig = stod(argv[64]);  /// - width
        gauss2.phi = {stod(argv[65]),stod(argv[66]),stod(argv[67])};  /// - phase shift
        if (stoi(argv[41])>2) {
            cerr<<"You can only use up to two Gaussian pulses." << endl;
            exit(1);
        }
        if (stoi(argv[41])>0) Gaussians1D.emplace_back(gauss1);
        if (stoi(argv[41])==2) Gaussians1D.emplace_back(gauss2);
    }
    //// Do not change this below ////
    int *interactions = &processOrder;
    Sim1D(CVodeTolerances,StencilOrder,physical_sidelength,latticepoints,
            periodic,interactions,simulationTime,numberOfSteps,
            outputDirectory,outputStep,outputStyle,
            planewaves,Gaussians1D);
    }
    ////////////////////////////////////////////////////


    if (stoi(argv[2]) == 2) {
    ///////////////////// -- 2D -- ///////////////////////
    /** A 2D simulation with specified */
    array<sunrealtype,2> CVodeTolerances={stod(argv[3]),stod(argv[4])};  /// - relative and absolute tolerances of the CVode solver
    int StencilOrder=stoi(argv[5]);  /// - accuracy order of the stencils in the range 1-13
    array<sunrealtype,2> physical_sidelengths={stod(argv[6]),stod(argv[7])};  /// - physical length of the lattice in the given dimensions in meters
    array<sunindextype,2> latticepoints_per_dim={stoi(argv[9]),stoi(argv[10])};  /// - number of lattice points per dimension
    array<int,2> patches_per_dim={stoi(argv[12]),stoi(argv[13])};  /// - slicing of discrete dimensions into patches
    constexpr bool periodic=true;  /// - periodic or vanishing boundary values
    int processOrder=stoi(argv[15]);  /// - included processes of the weak-field expansion, see README.md
    sunrealtype simulationTime=stod(argv[16]);  /// - physical total simulation time
    int numberOfSteps=stoi(argv[17]);  /// - discrete time steps
    int outputStep=stoi(argv[18]);  /// - output step multiples
    char outputStyle=*(argv[19]);  /// - output in csv (c) or binary (b) format

    /// Add electromagnetic waves.
    vector<planewave> planewaves;
    vector<gaussian2D> Gaussians2D;
    if (stoi(argv[20])) {  // if use of plane waves
        planewave plane1;  /// A plane wave with 
        plane1.k = {stod(argv[22]),stod(argv[23]),stod(argv[24])};  /// - wavevector (normalized to \f$ 1/\lambda \f$)
        plane1.p = {stod(argv[25]),stod(argv[26]),stod(argv[27])};  /// - amplitude/polarization
        plane1.phi = {stod(argv[28]),stod(argv[29]),stod(argv[30])};  /// - phase shift
        planewave plane2;  /// Another plane wave with
        plane2.k = {stod(argv[31]),stod(argv[32]),stod(argv[33])};  /// - wavevector
        plane2.p = {stod(argv[34]),stod(argv[35]),stod(argv[36])};  /// - amplitude/polarization
        plane2.phi = {stod(argv[37]),stod(argv[38]),stod(argv[39])};  /// - phase shift
        if (stoi(argv[21])>2) {
            cerr<<"You can only use up to two plane waves." << endl;
            exit(1);
        }
        if (stoi(argv[21])>0) planewaves.emplace_back(plane1);
        if (stoi(argv[21])==2) planewaves.emplace_back(plane2);
    }
    if (stoi(argv[68])) {  // if use of Gaussians
        gaussian2D gauss1;  /// A Gaussian wave with
        gauss1.x0 = {stod(argv[70]),stod(argv[71])};  /// - center it approaches
        gauss1.axis = {stod(argv[73]),stod(argv[74])};  /// - normalized direction _from_ which the wave approaches the center
        gauss1.amp = stod(argv[75]);  /// - amplitude
        gauss1.phip = atan(stod(argv[76]));  /// - polarization rotation from TE-mode (z-axis)
        gauss1.w0 = stod(argv[77]);  /// - taille
        gauss1.zr = stod(argv[78]);  /// - Rayleigh length
        /// the wavelength is determined by the relation \f$ \lambda = \pi*w_0^2/z_R \f$
        gauss1.ph0 = stod(argv[79]);  /// - beam center
        gauss1.phA = stod(argv[80]);  /// - beam length
        gaussian2D gauss2;  /// Another Gaussian wave with
        gauss2.x0 = {stod(argv[81]),stod(argv[82])};  /// - center it approaches
        gauss2.axis = {stod(argv[84]),stod(argv[85])};  /// - normalized direction from which the wave approaches the center
        gauss2.amp = stod(argv[86]);  /// - amplitude
        gauss2.phip = atan(stod(argv[87]));  /// - polarization rotation fom TE-mode (z-axis)
        gauss2.w0 = stod(argv[88]);  /// - taille
        gauss2.zr = stod(argv[89]);  /// - Rayleigh length
        gauss2.ph0 = stod(argv[90]);  /// - beam center
        gauss2.phA = stod(argv[91]);  /// - beam length
        if (stoi(argv[69])>2) {
            cerr<<"You can only use up to two Gaussian pulses." << endl;
            exit(1);
        }
        if (stoi(argv[69])>0) Gaussians2D.emplace_back(gauss1);
        if (stoi(argv[69])==2) Gaussians2D.emplace_back(gauss2);
    }
    //// Do not change this below ////
    int * interactions = &processOrder;
    Sim2D(CVodeTolerances,StencilOrder,physical_sidelengths,
            latticepoints_per_dim,patches_per_dim,periodic,interactions,
            simulationTime,numberOfSteps,outputDirectory,outputStep,
            outputStyle,planewaves,Gaussians2D);
    }
    ////////////////////////////////////////////////////


    if (stoi(argv[2]) == 3) {
    ///////////////////// -- 3D -- ///////////////////////
    /** A 3D simulation with specified */
    array<sunrealtype,2> CVodeTolerances={stod(argv[3]),stod(argv[4])};  /// - relative and absolute tolerances of the CVode solver
    int StencilOrder=stoi(argv[5]);  /// - accuracy order of the stencils in the range 1-13
    array<sunrealtype,3> physical_sidelengths={stod(argv[6]),stod(argv[7]),stod(argv[8])}; /// - physical dimensions in meters
    array<sunindextype,3> latticepoints_per_dim={stoi(argv[9]),stoi(argv[10]),stoi(argv[11])};  /// - number of lattice points in any dimension
    array<int,3> patches_per_dim={stoi(argv[12]),stoi(argv[13]),stoi(argv[14])};  /// - slicing of discrete dimensions into patches
    constexpr bool periodic=true;  /// - perodic or non-periodic boundaries
    int processOrder=stoi(argv[15]);  /// - processes of the weak-field expansion, see README.md
    sunrealtype simulationTime=stod(argv[16]);  /// - physical total simulation time
    int numberOfSteps=stoi(argv[17]);  /// - discrete time steps
    int outputStep=stoi(argv[18]);  /// - output step multiples
    char outputStyle=*(argv[19]);  /// - output in csv (c) or binary (b) format

    /// Add electromagnetic waves.
    vector<planewave> planewaves;
    vector<gaussian3D> Gaussians3D;
    if (stoi(argv[20])) {  // if use of plane waves
        planewave plane1;  /// A plane wave with 
        plane1.k = {stod(argv[22]),stod(argv[23]),stod(argv[24])};  /// - wavevector (normalized to \f$ 1/\lambda \f$)
        plane1.p = {stod(argv[25]),stod(argv[26]),stod(argv[27])};  /// - amplitude/polarization
        plane1.phi = {stod(argv[28]),stod(argv[29]),stod(argv[30])};  /// - phase shift
        planewave plane2;  /// Another plane wave with
        plane2.k = {stod(argv[31]),stod(argv[32]),stod(argv[33])};  /// - wavevector
        plane2.p = {stod(argv[34]),stod(argv[35]),stod(argv[36])};  /// - amplitude/polarization
        plane2.phi = {stod(argv[37]),stod(argv[38]),stod(argv[39])};  /// - phase shift
        if (stoi(argv[21])>2) {
            cerr<<"You can only use up to two plane waves." << endl;
            exit(1);
        }
        if (stoi(argv[21])>0) planewaves.emplace_back(plane1);
        if (stoi(argv[21])==2) planewaves.emplace_back(plane2);
    }
    if (stoi(argv[68])) {  // if use of Gaussians
        gaussian3D gauss1;  /// A Gaussian wave with
        gauss1.x0 = {stod(argv[70]),stod(argv[71]),stod(argv[72])};  /// - center it approaches
        gauss1.axis = {stod(argv[73]),stod(argv[74]),0};  /// - normalized direction _from_ which the wave approaches the center
        gauss1.amp = stod(argv[75]);  /// - amplitude
        gauss1.phip = atan(stod(argv[76]));  /// - polarization rotation from TE-mode (z-axis)
        gauss1.w0 = stod(argv[77]);  /// - taille
        gauss1.zr = stod(argv[78]);  /// - Rayleigh length
        /// the wavelength is determined by the relation \f$ \lambda = \pi*w_0^2/z_R \f$
        gauss1.ph0 = stod(argv[79]);  /// - beam center
        gauss1.phA = stod(argv[80]);  /// - beam length
        gaussian3D gauss2;  /// Another Gaussian wave with
        gauss2.x0 = {stod(argv[81]),stod(argv[82]),stod(argv[83])};  /// - center it approaches
        gauss2.axis = {stod(argv[84]),stod(argv[85]),0};  /// - normalized direction from which the wave approaches the center
        gauss2.amp = stod(argv[86]);  /// - amplitude
        gauss2.phip = atan(stod(argv[87]));  /// - polarization rotation from TE-mode (z-axis)
        gauss2.w0 = stod(argv[88]);  /// - taille
        gauss2.zr = stod(argv[89]);  /// - Rayleigh length
        gauss2.ph0 = stod(argv[90]);  /// - beam center
        gauss2.phA = stod(argv[91]);  /// - beam length
        if (stoi(argv[69])>2) {
            cerr<<"You can only use up to two Gaussian pulses." << endl;
            exit(1);
        }
        if (stoi(argv[69])>0) Gaussians3D.emplace_back(gauss1);
        if (stoi(argv[69])==2) Gaussians3D.emplace_back(gauss2);
    }
    //// Do not change this below ////
    int *interactions = &processOrder;
    Sim3D(CVodeTolerances,StencilOrder,physical_sidelengths,
            latticepoints_per_dim,patches_per_dim,periodic,interactions,
            simulationTime,numberOfSteps,outputDirectory,outputStep,
            outputStyle,planewaves,Gaussians3D);
    }
    ////////////////////////////////////////////////////

    //------------- END OF CONFIGURATION -------------//

    // Finalize MPI environment
    MPI_Finalize();

    return 0;
}
