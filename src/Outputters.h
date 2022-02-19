/////////////////////////////////////////////////////////////
/// @file Outputters.h
/// @brief OutputManager class to outstream simulation data
/////////////////////////////////////////////////////////////

// Include Guard
#ifndef OUTPUTTERS
#define OUTPUTTERS

// perform operations on the filesystem
#include <filesystem>

// output controlling with limits and timestep
#include <chrono>
#include <fstream>
#include <limits>

// project subfile header
#include "LatticePatch.h"

using namespace std;
namespace fs = std::filesystem;
namespace chrono = std::chrono;

/** @brief Output Manager class to generate and coordinate output writing to
 * disk */
class OutputManager {
private:
  /// function to create the Code of the Simulations
  static string SimCodeGenerator();
  /// varible to safe the SimCode generated at execution
  string simCode;
  /// variable for the path to the output folder
  string Path;
  /// process ID
  int myPrc;

public:
  /// default constructor
  OutputManager();
  /// function that creates folder to save simulation info
  void generateOutputFolder(const string &dir);
  /// output function for the whole lattice
  void outUState(const int &state, const LatticePatch &latticePatch);
  /// simCode getter function
  string getSimCode();
};

// End of Includeguard
#endif
