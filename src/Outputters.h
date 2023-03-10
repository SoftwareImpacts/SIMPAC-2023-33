/////////////////////////////////////////////////////////////
/// @file Outputters.h
/// @brief OutputManager class to outstream simulation data
/////////////////////////////////////////////////////////////

#pragma once

// perform operations on the filesystem
#include <filesystem>

// output controlling with limits and timestep
#include <chrono>
#include <fstream>
#include <limits>

// project subfile header
#include "LatticePatch.h"

/** @brief Output Manager class to generate and coordinate output writing to
 * disk */
class OutputManager {
private:
  /// function to create the Code of the Simulations
  static std::string SimCodeGenerator();
  /// varible to safe the SimCode generated at execution
  std::string simCode;
  /// variable for the path to the output folder
  std::string Path;
  /// output style; csv or binary
  char outputStyle;
public:
  /// default constructor
  OutputManager();
  /// function that creates folder to save simulation data
  void generateOutputFolder(const std::string &dir);
  /// set the output style
  void set_outputStyle(const char _outputStyle);
  /// function to write data to disk in specified way
  void outUState(const int &state, const Lattice &lattice,
          const LatticePatch &latticePatch);
  /// simCode getter function
  [[nodiscard]] const std::string &getSimCode() const { return simCode; }
};

